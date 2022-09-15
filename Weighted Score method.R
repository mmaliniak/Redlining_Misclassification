#Purpose: Classify area level data (blocks, block groups, tracts) using Meier approach (i.e., Weighted score method)
#         which assigns each area a historical redlining score

#Reference: https://www.openicpsr.org/openicpsr/project/141121/version/V2/view 

#Applying to Atlanta HOLC map but could run for any city with HOLC map

#Last updated: 9-15-2022


library(tidyverse)
library(sf)        # for handling simple features spatial data
library(tmap)      # for producing thematic maps
library(RColorBrewer)
library(tidycensus)
library(tigris)
library(purrr)


#read in ATL HOLC map
ATL_red <- read_sf('H:/Redlining and obesity/Data/Redlining/Atlanta/cartodb-query.shp') 


#function to calculate historical redlining score by area (blocks, block groups, tracts)
weighted_fun <- function(city, HOLCmap, st, counties, yr, cutoff){
  
  #read in HOLC map
  Red <- HOLCmap
  
  #pull in water polygons for counties that overlap HOLC areas
  water <- area_water(state=st, 
                      county = counties, 
                      year = yr) %>% st_as_sf()

  weighted_calc <- function(areadf){
      
      #GEOID column name is different depending on year so need to rename so always the same regardless of year pulling - this is clunky bc sf
      geoid_col <- colnames(areadf)[grepl("GEOID",colnames(areadf))]
      GEOID2 <- areadf[geoid_col]
      st_geometry(GEOID2) <- NULL
      GEOID2$GEOID_NEW <- GEOID2[,geoid_col]
      areadf <- left_join(areadf, GEOID2, by=geoid_col)

      #set crs to be the same - these are for GA so change for other areas
      areadf <- st_transform(areadf, crs = 26916)
      Red <- st_transform(Red, crs = 26916)
      water <- st_transform(water, crs = 26916)
      
      #restrict to census areas that overlap HOLC map (reduce the amount of time function takes to run)
      areadf2 <- areadf[Red,]
      water <- water %>% .[areadf2,] %>% st_union()
      
      # Subtract water polygons from total area to get land area for each census area
      landarea <- st_difference(areadf2, water)
      
      # Calculate land area (m^2) 
      landarea2 <- mutate(landarea, poly_area = st_area(landarea))

      # Calculate HOLC area within each census area and tidy up
      # Get warning message that attribute values are assumed to be spatially constant throughout all geometries
      intersect_pct <- st_intersection(landarea2, Red) %>% 
        mutate(intersect_area = st_area(.)) %>%   # create new column with shape area (m^2)
        dplyr::select(GEOID_NEW, holc_id, holc_grade, intersect_area) %>%   # only keep necessary columns
        st_drop_geometry()  # drop geometry as we don't need it
      
      # Merge HOLC areas with census areas (typically there are multiple HOLC areas within a census area)
      areadf3 <- merge(landarea2, intersect_pct, by = "GEOID_NEW", all.x = TRUE)
      
      # Calculate coverage of each HOLC area within each census area
      areadf3 <- areadf3 %>% 
        mutate(poly_prop = as.numeric(intersect_area/poly_area))
      
      # Sum up proportion HOLC-graded for each census area
      areadf3 <- areadf3 %>%
        group_by(GEOID_NEW) %>% #group by census area
        mutate(poly_prop=ifelse(!holc_grade %in% c("A","B","C","D"),NA,poly_prop), #sets holc grades other than A-D (e.g., E grade in Savannah to NA)
               sum_poly_prop = sum(poly_prop, na.rm=TRUE)) %>% 
        ungroup() #ungroup so still have total number of HOLC areas 
      
      #Exclude census areas where < cutoff of census area was graded (this is 20% in original Weighted Score approach)
      areadf3 <- mutate(areadf3, 
                          Exclude = ifelse(sum_poly_prop < cutoff, 1,0)) 
      
      #create weighted score
      
      #convert holc grades to numeric (A=1, B=2, C=3, D=4)
      #create weights
      areadf3 <- mutate(areadf3,
                          holc_num = ifelse(holc_grade=="A",1,
                                            ifelse(holc_grade=="B",2,
                                                   ifelse(holc_grade=="C",3,
                                                          ifelse(holc_grade=="D",4,NA)))), #if not A-D (e.g, E in Savannah) then set to NA (i.e., ungraded)
                          hrs_wt = (poly_prop/sum_poly_prop)*holc_num)
      
      #sum over weights to get final score
      areadf4 <- areadf3 %>%
        group_by(GEOID_NEW) %>% #group by census area
        mutate(hrs_score = round(sum(hrs_wt, na.rm=TRUE),2)) %>% 
        ungroup() #ungroup so still have total number of HOLC areas 

      #remove duplicates (those with multiple HOLC areas in them)
      areadf4 <- areadf4[!duplicated(areadf4$GEOID_NEW), ]
      
      #set hrs_score of excluded census areas to NA
      areadf4$hrs_score <- ifelse(areadf4$Exclude!=0,NA,areadf4$hrs_score)
      
      #categorize hrs score (rounded to the nearest grade)
      areadf4 <- mutate(areadf4,
                          hrs_cat = ifelse(is.na(hrs_score),9,
                                           ifelse(hrs_score < 1.5,1,
                                                  ifelse(hrs_score < 2.5,2,
                                                         ifelse(hrs_score < 3.5,3, 
                                                                ifelse(hrs_score>=3.5,4,9))))),
                          hrs_cat = factor(hrs_cat, c(1:4,9), c("1-<1.5","1.5-<2.5","2.5-<3.5","3.5-4","Excluded")))

      #keep columns that we need
      areadf4 <- select(areadf4, GEOID_NEW, sum_poly_prop, Exclude, hrs_score, hrs_cat)
      
      return(areadf4)
      
    }

  ################
  #Census blocks #
  ################
  
  #bring in census block data
  options(tigris_use_cache = T)
  blocksdf <- blocks(
    state="GA",
    county = counties, 
    class = 'sf',
    year = yr)
  
  
  blocks_weighted <- weighted_calc(blocksdf)
  
  
  ######################
  #Census block groups #
  ######################
  
  #bring in census block group data
  options(tigris_use_cache = T)
  bgs_df <- block_groups(
    state="GA",
    county = counties, 
    class = 'sf',
    year = yr)
  
  
  bgs_weighted <- weighted_calc(bgs_df)
  
  ######################
  #Census tracts       #
  ######################
  
  #bring in census tract data
  options(tigris_use_cache = T)
  tractsdf <- tracts(
    state="GA",
    county = counties, 
    class = 'sf',
    year = yr)
  
  
  tracts_weighted <- weighted_calc(tractsdf)
  
  
  ###############
  # Output     #
  ##############
  
  cutoff2 <- as.character(cutoff*100)
  
  weighted_data <- list("blocks_weighted"=blocks_weighted, "bgs_weighted"=bgs_weighted, "tracts_weighted"=tracts_weighted)
  saveRDS(weighted_data, file = paste0("H:/Redlining area data/Weighted Score/", city, "/", "Cutoff ", cutoff2, "/", city, "_weighted", yr,".RDS"))

}
  
#############################
# Function specifications
##############################

# City - specify city
# HOLCmap - specify HOLC map to be used
# st - specifiy state
# counties - specify counties that overlap with HOLC map to reduce number of areas pulled in 

        #get counties that intersect with Atlanta HOLC map
        counties(
          state="GA",
          class = 'sf',
          year = 2020) %>%
          st_as_sf() %>% 
          st_transform(4326) %>% #set crs to be the same as HOLC map
          .[ATL_red,] %>% #HOLC map
          select(COUNTYFP, NAME)

# year - specify year (2010, 2020) 
# cutoff - specify cutoff for exclusion based on proportion of area ungraded



#######
# ATL #
#######

#2020
ATLcounties <- c("063", "089","121") #Clayton, DeKalb, and Fulton
cutoff <- c(.10,.20,.30,.40,.50)
purrr::map(.x = cutoff, 
    .f = weighted_fun, 
    city = "Atlanta", 
    HOLCmap = ATL_red,
    st = "GA",
    counties=ATLcounties,
    yr="2020")


