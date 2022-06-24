#Purpose: Classify area level data (blocks, block groups, tracts, zipcodes) using highest HOLC approach 
#         which assigns highest HOLC grade within an area

#Reference(s): Li M, Yuan F. Race Soc Probl. 2021 Jun 18;1–16

#Applying to: Four MTAs with HOLC maps in Georgia but could run for any city with HOLC map

#Last updated: 4-24-2022


library(tidyverse)
library(sf)        # for handling simple features spatial data
library(tmap)      # for producing thematic maps
library(RColorBrewer)
library(tidycensus)
library(tigris)


#read in ATL HOLC map
ATL_red <- read_sf('H:/Redlining and obesity/Data/Redlining/Atlanta/cartodb-query.shp') 


#function to assign redlining status by area (blocks, block groups, tracts, zipcodes)
highest_fun <- function(city, HOLCmap, counties, yr, cutoff){
  
  #read in HOLC map
  Red <- HOLCmap
  
  
  highest_calc <- function(areadf){
    
    #GEOID column name is different depending on year so need to rename so always the same regardless of year pulling
    #NOTE: this is clunky AF but can't figure out better way to do this for sf (if do it normal way - ends up making it into a list with geometry)
    geoid_col <- colnames(areadf)[grepl("GEOID",colnames(areadf))]
    GEOID2 <- areadf[geoid_col]
    st_geometry(GEOID2) <- NULL
    GEOID2$GEOID_NEW <- GEOID2[,geoid_col]
    areadf <- left_join(areadf, GEOID2, by=geoid_col)
    
    #set crs to be the same 
    areadf2 <- st_transform(areadf, crs = 26916)
    Red <- st_transform(Red, crs = 26916)
    
    # Calculate area and tidy up
    intersect_pct <- st_intersection(areadf2, Red) %>% 
      mutate(intersect_area = st_area(.)) %>%   # create new column with shape area
      dplyr::select(GEOID_NEW, holc_id, holc_grade, intersect_area) %>%   # only select columns needed to merge
      st_drop_geometry()  # drop geometry as we don't need it
    
    # Create a fresh area variable (will give slightly different results than ALAND provided by CB (even for areas without water))
    areadf2 <- mutate(areadf2, poly_area = st_area(areadf2))
    
    # Merge by census area (block, block group, tract, or zcta)
    areadf3 <- merge(areadf2, intersect_pct, by = "GEOID_NEW", all.x = TRUE)
    
    # Calculate coverage
    areadf3 <- areadf3 %>% 
      mutate(poly_prop = as.numeric(intersect_area/poly_area),
             poly_prop = ifelse(!holc_grade %in% c("A","B","C","D"),NA,poly_prop)) #sets holc grades other than A-D (e.g., E grade in Savannah to NA))
    
    #sum up poly_prop
    areadf4 <- areadf3 %>%
      group_by(GEOID_NEW) %>% #group by census area
      mutate(sum_poly_prop = sum(poly_prop, na.rm=TRUE)) %>% 
      ungroup()
    
    #Exclude census areas where <=10% of area was graded
    areadf4 <- mutate(areadf4, 
                      Exclude = ifelse(sum_poly_prop > 0 & sum_poly_prop < cutoff, 1, #some overlap but < cutoff
                                       ifelse(sum_poly_prop==0, 2, 0))) #no overlap
    
    #sort by HOLC grade from highest to lowest
    areadf4 <- arrange(areadf4, holc_grade)

    #take first row per census area to get highest HOLC grade in area
    areadf5 <- areadf4 %>% group_by(GEOID_NEW) %>% filter(row_number(GEOID_NEW) == 1)
    
    #categorize
    areadf5 <- mutate(areadf5,
                        holc_grade_highest = ifelse(Exclude %in% c(1,2),"Exclude", holc_grade))
    
    #exclude areas with no overlap
    areadf5 <- filter(areadf5, Exclude %in% c(0,1))
    
    #keep columns that we need
    areadf5 <- select(areadf5, GEOID_NEW, sum_poly_prop, Exclude, holc_grade_highest)
    
    return(areadf5)
    
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
  
  blocks_highest <- highest_calc(blocksdf)
  
  
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
  
  
  bgs_highest <- highest_calc(bgs_df)
  
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
  
  
  tracts_highest <- highest_calc(tractsdf)
  
  
  ###############
  # Output     #
  ##############
  
  cutoff2 <- as.character(cutoff*100)
  
  highest_data <- list("blocks_highest"=blocks_highest, "bgs_highest"=bgs_highest, "tracts_highest"=tracts_highest)
  saveRDS(highest_data, file = paste0("H:/Redlining and obesity/Data/Redlining/Redlining area data/Highest HOLC/", city, "/", "Cutoff ", cutoff2, "/", city, "_highest", yr,".RDS"))
  
}

#############################
# Function specifications
##############################

# city - specifiy city 
# HOLCmap - specify HOLC map to be used
# counties - specify counties that overlap with HOLC map to reduce number of areas pulled in 
# year - specify year (2010, 2020) NOTE: zctas are always 2010 because 2020 data are unavailable 
# cutoff - specify cutoff for exclusion based on proportion of area ungraded

#Note: will get warnings about "attribute variables are assumed to be spatially constant throughout all geometries" - safe to ignore


#######
# ATL #
#######

#2020
ATLcounties <- c("063", "089","121") #Clayon, DeKalb, and Fulton
cutoff <- c(.10,.20,.30,.40,.50)
map(.x = cutoff, 
    .f = highest_fun, 
    city = "Atlanta", 
    HOLCmap = ATL_red, 
    counties=ATLcounties,
    yr="2020")



