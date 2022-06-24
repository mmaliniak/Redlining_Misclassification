#Purpose: Classify area level data (blocks, block groups, tracts) using Meier approach (i.e., Weighted score method)
#         which assigns each area a historical redlining score

#Reference to Meier method: https://www.openicpsr.org/openicpsr/project/141121/version/V2/view 

#Applying to: Four MTAs with HOLC maps in Georgia but could run for any city with HOLC map

#Last updated: 4-24-2022


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
meier_fun <- function(city, HOLCmap, counties, yr, cutoff){
  
  #read in HOLC map
  Red <- HOLCmap
  

  meier_calc <- function(areadf){
      
      #GEOID column name is different depending on year so need to rename so always the same regardless of year pulling
      #NOTE: this is clunky AF but can't figure out better way to do this for sf (if do it normal way - ends up making it into a list with geometry)
      geoid_col <- colnames(areadf)[grepl("GEOID",colnames(areadf))]
      GEOID2 <- areadf[geoid_col]
      st_geometry(GEOID2) <- NULL
      GEOID2$GEOID_NEW <- GEOID2[,geoid_col]
      areadf <- left_join(areadf, GEOID2, by=geoid_col)
      
      #set crs to be the same 
      areadf <- st_transform(areadf, crs = 26916)
      Red <- st_transform(Red, crs = 26916)
      
      # Calculate area and tidy up
      intersect_pct <- st_intersection(areadf, Red) %>% 
        mutate(intersect_area = st_area(.)) %>%   # create new column with shape area
        dplyr::select(GEOID_NEW, holc_id, holc_grade, intersect_area) %>%   # only select columns needed to merge
        st_drop_geometry()  # drop geometry as we don't need it
      
      # Create a fresh area variable (will give slightly different results than ALAND provided by CB (even for areas without water))
      areadf <- mutate(areadf, poly_area = st_area(areadf))
      
      # Merge by census area (block, block group, tract)
      areadf2 <- merge(areadf, intersect_pct, by = "GEOID_NEW", all.x = TRUE)
      
      # Calculate coverage
      areadf2 <- areadf2 %>% 
        mutate(poly_prop = as.numeric(intersect_area/poly_area))
      
      #Calculate historical redlining score 
      #sum up poly_prop
      areadf2 <- areadf2 %>%
        group_by(GEOID_NEW) %>% #group by census area
        mutate(poly_prop=ifelse(!holc_grade %in% c("A","B","C","D"),NA,poly_prop), #sets holc grades other than A-D (e.g., E grade in Savannah to NA)
               sum_poly_prop = sum(poly_prop, na.rm=TRUE)) %>% 
        ungroup()
      
      #Exclude census areas where < cutoff of area was graded
      areadf2 <- mutate(areadf2, 
                          Exclude = ifelse(sum_poly_prop > 0 & sum_poly_prop < cutoff, 1, #some overlap but < cutoff
                                           ifelse(sum_poly_prop==0, 2, 0))) #no overlap
      
      #create weighted score according to Meier method
      
      #convert holc grades to numeric (A=1, B=2, C=3, D=4)
      #create weights
      areadf2 <- mutate(areadf2,
                          holc_num = ifelse(holc_grade=="A",1,
                                            ifelse(holc_grade=="B",2,
                                                   ifelse(holc_grade=="C",3,
                                                          ifelse(holc_grade=="D",4,NA)))), #if not A-D (e.g, E in Savannah) then set to NA (i.e., ungraded)
                          hrs_wt = (poly_prop/sum_poly_prop)*holc_num)
      
      #sum over weights to get final score
      areadf2 <- areadf2 %>%
        group_by(GEOID_NEW) %>% #group by census area
        mutate(hrs_score = round(sum(hrs_wt, na.rm=TRUE),2)) %>% 
        ungroup()
      
      
      #remove duplicates (those with multiple HOLC areas in them)
      areadf2 <- areadf2[!duplicated(areadf2$GEOID_NEW), ]
      
      #set hrs_score of excluded census areas to NA
      areadf2$hrs_score <- ifelse(areadf2$Exclude!=0,NA,areadf2$hrs_score)
      
      #categorize hrs score (rounded to the nearest grade)
      areadf2 <- mutate(areadf2,
                          hrs_cat = ifelse(is.na(hrs_score),9,
                                           ifelse(hrs_score < 1.5,1,
                                                  ifelse(hrs_score < 2.5,2,
                                                         ifelse(hrs_score < 3.5,3, 
                                                                ifelse(hrs_score>=3.5,4,9))))),
                          hrs_cat = factor(hrs_cat, c(1:4,9), c("1-<1.5","1.5-<2.5","2.5-<3.5","3.5-4","Excluded")))
      
      #exclude areas with no overlap
      areadf2 <- filter(areadf2, Exclude %in% c(0,1))
      
      #keep columns that we need
      areadf2 <- select(areadf2, GEOID_NEW, sum_poly_prop, Exclude, hrs_score, hrs_cat)
      
      return(areadf2)
      
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
  
  
  blocks_meier <- meier_calc(blocksdf)
  
  
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
  
  
  bgs_meier <- meier_calc(bgs_df)
  
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
  
  
  tracts_meier <- meier_calc(tractsdf)
  
  
  ###############
  # Output     #
  ##############
  
  cutoff2 <- as.character(cutoff*100)
  
  meier_data <- list("blocks_meier"=blocks_meier, "bgs_meier"=bgs_meier, "tracts_meier"=tracts_meier)
  saveRDS(meier_data, file = paste0("H:/Redlining and obesity/Data/Redlining/Redlining area data/Meier/", city, "/", "Cutoff ", cutoff2, "/", city, "_meier", yr,".RDS"))

}
  
#############################
# Function specifications
##############################

# HOLCmap - specify HOLC map to be used
# counties - specify counties that overlap with HOLC map to reduce number of areas pulled in 
# year - specify year (2010, 2020) 
# cutoff - specify cutoff for exclusion based on proportion of area ungraded



#######
# ATL #
#######

#2020
ATLcounties <- c("063", "089","121") #Clayon, DeKalb, and Fulton
cutoff <- c(.10,.20,.30,.40,.50)
map(.x = cutoff, 
    .f = meier_fun, 
    city = "Atlanta", 
    HOLCmap = ATL_red, 
    counties=ATLcounties,
    yr="2020")


