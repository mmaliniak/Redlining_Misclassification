#Purpose: Classify area level data (blocks, block groups, tracts) using modified Majority Land Area method
#         using different cut-offs for % ungraded (10,20,30,40,50%)

#Reference(s): Krieger N, et al. Am J Epidemiol. 2020 Oct 1;189(10):1065–75.
#              Krieger N, et al. Am J Public Health. 2020 Jul;110(7):1046–53.


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


#function to assign redlining status by area (blocks, block groups, tracts, zipcodes)
modmajority_fun <- function(city, HOLCmap, st, counties, yr, cutoff){
  
  #read in HOLC map
  Red <- HOLCmap
  
  #pull in water polygons for counties that overlap HOLC areas
  water <- area_water(state=st, 
                      county = counties, 
                      year = yr) %>% st_as_sf()
  
  
  majority_calc <- function(areadf){
    
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
    
    # Merge HOLC areas with census areas
    areadf3 <- merge(landarea2, intersect_pct, by = "GEOID_NEW", all.x = TRUE)
    
    # Calculate coverage of each HOLC area within each census area
    areadf3 <- areadf3 %>% 
      mutate(poly_prop = as.numeric(intersect_area/poly_area),
             poly_prop = ifelse(!holc_grade %in% c("A","B","C","D"),NA,poly_prop)) #sets holc grades other than A-D (e.g., E grade in Savannah to NA))
    
    #sum up poly_prop by HOLC grade to get % HOLC grade by census area
    areadf4 <- areadf3 %>%
      group_by(GEOID_NEW, holc_grade) %>% #group by geoid and holc grade
      summarise(sum_poly_prop = sum(poly_prop, na.rm=TRUE)) %>% 
      st_set_geometry(NULL) %>%
      spread(holc_grade, sum_poly_prop, fill = 0) #this will give proportion of tract in A, B, C, and D grades as columns
    
    #Modified Majority Land Area method: 
    # % ungraded and threshold for assigning HOLC grade is the same but if more than one HOLC grade is
    #above threshold then assign the one with the greatest proportion within the area (e.g., if cutoff is 10% then 
    #the area must have at least 10% graded and one grade must make up at least 10% of the area - if more than 1 HOLC grade exceeds 10% then go with max)
    areadf4 <- mutate(areadf4,
                      
                      #sum across to get total proportion graded
                      sum_poly_prop = A + B + C + D,
                      
                      #Flag areas that don't have any overlap with HOLC map (=2) or areas with some overlap but < cutoff (=1)
                      Exclude = ifelse(sum_poly_prop > 0 & sum_poly_prop < cutoff, 1, #1=some overlap but < cutoff
                                       ifelse(sum_poly_prop==0, 2, 0)), #2=no overlap, 0=some overlap and meets cutoff
                      
                      #Flag areas that satisfy cutoff for assigning HOLC grades
                      HOLC_cut = ifelse(A >=cutoff | B >=cutoff | C >=cutoff | D >= cutoff, 1, 0),
                      
                      #For those that satisfy cutoffs then find HOLC grade that accounts for the highest proportion of the area
                      HOLC_max = ifelse(HOLC_cut==1, pmax(A,B,C,D), NA),
                      
                      #Assign HOLC grade - in the event that multiple grades were equal then would assign highest grade (A > B > C > D)
                      #don't think there are any instances where the proportions are exactly the same but just to have the code
                      HOLC_grade_majority = ifelse(is.na(HOLC_max),9, #No grade meets requirement
                                                  ifelse(HOLC_max==A,1, #A
                                                         ifelse(HOLC_max==B,2, #B
                                                                ifelse(HOLC_max==C,3, #C
                                                                       ifelse(HOLC_max==D,4,9))))), #D
                      #factor
                      HOLC_grade_majority=factor(HOLC_grade_majority,
                                                c(1:4,9),
                                                c("A","B","C","D","No grade")))

    
    #add geometry back for mapping
    areadf5 <- left_join(areadf4, areadf, by="GEOID_NEW")
    
    #keep columns that we need
    areadf5 <- select(areadf5, GEOID_NEW, A, B, C, D, sum_poly_prop, Exclude, HOLC_cut, HOLC_max, HOLC_grade_majority, geometry)
    
    #convert back to sf 
    areadf5 <- st_as_sf(areadf5)
    
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
  
  blocks_majority <- majority_calc(blocksdf)
  
  
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
  
  
  bgs_majority <- majority_calc(bgs_df)
  
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
  
  
  tracts_majority <- majority_calc(tractsdf)
  

  ###############
  # Output     #
  ##############
  
  cutoff2 <- as.character(cutoff*100)
  
  majority_data <- list("blocks_majority"=blocks_majority, "bgs_majority"=bgs_majority, "tracts_majority"=tracts_majority)
  
  saveRDS(majority_data, file = paste0("H:/Redlining area data/Majority Land Area/", city, "/", "Cutoff ", cutoff2, "/", city, "_modmajority", yr,".RDS"))
  
}

#############################
# Function specifications
##############################

# HOLCmap - specify HOLC map to be used
# st - specify state
# counties - specify counties that overlap with HOLC map to reduce number of areas pulled in 

        #get counties that intersect with HOLC map
                counties(
                    state="GA",
                    class = 'sf',
                    year = 2020) %>%
                    st_as_sf() %>% 
                    st_transform(4326) %>% 
                    .[ATL_red,] %>% 
                    select(COUNTYFP, NAME)

# year - specify year (2010, 2020) 


#######
# ATL #
#######


#2020
ATLcounties <- c("063", "089","121") #Clayon, DeKalb, and Fulton
cutoff <- c(.10,.20,.30,.40,.50)
purrr::map(.x = cutoff, 
    .f = modmajority_fun, 
    city = "Atlanta", 
    st = "GA",
    HOLCmap = ATL_red, 
    counties=ATLcounties,
    yr="2020")
