#Purpose: Classify area level data (blocks, block groups, tracts) using Majority Land Area method
#         which assigns Assign HOLC grade of majority

#Reference(s): Krieger N, et al. Am J Epidemiol. 2020 Oct 1;189(10):1065–75.
#              Krieger N, et al. Am J Public Health. 2020 Jul;110(7):1046–53.


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


#function to assign redlining status by area (blocks, block groups, tracts)
krieger_fun <- function(HOLCmap, counties, yr){
  
  #read in HOLC map
  Red <- HOLCmap
  
  
  krieger_calc <- function(areadf){
    
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
    
    #restrict to blocks that overlap HOLC map (reduce the amount of time function takes to run)
    areadf2 <- st_join(areadf, Red, join = st_intersects, left=FALSE) #only keep spatial matches
    
    #remove duplicate blocks (those with multiple HOLC areas in them)
    areadf2 <- areadf2[!duplicated(areadf2$GEOID_NEW), ]
    
    #remove columns from HOLC df
    areadf2 <- select(areadf2, -name, -holc_id, -holc_grade)
    
    # Calculate area and tidy up
    intersect_pct <- st_intersection(areadf2, Red) %>% 
      mutate(intersect_area = st_area(.)) %>%   # create new column with shape area
      dplyr::select(GEOID_NEW, holc_id, holc_grade, intersect_area) %>%   # only select columns needed to merge
      st_drop_geometry()  # drop geometry as we don't need it
    
    # Create a fresh area variable (will give slightly different results than ALAND provided by CB (even for areas without water))
    areadf2 <- mutate(areadf2, poly_area = st_area(areadf2))
    
    # Merge by census area (block, block group, tract)
    areadf3 <- merge(areadf2, intersect_pct, by = "GEOID_NEW", all.x = TRUE)
    
    # Calculate coverage
    areadf3 <- areadf3 %>% 
      mutate(poly_prop = as.numeric(intersect_area/poly_area),
             poly_prop = ifelse(!holc_grade %in% c("A","B","C","D"),NA,poly_prop)) #sets holc grades other than A-D (e.g., E grade in Savannah to NA))
    
    #sum up poly_prop
    areadf4 <- areadf3 %>%
      group_by(GEOID_NEW, holc_grade) %>% #group by geoid and holc grade
      summarise(sum_poly_prop = sum(poly_prop, na.rm=TRUE)) %>% 
      st_set_geometry(NULL) %>%
      spread(holc_grade, sum_poly_prop, fill = 0) #this will give proportion of tract in A, B, C, and D grades as columns
    
    
    #Identify tracts where 100% of tract is contained in 1 HOLC area
    areadf4 <- mutate(areadf4, 
                   Containment = ifelse(A==1 | B==1 | C==1 | D==1,1,
                                        ifelse( (A < 1 & A >=0.5) | 
                                                  (B < 1 & B >=0.5) | 
                                                  (C < 1 & C >=0.5) | 
                                                  (D < 1 & D >= 0.5),2,
                                                ifelse( (A + B + C + D) >=0.5,3,4))),
                   Containment = factor(Containment, c(1,2,3,4), c("Fully contained", "Majority contained", "Mixed", "No grade")),
                   HOLC_grade_Krieger = ifelse(A >=0.5,1,
                                               ifelse(B >=0.5,2,
                                                      ifelse(C >=0.5,3,
                                                             ifelse(D >= 0.5,4,
                                                                    ifelse(Containment=="Mixed",8,9))))),
                   HOLC_grade_Krieger=factor(HOLC_grade_Krieger,
                                             c(1:4,8,9),
                                             c("50+% A","50+% B","50+% C","50+% D","50+% Mixed", "<50%: No grade")),
                   
                   #collapse Mixed and no grade together
                   HOLC_grade_Krieger2 = ifelse(A>=0.5,1,
                                                ifelse(B>=0.5, 2,
                                                       ifelse(C >=0.5, 3,
                                                              ifelse(D >= 0.5,4,
                                                                     ifelse(Containment %in% c("Mixed","No grade"),9,9))))),
                   HOLC_grade_Krieger2=factor(HOLC_grade_Krieger2,
                                              c(1:4,9),
                                              c("50+% A","50+% B","50+% C","50+% D","Excluded")))
    
    
    #remove duplicates (those with multiple HOLC areas in them)
    # areadf2 <- areadf2[!duplicated(areadf2$GEOID_NEW), ]
    
    #add geometry back for mapping
    areadf5 <- left_join(areadf4, areadf, by="GEOID_NEW")
    
    #keep columns that we need
    areadf5 <- select(areadf5, GEOID_NEW, A, B, C, D, Containment, HOLC_grade_Krieger, HOLC_grade_Krieger2, geometry)
    
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
  
  blocks_krieger <- krieger_calc(blocksdf)
  
  
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
  
  
  bgs_krieger <- krieger_calc(bgs_df)
  
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
  
  
  tracts_krieger <- krieger_calc(tractsdf)
  

  ###############
  # Output     #
  ##############
  
  return(list("blocks_krieger"=blocks_krieger, "bgs_krieger"=bgs_krieger, "tracts_krieger"=tracts_krieger))
}

#############################
# Function specifications
##############################

# HOLCmap - specify HOLC map to be used
# counties - specify counties that overlap with HOLC map to reduce number of areas pulled in 
# year - specify year (2010, 2020) 


#######
# ATL #
#######

#2020
ATL_krieger20 <- krieger_fun(HOLCmap=ATL_red, #specify HOLC map
                         counties=c("063", #Clayton
                                    "089", #DeKalb
                                    "121"), #Fulton
                         yr="2020")

