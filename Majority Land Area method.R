#Purpose: Classify area level data (census blocks, block groups, tracts) using Majority Land Area method
#         which assigns Assign HOLC grade of majority

#Reference(s): Krieger N, et al. Am J Epidemiol. 2020 Oct 1;189(10):1065–75.
#              Krieger N, et al. Am J Public Health. 2020 Jul;110(7):1046–53.


#Applying to: Atlanta, GA HOLC map but could run for any city with HOLC map

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
majority_fun <- function(HOLCmap, st, counties, yr){
  
  #read in HOLC map
  Red <- HOLCmap
  
  #pull in water polygons for counties that overlap HOLC areas
  water <- area_water(state=st, 
                      county = counties, 
                      year = yr) %>% st_as_sf()
  
  
  majority_calc <- function(areadf){
    
    #GEOID column name is different depending on year so need to rename so always the same regardless of year pulling
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
    
    
    #HOLC grade assignment to each census area
    areadf4 <- mutate(areadf4, 
                   Containment = ifelse(A==1 | B==1 | C==1 | D==1,1, #Identify tracts where 100% of tract is contained in 1 HOLC area
                                        ifelse( (A < 1 & A >=0.5) |  #Identify tracts where 50+% of census area is contained in 1 HOLC grade
                                                  (B < 1 & B >=0.5) | 
                                                  (C < 1 & C >=0.5) | 
                                                  (D < 1 & D >= 0.5),2,
                                                ifelse( (A + B + C + D) >=0.5,3,4))),
                   Containment = factor(Containment, c(1,2,3,4), c("Fully contained", "Majority contained", "Mixed", "No grade")),
                   HOLC_grade_majority = ifelse(A >=0.5,1,
                                               ifelse(B >=0.5,2,
                                                      ifelse(C >=0.5,3,
                                                             ifelse(D >= 0.5,4,
                                                                    ifelse(Containment=="Mixed",8,9))))),
                   HOLC_grade_majority=factor(HOLC_grade_majority,
                                             c(1:4,8,9),
                                             c("50+% A","50+% B","50+% C","50+% D","50+% Mixed", "<50%: No grade")),
                   
                   #collapse Mixed and no grade together
                   HOLC_grade_majority2 = ifelse(A>=0.5,1,
                                                ifelse(B>=0.5, 2,
                                                       ifelse(C >=0.5, 3,
                                                              ifelse(D >= 0.5,4,
                                                                     ifelse(Containment %in% c("Mixed","No grade"),9,9))))),
                   HOLC_grade_majority2=factor(HOLC_grade_majority2,
                                              c(1:4,9),
                                              c("50+% A","50+% B","50+% C","50+% D","Excluded")))
    
  
    
    #add geometry back for mapping
    areadf5 <- left_join(areadf4, areadf, by="GEOID_NEW")
    
    #keep columns that we need
    areadf5 <- select(areadf5, GEOID_NEW, A, B, C, D, Containment, HOLC_grade_majority, HOLC_grade_majority2, geometry)
    
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
    state=st,
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
    state=st,
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
    state=st,
    county = counties, 
    class = 'sf',
    year = yr)
  
  
  tracts_majority <- majority_calc(tractsdf)
  

  ###############
  # Output     #
  ##############
  
  return(list("blocks_majority"=blocks_majority, "bgs_majority"=bgs_majority, "tracts_majority"=tracts_majority))
}

#############################
# Function specifications
##############################

# HOLCmap - specify HOLC map to be used (loaded above)

# st - specificy state

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

# year - specify census year for getting data (2010, 2020) 


#######
# ATL #
#######

#2020
ATL_majority20 <- majority_fun(HOLCmap=ATL_red, #specify HOLC map
                         st="GA",
                         counties=c("063", #Clayton
                                    "089", #DeKalb
                                    "121"), #Fulton
                         yr="2020")


