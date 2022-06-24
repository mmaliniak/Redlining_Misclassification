#Purpose: Classify area level data (blocks, block groups, tracts) using centroid approach 
#         which assigns HOLC grade of area's centroid using different cutoffs for % ungraded

#Reference(s): Nardone A, et al. Environmental Justice. 2020 Aug;13(4):109–19.
#              Nardone A, et al. The Lancet Planetary Health. 2020 Jan 1;4(1):e24–31

#Applying to: Four MTAs with HOLC maps in Georgia but could run for any city with HOLC map

#Last updated: 4-25-2022


library(tidyverse)
library(sf)        # for handling simple features spatial data
library(tmap)      # for producing thematic maps
library(RColorBrewer)
library(tidycensus)
library(tigris)
library(purrr)


#read in GA HOLC maps
ATL_red <- read_sf('H:/Redlining and obesity/Data/Redlining/Atlanta/cartodb-query.shp') 


#read in 2020 Population-weighted centroid data + create GEOID
#Note: block group and tract data are for GA only (downloaded from US Census Bureau: https://www.census.gov/geographies/reference-files/time-series/geo/centers-population.html) 
cent_bgs <- read.xlsx("H:/Redlining and obesity/Data/Centroid data/Census_2020_PopCentroids_BlockGrp.xlsx") %>% 
  mutate(GEOID=paste0(STATEFP, COUNTYFP,TRACTCE,BLKGRPCE))

cent_tracts <- read.xlsx("H:/Redlining and obesity/Data/Centroid data/Census_2020_PopCentroids_Tract.xlsx") %>% 
  mutate(GEOID=paste0(STATEFP, COUNTYFP,TRACTCE))



#centroid function


#function to assign redlining status by area (blocks, block groups, tracts)
centroid_fun <- function(city, HOLCmap, counties, yr, cutoff){
  
  #read in HOLC map
  Red <- HOLCmap
  
  HOLC_match_cut <- function(areadf){
    
    #GEOID column name is different depending on year so need to rename so always the same regardless of year pulling
    #NOTE: this is clunky AF but can't figure out better way to do this for sf (if do it normal way - ends up making it into a list with geometry)
    geoid_col <- colnames(areadf)[grepl("GEOID",colnames(areadf))]
    GEOID2 <- areadf[geoid_col]
    st_geometry(GEOID2) <- NULL
    GEOID2$GEOID_NEW <- GEOID2[,geoid_col]
    areadf <- left_join(areadf, GEOID2, by=geoid_col)
    
    #set crs to be the same 
    areadf <- st_transform(areadf, crs = 26916)
    Red1 <- st_transform(Red, crs = 26916)
    
    # Calculate area and tidy up
    intersect_pct <- st_intersection(areadf, Red1) %>% 
      mutate(intersect_area = st_area(.)) %>%   # create new column with shape area
      dplyr::select(GEOID_NEW, holc_id, holc_grade, intersect_area) %>%   # only select columns needed to merge
      st_drop_geometry()  # drop geometry as we don't need it
    
    # Create a fresh area variable (will give slightly different results than ALAND provided by CB (even for areas without water))
    areadf <- mutate(areadf, poly_area = st_area(areadf))
    
    # Merge by census area (block, block group, tract)
    areadf2 <- merge(areadf, intersect_pct, by = "GEOID_NEW", all.x = TRUE)
    
    # Calculate coverage (aka tract_prop if tract = % of tract within HOLC area)
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
                      Exclude = ifelse(sum_poly_prop > 0 & sum_poly_prop < cutoff, 1, #1=some overlap but < cutoff
                                       ifelse(sum_poly_prop==0, 2, 0))) #2=no overlap
    
    #remove duplicates (those with multiple HOLC areas in them)
    areadf2 <- areadf2[!duplicated(areadf2$GEOID_NEW), ]
    
    #exclude areas with no overlap
    areadf2 <- filter(areadf2, Exclude!=2) #0=overlap and meets ungraded requirements, 1=overlap but doesn't meet ungraded requirements
    
    #remove columns that we don't need
    areadf2 <- select(areadf2, -holc_id, -holc_grade, -intersect_area, -poly_area)
    
    
    return(areadf2)
  }
  
  #have to do blocks separately from other areas because no population-weighted centroid data for blocks bc block is lowest level with pop data
  centroid_calc_blocks <- function(areadf){
    
    #remove geometry so can reset geometry using coordinate for centroid (i.e., interior points for blocks)
    st_geometry(areadf) <- NULL
    
    #use lat and long provided by census bureau for each block
    proj4 <- st_crs(Red)$proj4string  ##make sure to use original HOLC map (for projection)
    areadf$INTPTLON20 <- as.numeric(areadf$INTPTLON20)
    areadf$INTPTLAT20 <- as.numeric(areadf$INTPTLAT20)
    
    areadf <- st_as_sf(areadf, coords = c('INTPTLON20', 'INTPTLAT20'), crs = proj4) #x=lon y=lat
    
    #determine HOLC grade of block centroid
    areadf2 <- areadf %>% mutate(
      intersection = as.integer(st_intersects(geometry, Red))
      , area = if_else(is.na(intersection), '', Red$holc_grade[intersection])
    ) 
    
    areadf2$holc_grade_centroid <- ifelse(areadf2$area=="A", 1,
                                          ifelse(areadf2$area=="B",2,
                                                 ifelse(areadf2$area=="C",3,
                                                        ifelse(areadf2$area=="D",4,9))))
    
    #recode those that don't meet cutoff requirement for % ungraded
    areadf2$holc_grade_centroid <- ifelse(areadf2$Exclude==1,9,areadf2$holc_grade_centroid)
    
    areadf2$holc_grade_centroid <- factor(areadf2$holc_grade_centroid, 
                                          c(1:4,9),
                                          c("A","B","C","D","Excluded"))
    
    return(areadf2)
  }
  
  centroid_calc <- function(centdf){
    
    #set projections so they match
    proj4 <- st_crs(Red)$proj4string  ##make sure to use original HOLC map (for projection)
    
    centdata <- st_as_sf(centdf, 
                         coords = c('LONGITUDE', 'LATITUDE'), #x=longitude y=latitude
                         crs = proj4)
    
    centdata2 <- centdata %>% mutate(
      intersection = as.integer(st_intersects(geometry, Red))
      , area = if_else(is.na(intersection), '', Red$holc_grade[intersection])
    ) 
    
    
    centdata2$holc_grade_centroid <- ifelse(centdata2$area=="A", 1,
                                            ifelse(centdata2$area=="B",2,
                                                   ifelse(centdata2$area=="C",3,
                                                          ifelse(centdata2$area=="D",4,9))))
    
    # centdata2$holc_grade_centroid <- factor(centdata2$holc_grade_centroid, 
    #                                         c(1:4,9),
    #                                         c("A","B","C","D","Excluded"))
    return(centdata2)
  }
  
  
  ################
  #Census blocks #
  ################
  
  #bring in census block data
  blocksdf <- blocks(
    state="GA",
    county = counties, 
    class = 'sf',
    year = yr)
  
  #restrict to blocks that overlap with HOLC areas for mapping and for function to run faster (otherwise returns all blocks within specified GA counties)
  HOLC_blocks <- HOLC_match_cut(blocksdf)
  
  #Get HOLC grades according to area's centroid (includes centroid geometry for mapping)
  blocks_centroid <- centroid_calc_blocks(HOLC_blocks)
  
  #remove centroid geometry for block only maps 
  blocks_centroid2 <- blocks_centroid
  st_geometry(blocks_centroid2) <- NULL
  
  #merge to get area geometry for mapping blocks (no centroids)
  blocks_centroid2 <- left_join(HOLC_blocks, blocks_centroid2)
  
  
  ######################
  #Census block groups #
  ######################
  
  #bring in census block group data
  bgs_df <- block_groups(
    state="GA",
    county = counties, 
    class = 'sf',
    year = yr)
  
  #restrict to block groups that overlap with HOLC areas for mapping and for function to run faster 
  HOLC_bgs <- HOLC_match_cut(bgs_df)
  
  #Get HOLC grades according to area's centroid (NOTE: this is for all of GA)
  bgs_centroid <- centroid_calc(cent_bgs)
  
  #subset to areas that overlap with HOLC map and recode holc grades per cutoff
  HOLC_bgs2 <- HOLC_bgs
  st_geometry(HOLC_bgs2) <- NULL
  bgs_centroid <- left_join(bgs_centroid, HOLC_bgs2) %>% filter(Exclude!=2)
  
  #recode those that don't meet cutoff requirement for % ungraded
  bgs_centroid$holc_grade_centroid <- ifelse(bgs_centroid$Exclude==1,9,bgs_centroid$holc_grade_centroid)
  
  bgs_centroid$holc_grade_centroid <- factor(bgs_centroid$holc_grade_centroid, 
                                          c(1:4,9),
                                          c("A","B","C","D","Excluded"))
  
  #remove centroid geometry for block group only maps 
  bgs_centroid2 <- bgs_centroid
  st_geometry(bgs_centroid2) <- NULL
  
  #merge to get area geometry for mapping block groups (no centroids)
  bgs_centroid2 <- left_join(HOLC_bgs, bgs_centroid2)
  
  
  ######################
  #Census tracts       #
  ######################
  
  tractsdf <- tracts(
    state="GA",
    county = counties, 
    class = 'sf',
    year = yr)
  
  
  #restrict to tracts that overlap with HOLC areas for mapping and for function to run faster 
  HOLC_tracts <- HOLC_match_cut(tractsdf)
  
  #Get HOLC grades according to area's centroid (NOTE: this is for all of GA)
  tracts_centroid <- centroid_calc(cent_tracts)
  
  #subset to areas that overlap with HOLC map and recode holc grades per cutoff
  HOLC_tracts2 <- HOLC_tracts
  st_geometry(HOLC_tracts2) <- NULL
  tracts_centroid <- left_join(tracts_centroid, HOLC_tracts2) %>% filter(Exclude!=2)
  
  #recode those that don't meet cutoff requirement for % ungraded
  tracts_centroid$holc_grade_centroid <- ifelse(tracts_centroid$Exclude==1,9,tracts_centroid$holc_grade_centroid)
  
  tracts_centroid$holc_grade_centroid <- factor(tracts_centroid$holc_grade_centroid, 
                                             c(1:4,9),
                                             c("A","B","C","D","Excluded"))
  
  #remove centroid geometry for tract only maps 
  tracts_centroid2 <- tracts_centroid
  st_geometry(tracts_centroid2) <- NULL
  
  #merge to get geometry for mapping
  tracts_centroid2 <- left_join(HOLC_tracts, tracts_centroid2)
  

  # ###############
  # # Output     #
  # ##############
  # 
  cutoff2 <- as.character(cutoff*100)

  centroid_data <- list("blocks_centroid"=blocks_centroid, "bgs_centroid"=bgs_centroid, "tracts_centroid"=tracts_centroid,
                     "blocks_centroid2"=blocks_centroid2, "bgs_centroid2"=bgs_centroid2, "tracts_centroid2"=tracts_centroid2)
  saveRDS(centroid_data, file = paste0("H:/Redlining and obesity/Data/Redlining/Redlining area data/Centroid/", city, "/", "Cutoff ", cutoff2, "/", city, "_Centroid", yr,".RDS"))
 
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
ATLcounties <- c("063", "089","121") #Clayton, DeKalb, and Fulton
cutoff <- c(.10,.20,.30,.40,.50)
map(.x = cutoff, 
    .f = centroid_fun, 
    city = "Atlanta", 
    HOLCmap = ATL_red, 
    counties=ATLcounties,
    yr="2020")



