#### Working directory and packages ####
# wd <- "C:/Users/user2380/Documents/Julie/macroclim/data/"
wd <- "/home/julie/macroclim/data/"
setwd(wd) 

pacman::p_load(data.table, dplyr, tidyverse, stringr, tibble,#data wrangling
               terra, tidyterra, #gis
               foreach, doParallel # to parallelize
)

select <- dplyr::select

#### Loading data and data wrangling ####

# boundary box for the shapefile
bbox_terra <- c(-10.5,32.702, 34.856, 71.31)

# Shapefile of selected rivers
rivers_selec <-terra::vect(paste0(wd,"gis/HydroRIVERS_hybas5_family_level.shp"),
                           extent=bbox_terra) #already cropped to the right extent

# Dataframe of river attributes
rivers_selec_df <- rivers_selec %>%
  as.data.frame() %>%
  select(HYRIV_ID, NEXT_DOWN, MAIN_RIV, LENGTH_KM, HYBAS_L12) %>%
  as.data.table()

# Dataframe with the ID of all downstream reaches for each river segment
# (IDs in column 'ALL_DOWN' are in the order from up to downstream)
river_downstream_segments  <- read.csv(paste0(wd,"hydroatlas/hydroshed_downstream_segments_v2.csv"),
                                       header=T, sep=",", stringsAsFactors = FALSE)

# # Spatial join between HydroRIVERS and Hybas_level3
# hybas_hyriv <- read.csv(paste0(wd,"hydroatlas/hydroriver_hybas3.csv"),
#                         header=T, sep=";", stringsAsFactors = FALSE) %>%
#   select(HYRIV_ID, HYBAS_ID)

# Data wrangling on river_downstream_segments to prepare the loop
river_downstream_segments <- river_downstream_segments %>%
  filter(HYRIV_ID %in% rivers_selec_df$HYRIV_ID) %>%  # keeping rivers from the selection
  filter(ALL_DOWN != "0") %>% # removing outlets, only if directional scenario = only up -> down distances
  mutate(ALL_DOWN = str_replace_all(ALL_DOWN, ', 0$','')) %>%
  as.data.table()
  # left_join(hybas_hyriv, by = "HYRIV_ID") %>% # adding hybas_level3 information
  # filter(!is.na(HYBAS_ID))

# # vector of all hybas_level3 identifiers
# hybas3_vec <- unique(river_downstream_segments$HYBAS_ID)
# hybas3_vec <- hybas3_vec[order(hybas3_vec)]

#### Parallelization parameters #### 
n.cores <- 10
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

# #### Calculating nb of reaches between outlets of each pairs of connected reaches - ASYMETRIC DISTANCE ####
# 
# print("Starting loop")
# Sys.time()
# 
# # Creation of a dataframe "length_down" giving the length between outlets of pairs of river reaches.
# # The foreach loop is adding rows/ blocks of rows to "length_down". One row/ block of rows ("bound" object) has the following characteristics:
# # id of one considered HydroRIVER reach, id of one of its downstream reaches,
# # the length between the outlets of the two reaches, and the "boundary" (Virgilio Hermoso's definition to use in Marxan)
# # If a given river has 3 downstream reaches, the associated "bound" object should have 3 rows.
# count_down <- foreach (i = 1:nrow(river_downstream_segments),
#                         .packages = 'dplyr',
#                         .combine = 'rbind') %dopar% {
#                           
#                           # for one giver HydroRIVER reach, I obtain a vector with the IDs of all downstream reaches
#                           all_down <- unlist(strsplit(river_downstream_segments[i,'ALL_DOWN'], split=", "))
#                           n = length(all_down)
#                           
#                           # for all downstream reaches, I'm gonna create a new row in the "bound" object 
#                           for (j in 1:n){
#                             
#                             if (j == 1){
#                               count <- 1
#                               boundary <- 1/(sqrt(count))
#                               bound <- c(river_downstream_segments[i, 'HYRIV_ID'],
#                                          all_down[j],
#                                          count,
#                                          boundary)
#                             }
#                             else if (count < 100 & j > 1){ # currently there is a limit set at 100km: pairs of sites that are further apart are not considered
#                               # adding the length of the j-th downstream reach to the "length" previously calculated (because the IDs are in order from up to downstream)
#                               count <- count + 1
#                               boundary <- 1/(sqrt(count))
#                               bound <- rbind(bound,
#                                              c(river_downstream_segments[i, 'HYRIV_ID'],
#                                                all_down[j],
#                                                count,
#                                                boundary
#                                              ))
#                               
#                             }
#                             
#                           }
#                           return(bound)
#                           
#                           rm(list=c('count','boundary','all_down'))
#                         }
# 
# Sys.time()
# print("End of loop")
# 
# count_down <- count_down %>%
#   as.data.frame() %>%
#   setNames(c("id1", "id2", "count", "boundary"))
# 
# write.csv(count_down,
#           paste0(wd,"marxan/boundary_limit.csv"), row.names = FALSE)
# 
# print("Saved")              
# 
# 
#### Calculating total length between outlets of each pairs of connected reaches - ASYMETRIC DISTANCE ####

print("Starting loop")
Sys.time()

# Creation of a dataframe "length_down" giving the length between outlets of pairs of river reaches.
# The foreach loop is adding rows/ blocks of rows to "length_down". One row/ block of rows ("bound" object) has the following characteristics:
# id of one considered HydroRIVER reach, id of one of its downstream reaches,
# the length between the outlets of the two reaches, and the "boundary" (Virgilio Hermoso's definition to use in Marxan)
# If a given river has 3 downstream reaches, the associated "bound" object should have 3 rows.
length_down <- foreach (i = 500:1000,
                   # .packages = 'dplyr',
                   .packages = 'data.table',
                   .combine = 'rbind') %dopar% {
                     
                     # for (i in 1:1000){

                     # for one giver HydroRIVER reach, I obtain a vector with the IDs of all downstream reaches
                     # all_down <- unlist(strsplit(river_downstream_segments[i,'ALL_DOWN'], split=", "))
                     all_down <- unlist(strsplit(river_downstream_segments[i,ALL_DOWN], split=", "))
                     n = length(all_down)

                     # for all downstream reaches, I'm gonna create a new row in the "bound" object
                     for (j in 1:n){

                       if (j == 1){
                         # length <- (rivers_selec_df %>% filter(HYRIV_ID ==  all_down[j]))$LENGTH_KM # obtaining the length of the downstream reach "j"
                         # length <-rivers_selec_df[which(rivers_selec_df$HYRIV_ID ==  all_down[j]),]$LENGTH_KM
                         length <-rivers_selec_df[HYRIV_ID ==  all_down[j],LENGTH_KM]
                         boundary <- 1/(sqrt(length))
                         bound <- c(river_downstream_segments[i, 'HYRIV_ID'],
                                    all_down[j],
                                    length,
                                    boundary)
                         # bound <- rbind(bound,  # to test with the simple for-loop
                         #                c(river_downstream_segments[i, 'HYRIV_ID'],
                         #                  all_down[j],
                         #                  length,
                         #                  boundary
                         #                ))
                       }
                       else if (length < 100 & j > 1){ # currently there is a limit set at 100km: pairs of sites that are further apart are not considered
                         # adding the length of the j-th downstream reach to the "length" previously calculated (because the IDs are in order from up to downstream)
                         # length <- length + (rivers_selec_df %>% filter(HYRIV_ID ==  all_down[j]))$LENGTH_KM
                         length <- length + rivers_selec_df[HYRIV_ID ==  all_down[j],LENGTH_KM]
                         boundary <- 1/(sqrt(length))
                         bound <- rbind(bound,
                                        c(river_downstream_segments[i, 'HYRIV_ID'],
                                          all_down[j],
                                          length,
                                          boundary
                                        ))

                       }

                     }
                     return(bound)

                     rm(list=c('length','boundary','all_down'))
                   }

Sys.time()
print("End of loop")

length_down <- length_down %>%
  as.data.frame() %>%
  setNames(c("id1", "id2", "length", "boundary"))

write.csv(length_down,
          paste0(wd,"marxan/boundary_limit_length10000.csv"), row.names = FALSE)

print("Saved")

# #### Calculating total length between outlets of each pairs of connected reaches - complementary down > up DISTANCE ####
# 
# Sys.time()
# 
# length_down <- foreach (i = 1:nrow(river_downstream_segments),
#                         .packages = 'dplyr',
#                         .combine = 'rbind') %dopar% {
#                           
#                           all_down <- unlist(strsplit(river_downstream_segments[i,'ALL_DOWN'], split=", "))
#                           all_down <- all_down[-length(all_down)]
#                           n = length(all_down)
#                           
#                           hyriv = river_downstream_segments$HYRIV_ID[i]
#                           mainriv = (rivers_selec_df %>% filter(HYRIV_ID == hyriv))$MAIN_RIV
#                           same_catch <- (rivers_selec_df %>% filter(MAIN_RIV == mainriv))$HYRIV_ID
#                           river_sub <- river_downstream_segments %>%
#                             filter(HYRIV_ID %in% same_catch) %>%
#                             filter(HYRIV_ID != hyriv) %>%
#                             filter(!(HYRIV_ID %in% all_down))
#                           
#                           bound <- vector()
#                           
#                           for (j in 1 :nrow(river_sub)){
#                             
#                             hyriv2 =  river_sub$HYRIV_ID[j]
#                             all_down2 <- unlist(strsplit(river_sub[j,'ALL_DOWN'], split=", "))
#                             all_down2 <- all_down2[-length(all_down2)]
#                             
#                             if (identical(all_down, all_down2)){
#                               length <- 0
#                               boundary <- 100
#                                 
#                             } else {
# 
#                               n1 = which(all_down %in% intersect(all_down, all_down2))[1]
#                               n2 = which(all_down2 %in% intersect(all_down, all_down2))[1]
#                               
#                               reaches_to_sum = union(all_down[-(n1:length(all_down))], all_down2[-(n2:length(all_down2))])
#                               
#                               length <- sum((rivers_selec_df %>% filter(HYRIV_ID %in% reaches_to_sum))$LENGTH_KM)
#                               boundary <- 1/(sqrt(length))
#                               
#                             }
#  
#                               bound <- rbind(bound,
#                                              c(hyriv,
#                                                hyriv2,
#                                                length,
#                                                boundary
#                                              ))
# 
#                           }
#                           return(bound)
# 
#                           rm(list=c('length','boundary','all_down'))
#                         }
# 
# Sys.time()
# 
# length_down <- length_down %>%
#   as.data.frame() %>%
#   setNames(c("id1", "id2", "length", "boundary"))
# 
# write.csv(length_down,
#           paste0(wd,"marxan/boundary_limit_asym.csv"), row.names = FALSE)
# 
# 
# 
# 
# 
