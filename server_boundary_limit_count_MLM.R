pacman::p_load(data.table, dplyr, igraph, tidyverse, stringr, tibble,#data wrangling
               rprojroot, RcppAlgos,
               terra, tidyterra, #gis
               foreach, doParallel # to parallelize
)
rootdir <- rprojroot::find_root(rprojroot::has_dir('src'))
wd <- file.path(rootdir, 'data')
select <- dplyr::select

#Paths
in_rivers_geom <- file.path(wd,"gis", "HydroRIVERS_hybas5_family_level.shp")
in_rivers_downstream <- file.path(wd, "hydroatlas", "hydroshed_downstream_segments_v2.csv")

#### Loading data and data wrangling ####
# boundary box for the shapefile
bbox_terra <- c(-10.5,32.702, 34.856, 71.31) #This cuts off a few basins in Cyprus and northeastern Europe: eastern cutoff would be better with 35E 

# Shapefile of selected rivers
rivers_selec <- terra::vect(in_rivers_geom,
                            extent=bbox_terra) #already cropped to the right extent
rivers_dt <- as.data.table(rivers_selec)

#Compute number of connected segments
rivers_dt[, nsegs_network := .N, by=MAIN_RIV]

#Distinguish single-seg networks
rivers_coastal <- rivers_dt[nsegs_network == 1,]

#Get length of downstream segment because interested in computing distance between reach outlets
rivers_inland <- merge(rivers_dt[nsegs_network > 1,], 
                       rivers_dt[nsegs_network > 1, .(HYRIV_ID, LENGTH_KM)],
                       by.x='NEXT_DOWN', by.y='HYRIV_ID', 
                       all.x=T,
                       suffixes=c("", "_down")) %>%
  .[, .(HYRIV_ID, NEXT_DOWN, MAIN_RIV, ORD_STRA,
        LENGTH_KM, LENGTH_KM_down, nsegs_network)] 

#Get unique multiseg networks IDs (based on outlet's HYRIV_ID)
unique_main <- rivers_inland[, unique(MAIN_RIV)]

#### IGRAPH method #############################################################
#### Calculating total length between outlets of each pairs of connected reaches - ASYMETRIC DISTANCE ####
#Pairwise distance computation is extremely fast in igraph but creates a large matrix
#so it is limited by memory.

#Function to compute downstream distance for all segments in a HydroRIVERS basin
downdist_hydrorivers_basin <- function(in_rivers_format, in_main_riv, 
                                       max_dist=NULL, quiet=T) {
  #Subset basin
  basin_sel <- in_rivers_format[MAIN_RIV == in_main_riv,] 
  
  #Print info
  if (!quiet) {
    print(in_main_riv)
    print(basin_sel[1, nsegs_network])
  }
  
  #Get first order streams
  basin_sel_k1 <- basin_sel[ORD_STRA==1,]
  
  #If more than one confluence in network
  if (basin_sel[!is.na(LENGTH_KM_down), length(unique(NEXT_DOWN))] > 1) {
    #Make directed graph
    #-- Removing segments whose downstream segment was truncated from the network
    basin_graph <- igraph::graph_from_data_frame(
      basin_sel[!is.na(LENGTH_KM_down) & ORD_STRA>1,],  
      directed = T); #plot(graph)
    
    #Assign distance from reference outlet to the outlet of downstream segment
    E(basin_graph)$weight <- basin_sel[!is.na(LENGTH_KM_down) & ORD_STRA>1, 
                                       LENGTH_KM_down]
    #Compute pairwise distances between all segments except 1st order streams
    dist_mat <- distances(basin_graph, mode='out') 
    #Remove diagonal
    diag(dist_mat) <- NA
    from_ids <- as.integer(rownames(dist_mat))
    #Convert to long form data.table
    dist_mat <- as.data.table(dist_mat) %>%
      .[, HYRIV_ID_from := from_ids] %>%
      melt(id.vars='HYRIV_ID_from', 
           variable.name='HYRIV_ID_to', value.name='dist_km') %>%
      .[!is.na(dist_km) & !is.infinite(dist_km),]
    
    #Compute distance between first- and second-order stream segments
    dist_k1_to_k2<- basin_sel_k1[, .(HYRIV_ID, NEXT_DOWN, LENGTH_KM_down)] %>%
      setnames(c('HYRIV_ID', 'NEXT_DOWN', 'LENGTH_KM_down'),
               c('HYRIV_ID_from', 'HYRIV_ID_to', 'dist_km'))
    
    #Compute distance between first- and all downstream segments
    dist_k1_to_rest <- merge(basin_sel_k1[, .(HYRIV_ID, NEXT_DOWN, LENGTH_KM_down)],
                             dist_mat[HYRIV_ID_from %in% unique(basin_sel_k1$NEXT_DOWN),],
                             by.x='NEXT_DOWN',
                             by.y='HYRIV_ID_from',
                             suffix=c('_from', '_to'),
                             allow.cartesian=TRUE) %>%
      .[, dist_km := dist_km + LENGTH_KM_down] %>%
      setnames('HYRIV_ID', 'HYRIV_ID_from') %>%
      .[, .(HYRIV_ID_from, HYRIV_ID_to, dist_km)]
    
    #Merge everything
    dist_mat_fast <- rbindlist(list(dist_k1_to_k2, dist_k1_to_rest, dist_mat))
  } else {
    #For networks with only one confluence, simply use the segment-wise data
    sel_cols <- c('HYRIV_ID', 'NEXT_DOWN', 'LENGTH_KM_down')
    dist_mat_fast <- basin_sel_k1[, ..sel_cols] %>%
      setnames(sel_cols,
               c('HYRIV_ID_from', 'HYRIVID_to', 'dist_km')
    )
  }
  return(dist_mat_fast)
}

#Run for all basins (except the one with the most segment as my laptop does not have much memory)
dist_dt_all <- lapply(unique_main, function(basin_id) {
  if (basin_id != 20498112) {
    downdist_hydrorivers_basin(in_rivers_format = rivers_inland,
                               in_main_riv = basin_id,
                               quiet=F)
  }
})

#### Strahler order method #####################################################
rivers_inland <- rivers_inland[MAIN_RIV == 20431728,]

plot(igraph::graph_from_data_frame(
  rivers_inland[!is.na(LENGTH_KM_down),],  
  directed = T))

sel_cols <- c('HYRIV_ID', 'NEXT_DOWN', 'LENGTH_KM_down', 'ORD_STRA')
new_cols <- c('HYRIV_ID_from', 'HYRIV_ID_to', 'dist_km', 'ORD_STRA_from')
to_k2 <- rivers_inland[ORD_STRA==1, ..sel_cols] %>% 
  setnames(sel_cols, new_cols) %>%
  .[, ORD_STRA_to := ORD_STRA_from + 1]

k2_to_k3 <- rivers_inland[ORD_STRA==2 & !is.na(LENGTH_KM_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols) %>%
  .[, ORD_STRA_to := ORD_STRA_from + 1]

acc_to_k3 <- merge(to_k2, k2_to_k3, 
                  by.x='HYRIV_ID_to', by.y='HYRIV_ID_from',
                  suffix=c('_upst', '_extra')) %>%
  .[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
  .[, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
        'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
        'dist_km_acc'), with=F] %>%
  setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
           c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))

to_k3 <- rbind(acc_to_k3, k2_to_k3)







k3_to_k4 <- rivers_inland[ORD_STRA==3 & !is.na(LENGTH_KM_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols) %>%
  .[, ORD_STRA_to := ORD_STRA_from + 1]

acc_to_k4 <- merge(to_k3, k3_to_k4, 
                   by.x='HYRIV_ID_to', by.y='HYRIV_ID_from',
                   suffix=c('_upst', '_extra'),
                   allow.cartesian=TRUE) %>%
  .[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
  .[, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
        'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
        'dist_km_acc'), with=F] %>%
  setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
           c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))

to_k4 <- rbind(acc_to_k4, k3_to_k4)

k4_to_k5 <- rivers_inland[ORD_STRA==4 & !is.na(LENGTH_KM_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols) %>%
  .[, ORD_STRA_to := ORD_STRA_from + 1]

acc_to_k5 <- merge(to_k4, k4_to_k5, 
                   by.x='HYRIV_ID_to', by.y='HYRIV_ID_from',
                   suffix=c('_upst', '_extra'),
                   allow.cartesian=TRUE) %>%
  .[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
  .[, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
        'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
        'dist_km_acc'), with=F] %>%
  setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
           c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))

to_k5 <- rbind(acc_to_k5, k4_to_k5)

k5_to_k6 <- rivers_inland[ORD_STRA==5 & !is.na(LENGTH_KM_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols) %>%
  .[, ORD_STRA_to := ORD_STRA_from + 1]

acc_to_k6 <- merge(to_k5, k5_to_k6, 
                   by.x='HYRIV_ID_to', by.y='HYRIV_ID_from',
                   suffix=c('_upst', '_extra'),
                   allow.cartesian=TRUE) %>%
  .[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
  .[, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
        'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
        'dist_km_acc'), with=F] %>%
  setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
           c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))

to_k6 <- rbind(acc_to_k6, k5_to_k6)
