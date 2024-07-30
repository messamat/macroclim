#wd <- "/home/julie/macroclim/data/"
#setwd(wd) 

pacman::p_load(data.table, dplyr, igraph, tidyverse, stringr, tibble,#data wrangling
               rprojroot, RcppAlgos,
               terra, tidyterra, #gis
               foreach, doParallel # to parallelize
)
rootdir <- rprojroot::find_root(rprojroot::has_dir('src'))
wd <- file.path(rootdir, 'data')
# wd <- "/home/julie/macroclim/data/"
# setwd(wd) 
# select <- dplyr::select

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
                       rivers_dt[nsegs_network > 1, .(HYRIV_ID, LENGTH_KM, ORD_STRA)],
                       by.x='NEXT_DOWN', by.y='HYRIV_ID', 
                       all.x=T,
                       suffixes=c("", "_down")) %>%
  .[, .(HYRIV_ID, NEXT_DOWN, MAIN_RIV, ORD_STRA, ORD_STRA_down,
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
    
    #Compute distance between first- and segments downstream
    dist_k1_to_k2<- basin_sel_k1[, .(HYRIV_ID, NEXT_DOWN, LENGTH_KM_down)] %>%
      setnames(c('HYRIV_ID', 'NEXT_DOWN', 'LENGTH_KM_down'),
               c('HYRIV_ID_from', 'HYRIV_ID_to', 'dist_km'))
    
    #Accumulate distance between first- and all downstream segments
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
               c('HYRIV_ID_from', 'HYRIV_ID_to', 'dist_km')
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
}) %>% 
  rbindlist(use.names=T) %>%
  merge(rivers_inland[, .(HYRIV_ID, MAIN_RIV)], 
        by.x='HYRIV_ID_from', by.y='HYRIV_ID')


#### Strahler order method #####################################################
rivers_sub <- rivers_inland[MAIN_RIV == 20268761,]

# plot(igraph::graph_from_data_frame(
#   rivers_sub[!is.na(LENGTH_KM_down),],
#   directed = T))

sel_cols <- c('HYRIV_ID', 'NEXT_DOWN', 'LENGTH_KM_down', 'ORD_STRA', 'ORD_STRA_down')
new_cols <- c('HYRIV_ID_from', 'HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to')

#Compute distance from outlet of first-order streams to outlet of second-order streams
k_max <- max(rivers_sub$ORD_STRA_down, na.rm=T)
k_dt <- expand.grid(1:k_max, 1:k_max) %>%
  setDT %>%
  setnames(c('k_i', "k_j")) %>%
  .[order(k_i, k_j)] %>%
  .[!(k_i==1 & k_j<=2) & k_j >= k_i,]

all_acc <- rivers_sub[ORD_STRA==1 & ORD_STRA_down==2, ..sel_cols] %>% 
  setnames(sel_cols, new_cols) 

downdist_hydrorivers_basin_manual_inner <- function(
    in_rivers_format, in_all_acc, strahler_order_from, strahler_order_to) {
  
  k_list_sub <- in_rivers_format[ORD_STRA==strahler_order_from & 
                                   ORD_STRA_down==strahler_order_to, 
                                 ..sel_cols]
  
  if (nrow(k_list_sub) > 0) {
    to_next <- k_list_sub %>% 
      copy %>%
      setnames(sel_cols, new_cols)
    
    acc_to_k <- merge(all_acc, k_list_sub, 
                      by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                      suffix=c('_upst', '_extra')) 
    if (nrow(acc_to_k) > 0) {
      while (nrow(acc_to_k) > 0) {
        acc_to_k_format <- acc_to_k[, dist_km_acc := dist_km + LENGTH_KM_down] %>%
          .[, c('HYRIV_ID_from', 'NEXT_DOWN', 
                'ORD_STRA_from', 'ORD_STRA_down',
                'dist_km_acc'), with=F] %>%
          setnames(c('NEXT_DOWN', 'dist_km_acc', 'ORD_STRA_from', 'ORD_STRA_down'),
                   c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))
        
        # acc_to_k <- acc_to_k[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
        #   . [, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
        #          'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
        #          'dist_km_acc'), with=F] %>%
        #   setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
        #            c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))
        
        all_acc <- rbindlist(list(all_acc, to_next, acc_to_k_format), use.names=T)
        
        
        ####REMOVE DUPLICATES
        acc_to_k <- merge(all_acc, k_list_sub, 
                          by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                          suffix=c('_upst', '_extra')) %>%
          .[!(paste(HYRIV_ID_from, NEXT_DOWN, sep='_') %in% 
                all_acc[, paste(HYRIV_ID_from, HYRIV_ID_to, sep='_')]),]
      }
    } else {
      all_acc <- rbindlist(list(all_acc, to_next), use.names=T)
    }
    remove(k_list_sub, to_next, acc_to_k)
  }  
  
  return(all_acc)
}

for (k_comb_index in 1:nrow(k_dt)) {
  k_i <- k_dt[k_comb_index, k_i]
  k_j <- k_dt[k_comb_index, k_j]
  print(k_i)
  print(k_j)
  
  all_acc <- downdist_hydrorivers_basin_manual_inner(
    in_rivers_format = rivers_sub,
    in_all_acc = all_acc,
    strahler_order_from = k_i,
    strahler_order_to = k_j
  )
}


# check<- all_acc[!duplicated(paste(HYRIV_ID_from, HYRIV_ID_to, sep='_')),]
# check2 <-dist_dt_all[MAIN_RIV==20268761,]
# dist_dt_all[MAIN_RIV==20268761,][
#   !(paste(HYRIV_ID_from, HYRIV_ID_to, sep='_') %in%
#     all_acc[,  paste(HYRIV_ID_from, HYRIV_ID_to, sep='_')]),]
