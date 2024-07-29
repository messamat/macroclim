#wd <- "/home/julie/macroclim/data/"
#setwd(wd) 

pacman::p_load(data.table, dplyr, igraph, tidyverse, stringr, tibble,#data wrangling
               rprojroot, RcppAlgos,
               terra, tidyterra, #gis
               foreach, doParallel # to parallelize
)
# rootdir <- rprojroot::find_root(rprojroot::has_dir('src'))
# wd <- file.path(rootdir, 'data')
wd <- "/home/julie/macroclim/data/"
setwd(wd) 
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

sel_cols <- c('HYRIV_ID', 'NEXT_DOWN', 'LENGTH_KM_down', 'ORD_STRA', 'ORD_STRA_down')
new_cols <- c('HYRIV_ID_from', 'HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to')

#Compute distance from outlet of first-order streams to outlet of second-order streams
k_max <- max(rivers_inland$ORD_STRA_down, na.rm=T)
k_dt <- expand.grid(1:k_max, 1:k_max) %>%
  setDT %>%
  setnames(c('k_i', "k_j")) %>%
  .[order(k_i, k_j)] %>%
  .[!(k_i==1 & k_j<=2) & k_j >= k_i,]

all_acc <- rivers_inland[ORD_STRA==1 & ORD_STRA_down==2, ..sel_cols] %>% 
  setnames(sel_cols, new_cols) 
k_comb_index = 1

##########
k_i <- k_dt[k_comb_index, k_i]
k_j <- k_dt[k_comb_index, k_j]

k_list_sub <- rivers_inland[ORD_STRA==k_i & ORD_STRA_down==k_j, ..sel_cols]

if (nrow(k_list_sub) > 0) {
  to_next <- k_list_sub %>% 
    copy %>%
    setnames(sel_cols, new_cols)
  
  acc_to_k <- merge(all_acc, k_list_sub, 
                    by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                    suffix=c('_upst', '_extra')) 
  
  if (nrow(acc_to_k) > 0) {
    acc_to_k <- acc_to_k[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
      . [, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
             'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
             'dist_km_acc'), with=F] %>%
      setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
               c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))
    
    all_acc <- rbindlist(list(all_acc, to_next, acc_to_k), use.names=T)
  } else {
    all_acc <- rbindlist(list(all_acc, to_next), use.names=T)
  }
  
  k_comb_index = k_comb_index + 1
}  else {
  #break
  k_comb_index = k_comb_index + 1
}

##OLD###########
####################################################################################
#################################################################################

if 
k_list <- rivers_inland[ORD_STRA==k & !is.na(LENGTH_KM_down), ..sel_cols]

#Deal with segments that do no increment in stream order
same_kupdown <- merge(all_acc, k2_list[(ORD_STRA == ORD_STRA_down),], 
                      by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                      suffix=c('_upst', '_extra'))


#Check second order streams -----------------------------------------------------
k2_list <- rivers_inland[ORD_STRA==2 & !is.na(LENGTH_KM_down), ..sel_cols]

#Deal with segments that do no increment in stream order
same_kdown <- merge(all_acc, k2_list[(ORD_STRA == ORD_STRA_down),], 
                by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                suffix=c('_upst', '_extra'))

same_kdown_format <- k2_list[(ORD_STRA == ORD_STRA_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols)

while (nrow(same_kdown)>0) {
  same_kdown_acc <- same_kdown[, dist_km_acc := dist_km + LENGTH_KM_down] %>%
    .[, c('HYRIV_ID_from', 'NEXT_DOWN', 
          'ORD_STRA_from', 'ORD_STRA_down',
          'dist_km_acc'), with=F] %>%
    setnames(c('NEXT_DOWN', 'dist_km_acc', 'ORD_STRA_from', 'ORD_STRA_down'),
             c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))
  
  all_acc <- rbindlist(
    list(
      all_acc,
      same_kdown_format,
      same_kdown_acc
    ),
    use.names=T
  )
    
  same_kdown <- merge(all_acc, k2_list, 
                      by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                      suffix=c('_upst', '_extra')) %>%
    .[(ORD_STRA == ORD_STRA_down) & 
        !(NEXT_DOWN %in% unique(all_acc$HYRIV_ID_to)),]
}

#Deal with segments that increment in stream order
k2_down <- rivers_inland[ORD_STRA==2 & !is.na(LENGTH_KM_down) &
                           (ORD_STRA != ORD_STRA_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols)

acc_to_k2_down <- merge(all_acc, k2_down, 
                   by.x='HYRIV_ID_to', by.y='HYRIV_ID_from',
                   suffix=c('_upst', '_extra')) %>%
  .[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
  .[, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
        'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
        'dist_km_acc'), with=F] %>%
  setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
           c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))

all_acc <- rbindlist(list(all_acc, k2_down, acc_to_k2_down), use.names=T)
 
#Check third order streams -----------------------------------------------------
k3_list <- rivers_inland[ORD_STRA==3 & !is.na(LENGTH_KM_down), ..sel_cols]

#Deal with segments that do no increment in stream order
same_kdown <- merge(all_acc, k3_list[(ORD_STRA == ORD_STRA_down),], 
                   by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                   suffix=c('_upst', '_extra'))

same_kdown_format <- k3_list[(ORD_STRA == ORD_STRA_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols)

while (nrow(same_kdown)>0) {
  same_kdown_acc <- same_kdown[, dist_km_acc := dist_km + LENGTH_KM_down] %>%
    .[, c('HYRIV_ID_from', 'NEXT_DOWN', 
          'ORD_STRA_from', 'ORD_STRA_down',
          'dist_km_acc'), with=F] %>%
    setnames(c('NEXT_DOWN', 'dist_km_acc', 'ORD_STRA_from', 'ORD_STRA_down'),
             c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))
  
  all_acc <- rbindlist(
    list(
      all_acc,
      same_kdown_format,
      same_kdown_acc
    ),
    use.names=T
  )
  
  same_kdown <- merge(all_acc, k3_list, 
                      by.x='HYRIV_ID_to', by.y='HYRIV_ID',
                      suffix=c('_upst', '_extra')) %>%
    .[(ORD_STRA == ORD_STRA_down) & 
        !(NEXT_DOWN %in% unique(all_acc$HYRIV_ID_to)),]
}

#Deal with segments that increment in stream order
k3_down <- rivers_inland[ORD_STRA==3 & !is.na(LENGTH_KM_down) &
                           (ORD_STRA != ORD_STRA_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols)

if (nrow(k3_down)) {
  acc_to_k3_down <- merge(all_acc, k3_down, 
                          by.x='HYRIV_ID_to', by.y='HYRIV_ID_from',
                          suffix=c('_upst', '_extra')) %>%
    .[, dist_km_acc := dist_km_upst + dist_km_extra] %>%
    .[, c('HYRIV_ID_from', 'HYRIV_ID_to_extra', 
          'ORD_STRA_from_upst', 'ORD_STRA_to_extra',
          'dist_km_acc'), with=F] %>%
    setnames(c('HYRIV_ID_to_extra', 'dist_km_acc', 'ORD_STRA_from_upst', 'ORD_STRA_to_extra'),
             c('HYRIV_ID_to', 'dist_km', 'ORD_STRA_from', 'ORD_STRA_to'))=
  all_acc <- rbindlist(list(all_acc, k3_down, acc_to_k3_down), use.names=T)
}








####################3

#Loop through each segment level, keeping only those downstream 
#of a segment that must be kept
while (length(new_ids) > 0) {
  #print(length(new_ids))
  #Keep all segments downstream of those already kept
  new_ids <- netdt[get(id_colname) %in% new_ids, 
                   get(nextid_colname)]
  ids_keep <- c(ids_keep, new_ids)
  #Only keep segments which still have another segment upstream
  netdt <- netdt[get(id_colname) %in% 
                   c(get(nextid_colname), ids_keep),]
}

k2_to_k3 <- rivers_inland[ORD_STRA==2 & !is.na(LENGTH_KM_down), ..sel_cols] %>% 
  setnames(sel_cols, new_cols)

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
