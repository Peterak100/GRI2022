### STEP-BY-STEP TESTING FOR POST-PROCESSING

## 'isolation by distance' Post-processing in R --------

# run function on isolation_by_distance_taxa
# for these taxa Circuitscape processing steps are not needed 

# run on 'isolation_by_distance_taxa'
##  '31' IN THE CODE BELOW IS JUST FOR TESTING..

taxonA <- isolation_taxa[9, ]
taxonA[c(1,2)]

# import preclusters preliminary dataframe for given taxon
pre_clusters31 <- st_read(file.path(taxon_path(taxonA, taxapath),
        paste0(gsub(" ","_", taxonA$ala_search_term),
        "_preclusters_prelim", ".shp"))) # 7 columns

# column names in saved '.shp' file are truncated to 7 characters so rename
pre_clusters31 <- pre_clusters31 |> 
  rename(midcluster = mdclstr, precluster = prclstr, pix_count = pix_cnt,
         recent_year = rcnt_yr)
# add a new first column as placeholder for cluster number
pre_clusters31 <- cbind(cluster = pre_clusters31$midcluster, pre_clusters31)

# add new columns, to be filled or updated subsequently
pre_clusters31 <- pre_clusters31 |> add_column(pop_name = NA,
          gene_div_special = NA, pix_ignore = 0, Ne_override = NA,
          gene_div_weight = 0, proximity = 0) # 14 columns

## update populations for 'isolation by distance' taxa --------

### If and when available, record expert opinion
## or info derived from genetic data on
## effective pop'n size or proportion of genetic diversity contained within
## a given local population, along with relevant
## local population (precluster) names in a supplementary file

# import expert opinion estimates file for given taxon if file exists
expert1_filename <- suppressWarnings(file.path(taxon_path(taxonA, taxapath), 
                paste0("expert_op1_preclust", ".csv")))
expert2_filename <- suppressWarnings(file.path(taxon_path(taxonA, taxapath), 
                paste0("expert_op2_midclust", ".csv")))

if (file.exists(expert2_filename)){
  expert_op <- read.csv(expert2_filename, header = TRUE)
  } else {
    if (file.exists(expert1_filename)){
      expert_op <- read.csv(expert1_filename, header = TRUE)
      }
}

if (exists("expert_op")) {
  if (nrow(expert_op) == nrow(pre_clusters31)) {
    pre_clusters31$pop_name <- expert_op$pop_name
    pre_clusters31$gene_div_special <- expert_op$gene_div_special
    pre_clusters31$Ne_override <- expert_op$Ne_override
  } else {
    cat("\nWrong number of populations in expert info file!\n")
  }
}    


# then for any clusters where 'gene_div_special' is not NA; 
# copy value of 'pix_count' column to 'pix_ignore' column
for (i in 1:nrow(pre_clusters31)) {
  if(!is.na(pre_clusters31$gene_div_special[i])) {
    pre_clusters31$pix_ignore[i] <- pre_clusters31$pix_count[i]
  }
}

# code or function to calculate 'gen_diverse_weight'..
for (i in 1:nrow(pre_clusters31)) {
  if(!is.na(pre_clusters31$gene_div_special[i])) {
    pre_clusters31$gene_div_weight[i] <- pre_clusters31$gene_div_special[i]
  }
  else {
    pre_clusters31$gene_div_weight[i] <-
    pre_clusters31$pix_count[i] / (sum(pre_clusters31$pix_count) - 
    sum(pre_clusters31$pix_ignore))*(1-sum(pre_clusters31$gene_div_special,
        na.rm = TRUE))
  }
}

# add another column for which is the nearest polygon
nearest <- sf::st_nearest_feature(pre_clusters31)
pre_clusters31 <- cbind(pre_clusters31, nearest) # 15 columns

# create a 'units' matrix of distances between polygons
# then convert to a normal numeric matrix & convert to kilometres
prox31 <- sf::st_distance(pre_clusters31) |> as.data.frame() |> as.matrix()
prox31 <- prox31/1000

## update proximity column with nearest distance
# slice needs to operate on a data frame
for (i in 1:nrow(pre_clusters31)) {
  pre_clusters31$proximity[i] <- as.data.frame(prox31[,i]) |> 
    dplyr::slice(nearest[i])
}
# The above makes $proximity a list. Need to convert back to numeric
pre_clusters31$proximity <- as.numeric(pre_clusters31$proximity)

# add 3 columns for weighted size of precluster, epsilon proximity, and
## effective (Victorian) population size (rounded to whole numbers)
## Vic_popn_size should be == 0 if extinct in Vic or unknown
pre_clusters31 <- pre_clusters31 |> 
  mutate(pix_pc = (pix_count / sum(pix_count))) |>
  mutate(eps_prox = proximity / taxonA$epsilon) |> 
  mutate(Ne_est = pix_pc * taxonA$vic_popn_size * 0.1)
pre_clusters31$Ne_est <- round(pre_clusters31$Ne_est, digits = 0)
# 18 columns

# update with Ne values manually entered for particular populations
for (i in 1:nrow(pre_clusters31)) {
  if(!is.na(pre_clusters31$Ne_override[i])) {
    pre_clusters31$Ne_est[i] <- pre_clusters31$Ne_override[i]
  }
}

# add risk score calculation columns
pre_clusters31 <- pre_clusters31 |> add_column(size_nearest = 0,
        gene_flow_factor = 0, div_risk = 0, inbreed_risk = 0,
        cluster_score = 0) # should now be 23 columns
# derive size of nearest precluster..
for (i in 1:nrow(pre_clusters31)){
  pre_clusters31$size_nearest[i] <- dplyr::slice(pre_clusters31,
        pre_clusters31$nearest[i]) |> 
    dplyr::select(pix_count) |> sf::st_drop_geometry()
}

# derive 'gene_flow_factor'
### TO DO: MOVE 2.6 and 2.2 INTO config.toml file
for (i in 1:nrow(pre_clusters31)){
  if(pre_clusters31$eps_prox[i] > 2.6){
    pre_clusters31$gene_flow_factor[i] <- 0
  }
  else {
    pre_clusters31$gene_flow_factor[i] <- 
      1/((pre_clusters31$eps_prox[i]+1)^2.2)
  }
}
# code to derive diversity risk based on 'Ne_est' & 'gene_flow_factor'
## TO DO: WHAT IF NEAREST PRECLUSTER IS LARGE (OR SMALL)?
### TO DO: MOVE 1000 INTO config.toml file
for (i in 1:nrow(pre_clusters31)) {
  if(pre_clusters31$Ne_est[i] < 1000) {
    (pre_clusters31$div_risk[i] <-
      (1 - (pre_clusters31$Ne_est[i]/1000)) * 100) *
      (1 - pre_clusters31$gene_flow_factor[i])
  }
  else {
    pre_clusters31$div_risk[i] <- 0
  }
}

# derive 'adj_inbreed_risk'
### To do: MOVE 110 INTO config.toml file
for (i in 1:nrow(pre_clusters31)) {
  if(pre_clusters31$Ne_est[i] < 110) {
    if(pre_clusters31$Ne_est[i] < 10) {
      pre_clusters31$inbreed_risk[i] <-
        (2000 * (1 - pre_clusters31$gene_flow_factor[i]))
    }
    else {
      pre_clusters31$inbreed_risk[i] <-
        ((110-pre_clusters31$Ne_est[i])^2)/5 *
        (1-pre_clusters31$gene_flow_factor[i])
    }
  }
  else {
    pre_clusters31$inbreed_risk[i] <- 0
  }
}

## TO DO: ADJUST 'precluster_score' BASED ON 'recent_year'?
## WHAT PROPORTION OF OBS / PIXELS ARE BASED ON QUITE OLD RECORDS?
## What is the likelihood that for preclusters where all records are old,
##  the taxon is now locally extinct??
## TO DO: Also derive an uncertainty score per precluster?
##  or just for the taxon overall?

# code to derive 'precluster_score'
for (i in 1:nrow(pre_clusters31)) {
  pre_clusters31$cluster_score[i] <-
    (pre_clusters31$div_risk[i] + pre_clusters31$inbreed_risk[i]) *
    pre_clusters31$gene_div_weight[i]
}
# fragmented risk for a given taxon will be sum(pre_clusters31$cluster_score)
# frag_risk31 <- sum(pre_clusters31$cluster_score)
## TO DO: DATA BASE WILL NEED A NEW COLUMN FOR THIS,
##  AND THEN ANOTHER NEW COLUMN FOR THE FINAL RISK SCORE AFTER ADJUSTING 
##  FOR VARIOUS OTHER FACTORS

# write preclusters with derived risk scores to a file
## use 'taxon' in place of 'taxonA' for main obs4.R script
pre_clusters31 |> sf::st_drop_geometry() |> 
  write_csv(file.path(taxonpath, paste0(gsub(" ","_",
  taxonA$ala_search_term), "_preclusters31_info", ".csv")))



## import 1st Circuitscape output matrix --------
# need to import the saved 'sf' object

taxonpath <- taxon_path(taxonA, taxapath)
pre_clusters32 <- sf::st_read(file.path(taxonpath,
          paste0(gsub(" ","_", taxon$ala_search_term),
          "_preclusters_prelim", ".shp")), quiet = TRUE) # 7 columns
# column names in saved '.shp' file are truncated to 7 characters so rename
pre_clusters32 <- pre_clusters32 |> 
  rename(midcluster = mdclstr, precluster = prclstr, pix_count = pix_cnt,
         recent_year = rcnt_yr)
# add a new first column as placeholder for cluster number
pre_clusters32 <- cbind(cluster = pre_clusters32$midcluster, pre_clusters32)
# 8 columns

# derive clusters by iteratively grouping preclusters based on
# relative resistances calculated with Circuitscape
# will be space delimited (" "), so use sep = " "
CScape_FILE1 <- paste0("CSpwise1_resistances", ".out")
CScape_PATH1 <- file.path(taxonpath, CScape_FILE1)
if (file.exists(CScape_PATH1)){
  pwise1 <- read.table(CScape_PATH1, header = TRUE, sep = " ",
                row.names = 1, check.names = FALSE)
  # rownames are "1", "2", etc
  # but colnames are "1.0", "2.0", "3.0" etc
  colnames(pwise1) <- rownames(pwise1)
  # use different resistance threshold for generic resistance layer?
  SMP_MODEL <- paste0("SMP_", gsub(" ","_", taxon$ala_search_term), ".tif")
  SMP_filename <- suppressWarnings(file.path(SMPpath, SMP_MODEL))
  HIM_MODEL <- paste0("HIM_", gsub(" ","_", taxon$ala_search_term), ".tif")
  HIM_filename <- suppressWarnings(file.path(HIMpath, HIM_MODEL))
  if(file.exists(SMP_filename) | file.exists(HIM_filename)) {
    resist_factor <- (taxon$disperse_log * 3) + 8.8
  } else {
    resist_factor <- (taxon$disperse_log * 3) + 6.3
  }
  # run functions from igraph package to group preclusters into clusters
  connect1 <- which(pwise1 < resist_factor, arr.ind = TRUE)
  grph1 <- igraph::graph_from_data_frame(connect1, directed = FALSE)
  groups1 <- split(unique(as.vector(connect1)),
            igraph::clusters(grph1)$membership)
  fcluster1 <- lapply(groups1,
            FUN = function(list.cor){rownames(pwise1)[list.cor]})
  for (j in 1:nrow(pwise1)) {
    pre_clusters32$cluster[j] <- which(lapply(fcluster1,
            function(x) grep(paste0("^",j,"$"),x))!=0) |> unname()
  }
  } else {
  cat("Circuitscape results NOT FOUND for",
           taxon$ala_search_term, "\n")
}
# create summary table based on cluster number
clusters_info31 <- pre_clusters32 |> dplyr::group_by(cluster) |> 
  dplyr::summarise(midcluster = paste(midcluster, collapse = ", "),
          precluster = paste(precluster, collapse = ", "),
          pix_count = sum(pix_count), recent_year = max(recent_year),
          latin = max(latin), num_obs = sum(num_obs),
          geometry = sf::st_union(geometry)) # also 8 columns
# also make a new patches .tif file for newly grouped clusters
shapevect31 <- terra::vect(clusters_info31) # plot(shapevect31)
obs_retest_rast <- terra::rasterize(shapevect31, mask_layer,
          field = "cluster") # plot(obs_retest_rast31)
obs_retest_filename <- file.path(taxonpath, paste0("clusters", ".tif"))
crop_final_rast <- obs_retest_rast |> padded_trim()
terra::crop(obs_retest_rast, crop_final_rast,
            filename = obs_retest_filename, overwrite = TRUE)
# output clusters_info31 to a variable
return(clusters_info31)

## TO DO: is the above a sufficiently robust method to determine
##  a 'resist_factor' (minimum resistance threshold for an initial pair
##  of preclusters to be considered fully connected)?
### MAKE THIS SOME MULTIPLE OF EPSILON?
###  how much variance to allow for between high & low dispersal taxa??


## Run from Julia:
# using Circuitscape


## import 2nd Circuitscape output matrix --------

## TO DO: make a separate function from this
pwise2 <- read.table(file.path(taxonpath,
        paste0("CSpwise2_resistances_3columns", ".out")), header = TRUE,
        sep = " ", col.names = c("cluster", "closest", "resistance"),
        check.names = FALSE)
# transform to find which is the next nearest cluster by resistance
flip <- as.data.frame(cbind(pwise2[,2], pwise2[,1], pwise2[,3]))
names(flip)<-names(pwise2)
# this returns only the 1st minima (just in case there is more than one)
pwise2 <- rbind(pwise2, flip) |> dplyr::group_by(cluster) |> 
  dplyr::slice(which.min(resistance))


## 'isolation by resistance' Post-processing in R --------

# start with 8 column version of clusters info file
# add new columns, including for next nearest cluster by resistance
clusters_info31 <- clusters_info31 |> add_column(pop_name = NA,
        gene_div_special = NA, pix_ignore = 0, Ne_override = NA,
        gene_div_weight = 0, lowest_resist = NA, closest = NA) # 15 columns

# update lowest_resist column with next lowest resistance
for (i in 1:nrow(clusters_info31)) {
  clusters_info31$lowest_resist[i] <- pwise2$resistance[i]
}

# update populations for 'isolation by resistance' taxa --------

### If and when available, record expert opinion or info derived from
# genetic data on effective size of, or proportion of genetic diversity
# contained within, a given local population, along with relevant local
# population (precluster) names in a supplementary file

# import expert opinion estimates file for given taxon if file exists
expert3_filename <- suppressWarnings(file.path(taxon_path(taxon, taxapath), 
            paste0("expert_op3_clusters", ".csv")))
if (file.exists(expert3_filename)){
  expert3 <- read.csv(expert3_filename, header = TRUE)
  if (nrow(expert3) == nrow(clusters_info31)) {
    clusters_info31$pop_name <- expert3$pop_name
    clusters_info31$gene_div_special <- expert3$gene_div_special
    clusters_info31$Ne_override <- expert3$Ne_override
  }
  else {cat("\nWrong number of populations in expert info file!\n")
  }
}

## then for any clusters where 'gene_div_special' is not NA; 
## copy value of 'pix_count' column to 'pix_ignore' column
for (i in 1:nrow(clusters_info31)) {
  if(!is.na(clusters_info31$gene_div_special[i])) {
    clusters_info31$pix_ignore[i] <- clusters_info31$pix_count[i]
  }
}

# code or function to calculate 'gen_diverse_weight'..
for (i in 1:nrow(clusters_info31)) {
  if(!is.na(clusters_info31$gene_div_special[i])) {
    clusters_info31$gene_div_weight[i] <- clusters_info31$gene_div_special[i]
  }
  else {
    clusters_info31$gene_div_weight[i] <-
    clusters_info31$pix_count[i] / (sum(clusters_info31$pix_count) - 
    sum(clusters_info31$pix_ignore))*(1-sum(clusters_info31$gene_div_special,
          na.rm = TRUE))
  }
}

# update closest column with next closest cluster by resistance
for (i in 1:nrow(clusters_info31)) {
  clusters_info31$closest[i] <- pwise2$closest[i]
}

## TO DO: NEED NEW METHOD TO DERIVE 'eps_prox' HERE ##
### using 'lowest_resist' instead of 'proximity'
### lowest_resist / taxonA$epsilon ?# WHAT RANGE OF 'lowest_resist TO EXPECT?
## potentially quite low???

# add columns for weighted size of precluster, epsilon factor (used to 
##  calculate 'gene_flow_factor'), and effective (Victorian) population
##  size (rounded to whole numbers)
## Vic_popn_size should be == 0 if extinct in Vic or unknown
clusters_info31 <- clusters_info31 |> 
  mutate(pix_pc = (pix_count / sum(pix_count))) |>
  mutate(Ne_est = pix_pc * taxon$vic_popn_size * 0.1)
clusters_info31$Ne_est <- round(clusters_info31$Ne_est, digits = 0)
# 15 columns

# update with Ne values manually entered for particular populations
for (i in 1:nrow(clusters_info31)) {
  if(!is.na(clusters_info31$Ne_override[i])) {
    clusters_info31$Ne_est[i] <- clusters_info31$Ne_override[i]
  }
}

# add risk score calculation columns
clusters_info31 <- clusters_info31 |> add_column(size_nearest = 0,
      gene_flow_factor = 0, div_risk = 0, inbreed_risk = 0,
      cluster_score = 0) # 22 columns

# derive size of nearest cluster..
for (i in 1:nrow(clusters_info31)){
  clusters_info31$size_nearest[i] <- dplyr::slice(clusters_info31,
      clusters_info31$closest[i]) |> 
    dplyr::select(pix_count) |> sf::st_drop_geometry()
}

# derive 'gene_flow_factor'
## TO FIX: Is this method appropriate?
### TO DO: move 2.2 INTO config.toml file?
res_factor2 <- (taxon$disperse_log * 3) + 6.3
for (i in 1:nrow(clusters_info31)){
  if(clusters_info31$lowest_resist[i] > (res_factor2 * 1.75)){
    clusters_info31$gene_flow_factor[i] <- 0
  } else if (clusters_info31$lowest_resist[i] < res_factor2) {
    clusters_info31$gene_flow_factor[i] <- 1
  } else {
    clusters_info31$gene_flow_factor <-
      (((res_factor2 * 2)-clusters_info31$lowest_resist) *
         (1/res_factor2))^2.2
    }
}

# code to derive diversity risk based on 'Ne_est' & 'gene_flow_factor'
## WHAT IF NEAREST PRECLUSTER IS LARGE (OR SMALL)?
### TO DO: MOVE 1000 INTO config.toml file
for (i in 1:nrow(clusters_info31)) {
  if(clusters_info31$Ne_est[i] < 1000) {
    (clusters_info31$div_risk[i] <-
       (1 - (clusters_info31$Ne_est[i]/1000)) * 100) *
      (1 - clusters_info31$gene_flow_factor[i])
  }
  else {
    clusters_info31$div_risk[i] <- 0
  }
}

# derive 'adj_inbreed_risk'
### To do: MOVE 110 INTO config.toml file
for (i in 1:nrow(clusters_info31)) {
  if(clusters_info31$Ne_est[i] < 110) {
    if(clusters_info31$Ne_est[i] < 10) {
      clusters_info31$inbreed_risk[i] <-
        (2000 * (1 - clusters_info31$gene_flow_factor[i]))
    }
    else {
      clusters_info31$inbreed_risk[i] <-
        ((110-clusters_info31$Ne_est[i])^2)/5 *
        (1-clusters_info31$gene_flow_factor[i])
    }
  }
  else {
    clusters_info31$inbreed_risk[i] <- 0
  }
}

# code to derive 'cluster_score'
for (i in 1:nrow(clusters_info31)) {
  clusters_info31$cluster_score[i] <-
    (clusters_info31$div_risk[i] + clusters_info31$inbreed_risk[i]) *
    clusters_info31$gene_div_weight[i]
}

# write clusters with derived risk scores to a file
## use 'taxon' in place of 'taxonA' for main obs4.R script
clusters_info31 |> sf::st_drop_geometry() |> 
  write_csv(file.path(taxonpath, paste0(gsub(" ","_",
    taxon$ala_search_term), "_clusters_info31", ".csv")))


# fragmented risk for a given taxon will be sum(clusters_info31$cluster_score)
## DATA BASE WILL NEED A NEW COLUMN FOR THIS,
##  AND THEN ANOTHER NEW COLUMN FOR THE FINAL RISK SCORE AFTER ADJUSTING 
##  FOR VARIOUS OTHER FACTORS


