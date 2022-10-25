### STEP-BY-STEP TESTING FOR POST-PROCESSING

## 'isolation by distance' Post-processing in R --------

# for these taxa Circuitscape processing steps are not needed 
## if taxonA$filter_category == "isolation by distance"
# remove the first column ('precluster')?
# clusters_info32 <- pclust_info31[,-1] # 5 columns

### NEED TO PERFORM ALL STEPS THAT RELY ON geometry FIRST,
###  THEN DO st_drop_geometry?

# TURN ALL THIS INTO A FUNCTION?
##  '31' IN THE CODE IS JUST FOR TESTING..

# add new columns, to be filled or updated subsequently
pclust_info31 <- pclust_info31 |> add_column(pop_name = NA,
          gene_div_special = NA, pix_ignore = 0, Ne_override = NA,
          gene_div_weight = 0, proximity = 0) # 13 columns

# update populations for 'isolation by distance' taxa --------

### If and when available, record expert opinion or info derived from
# genetic data on effective size of, or proportion of genetic diversity
# contained within, a given local population, along with relevant local
# population (precluster) names in a supplementary file

# import expert opinion estimates file for given taxon if file exists
expert1_filename <- suppressWarnings(file.path(taxon_path(taxonA, taxapath), 
                paste0("expert_op1_preclust", ".csv")))
if (file.exists(expert1_filename)){
  expert1 <- read.csv(expert1_filename, header = TRUE)
  if (nrow(expert1) == nrow(pclust_info31)) {
    pclust_info31$pop_name <- expert1$pop_name
    pclust_info31$gene_div_special <- expert1$gene_div_special
    pclust_info31$Ne_override <- expert1$Ne_override
  }
  else {cat("\nWrong number of populations in expert info file!\n")
    }
}

# then for any clusters where 'gene_div_special' is not NA; 
# copy value of 'pix_count' column to 'pix_ignore' column
for (i in 1:nrow(pclust_info31)) {
  if(!is.na(pclust_info31$gene_div_special[i])) {
    pclust_info31$pix_ignore[i] <- pclust_info31$pix_count[i]
  }
}

# code or function to calculate 'gen_diverse_weight'..
for (i in 1:nrow(pclust_info31)) {
  if(!is.na(pclust_info31$gene_div_special[i])) {
    pclust_info31$gene_div_weight[i] <- pclust_info31$gene_div_special[i]
  }
  else {
    pclust_info31$gene_div_weight[i] <-
    pclust_info31$pix_count[i] / (sum(pclust_info31$pix_count) - 
    sum(pclust_info31$pix_ignore))*(1-sum(pclust_info31$gene_div_special,
        na.rm = TRUE))
  }
}

# add another column for which is the nearest polygon
nearest <- sf::st_nearest_feature(pclust_info31)
pclust_info31 <- cbind(pclust_info31, nearest) # 14 columns

## update proximity column with nearest distance
# slice needs to operate on a data frame
for (i in 1:nrow(pclust_info31)) {
  pclust_info31$proximity[i] <- as.data.frame(prox31[,i]) |> 
    dplyr::slice(nearest[i])
}
# The above makes $proximity a list. Need to convert back to numeric
pclust_info31$proximity <- as.numeric(pclust_info31$proximity)

# add columns for weighted size of precluster, epsilon proximity, and
## effective (Victorian) population size (rounded to whole numbers)
## Vic_popn_size should be == 0 if extinct in Vic or unknown
pclust_info31 <- pclust_info31 |> 
  mutate(pix_pc = (pix_count / sum(pix_count))) |>
  mutate(eps_prox = proximity / taxonA$epsilon) |> 
  mutate(Ne_est = pix_pc * taxonA$vic_popn_size * 0.1)
pclust_info31$Ne_est <- round(pclust_info31$Ne_est, digits = 0)

# update with Ne values manually entered for particular populations
for (i in 1:nrow(pclust_info31)) {
  if(!is.na(pclust_info31$Ne_override[i])) {
    pclust_info31$Ne_est[i] <- pclust_info31$Ne_override[i]
  }
}

# add risk score calculation columns
pclust_info31 <- pclust_info31 |> add_column(size_nearest = 0,
        gene_flow_factor = 0, div_risk = 0, inbreed_risk = 0,
        cluster_score = 0) # should now be 22 columns
# derive size of nearest precluster..
for (i in 1:nrow(pclust_info31)){
  pclust_info31$size_nearest[i] <- dplyr::slice(pclust_info31,
        pclust_info31$nearest[i]) |> 
    dplyr::select(pix_count) |> sf::st_drop_geometry()
}

# derive 'gene_flow_factor'
### TO DO: MOVE 2.6 and 2.2 INTO config.toml file
for (i in 1:nrow(pclust_info31)){
  if(pclust_info31$eps_prox[i] > 2.6){
    pclust_info31$gene_flow_factor[i] <- 0
  }
  else {
    pclust_info31$gene_flow_factor[i] <- 
      1/((pclust_info31$eps_prox[i]+1)^2.2)
  }
}
# code to derive diversity risk based on 'Ne_est' & 'gene_flow_factor'
## TO DO: WHAT IF NEAREST PRECLUSTER IS LARGE (OR SMALL)?
### TO DO: MOVE 1000 INTO config.toml file
for (i in 1:nrow(pclust_info31)) {
  if(pclust_info31$Ne_est[i] < 1000) {
    (pclust_info31$div_risk[i] <-
      (1 - (pclust_info31$Ne_est[i]/1000)) * 100) *
      (1 - pclust_info31$gene_flow_factor[i])
  }
  else {
    pclust_info31$div_risk[i] <- 0
  }
}

# derive 'adj_inbreed_risk'
### To do: MOVE 110 INTO config.toml file
for (i in 1:nrow(pclust_info31)) {
  if(pclust_info31$Ne_est[i] < 110) {
    if(pclust_info31$Ne_est[i] < 10) {
      pclust_info31$inbreed_risk[i] <-
        (2000 * (1 - pclust_info31$gene_flow_factor[i]))
    }
    else {
      pclust_info31$inbreed_risk[i] <-
        ((110-pclust_info31$Ne_est[i])^2)/5 *
        (1-pclust_info31$gene_flow_factor[i])
    }
  }
  else {
    pclust_info31$inbreed_risk[i] <- 0
  }
}

## TO DO: ADJUST 'precluster_score' BASED ON 'recent_year'?
## WHAT PROPORTION OF OBS / PIXELS ARE BASED ON QUITE OLD RECORDS?
## What is the likelihood that for preclusters where all records are old,
##  the taxon is now locally extinct??
## TO DO: Also derive an uncertainty score per precluster?
##  or just for the taxon overall?

# code to derive 'precluster_score'
for (i in 1:nrow(pclust_info31)) {
  pclust_info31$cluster_score[i] <-
    (pclust_info31$div_risk[i] + pclust_info31$inbreed_risk[i]) *
    pclust_info31$gene_div_weight[i]
}
# fragmented risk for a given taxon will be sum(pclust_info31$cluster_score)
## TO DO: DATA BASE WILL NEED A NEW COLUMN FOR THIS,
##  AND THEN ANOTHER NEW COLUMN FOR THE FINAL RISK SCORE AFTER ADJUSTING 
##  FOR VARIOUS OTHER FACTORS

# write preclusters with derived risk scores to a file
## use 'taxon' in place of 'taxonA' for main obs4.R script
pclust_info31 |> sf::st_drop_geometry() |> 
  write_csv(file.path(taxonpath, paste0(gsub(" ","_",
  taxonA$ala_search_term), "_preclusters_info31", ".csv")))


## import first Circuitscape output matrix --------

# derive clusters by iteratively grouping preclusters based on
# relative resistances calculated with Circuitscape
# will be space delimited (" "), so use sep = " "
pwise1 <- read.table(file.path(taxonpath,
            paste0("CSpwise1_resistances",".out")), header = TRUE,
            sep = " ", row.names = 1, check.names = FALSE)
# rownames are "1", "2", etc
# but colnames are "1.0", "2.0", "3.0" etc
colnames(pwise1) <- rownames(pwise1)

## TO DO: What is a sufficiently robust method to determine an appropriate
##  minimum resistance thresh for an initial pair of preclusters
##  to be considered fully connected?
### MAKE THIS SOME MULTIPLE OF EPSILON?
###  how much variance to allow for between high & low dispersal taxa??
### Use a different value when using the generic resistance (habitat)
##  layer as opposed to when using a species-specific one?

# run functions from igraph package to group preclusters into clusters
## TO DO: What is an appropriate level for res_factor?
eps_factor2 <- (taxonA$disperse_log * 3) + 6.3
connect1 <- which(pwise1 < res_factor, arr.ind = TRUE)
grph1 <- igraph::graph_from_data_frame(connect1, directed = FALSE)
groups1 <- split(unique(as.vector(connect1)),
              igraph::clusters(grph1)$membership)
fcluster1 <- lapply(groups1,
              FUN = function(list.cor){rownames(pwise1)[list.cor]})

## This will have returned the index of the relevant list of 
##  grouped preclusters, so use the result as the new cluster number

for (i in 1:nrow(pwise1)) {
  pclust_info31$cluster[i] <- which(lapply(fcluster1,
            function(x) grep(paste0("^",i,"$"),x))!=0) |> unname()
}

# also make a new patches .tif file for newly grouped clusters
clust_retest31 <- preclustered_obs31 |> add_column(cluster = NA)
for (i in 1:nrow(pwise1)) {
  clust_retest31$cluster[i] <- which(lapply(fcluster1,
            function(x) grep(paste0("^",i,"$"),x))!=0) |> unname()
}

clust_retest31 <- clust_retest31 |> dplyr::group_by(cluster) |> 
  dplyr::summarize(geometry = sf::st_union(geometry))

shapevect31 <- terra::vect(clust_retest31)
obs_retest_rast31 <- terra::rasterize(shapevect31, mask_layer,
            field = "cluster")
obs_test_filename <- file.path(taxonpath, paste0("clusters31", ".tif"))
terra::crop(obs_retest_rast31, crop_rast31,
            filename = obs_test_filename, overwrite = TRUE)

## run Circuitscape again --------

## NOW USE "clusters31.tif" as a new patches layer for a 2nd Circuitscape run
# read in the vanilla version of the Circuitscape.ini file with standard
# settings for a Circuitscape pairwise run
Cscape1path <- file.path(datapath, paste0("circuitscape_pwise.ini"))
CSrun <- read.ini(Cscape1path, encoding = getOption("encoding"))
# produces a list of 11

# update for source files relevant to the particular taxon:
CSrun$`Habitat raster or graph`$habitat_file <- 
  file.path(taxonpath, paste0("resistance.tif"))

CSrun$`Options for pairwise and one-to-all and all-to-one modes`$point_file <-
  file.path(taxonpath, paste0("clusters31.tif"))

CSrun$`Output options`$output_file <-
  file.path(taxonpath, paste0("CSpwise2"))

# save as a customized .ini file in the relevant directory
write.ini(CSrun, file.path(taxonpath, paste0("Circuitscape_custom2.ini")),
          encoding = getOption("encoding"))


## Run from Julia:
# using Circuitscape
# compute("/home/peter/data/taxa/Varanus_varius/Circuitscape_custom2.ini")

### both all-to-one and one-to-all modes fail(!?), but pairwise mode ok


## import second Circuitscape output matrix --------
pwise2 <- read.table(file.path(taxonpath,
            paste0("CSpwise2_resistances_3columns", ".out")), header = TRUE,
            sep = " ", col.names = c("cluster", "closest", "resistance"),
            check.names = FALSE)
# transform to find which is the next nearest cluster by resistance
flip <- as.data.frame(cbind(pwise2[,2], pwise2[,1], pwise2[,3]))
names(flip)<-names(pwise2)
pwise2 <- rbind(pwise2, flip) |> dplyr::group_by(cluster) |> 
  dplyr::slice(which.min(resistance)) # returns only the 1st minima


## 'isolation by resistance' Post-processing in R --------

## if taxonA$filter_category == "isolation by resistance"

# start with 7 column version of preclusters info file
## TO DO: store as a variable in obs4.R that can be retrieved here

# update the cluster number based on initial Circuitscape processing
for (i in 1:nrow(pwise1)) {
  pclust_info31$cluster[i] <- which(lapply(fcluster1,
      function(x) grep(paste0("^",i,"$"),x))!=0) |> unname()
}
# create summary table based on cluster number
clusters_info31 <- pclust_info31 |> dplyr::group_by(cluster) |> 
  dplyr::summarize(pix_count = sum(pix_count), num_obs = sum(num_obs),
        recent_year = max(recent_year),
        latin = max(latin), geometry = sf::st_union(geometry)) # 6 columns

# add new columns, including for next nearest cluster by resistance
clusters_info31 <- clusters_info31 |> add_column(pop_name = NA,
        gene_div_special = NA, pix_ignore = 0, Ne_override = NA,
        gene_div_weight = 0, lowest_resist = NA, closest = NA) # 13 columns

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
expert2_filename <- suppressWarnings(file.path(taxon_path(taxonA, taxapath), 
            paste0("expert_op2_clusters", ".csv")))
if (file.exists(expert2_filename)){
  expert2 <- read.csv(expert2_filename, header = TRUE)
  if (nrow(expert2) == nrow(clusters_info31)) {
    clusters_info31$pop_name <- expert2$pop_name
    clusters_info31$gene_div_special <- expert2$gene_div_special
    clusters_info31$Ne_override <- expert2$Ne_override
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
  mutate(Ne_est = pix_pc * taxonA$vic_popn_size * 0.1)
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
      cluster_score = 0) # 20 columns

# derive size of nearest cluster..
for (i in 1:nrow(clusters_info31)){
  clusters_info31$size_nearest[i] <- dplyr::slice(clusters_info31,
      clusters_info31$closest[i]) |> 
    dplyr::select(pix_count) |> sf::st_drop_geometry()
}

# derive 'gene_flow_factor'
## TO FIX: Is this method appropriate?
### TO DO: 2.2 INTO config.toml file?
for (i in 1:nrow(clusters_info31)){
  if(clusters_info31$lowest_resist[i] > (res_factor * 1.75)){
    clusters_info31$gene_flow_factor[i] <- 0
  } else if (clusters_info31$lowest_resist[i] < res_factor) {
    clusters_info31$gene_flow_factor[i] <- 1
  } else {
    clusters_info31$gene_flow_factor <-
      (((res_factor * 2)-clusters_info31$lowest_resist)*(1/res_factor))^2.2
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
    taxonA$ala_search_term), "_clusters_info31", ".csv")))


# fragmented risk for a given taxon will be sum(clusters_info31$cluster_score)
## DATA BASE WILL NEED A NEW COLUMN FOR THIS,
##  AND THEN ANOTHER NEW COLUMN FOR THE FINAL RISK SCORE AFTER ADJUSTING 
##  FOR VARIOUS OTHER FACTORS


