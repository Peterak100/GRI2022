# i.e. derive which preclusters to group into clusters first
# based on relative resistances calculated with Circuitscape

## import Circuitscape output matrix:
# will be space delimited (" "), so use sep = " "
# examples:
pwise2 <- read.table("Dusky1_resistances.out",header=TRUE,
      sep = " ", row.names = 1, check.names = FALSE)
# for P_cryodrama
pwise1 <- read.table("cryo1_resistances.out",header=TRUE,
      sep = " ", row.names = 1, check.names = FALSE)
# rownames are "1", "2", etc
# but colnames are "1.0", "2.0", "3.0" etc
colnames(pwise1) <- rownames(pwise1)

### NEED TO DETERMINE A RESISTANCE THRESHOLD FOR PRECLUSTERS TO BE
###  CONSIDERED FULLY CONNECTED
### MAKE THIS SOME MULTIPLE OF EPSILON? ## is 1.5 times good???
### DIFFERENT VALUES WHEN USING GENERIC RESISTANCE LAYER?
eps_factor1 <- taxonA$epsilon * 1.5 ## use 1.5 here?
# make a function out of this?
connect1 <- which(pwise1 < eps_factor1, arr.ind = TRUE)
grph1 <- igraph::graph_from_data_frame(connect1, directed = FALSE)
groups1 <- split(unique(as.vector(connect1)),
          igraph::clusters(grph1)$membership)
fcluster1 <- lapply(groups1,
          FUN = function(list.cor){rownames(pwise1)[list.cor]})

### THIS RETURNS THE INDEX OF THE RELEVANT LIST OF PRECLUSTERS
## SO USE THE RESULT AS THE NEW CLUSTER NUMBER
## Make a function out of this

for (i in 1:nrow(pwise1)) {
  pclust_info31$cluster[i] <- which(lapply(fcluster1,
          function(x) grep(paste0("^",i,"$"),x))!=0) |> unname()
}

## Post-processing in R --------

# for 'isolation by distance' taxa: post-processing steps should be ok
# for 'isolation by resistance' taxa:
# next need to summarize for each cluster - using dplyr::group_by()?
## sum(pix_count), max(recent_year), max(latin)
###  NEED TO GROUP THE geometry ?? OKAY TO BECOME DISCONTINUOUS SHAPES?
###   dplyr::summarize(geometry = st_union(geometry))
###  run circuitscape again (in all-to-one mode?) 
###  to get new relative resistances?
# https://gis.stackexchange.com/questions/316181/how-do-i-combine-geometries-in-a-shapefile-based-on-a-grouping-variable

### NEED TO PERFORM ALL STEPS THAT RELY ON geometry FIRST,
###  THEN DO st_drop_geometry?

mtcars |> group_by(cyl) |> summarize(cyl_sum = sum(cyl), mpg = max(mpg))


## TURN ALL THIS INTO A FUNCTION?
##  '31' IN THE CODE IS JUST FOR TESTING..

# manually update 'pop_name','gene_div_special' & 'Ne_override' columns;
# e.g. pclust_info31$pop_name[7] <- "Gibbo"
# e.g. pclust_info31$gene_div_special[7] <- 0.35
# e.g. pclust_info31$Ne_override[7] <- 300

## then for any clusters where 'gene_div_special' > 0; 
## copy value of 'pix_count' column to 'pix_ignore' column
for (i in 1:nrow(pclust_info31)) {
  if(pclust_info31$gene_div_special[i] > 0) {
    pclust_info31$pix_ignore[i] <- pclust_info31$pix_count[i]
  }
}

# code or function to calculate 'gen_diverse_weight'..
for (i in 1:nrow(pclust_info31)) {
  if(pclust_info31$gene_div_special[i] != 0) {
    pclust_info31$gene_div_weight[i] <- pclust_info31$gene_div_special[i]
  }
  else {
    pclust_info31$gene_div_weight[i] <-
      pclust_info31$pix_count[i] / (sum(pclust_info31$pix_count) - 
        sum(pclust_info31$pix_ignore))*(1-sum(pclust_info31$gene_div_special))
  }
}

# add count of observations per precluster
pclust_recs31 <- dplyr::filter(pclust_counts31, precluster !=0) |> 
  dplyr::count(precluster) |> sf::st_drop_geometry() |> 
  dplyr::rename(num_obs = n)
pclust_info31 <- left_join(pclust_info31, pclust_recs31, by = "precluster")
# add another column for which is the nearest polygon
nearest <- sf::st_nearest_feature(pclust_info31)
pclust_info31 <- cbind(pclust_info31, nearest)

## update proximity column with nearest distance
for (i in 1:nrow(pclust_info31)) {
  pclust_info31$proximity[i] <- prox31[,i] |> 
    dplyr::slice(nearest[i])
}

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
  if(pclust_info31$Ne_override[i] != 0) {
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
## NEED DIFFERENT CODE FOR 'isolation by resistance' TAXA
## A (linear?) FUNCTION OF LANDSCAPE RESISTANCE COMPUTED BY CIRCUITSCAPE
## e.g. below LSCAPE_RES_MIN; not applicable (connect relevant preclusters)
##  above LSCAPE_RES_MAX; = 0 (negligible effective gene flow)
### To do: MOVE 2.6 and 2.2 INTO config.toml file
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
## WHAT IF NEAREST PRECLUSTER IS LARGE (OR SMALL)?
### To do: MOVE 1000 INTO config.toml file
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

# code to derive 'cluster_score'
for (i in 1:nrow(pclust_info31)) {
  pclust_info31$cluster_score[i] <-
    (pclust_info31$div_risk[i] + pclust_info31$inbreed_risk[i]) *
    pclust_info31$gene_div_weight[i]
}
# fragmented risk for a given taxon will be sum(pclust_info31$cluster_score)
## DATA BASE WILL NEED A NEW COLUMN FOR THIS,
##  AND THEN ANOTHER NEW COLUMN FOR THE FINAL RISK SCORE AFTER ADJUSTING 
##  FOR VARIOUS OTHER FACTORS

