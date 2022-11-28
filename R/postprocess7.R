## see test5.R for step-by-step testing

## Score 'isolation_by_distance' taxa --------
# To be run on 'isolation_by_distance_taxa'
score_distance_taxa <- function(taxa, taxapath) {
  for (i in 1:nrow(taxa)) {
    taxon <- taxa[i, ]
    # read in previously saved ESRI shapefile
    pre_clusters <- st_read(file.path(taxonpath,
          paste0(gsub(" ","_", taxon$ala_search_term),
          "_preclusters_prelim", ".shp"))) # 7 columns
    # column names in saved '.shp' file are truncated
    # to 7 characters max so rename where appropriate
    pre_clusters <- pre_clusters |> 
      rename(midcluster = mdclstr, precluster = prclstr,
          pix_count = pix_cnt, recent_year = rcnt_yr)
    # add a new first column as placeholder for cluster number
    pre_clusters <- cbind(cluster = pre_clusters$midcluster, pre_clusters)
    # add new columns, to be subsequently filled or updated
    pre_clusters <- pre_clusters |> add_column(pop_name = NA,
          gene_div_special = NA, pix_ignore = 0, Ne_override = NA,
          gene_div_weight = 0, proximity = 0) # 14 columns
    # import expert opinion file for given taxon if one exists
    expert1_filename <- suppressWarnings(file.path(taxon_path(taxon,
          taxapath), paste0("expert_op1_preclust", ".csv")))
    if (file.exists(expert1_filename)){
      expert1 <- read.csv(expert1_filename, header = TRUE)
      if (nrow(expert1) == nrow(pre_clusters31)) {
        pre_clusters31$pop_name <- expert1$pop_name
        pre_clusters31$gene_div_special <- expert1$gene_div_special
        pre_clusters31$Ne_override <- expert1$Ne_override
      }
      else {cat("\nWrong number of populations in expert info file!\n")
      }
    }
    # for any clusters where 'gene_div_special' is not = NA; 
    # copy value of 'pix_count' column to 'pix_ignore' column
    for (ja in 1:nrow(pre_clusters)) {
      if(!is.na(pre_clusters$gene_div_special[ja])) {
        pre_clusters$pix_ignore[ja] <- pre_clusters$pix_count[ja]
      }
    }
    # calculate 'gen_diverse_weight'..
    for (jb in 1:nrow(pre_clusters)) {
      if(!is.na(pre_clusters$gene_div_special[jb])) {
        pre_clusters$gene_div_weight[jb] <- pre_clusters$gene_div_special[jb]
      }
      else {
        pre_clusters$gene_div_weight[jb] <-
          pre_clusters$pix_count[jb] / (sum(pre_clusters$pix_count) - 
          sum(pre_clusters$pix_ignore))*(1-sum(pre_clusters$gene_div_special,
            na.rm = TRUE))
      }
    }
    # compute which is the nearest polygon
    nearest <- sf::st_nearest_feature(pre_clusters)
    pre_clusters <- cbind(pre_clusters, nearest) # 15 columns
    # create a 'units' matrix of distances between polygons
    # then convert to a normal numeric matrix & convert to kilometres
    prox <- sf::st_distance(pre_clusters) |> as.data.frame() |> as.matrix()
    prox <- prox/1000
    ## update proximity column with nearest distance
    # slice needs to operate on a data frame
    for (jc in 1:nrow(pre_clusters)) {
      pre_clusters$proximity[jc] <- as.data.frame(prox[,jc]) |> 
        dplyr::slice(nearest[jc])
    }
    # The above makes $proximity a list. Need to convert back to numeric
    pre_clusters$proximity <- as.numeric(pre_clusters$proximity)
    # add 3 columns for weighted size of precluster, epsilon proximity, and
    ## effective (Victorian) population size (rounded to whole numbers)
    pre_clusters <- pre_clusters |> 
      mutate(pix_pc = (pix_count / sum(pix_count))) |>
      mutate(eps_prox = proximity / taxon$epsilon) |> 
      mutate(Ne_est = pix_pc * taxon$vic_popn_size * 0.1)
    pre_clusters$Ne_est <- round(pre_clusters$Ne_est, digits = 0)
    # 18 columns
    # update with Ne values manually entered for particular populations
    for (jd in 1:nrow(pre_clusters)) {
      if(!is.na(pre_clusters$Ne_override[jd])) {
        pre_clusters$Ne_est[jd] <- pre_clusters$Ne_override[jd]
      }
    }
    # add columns for risk score calculations
    pre_clusters <- pre_clusters |> add_column(size_nearest = 0,
        gene_flow_factor = 0, div_risk = 0, inbreed_risk = 0,
        cluster_score = 0) # should now be 23 columns
    # derive size of nearest precluster..
    for (je in 1:nrow(pre_clusters)){
      pre_clusters$size_nearest[je] <- dplyr::slice(pre_clusters,
          pre_clusters$nearest[je]) |> 
        dplyr::select(pix_count) |> sf::st_drop_geometry()
    }
    # derive 'gene_flow_factor'
    ### TO DO: MOVE 2.6 and 2.2 values INTO config.toml file?
    for (jf in 1:nrow(pre_clusters)){
      if(pre_clusters$eps_prox[jf] > 2.6){
        pre_clusters$gene_flow_factor[jf] <- 0
      }
      else {
        pre_clusters$gene_flow_factor[jf] <- 
          1/((pre_clusters$eps_prox[jf]+1)^2.2)
      }
    }
    # code to derive diversity risk based on 'Ne_est' & 'gene_flow_factor'
    ## TO DO: WHAT IF NEAREST PRECLUSTER IS LARGE (OR SMALL)?
    ### TO DO: MOVE 1000 value INTO config.toml file
    for (jg in 1:nrow(pre_clusters)) {
      if(pre_clusters$Ne_est[jg] < 1000) {
        (pre_clusters$div_risk[jg] <-
           (1 - (pre_clusters$Ne_est[jg]/1000)) * 100) *
          (1 - pre_clusters$gene_flow_factor[jg])
      }
      else {
        pre_clusters$div_risk[jg] <- 0
      }
    }
    # derive 'adj_inbreed_risk'
    ### To do: MOVE 110 value INTO config.toml file
    for (jh in 1:nrow(pre_clusters)) {
      if(pre_clusters$Ne_est[jh] < 110) {
        if(pre_clusters$Ne_est[jh] < 10) {
          pre_clusters$inbreed_risk[jh] <-
            (2000 * (1 - pre_clusters$gene_flow_factor[jh]))
        }
        else {
          pre_clusters$inbreed_risk[jh] <-
            ((110-pre_clusters$Ne_est[jh])^2)/5 *
            (1-pre_clusters$gene_flow_factor[jh])
        }
      }
      else {
        pre_clusters$inbreed_risk[jh] <- 0
      }
    }
    # code to derive 'precluster_score'
    for (ji in 1:nrow(pre_clusters)) {
      pre_clusters$cluster_score[ji] <-
        (pre_clusters$div_risk[ji] + pre_clusters$inbreed_risk[ji]) *
        pre_clusters$gene_div_weight[ji]
    }
    # overall risk due to fragmentation for the given taxon
    taxon <- taxon |> add_column(frag_risk = sum(pre_clusters$cluster_score))
    ## TO DO: Where does frag_risk column get returned???
    ## TO DO: Add another new column for the final risk score after 
    ##  adjusting for various other factors
    
    # write preclusters with derived risk scores to a file
    pre_clusters |> sf::st_drop_geometry() |> 
      write_csv(file.path(taxonpath, paste0(gsub(" ","_",
        taxon$ala_search_term), "_clusters_info", ".csv")))
  }
}


## TO DO: ADJUST 'precluster_score' BASED ON 'recent_year'?
## WHAT PROPORTION OF OBS / PIXELS ARE BASED ON QUITE OLD RECORDS?
## What is the likelihood that for preclusters where all records are old,
##  the taxon is now locally extinct??
## TO DO: Also derive an uncertainty score per precluster?
##  or just for the taxon overall?

## Import first Circuitscape output matrix --------

# derive clusters by iteratively grouping preclusters based on
# relative resistances calculated with Circuitscape
# To be run on 'isolation_by_resistance_taxa'
score_resistance_taxa1 <- function(taxa, taxapath) {
  for (i in 1:nrow(taxa)) {
    taxon <- taxa[i, ]
    # read in previously saved ESRI shapefile
    pre_clusters <- st_read(file.path(taxonpath,
          paste0(gsub(" ","_", taxon$ala_search_term),
          "_preclusters_prelim", ".shp"))) # 7 columns
    # column names in saved '.shp' file are truncated to 7 characters max
    # add a new first column as placeholder for cluster number
    pre_clusters <- cbind(cluster = pre_clusters$mdclstr, pre_clusters)
    # 8 columns
    # derive clusters by iteratively grouping preclusters based on
    # relative resistances calculated with Circuitscape
    CScape_FILE1 <- paste0("CSpwise1_resistances", ".out")
    CScape_PATH1 <- file.path(taxonpath, CScape_FILE1)
    if (file.exists(CScape_PATH1)){
      pwise1 <- read.table(CScape_PATH1, header = TRUE, sep = " ",
            row.names = 1, check.names = FALSE)
      # rownames are "1", "2", etc
      # but colnames are "1.0", "2.0", "3.0" etc
      colnames(pwise1) <- rownames(pwise1)
      # use a different resistance threshold for generic resistance layer?
      SMP_MODEL <- paste0("SMP_", gsub(" ","_",
            taxon$ala_search_term), ".tif")
      SMP_filename <- suppressWarnings(file.path(SMPpath, SMP_MODEL))
      HIM_MODEL <- paste0("HIM_", gsub(" ","_",
            taxon$ala_search_term), ".tif")
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
        pre_clusters$cluster[j] <- which(lapply(fcluster1,
            function(x) grep(paste0("^",j,"$"),x))!=0) |> unname()
      }
    } else {
      cat("Circuitscape results NOT FOUND for",
          taxon$ala_search_term, "\n")
    }
    # create summary table based on cluster number
    clusters_info <- pre_clusters |> dplyr::group_by(cluster) |> 
      dplyr::summarise(mdclstr = paste(mdclstr, collapse = ", "),
                prclstr = paste(prclstr, collapse = ", "),
                pix_cnt = sum(pix_cnt), rcnt_yr = max(rct_yr),
                latin = max(latin), num_obs = sum(num_obs),
                geometry = sf::st_union(geometry)) # also 8 columns
    # also make a new patches .tif file for newly grouped clusters
    shapevect <- terra::vect(clusters_info)
    obs_retest_rast <- terra::rasterize(shapevect, mask_layer,
                field = "cluster")
    obs_retest_filename <- file.path(taxonpath, paste0("clusters", ".tif"))
    crop_final_rast <- obs_retest_rast |> padded_trim()
    terra::crop(obs_retest_rast, crop_final_rast,
                filename = obs_retest_filename, overwrite = TRUE)
    # overwrite existing sf object file to retrieve for next step
    clusters_info |>  suppressWarnings(sf::write_sf(file.path(taxonpath,
                paste0(gsub(" ","_", taxon$ala_search_term),
                "_preclusters_prelim", ".shp")), append = FALSE))
  }
}

# after running Circuitscape a 2nd time
## Import second Circuitscape output matrix --------
score_resistance_taxa2 <- function(taxa, taxapath) {
  for (i in 1:nrow(taxa)) {
    taxon <- taxa[i, ]
    # read in previously saved ESRI shapefile a 2nd time
    final_clusters <- st_read(file.path(taxonpath,
            paste0(gsub(" ","_", taxon$ala_search_term),
            "_preclusters_prelim", ".shp"))) # 7 columns
    # column names in saved '.shp' file are truncated
    # to 7 characters max so rename where appropriate
    final_clusters <- final_clusters |> 
      rename(midcluster = mdclstr, precluster = prclstr,
             pix_count = pix_cnt, recent_year = rcnt_yr)
    pwise2 <- read.table(file.path(taxonpath,
        paste0("CSpwise2_resistances_3columns", ".out")), header = TRUE,
        sep = " ", col.names = c("cluster", "closest", "resistance"),
        check.names = FALSE)
    # transform to find which is the next nearest cluster by resistance
    flip <- as.data.frame(cbind(pwise2[,2], pwise2[,1], pwise2[,3]))
    names(flip)<-names(pwise2)
    # this returns only the 1st minima (in case there is more than one)
    pwise2 <- rbind(pwise2, flip) |> dplyr::group_by(cluster) |> 
      dplyr::slice(which.min(resistance))
    # add new columns, including for next nearest cluster by resistance
    final_clusters <- final_clusters |> add_column(pop_name = NA,
        gene_div_special = NA, pix_ignore = 0, Ne_override = NA,
        gene_div_weight = 0, lowest_resist = NA, closest = NA) # 15 columns
    # update lowest_resist column with next lowest resistance
    for (i in 1:nrow(final_clusters)) {
      final_clusters$lowest_resist[i] <- pwise2$resistance[i]
    }
    ## If and when available, record expert opinion or info derived from
    # genetic data on effective size of, or proportion of genetic diversity
    # contained within, a given local population, along with relevant local
    # population (precluster) names in a supplementary file
    
    # if available, import expert opinion file for clusters for given taxon
    expert3_filename <- suppressWarnings(file.path(taxon_path(taxon,
              taxapath), paste0("expert_op3_clusters", ".csv")))
    if (file.exists(expert3_filename)){
      expert3 <- read.csv(expert3_filename, header = TRUE)
      if (nrow(expert3) == nrow(final_clusters)) {
        final_clusters$pop_name <- expert3$pop_name
        final_clusters$gene_div_special <- expert3$gene_div_special
        final_clusters$Ne_override <- expert3$Ne_override
      }
      else {cat("\nWrong number of populations in expert info file!\n")
      }
    }
    # for any clusters where 'gene_div_special' is not NA; 
    # copy value of 'pix_count' column to 'pix_ignore' column
    for (i in 1:nrow(final_clusters)) {
      if(!is.na(final_clusters$gene_div_special[i])) {
        final_clusters$pix_ignore[i] <- final_clusters$pix_count[i]
      }
    }
    # calculate 'gen_diverse_weight'..
    for (i in 1:nrow(final_clusters)) {
      if(!is.na(final_clusters$gene_div_special[i])) {
        final_clusters$gene_div_weight[i] <- 
          final_clusters$gene_div_special[i]
      }
      else {
        final_clusters$gene_div_weight[i] <-
          final_clusters$pix_count[i] / (sum(final_clusters$pix_count) - 
            sum(final_clusters$pix_ignore)) *
          (1-sum(final_clusters$gene_div_special, na.rm = TRUE))
      }
    }
    # update 'closest' column with next closest cluster by resistance
    for (i in 1:nrow(final_clusters)) {
      final_clusters$closest[i] <- pwise2$closest[i]
    }
    # add columns for weighted size of clusters and effective (Victorian)
    # population size (rounded to whole numbers)
    final_clusters <- final_clusters |> 
      mutate(pix_pc = (pix_count / sum(pix_count))) |>
      mutate(Ne_est = pix_pc * taxon$vic_popn_size * 0.1)
    final_clusters$Ne_est <- round(final_clusters$Ne_est, digits = 0)
    # 15 columns
    # update with Ne values manually entered for particular populations
    for (i in 1:nrow(final_clusters)) {
      if(!is.na(final_clusters$Ne_override[i])) {
        final_clusters$Ne_est[i] <- final_clusters$Ne_override[i]
      }
    }
    # add risk score calculation columns
    final_clusters <- final_clusters |> add_column(size_nearest = 0,
        gene_flow_factor = 0, div_risk = 0, inbreed_risk = 0,
        cluster_score = 0) # 22 columns
    # derive size of nearest cluster
    for (i in 1:nrow(final_clusters)){
      final_clusters$size_nearest[i] <- dplyr::slice(final_clusters,
          final_clusters$closest[i]) |> 
        dplyr::select(pix_count) |> sf::st_drop_geometry()
    }
    # derive 'gene_flow_factor'
    ## TO DO: move 2.2 INTO config.toml file? Is this method appropriate?
    res_factor2 <- (taxon$disperse_log * 3) + 6.3
    for (i in 1:nrow(final_clusters)){
      if(final_clusters$lowest_resist[i] > (res_factor2 * 1.75)){
        final_clusters$gene_flow_factor[i] <- 0
      } else if (final_clusters$lowest_resist[i] < res_factor2) {
        final_clusters$gene_flow_factor[i] <- 1
      } else {
        final_clusters$gene_flow_factor <-
          (((res_factor2 * 2) - final_clusters$lowest_resist) *
             (1/res_factor2))^2.2
      }
    }
    # code to derive diversity risk based on 'Ne_est' & 'gene_flow_factor'
    ## TO DO: Move 1000 to config.toml file?
    ## TO DO: What if nearest cluster is small (or large)?
    for (i in 1:nrow(final_clusters)) {
      if(final_clusters$Ne_est[i] < 1000) {
        (final_clusters$div_risk[i] <-
           (1 - (final_clusters$Ne_est[i]/1000)) * 100) *
          (1 - final_clusters$gene_flow_factor[i])
      }
      else {
        final_clusters$div_risk[i] <- 0
      }
    }
    # derive 'adj_inbreed_risk'
    ### TO DO: Move 110 value to config.toml file?
    for (i in 1:nrow(final_clusters)) {
      if(final_clusters$Ne_est[i] < 110) {
        if(final_clusters$Ne_est[i] < 10) {
          final_clusters$inbreed_risk[i] <-
            (2000 * (1 - final_clusters$gene_flow_factor[i]))
        }
        else {
          final_clusters$inbreed_risk[i] <-
            ((110 - final_clusters$Ne_est[i])^2)/5 *
            (1 - final_clusters$gene_flow_factor[i])
        }
      }
      else {
        final_clusters$inbreed_risk[i] <- 0
      }
    }
    # code to derive 'cluster_score'
    for (i in 1:nrow(final_clusters)) {
      final_clusters$cluster_score[i] <-
        (final_clusters$div_risk[i] + final_clusters$inbreed_risk[i]) *
        final_clusters$gene_div_weight[i]
    }
    # overall risk due to fragmentation for the given taxon
    taxon <- taxon |> 
      add_column(frag_risk = sum(final_clusters$cluster_score))
    ## TO DO: Add another new column for the final risk score after 
    ##  adjusting for various other factors
    # write clusters with derived risk scores to a file
    final_clusters |> sf::st_drop_geometry() |> 
      write_csv(file.path(taxonpath, paste0(gsub(" ","_",
        taxon$ala_search_term), "_clusters_info", ".csv")))
  }
}
