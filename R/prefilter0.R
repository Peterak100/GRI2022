# reworking of prefilter.R (main script to run)
GRI_dir <- getwd()
setwd(paste0(GRI_dir,"/R"))
source("common1.R")
# source("verify2.R")
source("categorize3.R")
source("obs4.R")
source("resistance5.R")
source("postprocess7.R")

# Precategorize taxa based on record counts --------

# Precategorize using functions from g_categorize2.R
# Categorize risk using queries to ALA: slow.
## adds 4 columns to BATCH_TAXA_CSV:
## "state_count", "count", "risk", and "filter_category"
start_time1 <- Sys.time()
precategorized_taxa <- precategorize_risk(taxa)
end_time1 <- Sys.time()
end_time1 - start_time1 # ~1 to 1.3 mins for first 30 taxa records
# display selected columns for first 14 taxa..
# head(precategorized_taxa[,c(1,2,39, 94:97)],14)

write_csv(precategorized_taxa,
          file.path(groupingspath, "precategorized_taxa.csv"))

# Run processing of observations --------

## Filter for just those taxa that obs records are needed for
isolation_taxa <- dplyr::filter(precategorized_taxa,
        filter_category == "isolation by distance" | 
        filter_category == "isolation by resistance")

# Initial bulk download of all observations and write to individual
## .csv files, no further processing
start_time1 <- Sys.time()
retrieve_bulk_obs(isolation_taxa)
end_time1 <- Sys.time()
end_time1 - start_time1

# adds another 5 columns to isolation_taxa
# omitting 'throw_errors=TRUE' as 4th argument to function
start_time1 <- Sys.time()
preclustered_isolation_taxa <- process_observations(isolation_taxa,
        mask_layer, taxapath)
end_time1 <- Sys.time()
end_time1 - start_time1 # ~4.41 mins for first 30 taxon records


## WHAT DOES / WOULD 'throw_errors=TRUE' DO IN THE ABOVE FUNCTION?

## To examine first (10) rows of new dataframe..
# head(preclustered_isolation_taxa[,c(1,2,39,43,94:102)], 10)


## Circuitscape/isolation by resistance output --------

# Save separate lists of isolation by distance & resistance taxa to variables
isolation_by_distance_taxa <- dplyr::filter(preclustered_isolation_taxa,
        is.na(risk), filter_category == "isolation by distance")
isolation_by_resistance_taxa <- dplyr::filter(preclustered_isolation_taxa,
        is.na(risk), filter_category == "isolation by resistance")

# Write list of taxa to process in Circuitscape
#  as a single column job-list text file
job_file <- file(file.path(datapath, "batch_jobs.txt"))
underscored <- gsub(" ", "_", isolation_by_resistance_taxa$ala_search_term)
writeLines(underscored, job_file)
close(job_file)

# Write raster files for Circuitscape resistance models
prepare_resistance_files(isolation_by_resistance_taxa, taxapath)

## First Circuitscape run from Julia using 'batch_jobs.txt'
# see circuitscape6.R for instructions

## Post processing --------
score_distance_taxa(isolation_by_distance_taxa, taxapath)
score_resistance_taxa1(isolation_by_resistance_taxa, taxapath)

## Second Circuitscape run from Julia using 'batch_jobs.txt'
# see circuitscape6.R for instructions

score_resistance_taxa2(isolation_by_resistance_taxa, taxapath)

# Combine updated isolation taxa with non-processed taxa
other_taxa <- dplyr::filter(precategorized_taxa,
        filter_category != "isolation by distance" &
        filter_category != "isolation by resistance")
## more columns to add??
other_taxa <- add_column(other_taxa, num_preclusters = 0, num_orphans = 0, 
        precluster_cellcount = 0, orphan_cellcount = 0, error = NA,
        frag_risk = 0) 
categorized_taxa <- rbind(preclustered_isolation_taxa, other_taxa) # sort?

write_csv(categorized_taxa, file.path(groupingspath, "categorized_taxa.csv"))

