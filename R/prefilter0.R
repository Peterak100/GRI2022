# reworking of prefilter.R (main script to run)
# setwd("/home/peter/junk/Galah1")
setwd("/home/peter/GRI2022/R")
source("common1.R")
# source("verify2.R")
source("categorize3.R")
source("obs4.R")
source("resistance5.R") # working on this 7th Oct 2022


# Precategorize taxa based on record counts --------

# Precategorize using functions from g_categorize2.R
# Categorize risk using queries to ALA: slow.
## adds 4 columns to BATCH_TAXA_CSV:
## "state_count", "count", "risk", "filter_category"
start_time1 <- Sys.time()
precategorized_taxa <- precategorize_risk(taxa)
end_time1 <- Sys.time()
end_time1 - start_time1 # 59 secs for first 30 taxa records
# display selected columns for first 14 taxa..
# head(precategorized_taxa[,c(1,2,39, 94:97)],14)

write_csv(precategorized_taxa,
          file.path(groupingspath, "precategorized_taxa.csv"))

# Run processing of observations --------

## Filter for just those taxa that obs records are needed for
isolation_taxa <- dplyr::filter(precategorized_taxa,
        filter_category == "isolation by distance" | 
        filter_category == "isolation by resistance")

# adds another 5 columns to isolation_taxa
# omitting 'throw_errors=TRUE' as 4th argument to function
start_time1 <- Sys.time()
preclustered_isolation_taxa <- process_observations(isolation_taxa,
        mask_layer, taxapath)
end_time1 <- Sys.time()
end_time1 - start_time1 # ~4.41 mins for first 30 taxon records

## warnings that dir.create() already exists should now be suppressed


## WHAT DOES / WOULD 'throw_errors=TRUE' DO IN THE ABOVE FUNCTION?

## To examine first (10) rows of new dataframe..
# head(preclustered_isolation_taxa[,c(1,2,39,43,94:102)], 10)


# Combine updated isolation taxa with non-processed taxa
other_taxa <- dplyr::filter(precategorized_taxa,
        filter_category != "isolation by distance" &
        filter_category != "isolation by resistance")
other_taxa <- add_column(other_taxa, num_preclusters = 0, num_orphans = 0, 
        precluster_cellcount = 0, orphan_cellcount = 0, error = NA) 
categorized_taxa <- rbind(preclustered_isolation_taxa, other_taxa) # sort?

write_csv(categorized_taxa, file.path(groupingspath, "categorized_taxa.csv"))


# Circuitscape/isolation by resistance output --------

# Save separate lists of isolation by distance & resistance taxa to variables
isolation_by_distance_taxa <- dplyr::filter(preclustered_isolation_taxa,
          is.na(risk), filter_category == "isolation by distance")
isolation_by_resistance_taxa <- dplyr::filter(preclustered_isolation_taxa,
          is.na(risk), filter_category == "isolation by resistance")


# Write csv for taxa that we need to process with Circuitscape
write_csv(isolation_by_resistance_taxa, file.path(groupingspath,
          "isolation_by_resistance_taxa.csv"))

# Write list of taxa to process in Circuitscape
#  as a single column job-list text file
job_file <- file(file.path(datapath, "batch_jobs.txt"))
underscored <- gsub(" ", "_", isolation_by_resistance_taxa$ala_search_term)
writeLines(underscored, job_file)
close(job_file)




# Download and write raster files for Circuitscape resistance models
prepare_resistance_files(isolation_by_resistance_taxa, taxapath)

## Run from Julia:
# using Circuitscape
# compute("/home/peter/data/taxa/Varanus_varius/Circuitscape_custom1.ini")
## TO DO: make a small .jl script to run the above
###  for all taxa in "batch_jobs.txt"



