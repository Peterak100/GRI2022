suppressPackageStartupMessages(source("R/common1.R"))
source("categorize3.R")
source("obs4.R")
source("resistance5.R")
source("postprocess7.R")
precategorized_taxa <- read_csv(file.path(groupingspath,
          "precategorized_taxa.csv"), show_col_types = FALSE)
isolation_taxa <- dplyr::filter(precategorized_taxa,
          filter_category == 
            "isolation by distance" | filter_category == 
            "isolation by resistance")
preclustered_isolation_taxa <- process_observations(isolation_taxa,
          mask_layer, taxapath)
isolation_by_distance_taxa <- dplyr::filter(preclustered_isolation_taxa,
          is.na(risk), filter_category == "isolation by distance")
isolation_by_resistance_taxa <- dplyr::filter(preclustered_isolation_taxa,
          is.na(risk), filter_category == "isolation by resistance")
job_file <- file(file.path(datapath, "batch_jobs.txt"))
underscored <- gsub(" ", "_", isolation_by_resistance_taxa$ala_search_term)
writeLines(underscored, job_file)
cat(paste0("batch_jobs",".txt"," for Circuitscape written to ",
           file.path(datapath)))
close(job_file)
prepare_resistance_files(isolation_by_resistance_taxa, taxapath)
other_taxa <- dplyr::filter(precategorized_taxa,
          filter_category != "isolation by distance" &
          filter_category != "isolation by resistance")
other_taxa <- add_column(other_taxa, num_preclusters = 0, num_orphans = 0, 
          precluster_cellcount = 0, orphan_cellcount = 0, error = NA) 
categorized_taxa <- rbind(preclustered_isolation_taxa, other_taxa)
write_csv(categorized_taxa, file.path(groupingspath, "categorized_taxa.csv"))
