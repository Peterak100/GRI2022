suppressMessages(source("R/common1.R"))
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
retrieve_bulk_obs(isolation_taxa)