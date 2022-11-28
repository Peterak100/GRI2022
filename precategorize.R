source("R/common1.R")
source("categorize3.R")
source("obs4.R")
source("resistance5.R")
source("postprocess7.R")
precategorized_taxa <- precategorize_risk(taxa)
write_csv(precategorized_taxa, file.path(groupingspath,
                  "precategorized_taxa.csv"))