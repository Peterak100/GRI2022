setwd("/home/peter/GeneticRiskIndex/R")
galah_config(email = "peterkriesner@yahoo.com", download_reason_id = 0,
             verbose = TRUE)

# Retrieve state counts from ALA
### use galah_identify() since v1.4 instead of select_taxa()
### use atlas_counts() instead of ala_counts()
# taxa or taxa$ala_search_term

## retrieve state counts from ALA
## this yields count for all items in 'taxa' (25212 for test set of 30 reptiles)
## Need to group them...
# read in test set (if not already loaded)
taxa <- read.csv(BATCH_TAXA_CSV_PATH, header = TRUE)
# for an individual species / taxon:
taxon <- taxa[7,]
taxon$ala_search_term

# test function for raw count for Victoria (from categorize.R)
test1 <- get_state_counts(taxon)

# test function for all counts raw from ALA (from categorize.R)
test2 <- get_all_counts(taxon)


### original version of function
# get_state_counts <- function(taxa) {
#  ala_counts(
#    taxa = select_taxa(taxa$ala_search_term), 
#    filters = ALA_FILTERS,
#    group_by = "species",
#    type = "record",
#    limit = NULL
#  )
# }

# "Varanus rosenbergi"
# "Tiliqua"
# cl1048 is IBRA7 Region
# cl10930 is NRM Region 2017
# for all records for a given taxon:
test3 <- galah_call() |> galah_identify(taxon$ala_search_term) |>
  galah_filter(year >= 1960, year <= 2021, basisOfRecord == BASIS3)|>
  galah_select(basisOfRecord, species, cl1048, cl10930, cl22, year,
               group = "basic") |>
  atlas_occurrences()
write_csv(test3, "/media/sf_Ubuntu_102/temp1/test3.csv")

check1a <- galah_call() |>
  galah_identify("Tiliqua") |>
  galah_filter(year >= 1960, year <= 2021, basisOfRecord == BASIS2,
               stateProvince == STATE) |>
  galah_group_by("species") |>
  atlas_counts(type = "record", limit = NULL)

# for actual records, need to include a call to galah_select()
# for Victoria records only:
test4 <- galah_call() |>
  galah_identify(taxon$ala_search_term) |>
  galah_filter(year >= 1960, year <= 2021, basisOfRecord == BASIS2,
               stateProvince == STATE) |>
  galah_select(basisOfRecord, species, cl1048, cl10930, cl22, year,
               group = "basic") |>
  atlas_occurrences()


### this should work:
check3 <- galah_call() |>
  galah_identify("Varanus rosenbergi") |>
  galah_filter(stateProvince == STATE) |>
  galah_select(basisOfRecord, species, group = "basic") |>
  atlas_occurrences()
write_csv(check3, "/media/sf_Ubuntu_102/temp1/check3.csv")

# 89 possibilities for IBRA7 regions
check4 <- search_field_values("cl1048", limit = 100)
# not clear why this only returns 48 fields:
check2 <- seach_fields("all")
# but this returns all 791 fields:
check2a <- show_all_fields()

## DIFFERENCE BETWEEN precategorized_taxa and preclustered_taxa??

write_csv(precategorized_taxa, "/media/sf_Ubuntu_102/temp1/test7.csv")


## a test plot...
library(yarrr)
# set colour scheme
YarrCols <- c("coral2","steelblue3")

png(file.path(plotpath, paste0(gsub(" ", "_", taxon$ala_search_term),
            "_test59", ".png")), width = 2160, height = 1440, pointsize = 30)
# sets dev = png 2
# Create plot (default settings are too low res)
plot(x = pirates$height, y = pirates$weight,
     xlim = c(120, 225), ylim = c(30, 110), pch = 16, cex = 3,
     col = transparent(YarrCols[factor(pirates$sex)], trans.val = 0.8),
     xlab = "Height", ylab = "Weight",
     main = paste0(taxon$ala_search_term, " Orphans"), sub = NULL)
# cex.lab = 2? cex.axis = 2? cex.main = 2?

cat("plot for", taxon$ala_search_term, "saved to :", plotpath,"\n")
while (!is.null(dev.list())) dev.off()


zz2 <- matrix(letters[1:32], nrow = 4, ncol = 8, byrow = T) |> as.data.frame()
# writes csv file in plotpath directory...
write_csv(zz2, file.path(plotpath, "bing.csv")) # must be a data frame





