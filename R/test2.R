# Single taxon observations & preclustering for testing --------
# test80 <- isolation_taxa[1, ]
# substitute whatever taxon number is to be tested for "1" in above
# test80$ala_search_term

# test81 <- process_obs_single(test80, mask_layer, taxapath)
# head(test81[,c(1:2,99:107)])
process_obs_single <- function(taxon, mask_layer, taxapath) {
  # Add new columns to taxa dataframe
  preclustered_taxa <- add_column(taxon, 
          num_preclusters = 0, 
          num_orphans = 0, 
          precluster_cellcount = 0, 
          orphan_cellcount = 0,
          error = NA
  ) 
  num_preclusters <- 0 # ??
  cat("\nTaxon: ", taxon$ala_search_term, "\n")
  obs <- load_and_filter(taxon, taxapath)
  if (max(obs$precluster) != 0) {
    cell_counts <- write_precluster(obs, taxon, mask_layer, taxapath)
  } else {
    cell_counts <- c(0, 0)
  }
  preclustered_taxa$num_preclusters <- max(obs$precluster)
  preclustered_taxa$num_orphans <- sum(obs$precluster == 0)
  preclustered_taxa$precluster_cellcount <- cell_counts[1]
  preclustered_taxa$orphan_cellcount <- cell_counts[2] |> 
    label_by_clusters()
  return(preclustered_taxa)
}

# head(test81[,c(1:2,99:107)])

### to view file permissions in Bash terminal: ls -l
### possible file permissions are: -rwxrwxrwx
###  where r = read, w = write, x = eXecute
###  1st group of 3 = owner/user (u), 2nd group = group owner (g),
###  3rd group = other users (o)
### to change file permissions: chmod u+rw filename


# subfunction to obtain species counts (SHOULD BE 1 FOR ALL OF THEM)
## and alternative (updated) scientific_name / species name

check_unique_sp <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= as.character(TIME_START)) |>
    galah_group_by("species") |>
    atlas_counts(type = "species", limit = NULL)
}

add_sp_column <- function(taxa) {
  count_taxa <- filter(taxa, assess == "ALA") 
  all_species <- check_unique_sp(count_taxa) |> rename(ala_search_term = species)
  taxa <- left_join(taxa, all_species, by = "ala_search_term", all.x=TRUE)
  return(taxa)
}

test1112 <- add_sp_column(taxa)
test1113 <- dplyr::filter(test1112, is.na(count) | count > 1)
test1113[,c(1,2,99)]
test1114 <- search_taxa(test1113$ala_search_term)
test1115 <- test1114[,c(1,2,11:14)]


## NENAD'S IDEA
test1 <- taxa
results<-list()
recheck<-list()
for(i in nrow(test1)){
  try(taxonomyDF<-ala_species(taxa=test1))
  searchTerm<-test1$ala_search_term[i]
  if(nrow(taxonomyDF)>0){
    results[[test1]]<-taxonomyDF$genus
  } else {
    recheck[[test1]]<-test1
  }
}

