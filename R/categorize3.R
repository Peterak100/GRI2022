# Functions to modify taxa --------

# subfunction to obtain record counts for Victoria only
get_state_counts <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= as.character(TIME_START),
                 year <= as.character(TIME_END),
                 basisOfRecord == BASIS3,
                 stateProvince == STATE, profile = "ALA") |>
    atlas_counts(type = "record", limit = NULL)
}

# subfunction to obtain record counts for all of Australia
## DOES NOT WORK TO USE galah_group_by("species") as some search terms
## are subspecies, etc.
get_all_counts <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= as.character(TIME_START),
                 year <= as.character(TIME_END),
                 basisOfRecord == BASIS3, profile = "ALA") |>
    atlas_counts(type = "record", limit = NULL)
}

## function to append two record count columns to end of taxa
## this is quite slow, but returns results for species & subspecies
add_count_cols <- function(taxa){
  cat("\nGetting taxon counts for initial filtering...\n")
  count_taxa <- data.frame(matrix(nrow = nrow(taxa), ncol = 3))
  ALA1 <- colnames(taxa)[1]
  colnames(count_taxa) <- c(ALA1, "state_rec_count", "all_rec_count")
  for (i in 1:nrow(taxa)) {
    state_rec_count <- get_state_counts(taxa[i,])
    all_rec_count <- get_all_counts(taxa[i,])
    count_taxa[i, 1] <- taxa$ala_search_term[i]
    count_taxa[i, 2] <- state_rec_count
    count_taxa[i, 3] <- all_rec_count
    dplyr::bind_rows(count_taxa, .id = ALA1)
  }
  joined <- left_join(taxa, count_taxa, by = ALA1, all.x = TRUE) # new
  cat("Individual taxon counts retrieved. \n\n")
  return(joined)
}

## to test the above functions:
# test_A <- add_count_cols(taxa)
# test_A[1:9,c(1,2,97:100)]

## Add a risk column that classifies risk. Initially "unknown" so NA.
add_risk_col <- function(taxa) {
  taxa$risk <- rep(NA, length(taxa$ala_search_term))
  return(taxa)
}
# test_B <- add_risk_col(test_A)

## Label non ALA taxa (adds 'filter_category' column)
label_not_assessed <- function(taxa) {
  ids <- taxa$assess != "ALA"
  taxa$risk[ids] <- "not_assessed"
  taxa$filter_category[ids] <- "not_ALA_taxon"
  return(taxa)
}
# test_C <- label_not_assessed(test_B)

## Update filter_category column based on disperse_model column (#40)
label_isolation_by_resistance <- function(taxa) {
  ids <- is.na(taxa$filter_category) & taxa$disperse_model == "Habitat"
  taxa$filter_category[ids] <- "isolation by resistance"
  return(taxa)
}
# test_D <- label_isolation_by_resistance(test_C)
# test_D[1:25,c(1,2,40,97:102)]

label_isolation_by_distance <- function(taxa) {
  ids <- is.na(taxa$filter_category) & taxa$disperse_model == "Distance"
  taxa$filter_category[ids] <- "isolation by distance"
  return(taxa)
}
# test_E <- label_isolation_by_distance(test_D)
# test_E[1:25,c(1,2,40,97:102)]

## Label species with little relevance to Victoria
###  current threshold is < 5% using MIN_PROP_IN_STATE from config.toml

label_low_regional_relevance <- function(taxa) {
  ids <- (taxa$state_rec_count / taxa$all_rec_count) < MIN_PROP_IN_STATE
  taxa$risk[ids] <- "widespread"
  taxa$filter_category[ids] <- "low_proportion_in_state"
  return(taxa)
}
# test_F <- label_low_regional_relevance(test_E)
# test_F[,c(1,2,40,97:102)]

## Label very common species as "abundant" in risk column
## Add a filter_category column
## current thresholds: > 100K records for Vic, or > 300K records for ALA total
label_many_observations <- function(taxa) {
  ids <- taxa$state_rec_count > MAX_OBSERVATIONS | 
    taxa$all_rec_count > MAX_OBS_TOTAL
  taxa$risk[ids] <- "abundant"
  taxa$filter_category[ids] <- "many_observations"
  return(taxa)
}
# test_G <- label_many_observations(test_F)
## test_G[,c(1,2,97:102)]

## Label very rare species as "rare", current threshold is < 50
label_few_observations <- function(taxa) {
  ids <- taxa$all_rec_count < MIN_OBSERVATIONS
  taxa$risk[ids] <- "rare"
  taxa$filter_category[ids] <- "few_observations"
  return(taxa)
}
# test_H <- label_few_observations(test_G)
## test_H[,c(1,2,97:102)]

# need another filter for species that are 'data deficient'?
# combining all above subfunctions..
precategorize_chunk <- function(taxa) {
  taxa |> 
    add_count_cols() |> 
    add_risk_col() |> 
    label_not_assessed() |> 
    label_isolation_by_resistance() |> 
    label_isolation_by_distance() |> 
    label_low_regional_relevance() |> 
    label_many_observations() |> 
    label_few_observations() %>% # errors with built-in pipe here?!
    identity
}
# test_J <- precategorize_chunk(taxa)
## test_J[,c(1,2,97:102)] ## should be same result as test_H

# Main function --------

precategorize_risk <- function(taxa) {
  start_time1 <- Sys.time()
  n <- nrow(taxa)
  r <- rep(1:ceiling(n/GALAH_MAXROWS), each=GALAH_MAXROWS)[1:n]
  s <- lapply(split(taxa, r), precategorize_chunk)
  # Split apply combine chunks
  do.call(rbind, s)
  end_time1 <- Sys.time()
  precategorize_time <- end_time1 - start_time1
  cat(paste0("Time to precategorize taxa was ", precategorize_time,"\n"))
}
# test_K <- precategorize_risk(taxa)
## test_K[,c(1,2,38,98:102)] ## should also be same result as test_H
## rm(test_A,test_B,test_C,test_D,test_E,test_F,test_G,test_H,test_J,test_K)
