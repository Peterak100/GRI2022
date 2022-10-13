# verify that database records correspond uniquely to a single species
## i.e. performs count of 'species', not of 'records'
check_unique_sp <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= as.character(TIME_START)) |>
    galah_group_by("species") |>
    atlas_counts(type = "species", limit = NULL)
}

# subfunction to add an natra column to taxa
add_sp_column <- function(taxa) {
  count_taxa <- filter(taxa, assess == "ALA") 
  all_species <- check_unique_sp(count_taxa) |>
    rename(ala_search_term = species)
  taxa <- left_join(taxa, all_species, by = "ala_search_term", all.x=TRUE)
  return(taxa)
}
# test35 <- add_sp_column(taxa)
# test35[,c(1,2,99)]

# combined function
unverified_taxa <- function(taxa){
  sub1 <- add_sp_column(taxa) |>
    dplyr::filter(is.na(count) | count > 1) |>
    dplyr::select(ala_search_term, common_name, count)
  sub1 <- search_taxa(sub1$ala_search_term)
  sub1 <- sub1[,c(1,2,5,6,11:15)]
  return(sub1)
}

# splitting the above function into either taxa with count of NA
unverified_taxa_na <- function(taxa){
  sub1 <- add_sp_column(taxa) |>
    dplyr::filter(is.na(count)) |>
    dplyr::select(ala_search_term, common_name, count)
  sub1 <- search_taxa(sub1$ala_search_term)
  sub1 <- sub1[,c(1,2,5,6,11:15)]
  return(sub1)
}
# or taxa with count > 1 (i.e. including subspecies)
##  these are not necessarily wrong
unverified_taxa_ex <- function(taxa){
  sub1 <- add_sp_column(taxa) |>
    dplyr::filter(count > 1) |>
    dplyr::select(ala_search_term, common_name, count)
  sub1 <- search_taxa(sub1$ala_search_term)
  sub1 <- sub1[,c(1,2,5,6,11:15)]
  return(sub1)
}


# import a version of the taxa database
OCT4_22_TAXA_CSV <- "Vic_data_4Oct22a.csv"
OCT4_22_PATH <- file.path(datapath, OCT4_22_TAXA_CSV)

taxa2 <- read.csv(OCT4_22_PATH, header = TRUE)
# head(taxa2[,c(1:3,38)])


# omit non-ALA taxa (83 taxa as at 4th Oct 2022)
taxa2 <- taxa2 |> dplyr::filter(assess == "ALA")
# divide into 4 blocks as this fails for the full data set all at once
div4 <- nrow(taxa2)/4
taxa2a <- taxa2[1:div4,] # 1.8 mins, 53 records
taxa2b <- taxa2[(div4+1):(div4*2),] # 2.244 mins, 42 records
taxa2c <- taxa2[((div4*2)+1):(div4*3),] # 2.17 mins, 27 records
taxa2d <- taxa2[((div4*3)+1):(div4*4),] # 2.024 mins, 39 records

start_time1 <- Sys.time()
test35a <- unverified_taxa(taxa2a) # 53 obs
test35b <- unverified_taxa(taxa2b) # 42 obs
test35c <- unverified_taxa(taxa2c) # 27 obs
test35d <- unverified_taxa(taxa2d) # 39 obs
test36 <- rbind(test35a, test35b, test35c, test35d) # 161 obs total
end_time1 <- Sys.time()
end_time1 - start_time1
# test35[,c(1,2,99)]

test36_probs <- subset(test36, test36$match_type != "exactMatch" |
          test36$issues != "noIssue") ## 44 taxa with some sort of issue

start_time1 <- Sys.time()
test35a_na <- unverified_taxa_na(taxa2a) # 38 obs
test35b_na <- unverified_taxa_na(taxa2b) # 15 obs
test35c_na <- unverified_taxa_na(taxa2c) # 12 obs
test35d_na <- unverified_taxa_na(taxa2d) # 20 obs
test36_na <- rbind(test35a_na, test35b_na, test35c_na, test35d_na)
end_time1 <- Sys.time()
end_time1 - start_time1 # 6.33 mins, 85 obs total
# test36_na[,c(1,2,99)]

start_time1 <- Sys.time()
test35a_ex <- unverified_taxa_ex(taxa2a) # 15 obs
test35b_ex <- unverified_taxa_ex(taxa2b) # 27 obs
test35c_ex <- unverified_taxa_ex(taxa2c) # 15 obs
test35d_ex <- unverified_taxa_ex(taxa2d) # 19 obs
test36_ex <- rbind(test35a_ex, test35b_ex, test35c_ex, test35d_ex)
end_time1 <- Sys.time()
end_time1 - start_time1
# test36_ex[,c(1,2,99)] # 6.37 mins, 76 obs total






#### IS THIS ANY USE???
verified_taxa3 <- function(taxa) {
  for (i in 1:nrow(taxa)) {
    if(taxa$ala_search_term == unverified_taxa$scientific_name){
      # don't do anything
    }
    else if(taxa$ala_search_term != unverified_taxa$scientific_name){
      # taxa$ala_search_term < unverified_taxa$$scientific_name)
    }
    else{
      # do something else - use gsub, or remove from taxa to process
      #  and save as separate list to verify manually (and naport to .csv?)
    }
  }
}





