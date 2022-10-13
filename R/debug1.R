# Can't rename columns that don't exist
# Column 'species' doesn't exist
rlang::last_error() # to see where the error occurred
rlang::last_trace() # to see full context

precategorized_taxa <- precategorize_risk(taxa)

# precategorize_risk(taxa) is from categorize.R
## is this function even needed?
##   just split main database into chunks manually from the start?

precategorized_taxa <- precategorize_chunk(taxa)

## where is load_and_filter defined -> line 76 observations.R
# setting force_download=FALSE prior to running this does the same thing
obs2 <- load_and_filter(taxon, taxapath, force_download=FALSE)


gsub(character, new_character, string) 
zz1 <- "donkey test"
zz1a <- gsub(" ", "_", zz1)
