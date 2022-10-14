# for testing of these functions see: test4.R

# Invert percentage from % habitat quality to % dispersal resistance
habitat_to_resistance <- function(habitat_raster) {
  101 - habitat_raster
}

crop_resistance <- function(resistance_raster, crop_filename) {
  crop_template <- terra::rast(crop_filename)
  # Crop the resistance_raster to match "preclusters.tif"
  terra::crop(resistance_raster, crop_template)
}

# if (taxon$resist_model_type[[1]] == "Species")
## -> run  download_hdm(taxon, taxapath, crop_filename)
## urls to enable download of hdms seems to have changed
## otherwise run use_generic_hdm

use_generic_hdm <- function(taxon, taxapath, crop_filename) {
  cat("Using generic resistance HDM for", taxon$ala_search_term, "\n")
  resistance_filename <- suppressWarnings(file.path(taxon_path(taxon, taxapath),
          RESISTANCE_RASTER))
  terra::rast(HABITAT_RASTER_PATH) |> 
    habitat_to_resistance() |> 
    crop_resistance(crop_filename) |> 
    terra::writeRaster(filename = resistance_filename, overwrite=TRUE)
}

## this is a minimal version of this function for now without a tryCatch
##  and only using the generic resistance layer
prepare_resistance_files <- function(taxa, taxapath) {
  if (nrow(taxa) == 0) {
    return()
  }
  for (i in 1:nrow(taxa)) {
    taxon <- taxa[i, ]
    crop_filename <- suppressWarnings(file.path(taxon_path(taxon, taxapath),
            paste0("preclusters", ".tif")))
    if (file.exists(crop_filename)) {
      use_generic_hdm(taxon, taxapath, crop_filename)
    }
  }
  return()
}

