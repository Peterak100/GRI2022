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


use_generic_hdm <- function(taxon, taxapath, crop_filename) {
  cat("Using generic resistance HDM for", taxon$ala_search_term, "\n")
  resistance_filename <- file.path(taxon_path(taxon, taxapath),
          RESISTANCE_RASTER)
  terra::rast(HABITAT_RASTER_PATH) |> 
    habitat_to_resistance() |> 
    crop_resistance(crop_filename) |> 
    terra::writeRaster(filename = resistance_filename, overwrite=TRUE)
}