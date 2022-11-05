# for testing of these functions see: test4.R

# Invert percentage from % habitat quality to % dispersal resistance
habitat_to_resistance <- function(habitat_raster) {
  101 - habitat_raster
}

crop_resistance <- function(resistance_raster, crop_filename) {
  crop_template <- terra::rast(crop_filename)
  # Crop the resistance_raster to match "preclusters.tif"
  ## this is actually midclusters if there is an AoO model
  terra::crop(resistance_raster, crop_template)
}

# if (taxon$resist_model_type[[1]] == "Species")
## -> run  download_hdm(taxon, taxapath, crop_filename)
## urls to enable download of hdms seems to have changed
## otherwise run use_generic_hdm

use_generic_hdm <- function(taxon, taxapath, crop_filename) {
  cat("Using generic resistance model for", taxon$ala_search_term, "\n")
  resistance_filename <- suppressWarnings(file.path(taxon_path(taxon,
          taxapath), RESISTANCE_RASTER))
  terra::rast(HABITAT_RASTER_PATH) |> 
    habitat_to_resistance() |> 
    crop_resistance(crop_filename) |> 
    terra::writeRaster(filename = resistance_filename, overwrite = TRUE)
}

# should be okay just to overwrite "resistance.tif" with either model?
## TO DO: change 'resistance2_filename' back to 'resistance_filename'
use_SMP_layer <- function(taxon, taxapath, crop_filename) {
  cat("Using HDM-derived resistance layer for", taxon$ala_search_term, "\n")
  resistance2_filename <- suppressWarnings(file.path(taxon_path(taxon,
          taxapath), "resistance2.tif"))
  SMP_MODEL <- paste0("SMP_",gsub(" ","_",
              taxon$ala_search_term), ".tif")
  SMP_filename31 <- suppressWarnings(file.path(taxon_path(taxonB,
              taxapath), SMP_MODEL31))
  terra::rast(SMP_filename31) |> 
    habitat_to_resistance() |> 
    crop_resistance(crop_filename31) |> 
    terra::writeRaster(filename = resistance2_filename, overwrite = TRUE)
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
      # read in a vanilla version of the Circuitscape.ini file with standard
      # settings for a Circuitscape pairwise run
      Cscape1path <- file.path(datapath, paste0("circuitscape_pwise.ini"))
      CSrun <- read.ini(Cscape1path, encoding = getOption("encoding"))
      # produces a list of 11
      # update the .ini file with source files relevant to the particular taxon:
      CSrun$`Habitat raster or graph`$habitat_file <- 
        file.path(taxonpath, paste0("resistance.tif"))
      CSrun$`Options for pairwise and one-to-all and all-to-one modes`$point_file <-
        file.path(taxonpath, paste0("preclusters31.tif"))
      CSrun$`Output options`$output_file <-
        file.path(taxonpath, paste0("CSpwise1"))
      # save as a customized .ini file in the relevant directory
      write.ini(CSrun, file.path(taxonpath,
            paste0("Circuitscape_custom1.ini")),
            encoding = getOption("encoding"))
    }
  }
  return()
}




