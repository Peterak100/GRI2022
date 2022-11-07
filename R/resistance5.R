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

# subfunction to choose the correct model for the resistance layer
## priority is to use the species-specific Habitat Distribution Model (SMP),
## then the Habitat Importance Model (HIM) if there is no SMP file,
## otherwise use the generic Strategic Biodiversity Values (SBV) layer
choose_resistance_model <- function(taxon, taxapath, crop_filename) {
  SMP_MODEL <- paste0("SMP_", gsub(" ","_", taxon$ala_search_term), ".tif")
  SMP_filename <- suppressWarnings(file.path(SMPpath, SMP_MODEL))
  HIM_MODEL <- paste0("HIM_", gsub(" ","_", taxon$ala_search_term), ".tif")
  HIM_filename <- suppressWarnings(file.path(HIMpath, HIM_MODEL))
  if(file.exists(SMP_filename)) {
    res_model <- SMP_filename
    cat("Using HDM-derived resistance layer for",
        taxon$ala_search_term, "\n")
  } else if (file.exists(HIM_filename)) {
    res_model <- HIM_filename
    cat("Using HIM-derived resistance layer for",
        taxon$ala_seach_term, "\n")
    } else {
        res_model <- HABITAT_RASTER_PATH
        cat("Using generic resistance model for",
            taxon$ala_search_term, "\n")
    }
  return(res_model)
}

## this function is without a tryCatch for now
prepare_resistance_files <- function(taxa, taxapath) {
  if (nrow(taxa) == 0) {
    return()
  } else {
    for (i in 1:nrow(taxa)) {
      taxon <- taxa[i, ]
      crop_filename <- suppressWarnings(file.path(taxon_path(taxon,
                taxapath), paste0("preclusters", ".tif")))
      resistance_filename <- suppressWarnings(file.path(taxon_path(taxon,
                taxapath), RESISTANCE_RASTER))
    if (file.exists(crop_filename)) {
      choose_resistance_model(taxon, taxapath, crop_filename)
      terra::rast(res_model) |> 
        habitat_to_resistance() |> 
        crop_resistance(crop_filename) |> 
        terra::writeRaster(filename = resistance_filename, overwrite = TRUE)
      # read in a vanilla version of the Circuitscape.ini file with
      # standard settings for a Circuitscape pairwise run
      Cscape1path <- file.path(datapath,
                paste0("circuitscape_pwise", ".ini"))
      CSrun <- ini::read.ini(Cscape1path, encoding = getOption("encoding"))
      # produces a list of 11
      # update .ini file with source files relevant to particular taxon:
      CSrun$`Habitat raster or graph`$habitat_file <-
        file.path(taxonpath, paste0("resistance", ".tif"))
      CSrun$`Options for pairwise and one-to-all and
      all-to-one modes`$point_file <-
        file.path(taxonpath, paste0("preclusters.tif"))
      CSrun$`Output options`$output_file <-
        file.path(taxonpath, paste0("CSpwise1"))
      # save as a customized .ini file in the relevant directory
      write.ini(CSrun, file.path(taxonpath,
        paste0("Circuitscape_custom1", ".ini")))
    }
  }  
  }
}

