# test preparation of single resistance files for Circuitscape..
## Habitat Importance Models (HIMs) are single geotiffs, 75m2 / pixel
## Habitat Distribution Models (HDMs) are zip files, containing
### a single geotiff, 225m2 / pixel

# isolation_by_resistance_taxa[, c(1,2,39,42,94:102)]
taxonB <- isolation_by_resistance_taxa[7,]
taxonB[c(1,2)]

taxonpath <- suppressWarnings(taxon_path(taxonB, taxapath))
paste0(gsub(" ", "_", taxon$ala_search_term), "_orphans", ".csv")
## path to relevant 'preclusters.tif' file
crop_filename31 <- suppressWarnings(file.path(taxon_path(taxonB, taxapath),
            paste0("preclusters31", ".tif")))

HIM_FILE31 <- paste0("HIM_",gsub(" ","_", taxonB$ala_search_term), ".tif")
HIM_PATH31 <- file.path(taxonpath, HIM_FILE31)
SMP_FILE31 <- paste0("SMP_",gsub(" ","_", taxonB$ala_search_term), ".tif")
SMP_PATH31 <- file.path(taxonpath, SMP_FILE31)

## TO DO: add if (file.exists(HIM_PATH31)) {} as a subsequent step..
if (file.exists(SMP_PATH31)) {
  use_SMP_layer(taxonB, taxapath, crop_filename31)
  } else {
    use_generic_hdm(taxonB, taxapath, crop_filename31)
}


# this generates and saves a particular version of the 'resistance.tif' file
## derived from the generic 'habitat.tif' but converted to resistance
## values and cropped to match the particular 'preclusters.tif' file
use_generic_hdm(taxonB, taxapath, crop_filename31)


use_generic_hdm <- function(taxon, taxapath, crop_filename) {
  cat("Using generic resistance model for", taxon$ala_search_term, "\n")
  resistance_filename <- suppressWarnings(file.path(taxon_path(taxon,
                  taxapath), RESISTANCE_RASTER))
  terra::rast(HABITAT_RASTER_PATH) |> 
    habitat_to_resistance() |> 
    crop_resistance(crop_filename) |> 
    terra::writeRaster(filename = resistance_filename, overwrite = TRUE)
}

use_SMP_layer <- function(taxon, taxapath, crop_filename) {
  cat("Using HDM-derived resistance layer for", taxon$ala_search_term, "\n")
  terra::rast(SMP_filename) |> 
    habitat_to_resistance() |> 
    crop_resistance(crop_filename) |> 
    terra::writeRaster(filename = resistance_filename, overwrite = TRUE)
}





habitat_filename <- Sys.glob(file.path(download_dir, "*.tif"))[1]
resistance_filename <- file.path(taxon_path(taxon, taxapath),
                RESISTANCE_RASTER)
terra::rast(habitat_filename) %>%
habitat_to_resistance() %>%
crop_resistance(crop_filename) %>%
terra::writeRaster(filename=resistance_filename, overwrite=TRUE)


### MANUALLY TEST HDM .tif file for Varanus varius (as a test)
###  because automated downloads not working (yet)..
###  manually copied to datapath and renamed as 'habitat_varanus.tif'
### need to import this as "Varanus_resistance.tif"


########## NOT YET TRIED.....
VARANUS_RASTER <- "habitat_varanus.tif"
VARANUS_RES_RASTER <- "Varanus_resistance.tif"
VARANUS_RASTER_PATH <- file.path(datapath, VARANUS_RASTER)
resistance2_filename <- file.path(taxon_path(taxon, taxapath),
              VARANUS_RES_RASTER)
terra::rast(VARANUS_RASTER_PATH) |> 
  habitat_to_resistance() |> 
  crop_resistance(crop_filename) |> 
  terra::writeRaster(filename = resistance_filename, overwrite=TRUE)








# Use the generic resistance layer --------
use_generic_hdm <- function(taxon, taxapath, crop_filename)

# cat("Using generic resistance HDM for", taxonB$ala_search_term, "\n")

resistance_filename <- file.path(taxon_path(taxonB, taxapath),
              RESISTANCE_RASTER)
terra::rast(HABITAT_RASTER_PATH) |> 
      habitat_to_resistance() |> 
      crop_resistance(crop_filename23) |> 
      terra::writeRaster(filename = resistance_filename, overwrite=TRUE)







