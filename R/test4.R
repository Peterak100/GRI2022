# test preparation of single resistance files for Circuitscape..
# isolation_by_resistance_taxa[, c(1,2,39,42,94:102)]
taxonB <- isolation_by_resistance_taxa[14,]
taxonB[c(1,2)]


## path to relevant 'preclusters.tif' file
crop_filename23 <- suppressWarnings(file.path(taxon_path(taxonB, taxapath),
            paste0("preclusters", ".tif")))

        if (taxonB$resist_model_type[[1]] == "Species") {
          download_hdm(taxonB, taxapath, crop_filename)
        } else {
          use_generic_hdm(taxonB, taxapath, crop_filename)
        }



# this generates and saves a particular version of the 'resistance.tif' file
## derived from the generic 'habitat.tif' but converted to resistance
## values and cropped to match the particular 'preclusters.tif' file
use_generic_hdm(taxonB, taxapath, crop_filename23)



# Download the HDM layer for this taxon
# Formatted as:
# "https://maps2.biodiversity.vic.gov.au/Models/
##     SMP_Dromaius%20novaehollandiae_Emu_10001.zip"

# download_hdm <- function(taxon, taxapath, crop_filename) 

# cat("Downloading specific habitat layer for", taxonB$ala_search_term, "\n")

taxonB_id <- taxonB$vic_taxon_id[[1]] # col 4 (a 5-digit number)

taxonB_escaped <- gsub(" ", "%20", taxonB$delwp_taxon)[[1]]
common_nameB <- gsub(" ", "%20", taxonB$delwp_common_name)[[1]]
url <- paste0("https://maps2.biodiversity.vic.gov.au/Models/SMP_",
            taxonB_escaped, "_", common_nameB, "_", taxonB_id, ".zip")
taxonB_dir <- taxon_path(taxonB, taxapath)
downloadB_dir <- file.path(taxonB_dir, "download")
dir.create(downloadB_dir, recursive = TRUE)
zippathB <- file.path(downloadB_dir, "hdm.zip")
download.file(url, zippathB) ### DOESN'T WORK changed url?


# If the download doesn't exist, download it
if (!dir.exists(download_dir) || is.na(Sys.glob(file.path(download_dir,
                  "*.tif"))[1])) { 
    download.file(url, zippath)
    unzip(zippath, exdir=download_dir)
    file.remove(zippath)
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







