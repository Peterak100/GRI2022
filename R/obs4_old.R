## previous code for obs4.R

pixel_freq |> 
  write_csv(file.path(taxonpath, paste0(gsub(" ","_", taxon$ala_search_term),
                  "_precluster_sizes", ".csv")))
pclust_counts <- left_join(shapes, pixel_freq, copy = TRUE,
                  by = c("precluster" = "value")) |> 
  write_csv(file.path(taxonpath, paste0(gsub(" ","_", taxon$ala_search_term),
                  "_precluster_counts", ".csv")))

# also save an info file with most recent year, pixel count & latin name
## need to retrieve this for post-processing
pclust_info <- dplyr::filter(pclust_counts, precluster != 0) |> 
  dplyr::group_by(precluster) |> 
  summarize(pix_count = max(count, na.rm = TRUE),
            recent_year = max(year, na.rm = TRUE),
            latin = max(scientificName))
# add new column for cluster number 
pclust_info <- pclust_info |> 
  add_column(cluster = pclust_info$precluster, pop_name = NA)

# add count of observations per precluster
pclust_recs <- dplyr::filter(pclust_counts, precluster !=0) |> 
  dplyr::count(precluster) |> sf::st_drop_geometry() |> 
  dplyr::rename(num_obs = n)
pclust_info <- left_join(pclust_info, pclust_recs, by = "precluster")

# create a 'units' matrix of distances between polygons
# then convert to normal numeric matrix & convert to kilometres
prox <- sf::st_distance(pclust_info) |> as.data.frame() |> 
  as.matrix()/1000
# add two rows and save for use as mask_file in Circuitscape?
# Not sure this even works, or is particularly useful
# mask <- prox[1:2,]
# mask[1:2,] <- "" # or NA ?
# mask[1,1] <- "min"
# mask[1,2] <- 0
# mask[2,1] <- "max"
# mask[2,2] <- taxonA$epsilon * 40
# mask <- rbind(mask, prox)
## NEED TO WRITE TO CORRECT FILE PATH (i.e. taxapath) & FILE NAME
# write.table(mask, "mask.txt", row.names = FALSE, col.names = FALSE)

# write as a separate csv file for subsequent post processing
pclust_info |> write_csv(file.path(taxonpath, paste0(gsub(" ","_",
          taxon$ala_search_term), "_preclusters_prelim", ".csv"))) 

# Create a dataframe and raster for orphans
orphan_obs <- buffer_orphans(shapes, scaled_eps)
orphan_rast <- shapes_to_raster(orphan_obs, taxon, mask_layer, taxonpath)

# Save a plot of orphans
png(file = file.path(taxonpath, paste0(gsub(" ","_", taxon$ala_search_term),
          "_orphans", ".png")), width = 2160, height = 1440, pointsize = 30)
plot(orphan_rast, main = paste(taxon$ala_search_term, "orphans"))
cat("orphans plot for", taxon$ala_search_term, "saved to:",
    taxonpath, "\n")
while (!is.null(dev.list())) dev.off()

# Write an orphans csv # IS THIS NEEDED?
# need to add extra step to convert 'x' & 'y' back to lon/lat?

dplyr::filter(shapes, precluster == 0) |>
  mut_euclidean_coords() |> write_csv(file.path(taxonpath,
          paste0(gsub(" ", "_", taxon$ala_search_term), "_orphans", ".csv")))

# Make a crop template by trimming the empty values from a
# combined precluster/orphan raster, with some added padding.
crop_rast <- terra::merge(precluster_rast, orphan_rast) |>  
  padded_trim()

# Crop and write rasters as .tif files to relevant taxon folder
precluster_filename <- file.path(taxonpath, paste0("preclusters", ".tif"))
orphan_filename <- file.path(taxonpath, paste0("orphans", ".tif"))
terra::crop(precluster_rast, crop_rast,
          filename = precluster_filename, overwrite = TRUE)
terra::crop(orphan_rast, crop_rast,
          filename = orphan_filename, overwrite = TRUE)

precluster_cellcount <- sum(terra::freq(precluster_rast))
orphan_cellcount <- sum(terra::freq(orphan_rast)) # sum makes this numeric
# this is returned as cell_counts in try_taxon_observations()
## ALSO NEED TO RETURN pclust_info FOR POST-PROCESSING
return(c(precluster_cellcount, orphan_cellcount))
}

# Try taxon observations function --------
### Oct 2022: this function not currently used due to assumption that
###  everything will run without errors!
## Block to run either with or without a try/catch block
## downloads, filters & preclusters observation records
## include 'force_download' as 3rd argument?
## also as 3rd argument for 2nd line?
try_taxon_observations <- function(taxon, taxapath) {
  obs <- load_and_filter(taxon, taxapath)
  # Create rasters with numbered preclustered observations
  # If there are any preclusters
  if (max(obs$precluster) != 0) {
    cell_counts <- write_precluster(obs, taxon, mask_layer, taxapath)
  } else {
    cell_counts <- c(0, 0)
  }
  taxon$error_string <- NA # ?? for errors?
  taxon$filter_category <- taxon$filter_category # ?? for errors?
  taxon$num_preclusters <- max(obs$precluster)
  taxon$num_orphans <- sum(obs$precluster == 0)
  taxon$precluster_cellcount <- cell_counts[1]
  taxon$orphan_cellcount <- cell_counts[2]
  list(taxon$error_string, taxon$filter_category, taxon$num_preclusters,
       taxon$num_orphans, taxon$precluster_cellcount, taxon$orphan_cellcount)
}

