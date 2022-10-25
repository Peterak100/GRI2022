## step through the 'process_observation' function to help find any errors..
# check that 'isolation_taxa' looks as expected
isolation_taxa[, c(1,2,39,43,94:97)]
# test a single taxon (substitute row number for taxon to be tested)

taxonA <- isolation_taxa[18, ]
taxonA[c(1,2)]

  # Add new columns to single 'taxon' dataframe
taxonA <- add_column(taxonA, 
            num_preclusters = 0, 
            num_orphans = 0, 
            precluster_cellcount = 0, 
            orphan_cellcount = 0,
            error = NA)

# num_preclusters <- 0 # still not sure if this is needed??
## only used inside the function?

# 'test79' here replaces 'obs'
test79 <- load_and_filter(taxonA, taxapath) ## up to testing this line..


# load_and_filter combines:
###  load_or_download_obs   # returns obs
# test73 <- load_or_download_obs(taxonA, taxapath)
##  filter_observations  # removes missing coords & location duplicates
# test74 <- filter_observations(test73, taxonA)
##  precluster_observations  # adds precluster number based on epsilon
# test77 <- precluster_observations(test74, taxonA)
# to count orphans:
# subset(test77, test77[9]==0) |> count()
# to count preclusters:
# max(test77[9])


# if there is at least one (or more) precluster..
# # # test41 <- write_precluster(test79, taxonA, mask_layer, taxapath)
## displays Num preclusters
## writes precluster and orphans plots to taxon directory

### breaking down the 'write_precluster' function..
taxonpath <- suppressWarnings(taxon_path(taxonA, taxapath))
shapes31 <- sf::st_as_sf(test79, coords = c("x", "y"), crs = METRIC_EPSG)
scaled_eps31 <- taxonA$epsilon * 1000 * EPSILON_SENSITIVITY_SCALAR
preclustered_obs31 <- buffer_preclustered(shapes31, scaled_eps31)
precluster_rast31 <- shapes_to_raster(preclustered_obs31, taxonA,
          mask_layer, taxonpath)

png(file = file.path(taxonpath, paste0(gsub(" ", "_",
    taxonA$ala_search_term), "_preclusters31", ".png")),
    width = 2160, height = 1440, pointsize = 30)
plot(precluster_rast31, main = paste(taxonA$ala_search_term, "preclusters31"))
while (!is.null(dev.list())) dev.off()

pixel_freq31 <- terra::freq(precluster_rast31)
# pixel_freq31 should have pixel frequencies for each precluster
pclust_counts31 <- left_join(shapes31, pixel_freq31, copy = TRUE,
      by = c("precluster" = "value")) |> 
  write_csv(file.path(taxonpath, paste0(gsub(" ","_",
        taxonA$ala_search_term), "_precluster_counts31", ".csv")))

# also save an info file with most recent year, pixel count & latin name
## need to retrieve this for post-processing
pclust_info31 <- dplyr::filter(pclust_counts31, precluster != 0) |> 
  dplyr::group_by(precluster) |> 
  dplyr::summarise(pix_count = max(count, na.rm = TRUE),
        recent_year = max(year, na.rm = TRUE), latin = max(scientificName))
pclust_info31 <- pclust_info31 |> 
  add_column(cluster = pclust_info31$precluster) # 6 columns

# add count of observations per precluster
pclust_recs31 <- dplyr::filter(pclust_counts31, precluster !=0) |> 
  dplyr::count(precluster) |> sf::st_drop_geometry() |> 
  dplyr::rename(num_obs = n)
pclust_info31 <- left_join(pclust_info31, pclust_recs31, by = "precluster")
# 7 columns

# save a .csv file to retrieve for postprocessing
pclust_info31 |> write_csv(file.path(taxonpath, paste0(gsub(" ","_", 
          taxonA$ala_search_term), "_preclusters_prelim", ".csv")))


# create a 'units' matrix of distances between polygons
# then convert to normal numeric matrix & convert to kilometres
prox31 <- sf::st_distance(pclust_info31) |> as.data.frame() |> as.matrix()
prox31 <- prox31/1000

# add two rows and save for use as mask_file in Circuitscape?
## not sure this either works, or is particularly useful
# mask31 <- prox31[1:2,]
# mask31[1:2,] <- "" # or NA ?
# mask31[1,1] <- "min"
# mask31[1,2] <- 0
# mask31[2,1] <- "max"
# mask31[2,2] <- taxonA$epsilon * 40 # IS 40 A GOOD VALUE HERE?
# mask31 <- rbind(mask31, prox31)
## NEED TO WRITE TO CORRECT FILE PATH (i.e. taxapath)
# write.table(mask31, "mask.txt", row.names = FALSE, col.names = FALSE)

orphan_obs31 <- buffer_orphans(shapes31, scaled_eps31)
orphan_rast31 <- shapes_to_raster(orphan_obs31, taxon, mask_layer, taxonpath)

png(file = file.path(taxonpath, paste0(gsub(" ","_", taxonA$ala_search_term),
      "_orphans31", ".png")), width = 2160, height = 1440, pointsize = 30)
plot(orphan_rast31, main = paste(taxonA$ala_search_term, "orphans31"))
while (!is.null(dev.list())) dev.off()

dplyr::filter(shapes31, precluster == 0) |>
  mut_euclidean_coords() |> write_csv(file.path(taxonpath,
      paste0(gsub(" ", "_", taxonA$ala_search_term), "_orphans31", ".csv")))

crop_rast31 <- terra::merge(precluster_rast31, orphan_rast31) |>  
  padded_trim()

## crop and write rasters as .tif files
precluster_filename <- file.path(taxonpath, paste0("preclusters31", ".tif"))
orphan_filename <- file.path(taxonpath, paste0("orphans31", ".tif"))
terra::crop(precluster_rast31, crop_rast31,
            filename = precluster_filename, overwrite = TRUE)
terra::crop(orphan_rast31, crop_rast31,
            filename = orphan_filename, overwrite = TRUE)

precluster_cellcount31 <- sum(terra::freq(precluster_rast31))
orphan_cellcount31 <- sum(terra::freq(orphan_rast31)) # sum makes it numeric
# this is returned as cell_counts in try_taxon_observations()
output31 <- (c(precluster_cellcount31, orphan_cellcount31))

# taxon$filter_category <- taxon$filter_category ## redundant?
taxonA$num_preclusters <- max(test79$precluster)
taxonA$num_orphans <- sum(test79$precluster == 0)
taxonA$precluster_cellcount <- output31[1]
taxonA$orphan_cellcount <- output31[2]

# check that added columns contain correct (plausible) data..
taxonA[c(1,2,94:102)]
taxonB <- taxonA
taxonB$orphan_cellcount <- 0

# check that label function for 'risk' & 'filter_category' works as expected
taxon_processed79 <- label_by_clusters(taxonA)
taxon_processed79[c(1,2,94:102)]


# as above, but broken down by individual label function..
## if further testing is needed
# taxon_processed79 <- label_high_orphan_area(taxonA)
# taxon_processed79 <- label_many_clusters(taxon_processed79)
# taxon_processed79 <- label_few_clusters(taxon_processed79)
# taxon_processed79 <- label_no_clusters(taxon_processed79)



# Dispersal categories --------
### increments are based on a 10^x log scale
###  where x increases in increments of 0.2

# 0.63 km : 'Minimal' 10^-0.2
# 1.00 km : 'Slight' 10^0
# 1.58 km : 'Very low'  10^0.2
# 2.51 km : 'Low'  10^0.4
# 3.16 km : 'Low+'  ## EXTRA CATEGORY  10^0.5
# 3.98 km : 'Modest' 10^0.6
# 6.31 km : 'Modest+' 10^0.8
# 10.0 km : 'Medium' 10^1
# 15.8 km : 'Medium+' 10^1.2
# 25.1 km : 'High' 10^1.4
# 39.8 km : 'High+' 10^1.6
# 63.1 km : 'Very high' 10^1.8
# 100 km : 'Extensive' 10^2

