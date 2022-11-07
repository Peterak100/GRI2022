## step through the 'process_observations' function to help find any errors..
# check that 'isolation_taxa' looks as expected
isolation_taxa[, c(1,2,39,43,94:97)]
# test a single taxon (substitute row number for taxon to be tested)

taxonA <- isolation_taxa[9, ]
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
test76 <- read_cached_observations(taxonA, taxapath)
test77 <- download_observations(taxonA)
test78 <- load_or_download_obs(taxonA, taxapath)
test79 <- load_and_filter(taxonA, taxapath)


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
preclustered_obs31 <- buffer_preclustered(shapes31, scaled_eps31) # 2 columns
precluster_rast31 <- shapes_to_raster(preclustered_obs31, taxonA,
          mask_layer, taxonpath)

png(file = file.path(taxonpath, paste0(gsub(" ", "_",
    taxonA$ala_search_term), "_preclusters31", ".png")),
    width = 2160, height = 1440, pointsize = 30)
plot(precluster_rast31, main = paste(taxonA$ala_search_term, "preclusters31"))
while (!is.null(dev.list())) dev.off()

pixel_freq31 <- terra::freq(precluster_rast31) |> 
  dplyr::rename(pix_count = count)
# pixel_freq31 should have pixel frequencies for each precluster
pclust_info31 <- left_join(preclustered_obs31, pixel_freq31[-1],
        by = c("precluster" = "value")) # don't need 'layer' # 3 columns
# this messes up geometry, so drop it from this dataframe
pclust_summary31 <- dplyr::filter(shapes31, precluster != 0) |> 
  dplyr::group_by(precluster) |> 
  dplyr::summarise(recent_year = max(year, na.rm = TRUE),
        latin = max(scientificName)) |> st_drop_geometry()

pclust_info31 <- left_join(pclust_info31, pclust_summary31, by = "precluster")
# 5 columns

pclust_num31 <- dplyr::filter(shapes31, precluster !=0) |> 
  dplyr::count(precluster) |> sf::st_drop_geometry() |> 
  dplyr::rename(num_obs = n)
pclust_info31 <- left_join(pclust_info31, pclust_num31, by = "precluster")
# 6 columns
midcluster_cellcount31 <- 0

orphan_obs31 <- buffer_orphans(shapes31, scaled_eps31)
orphan_rast31 <- shapes_to_raster(orphan_obs31, taxon, mask_layer, taxonpath)

png(file = file.path(taxonpath, paste0(gsub(" ","_", taxonA$ala_search_term),
      "_orphans31", ".png")), width = 2160, height = 1440, pointsize = 30)
plot(orphan_rast31, main = paste(taxonA$ala_search_term, "orphans31"))
while (!is.null(dev.list())) dev.off()

dplyr::filter(shapes31, precluster == 0) |>
  mut_euclidean_coords() |> write_csv(file.path(taxonpath,
      paste0(gsub(" ", "_", taxonA$ala_search_term), "_orphans31", ".csv")))

# crop and write orphan rasters as .tif file
orphan_filename <- file.path(taxonpath, paste0("orphans31", ".tif"))
crop_orph_rast31 <- orphan_rast31 |>  padded_trim()
terra::crop(orphan_rast31, crop_orph_rast31,
      filename = orphan_filename, overwrite = TRUE)


## if AoO file exists --------
# Note the 7 supporting files in addition to .shp must also be present
# this creates and returns an new version of pclust_info31
AoO_FILE <- paste0("AoO_",gsub(" ","_", taxonA$ala_search_term), ".shp")
AoO_PATH <- file.path(AoOpath, AoO_FILE)
if (file.exists(AoO_PATH)){
  AoO1 <- sf::st_read(AoO_PATH)
  AoO1 <- AoO1[,c(2,13)] # a single multipolygon
  # separates single multipolygon into separate polygons, then buffers
  ## buffer distance here is < buffer for individual observations
  AoO1 <- suppressWarnings(sf::st_cast(AoO1, "POLYGON")) |> 
    sf::st_buffer(dist = (taxonA$epsilon * 1000) / 3) |> 
    dplyr::rename(latin = SCIENTIFIC)
  # number new preclusters
  # numbering to following on from preclusters derived from obs data
  AoO1 <- AoO1 |> 
    add_column(precluster = 
                 1+nrow(pclust_info31):(nrow(AoO1)+nrow(pclust_info31)-1),
               pix_count = 0, recent_year = 2018, num_obs = 0)
  pclust_sort <- c("precluster","geometry","pix_count",
                   "recent_year","latin","num_obs")
  # why does position of 'geometry' change?
  AoO1 <- rbind(pclust_info31, AoO1[,pclust_sort])
  AoO1_parts <- sf::st_cast(sf::st_union(AoO1), "POLYGON")
  midcluster <- unlist(sf::st_intersects(AoO1, AoO1_parts))
  mclust_info31 <- cbind(AoO1, midcluster) |> dplyr::group_by(midcluster) |> 
    dplyr::summarise(precluster = paste(precluster, collapse = ", "),
        pix_count = max(pix_count, na.rm = TRUE),
        recent_year = max(recent_year, na.rm = TRUE),
        latin = max(latin, na.rm = TRUE),
        num_obs = max(num_obs, na.rm = TRUE))
  pclust_info31 <- mclust_info31
  
#  getting rid of this subfunction - too confusing  
#  pclust_info31 <- merge_AoO_regions(taxonA)
# recalculate 'pix_count' for midclusters
  mclust_retain <- c("midcluster","geometry")
  midcluster_obs31 <- pclust_info31[, mclust_retain]
  # need to use a different function here because 'shapes_to_raster()'
  #  uses field = "precluster"
  midcluster_rast31 <- mid_shapes_to_raster(midcluster_obs31, taxonA,
          mask_layer, taxonpath) # plot(midcluster_rast31)
  # also need to save a .tif file for Circuitscape patches layer
  crop_mid_rast31 <- terra::merge(midcluster_rast31, orphan_rast31) |> 
    padded_trim()
  # plot(crop_mid_rast31)
  preclust_filename <- file.path(taxonpath, paste0("preclusters31", ".tif"))
  # crop and write midcluster raster as .tif file
  terra::crop(midcluster_rast31, crop_mid_rast31,
              filename = preclust_filename, overwrite = TRUE)
  # recalculate 'pix_count'
  pixel_mid_freq31 <- terra::freq(midcluster_rast31)
  pclust_info31$pix_count <- pixel_mid_freq31$count
  midcluster_cellcount31 <- sum(terra::freq(midcluster_rast31))
} else {
  mclust_order <- c("midcluster", "precluster", "pix_count", "recent_year",
              "latin", "num_obs", "geometry")
  pclust_info31 <- pclust_info31 |> 
    add_column(midcluster = pclust_info31$precluster)
  pclust_info31 <- pclust_info31[, mclust_order]
  crop_rast31 <- terra::merge(precluster_rast31, orphan_rast31) |>  
    padded_trim()
  # plot(crop_rast31)
  # crop and write precluster raster as .tif file
  preclust_filename <- file.path(taxonpath, paste0("preclusters31", ".tif"))
}

# save mid / pre clusters as sf object file to retrieve for post processing
# see: https://r-spatial.github.io/sf/articles/sf2.html
pclust_info31 |> sf::write_sf(file.path(taxonpath, paste0(gsub(" ","_",
          taxonA$ala_search_term), "_preclusters_prelim", ".shp")))
# throws a warning: Field names abbreviated for ESRI Shapefile driver

### Make a mask layer???
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


precluster_cellcount31 <- sum(terra::freq(precluster_rast31))
grouped_cellcount31 <- max(precluster_cellcount31, midcluster_cellcount31)
orphan_cellcount31 <- sum(terra::freq(orphan_rast31)) # sum makes it numeric
# this is returned as cell_counts in try_taxon_observations()
output31 <- (c(grouped_cellcount31, orphan_cellcount31))

# taxon$filter_category <- taxon$filter_category ## redundant?
taxonA$num_preclusters <- max(pclust_info31$midcluster)
taxonA$num_orphans <- sum(test79$precluster == 0)
taxonA$precluster_cellcount <- output31[1]
taxonA$orphan_cellcount <- output31[2]

# check that added columns contain correct (plausible) data..
taxonA[c(1,2,94:102)]

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

