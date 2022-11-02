# To get a list of all files in a given directory sent to a .csv file
## use ls -1 > file.csv

## 1565 taxa in the database have AoO models (2462 do not)

library(RColorBrewer)

# sf object has 13 columns including TAXON_ID, SCIENTIFIC, and geometry
# geometry is a single multipolygon
# xmin: 2170810 ymin: 2539457 xmax: 2444810 ymax: 2561457
# x-extent: 274000 y-extent: 22000

# st_buffer creates rounded corners for AoO regions by default
## buffering for AoO (Area of Occupancy) should be less than
## buffering for point observation data because AoO is expected to cover
## the limits of occupied territory, but individual observations are
## unlikely to made at the range limit, just somewhere inside it so epsilon
## should be generous to allow for dispersal to, and then beyond limit
## of actual AoO, not just dispersal beyond where individual was observed
### USING taxa$epsilon / 3 FOR AoO BUFFERING FOR NOW

# new function to rasterize midclusters
mid_shapes_to_raster <- function(shapes, taxon, mask_layer, taxonpath) {
  print(shapes)
  shapevect <- terra::vect(shapes)
  print(shapevect)
  if (length(shapevect) > 0) {
    obs_raster <- terra::rasterize(shapevect, mask_layer,
                            field = "midcluster")
  } else {
    mask_layer * 0
  }
}

## AoO function --------

AoO_FILE <- paste0("AoO_",gsub(" ","_", taxonA$ala_search_term), ".shp")
# AoO_FILE <- "12992.shp" # for now
AoO_PATH <- file.path(taxonpath, AoO_FILE)
AoO1 <- st_read(AoO_PATH)
AoO1 <- AoO1[,c(2,13)] # a single multipolygon
# latin <- AoO1[1,1] |> st_drop_geometry() |> as.character()

# separates single multipolygon into separate polygons, then buffers
## buffer distance here is < buffer for individual observations
AoO1 <- suppressWarnings(st_cast(AoO1, "POLYGON")) |> 
  st_buffer(dist = (taxonA$epsilon * 1000) / 3) |> rename(latin = SCIENTIFIC)
# plot(AoO1) # for just the buffered AoO regions

# number new preclusters, following on from preclusters derived from obs data
AoO1 <- AoO1 |> 
  add_column(precluster = 
      1+nrow(pclust_info31):(nrow(AoO1)+nrow(pclust_info31)-1))

AoO1 <- AoO1 |> add_column(pix_count = 0, recent_year = 2019, num_obs = 0)
# sanity check for valid geometries: st_is_valid(AoO1)
AoO1 <- rbind(pclust_info31, AoO1[,c(3,2,4,5,1,6)])
# plot(AoO1[1])
AoO1_parts <- st_cast(st_union(AoO1), "POLYGON") # list of 18?
# plot(AoO1_parts)
midcluster <- unlist(st_intersects(AoO1, AoO1_parts)) # a vector of clusters
mclust_info31 <- cbind(AoO1, midcluster) |> group_by(midcluster) |> 
  summarise(precluster = paste(precluster, collapse = ", "),
        pix_count = max(pix_count, na.rm = TRUE),
        recent_year = max(recent_year, na.rm = TRUE),
        latin = max(latin, na.rm = TRUE),
        num_obs = max(num_obs, na.rm = TRUE)) # 7 columns
# plot(mclust_info31[1], col = rainbow(nrow(mclust_info31)))
mclust_retain <- c("midcluster","geometry")

midcluster_obs31 <- mclust_info31[, mclust_retain]

## NEED TO USE A DIFFERENT FUNCTION HERE BECAUSE 'shapes_to_raster()'
##  (line 196 of obs4.R) uses field = "precluster"
midcluster_rast31 <- mid_shapes_to_raster(midcluster_obs31, taxonA,
        mask_layer, taxonpath)
# plot(midcluster_rast31)
# recalculate 'pix_count'
pixel_mid_freq31 <- terra::freq(midcluster_rast31)
mclust_info31$pix_count <- pixel_mid_freq31$count

# add column for cluster number
mclust_info31 <- mclust_info31 |> 
  add_column(cluster = mclust_info31$precluster) # 8 columns

### NEED TO RETURN 'mclust_info31'
