# To get a list of all files in a given directory sent to a .csv file
## use ls -1 > file.csv
## 1565 taxa in the database have AoO models (2462 do not)
# library(RColorBrewer) # just for plotting colour options

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


# main AoO function
merge_AoO_regions <- function(taxon, taxonpath) {
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
  return(mclust_info31)
}


