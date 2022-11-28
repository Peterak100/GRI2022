# re-write of the observations.R script...

# Download observation data ------

## for testing a single taxon see test3.R

# subfunction that downloads observation records from ALA
## and saves to the variable 'obs'
download_observations <- function(taxon) {
  obs_csv_path <- suppressWarnings(file.path(taxon_path(taxon, taxapath),
            "observations.csv"))
  cat(" Retrieving observations from ALA for", taxon$ala_search_term,
      "...\n")
  obs <- galah_call() |> galah_identify(taxon$ala_search_term) |>
    galah_filter(year >= as.character(TIME_START),
            year <= as.character(TIME_END),
            basisOfRecord == BASIS3, stateProvince == STATE,
            profile = "ALA") |>
    galah_select(group = "basic","year") |>
    atlas_occurrences()
  cat("  Observations retrieved successfully\n")
  write_csv(obs, obs_csv_path)
  return(obs)
}

# subfunction that retrieves taxon observations already saved locally
read_cached_observations <- function(taxon, obs_csv_path) {
  obs_csv_path <- suppressWarnings(file.path(taxon_path(taxon, taxapath),
            "observations.csv"))
  cat("  Loading cached observations for ", taxon$ala_search_term, "\n")
  obs <- read.csv(obs_csv_path, header = TRUE)
  return(obs)
}

#  combined function: either reads previously saved records, or
#  downloads records from ALA and saves them in relevant folder
##  overwrites existing data
load_or_download_obs <- function(taxon, taxapath, force_download=FALSE) {
  obs_csv_path <- suppressWarnings(file.path(taxon_path(taxon, taxapath),
                  "observations.csv"))
  if (!force_download && file_exists(obs_csv_path)) {
    obs <- read_cached_observations(taxon, obs_csv_path)
    return(obs)
  } else {
    obs <- download_observations(taxon)
    cat(" Writing", obs_csv_path, "\n")
    write_csv(obs, obs_csv_path)
    return(obs)
  }
}

# Clean observation data ------

## subfunctions:
# Remove coordinates with NA values
remove_missing_coords <- function(obs) {
  drop_na(obs, any_of(c("decimalLatitude", "decimalLongitude")))
}

# Remove duplicate observations from the same location
## sorts by eventDate first so that the newest record is retained
remove_location_duplicates <- function(obs) {
  obs |> arrange(desc(eventDate)) |>
    distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE)
}
## TO DO: Include 'year' in distinct() above to allow for same location
##  but different year? Will not change size of clusters, but
##  will mean more records per cluster

## TO DO: further filter for nonsensical locations using galah_geolocate()
##  with a predefined polygon that covers geographic limits of Australia??

# Filter observations data
## also include 'maybe_remove_subspecies(taxon)' # not needed?
## also include 'filter_by_fire_severity(taxon)' # from fire_severity.R
filter_observations <- function(obs, taxon) {
  obs  |> 
    remove_missing_coords() |>
    remove_location_duplicates()
}

# Cluster observation data ------

## three subfunctions.. 
# converts to an 'sf' object with a geometry column
# Add transformed coordinates "x" and "y" for accurate distance calculations
## removes decimalLongitude and decimalLatitude columns on the left (default)
## and converts them to a geometry column on the right
add_euclidean_coords <- function(obs) {
  sf::st_as_sf(obs, coords = c("decimalLongitude", "decimalLatitude"),
            crs = LATLON_EPSG)  |> 
    sf::st_transform(crs = METRIC_EPSG)
}

## had to split previous subfunction into two parts
##  converts geometry column to "x" and "y" columns on the right
## then drops the geometry column (so no longer an 'sf' object)
mut_euclidean_coords <- function(obs){
  mutate(obs, x = sf::st_coordinates(obs)[,1],
        y = sf::st_coordinates(obs)[,2]) |>
    sf::st_drop_geometry()
}

# adds a new column for 'precluster' number
## try changing the epsilon value to test different scenarios
## add 'method = "hybrid" to arguments of this function?
scan_clusters <- function(obs, eps) {
  preclusters <- fpc::dbscan(obs[, c("x", "y")],
      eps = eps * 1000,
      MinPts = 3)
  mutate(obs, precluster = preclusters$cluster)
}

## put 'eps = taxon$epsilon * 1000' as 3rd line for testing
# scan_clusters2 <- function(obs, eps) {
#   preclusters <- fpc::dbscan(obs[, c("x", "y")],
#     eps = taxon$epsilon * 1000,
#     MinPts = 3)
# mutate(obs, precluster = preclusters$cluster)
# }

### adjust epsilon to test relative effect of dispersal
## taxon$epsilon <- taxon$epsilon*2

# Categorise preclusters and add preclusters column to dataframe
precluster_observations <- function(obs, taxon) {
  obs |>
    add_euclidean_coords() |>
    mut_euclidean_coords() |> 
    scan_clusters(taxon$epsilon)
}

# Load and filter observations --------

### include 'force_download' as 3rd argument to this function?
###  also include 'force_download as 3rd argument of 2nd line?
load_and_filter <- function(taxon, taxapath) {
  load_or_download_obs(taxon, taxapath) |> 
    filter_observations(taxon) |> 
    precluster_observations(taxon)
}

# Buffer preclusters & Write observation data ------

# buffers shapes by scaled_eps (write_precluster converts obs to shapes)
# allows for individual males & individual females to 'meet-in-the-middle'
buffer_obs <- function(shapes, scaled_eps) {
  shapes |> sf::st_buffer(dist = scaled_eps) |>  
    dplyr::group_by(precluster) |>  
    dplyr::summarise(precluster = unique(precluster))
}

# Add buffer around preclusters (dataframe with x and y)
buffer_preclustered <- function(shapes, scaled_eps) {
  dplyr::filter(shapes, precluster != 0) |>  
    buffer_obs(scaled_eps)
}

# Add buffer around orphans
## not sure why this is a function of 'shapes' (sf object with geometry)
buffer_orphans <- function(shapes, scaled_eps) {
  dplyr::filter(shapes, precluster == 0) |>  
    buffer_obs(scaled_eps)
}

# Convert points to a raster file matching mask_layer
## runs with previously buffered sf object
## this function errors if mask_layer is not reset during current session
##  because it is a non-exportable object?
## why print(shapes) & print(shapevect)? (lines 2 & 4 of this function)
pre_shapes_to_raster <- function(shapes, taxon, mask_layer, taxonpath) {
  # print(shapes)
  shapevect <- terra::vect(shapes)
  # print(shapevect)
  if (length(shapevect) > 0) {
    obs_raster <- terra::rasterize(shapevect, mask_layer,
                    field = "precluster")
  } else {
    mask_layer * 0
  }
}

# new function to rasterize midclusters
mid_shapes_to_raster <- function(shapes, taxon, mask_layer, taxonpath) {
  # print(shapes)
  shapevect <- terra::vect(shapes)
  # print(shapevect)
  if (length(shapevect) > 0) {
    obs_raster <- terra::rasterize(shapevect, mask_layer,
                    field = "midcluster")
  } else {
    mask_layer * 0
  }
}

# Trim removes all NA values, then pads by up to ten pixels each side
## but so that the padded version does not extend beyond original raster!
## TO DO: move standard pad_by value to config.toml
## TO DO: make overall pad_by value a function of epsilon for given taxon?
padded_trim <- function(gen_raster) {
  xresolution <- terra::xres(gen_raster)
  yresolution <- terra::yres(gen_raster)
  pad_by = 10
  trimmed <- terra::trim(gen_raster)
  if (terra::xmax(gen_raster) - (xresolution * pad_by) >
      terra::xmax(trimmed)){
    x1 =  terra::xmax(trimmed) + (xresolution * pad_by)
  } else {
    x1 = terra::xmax(trimmed) +
          abs((terra::xmax(gen_raster) - terra::xmax(trimmed)))
  }
  if (terra::xmin(gen_raster) + (xresolution * pad_by) <
      terra::xmin(trimmed)){
    x2 = terra::xmin(trimmed) - (xresolution * pad_by)
  } else {
    x2 = terra::xmin(trimmed) -
          abs((terra::xmin(trimmed) - terra::xmin(gen_raster)))
  }
  if (terra::ymax(gen_raster) - (yresolution * pad_by) >
      terra::ymax(trimmed)){
    y1 =  terra::ymax(trimmed) + (xresolution * pad_by)
  } else {
    y1 = terra::ymax(trimmed) +
      abs((terra::ymax(gen_raster) - terra::ymax(trimmed)))
  }
  if (terra::ymin(gen_raster) + (yresolution * pad_by) <
      terra::ymin(trimmed)){
    y2 = terra::ymin(trimmed) - (yresolution * pad_by)
  } else {
    y2 = terra::ymin(trimmed) -
      abs((terra::ymin(trimmed) - terra::ymin(gen_raster)))
  }
  padding <- ext(x2, x1, y2, y1)
  padded <- terra::extend(trimmed, padding)
  return(padded)
}

# Rename AoO files --------

## make a translation table in Excel and save as a .csv file
## needs old name ($1), new name ($2), and dummy info ($3) columns
## cd to relevant Ubuntu directory containing AoO files
## copy the .csv translation table there and then run this in terminal:
# awk -F',' 'system("mv " $1 " " $2)' convert1.csv
## -F',' indicates the translation file is comma separated
## having a 3rd column prevents $'\r' being added to new file names


# Write Precluster function --------

## Write the preclustered and orphan observations to raster files
write_precluster <- function(obs, taxon, mask_layer, taxapath) {
  suppressWarnings(taxonpath <- taxon_path(taxon, taxapath))
  shapes <- sf::st_as_sf(obs, coords = c("x", "y"), crs = METRIC_EPSG)
  scaled_eps <- taxon$epsilon * 1000 * EPSILON_SENSITIVITY_SCALAR
  preclustered_obs <- buffer_preclustered(shapes, scaled_eps)
  cat("Num preclusters for", taxon$ala_search_term,":",
      nrow(preclustered_obs), "\n")
  precluster_rast <- pre_shapes_to_raster(preclustered_obs, taxon,
        mask_layer, taxonpath)
  # save a plot of preclusters
  png(file = file.path(taxonpath, paste0(gsub(" ", "_", taxon$ala_search_term),
      "_preclusters", ".png")), width = 2160, height = 1440, pointsize = 30)
  plot(precluster_rast, main = paste(taxon$ala_search_term, "preclusters"))
  cat("precluster plot for", taxon$ala_search_term, "saved to:",
      taxonpath,"\n")
  while (!is.null(dev.list())) dev.off()
  # add columns for layer and value to right of shapes (NA for orphans)
  # then save as a csv (pixel count is repeated for each precluster)
  pixel_freq <- terra::freq(precluster_rast) |> 
    dplyr::rename(pix_count = count)
  # pixel_freq should have pixel frequencies for each precluster
  pclust_info <- left_join(preclustered_obs, pixel_freq[-1],
          by = c("precluster" = "value")) # don't need 'layer' # 3 columns
  # this messes up geometry, so drop it from this dataframe
  pclust_summary <- dplyr::filter(shapes, precluster != 0) |> 
    dplyr::group_by(precluster) |> 
    dplyr::summarise(recent_year = max(year, na.rm = TRUE),
          latin = max(scientificName)) |> st_drop_geometry()
  pclust_info <- left_join(pclust_info, pclust_summary, by = "precluster")
  # 5 columns
  pclust_num <- dplyr::filter(shapes, precluster !=0) |> 
    dplyr::count(precluster) |> sf::st_drop_geometry() |> 
    dplyr::rename(num_obs = n)
  pclust_info <- left_join(pclust_info, pclust_num, by = "precluster")
  # 6 columns
  midcluster_cellcount <- 0
  orphan_obs <- buffer_orphans(shapes, scaled_eps)
  orphan_rast <- pre_shapes_to_raster(orphan_obs, taxon, mask_layer,
          taxonpath)
  png(file = file.path(taxonpath, paste0(gsub(" ","_", taxon$ala_search_term),
          "_orphans", ".png")), width = 2160, height = 1440, pointsize = 30)
  plot(orphan_rast, main = paste(taxon$ala_search_term, "orphans"))
  while (!is.null(dev.list())) dev.off()
  dplyr::filter(shapes, precluster == 0) |>
    mut_euclidean_coords() |> write_csv(file.path(taxonpath,
          paste0(gsub(" ", "_", taxon$ala_search_term), "_orphans", ".csv")))
  # crop and write orphan rasters as .tif file
  orphan_filename <- file.path(taxonpath, paste0("orphans", ".tif"))
  crop_orph_rast <- orphan_rast |>  padded_trim()
  terra::crop(orphan_rast, crop_orph_rast,
          filename = orphan_filename, overwrite = TRUE)
  # use AoO (Area of Occupancy) file if there is one
  # the (7?) supporting files in addition to *.shp must also be present
  # this creates and returns a totally new version of pclust_info
  AoO_FILE <- paste0("AoO_", gsub(" ","_", taxon$ala_search_term), ".shp")
  AoO_PATH <- file.path(AoOpath, AoO_FILE)
  if (file.exists(AoO_PATH)){
    AoO1 <- sf::st_read(AoO_PATH, quiet = TRUE)
    AoO1 <- AoO1[,c(2,13)] # a single multipolygon
    # separates single multipolygon into separate polygons, then buffers
    ## buffer distance here is < buffer for individual observations
    AoO1 <- suppressWarnings(sf::st_cast(AoO1, "POLYGON")) |> 
      sf::st_buffer(dist = (taxon$epsilon * 1000) / 3) |> 
      dplyr::rename(latin = SCIENTIFIC)
    # number new preclusters
    # numbering to following on from preclusters derived from obs data
    AoO1 <- AoO1 |> 
      add_column(precluster = 
              1+nrow(pclust_info):(nrow(AoO1)+nrow(pclust_info)-1),
              pix_count = 0, recent_year = 2018, num_obs = 0)
    pclust_sort <- c("precluster","geometry","pix_count",
                     "recent_year","latin","num_obs")
    # why does position of 'geometry' change?
    AoO1 <- rbind(pclust_info, AoO1[,pclust_sort])
    AoO1_parts <- sf::st_cast(sf::st_union(AoO1), "POLYGON")
    midcluster <- unlist(sf::st_intersects(AoO1, AoO1_parts))
    mclust_info <- cbind(AoO1, midcluster) |> 
      dplyr::group_by(midcluster) |> 
      dplyr::summarise(precluster = paste(precluster, collapse = ", "),
          pix_count = max(pix_count, na.rm = TRUE),
          recent_year = max(recent_year, na.rm = TRUE),
          latin = max(latin, na.rm = TRUE),
          num_obs = max(num_obs, na.rm = TRUE))
    pclust_info <- mclust_info
    # pclust_info <- merge_AoO_regions(taxon)
    # recalculate 'pix_count' for midclusters
    mclust_retain <- c("midcluster","geometry")
    midcluster_obs <- pclust_info[, mclust_retain]
    # need to use a different function here because 'pre_shapes_to_raster()'
    #  uses field = "precluster"
    midcluster_rast <- mid_shapes_to_raster(midcluster_obs, taxon,
          mask_layer, taxonpath)
    # also need to save a .tif file for Circuitscape patches layer
    crop_mid_rast <- terra::merge(midcluster_rast, orphan_rast) |> 
      padded_trim()
    preclust_filename <- file.path(taxonpath, paste0("preclusters", ".tif"))
    # crop and write midcluster raster as .tif file
    terra::crop(midcluster_rast, crop_mid_rast,
          filename = preclust_filename, overwrite = TRUE)
    # recalculate 'pix_count'
    pixel_mid_freq <- terra::freq(midcluster_rast)
    pclust_info$pix_count <- pixel_mid_freq$count
    midcluster_cellcount <- sum(terra::freq(midcluster_rast))
  } else {
    mclust_order <- c("midcluster", "precluster", "pix_count",
                    "recent_year", "latin", "num_obs", "geometry")
    pclust_info <- pclust_info |> 
      add_column(midcluster = pclust_info$precluster)
    pclust_info <- pclust_info[, mclust_order]
    crop_rast <- terra::merge(precluster_rast, orphan_rast) |>  
      padded_trim()
    # crop and write precluster raster as .tif file
    preclust_filename <- file.path(taxonpath, paste0("preclusters", ".tif"))
  }
  # save mid/pre clusters as sf object file to retrieve for post processing
  # see: https://r-spatial.github.io/sf/articles/sf2.html
  pclust_info |>  suppressWarnings(sf::write_sf(file.path(taxonpath,
          paste0(gsub(" ","_", taxon$ala_search_term),
          "_preclusters_prelim", ".shp"))))
  # suppressed warning: Field names abbreviated for ESRI Shapefile driver
  precluster_cellcount <- sum(terra::freq(precluster_rast))
  grouped_cellcount <- max(precluster_cellcount, midcluster_cellcount)
  orphan_cellcount <- sum(terra::freq(orphan_rast)) # sum makes it numeric
  # this is returned as cell_counts in try_taxon_observations()
  output <- (c(grouped_cellcount, orphan_cellcount))
  taxon$num_preclusters <- max(pclust_info$midcluster) # overwritten later?
  taxon$num_orphans <- sum(obs$precluster == 0)
  taxon$precluster_cellcount <- output[1]
  taxon$orphan_cellcount <- output[2]
  return(c(grouped_cellcount, orphan_cellcount))
}

# Label taxa that don't need to be processed due to cluster numbers -------

# - That have too many orphan cells compared to precluster cells
label_high_orphan_area <- function(taxa) {
  id <- (taxa$orphan_cellcount/MAX_ORPHAN_PRECLUSTER_RATIO) > 
    taxa$precluster_cellcount
  taxa$filter_category[id] <- "high_ratio_orphan_cells"
  return(taxa)
}

# - That have too many preclusters
label_many_clusters <- function(taxa) {
  id <- taxa$num_preclusters >= MAX_CLUSTERS
  taxa$risk[id] <- "highly fragmented"
  taxa$filter_category[id] <- "many_clusters"  
  return(taxa)
}

# - That have too few preclusters
label_few_clusters <- function(taxa) {
  id <- taxa$num_preclusters <= MIN_CLUSTERS & taxa$num_preclusters > 0
  taxa$filter_category[id] <- "single_cluster"  
  return(taxa)
}

# - That have no preclusters
label_no_clusters <- function(taxa) {
  id <- taxa$num_preclusters == 0
  taxa$filter_category[id] <- "no_clusters"  
  return(taxa)
}

label_by_clusters <- function(taxa) {
  taxa |> 
    label_high_orphan_area() |> 
    label_many_clusters() |> 
    label_few_clusters() |> 
    label_no_clusters()
}

# Main Script to Process Observations --------
## Gets observation data, preclusters it into numbered groups,
##  and writes to both csv and raster files

process_observations <- function(taxa, mask_layer, taxapath) {
  start_time1 <- Sys.time()
  # Add new columns to taxa dataframe
  preclustered_taxa <- add_column(taxa, num_preclusters = 0,
              num_orphans = 0, precluster_cellcount = 0,
              orphan_cellcount = 0, error = NA) 
  num_preclusters <- 0
  # Loop over each taxon
  for (i in 1:nrow(preclustered_taxa)) {
    taxon <- preclustered_taxa[i, ] 
    cat("\nTaxon: ", taxon$ala_search_term, "\n")
    # Download, filter and precluster observation records
    obs <- load_and_filter(taxon, taxapath)
    
    # Create rasters with numbered preclustered observations
    # If there are any preclusters
    if (max(obs$precluster) != 0) {
      cell_counts <- write_precluster(obs, taxon, mask_layer, taxapath)
    } else {
      cell_counts <- c(0, 0)
    }
    error_string <- NA
    preclustered_taxa[i, ]$filter_category <- taxon$filter_category 
    preclustered_taxa[i, ]$num_preclusters <- max(obs$precluster)
    preclustered_taxa[i, ]$num_orphans <- sum(obs$precluster == 0)
    preclustered_taxa[i, ]$precluster_cellcount <- cell_counts[1]
    preclustered_taxa[i, ]$orphan_cellcount <- cell_counts[2]
    }
  end_time1 <- Sys.time()
  process_obs_time <- end_time1 - start_time1
  print(process_obs_time)
  return(label_by_clusters(preclustered_taxa))
}

## is this needed to download obs in chunks?
# retrieve_chunk_obs <- function(taxa) {
#  n1 <- nrow(taxa)
#  r1 <- rep(1:ceiling(n1/GALAH_MAXROWS), each = GALAH_MAXROWS)[1:n1]
#  s1 <- lapply(split(taxa, r1), retrieve_bulk_obs)
  # Split apply combine chunks
#  do.call(rbind, s1)
# }

retrieve_bulk_obs <- function(taxa) {
  start_time1 <- Sys.time()
  # Loop over each taxon
  for (i in 1:nrow(taxa)) {
    taxon <- taxa[i, ]
    # Download observation records and write to csv files
    obs <- download_observations(taxon)
  }
  end_time1 <- Sys.time()
  obs_retrieve_time <- end_time1 - start_time1
  print(obs_retrieve_time)
}
