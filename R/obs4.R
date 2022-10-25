# re-write of the observations.R script...

# Download observation data ------

## for testing a single taxon see test3.R

# subfunction that downloads observation records from ALA
## and saves to variable 'obs'
download_observations <- function(taxon) {
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
  return(obs)
}

# test71 <- download_observations(taxon)

# subfunction that retrieves taxon observations already saved locally
read_cached_observations <- function(taxon, obs_csv_path) {
  cat("  Loading cached observations for ", taxon$ala_search_term, "\n")
  obs <- read.csv(obs_csv_path, header = TRUE)
  return(obs)
}

### This works, but not needed for now? ##
#  combined function: either reads previously saved records, or
#  downloads records from ALA and saves them in relevant folder
#load_or_download_obs2 <- function(taxon, taxapath, force_download=FALSE) {
#  obs_csv_path <- file.path(taxon_path(taxon, taxapath), "observations.csv")
#  if (!force_download && file.exists(obs_csv_path)) {
#    obs <- read_cached_observations(taxon, obs_csv_path)
#    return(obs)
#  } else {
#    obs <- download_observations(taxon)
#    cat("  Writing ", obs_csv_path, "\n")
#    write_csv(obs, obs_csv_path)
#    return(obs)
#  }
#}
# test72 <- load_or_download_obs2(taxon, taxapath)

## omitting option to load from previously downloaded records..
##  overwrites existing data
##  is 'force_download=FALSE' needed?
load_or_download_obs <- function(taxon, taxapath, force_download=FALSE) {
  obs_csv_path <- suppressWarnings(file.path(taxon_path(taxon, taxapath),
                  "observations.csv"))
  obs <- download_observations(taxon)
  cat("  Writing", obs_csv_path, "\n")
  write_csv(obs, obs_csv_path)
  return(obs)
}
# test73 <- load_or_download_obs(taxonA, taxapath)

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
## INCLUDE 'year' in distinct() TO ALLOW FOR SAME PLACE BUT DIFFERENT TIME?
##  does not change size of clusters, but means more records per cluster

###  further filter for nonsensical locations using galah_geolocate()
###  with a predefined polygon that covers geographic limits of Australia??

# Filter observations data
## also include 'maybe_remove_subspecies(taxon)' # not needed?
## also include 'filter_by_fire_severity(taxon)' # from fire_severity.R
filter_observations <- function(obs, taxon) {
  obs  |> 
    remove_missing_coords() |>
    remove_location_duplicates()
}
# test74 <- filter_observations(test73, taxon)

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
# test75 <- add_euclidean_coords(test74) |> mut_euclidean_coords()

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

# test76 <- scan_clusters2(test75)
# plot(test76$x, test76$y, pch = 19, cex = 0.4,
#   xlim = c(2100000, 3000000), ylim = c(2250000, 2850000))
## title(main = paste("All observations for", taxon$ala_search_term))

## for no orphans (orphans are not coloured)
# plot(test76$x, test76$y, pch = 19, cex = 0.4, col = test76$precluster,
#   xlim = c(2100000, 3000000), ylim = c(2250000, 2850000))
## title(main = paste(taxon$ala_search_term, "- orphans not shown"))

### adjust epsilon to test relative effect of dispersal
## taxon$epsilon <- taxon$epsilon*2

# Categorise preclusters and add preclusters column to dataframe
precluster_observations <- function(obs, taxon) {
  obs |>
    add_euclidean_coords() |>
    mut_euclidean_coords() |> 
    scan_clusters(taxon$epsilon)
}
# test77 <- precluster_observations(test74, taxon)

# Load and filter observations --------

### include 'force_download' as 3rd argument to this function?
###  also include 'force_download as 3rd argument of 2nd line?
load_and_filter <- function(taxon, taxapath) {
  load_or_download_obs(taxon, taxapath) |> 
    filter_observations(taxon) |> 
    precluster_observations(taxon)
}
# test79 <- load_and_filter(taxon, taxapath) # should be same as test77

# Buffer preclusters & Write observation data ------

# buffers shapes by scaled_eps (write_precluster converts obs to shapes)
# allows for individual males & individual females to 'meet-in-the-middle'
buffer_obs <- function(shapes, scaled_eps) {
  shapes |> sf::st_buffer(dist = scaled_eps) |>  
    dplyr::group_by(precluster) |>  
    dplyr::summarise(precluster = unique(precluster))
}
# test31 <- buffer_obs(shapes, scaled_eps)

# Add buffer around preclusters (dataframe with x and y)
buffer_preclustered <- function(shapes, scaled_eps) {
  dplyr::filter(shapes, precluster != 0) |>  
    buffer_obs(scaled_eps)
}
# test31preclus <- buffer_preclustered(shapes, scaled_eps)

# Add buffer around orphans WHY FUNCTION OF SHAPES?? (sf object with geometry)
buffer_orphans <- function(shapes, scaled_eps) {
  dplyr::filter(shapes, precluster == 0) |>  
    buffer_obs(scaled_eps)
}
# test31orph <- buffer_orphans(shapes, scaled_eps)

# Convert points to raster file matching mask_layer
### runs with previously buffered sf object
### this function errors if mask_layer is not reset during current session
###  because it is a non-exportable object?
## why print(shapes) & print(shapevect)? (lines 2 & 4 of this function)
shapes_to_raster <- function(shapes, taxon, mask_layer, taxonpath) {
  print(shapes)
  shapevect <- terra::vect(shapes)
  print(shapevect)
  if (length(shapevect) > 0) {
    obs_raster <- terra::rasterize(shapevect, mask_layer,
                    field = "precluster")
  } else {
    mask_layer * 0
  }
}

# Trim removes all NA values, then pad by up to ten pixels each side
## but so that the padded version does not extend beyond original raster!
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

# pad3test <- padded_trim(crop_rast)

# Write Precluster function --------

## Write the preclustered and orphan observations to raster files
write_precluster <- function(obs, taxon, mask_layer, taxapath) {
  suppressWarnings(taxonpath <- taxon_path(taxon, taxapath))
  shapes <- sf::st_as_sf(obs, coords = c("x", "y"), crs = METRIC_EPSG)
  scaled_eps <- taxon$epsilon * 1000 * EPSILON_SENSITIVITY_SCALAR
  
  preclustered_obs <- buffer_preclustered(shapes, scaled_eps)
  cat("Num preclusters for", taxon$ala_search_term,":",
      nrow(preclustered_obs), "\n")
  precluster_rast <- shapes_to_raster(preclustered_obs, taxon, mask_layer,
                  taxonpath)
  # save a plot of preclusters
  png(file = file.path(taxonpath, paste0(gsub(" ", "_", taxon$ala_search_term),
      "_preclusters", ".png")), width = 2160, height = 1440, pointsize = 30)
  plot(precluster_rast, main = paste(taxon$ala_search_term, "preclusters"))
  cat("precluster plot for", taxon$ala_search_term, "saved to:",
      taxonpath,"\n")
  while (!is.null(dev.list())) dev.off()
  
  # add columns for layer and value to right of shapes (NA for orphans)
  # then save as a csv (pixel count is repeated for each precluster)
  pixel_freq <- terra::freq(precluster_rast)
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

## removed 'force_download = FALSE' and 'throw_errors = FALSE'
##  as 4th & 5th arguments to this function
process_observations <- function(taxa, mask_layer, taxapath) {
  # Add new columns to taxa dataframe
  preclustered_taxa <- add_column(taxa, num_preclusters = 0, num_orphans = 0, 
          precluster_cellcount = 0, orphan_cellcount = 0, error = NA) 
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
    #what does list do here?
    }
  return(label_by_clusters(preclustered_taxa))
}
