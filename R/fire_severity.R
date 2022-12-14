# Remove observations according to fire severity at the location
# obs are observations returned by galah,
# taxon is a single row of the taxa dataframe

# See https://www.environment.vic.gov.au/biodiversity/naturekit/nk-datalists
## no longer available online (March 2022)
## Two fire severity maps were available:
### the initial map layer had 5 severity classes and is a GeoTiff of ~9.6 Mb
### the final map layer had 6 severity classes and is a GeoTiff of ~348 Mb

# The Fire Severity GeoTiff layer shows the Fire severity map of the major fires
# in Gippsland and north east Victoria in 2019/20.

# The 6 fire severity classes are:
#  Class 6: Canopy fully burnt
#  Class 5: High canopy scorch (>80% of canopy foliage is scorched)
#  Class 4: Medium canopy scorch
#    (Canopy is a mosaic of both unburnt & scorched foliage, 20 - 80%)
#  Class 3: Low canopy scorch
#    (Canopy foliage largely unaffected (<20% scorched),
#      but understorey has been burnt)
#  Class 2: Unburnt (Canopy and understorey foliage are largely (>90%) unburnt)
#  Class 1: Non-woody vegetation (unclassified)
#  Class 0: No Data
#     (e.g. due to obscuration by cloud, cloud-shadow and/or smoke and haze)

## The 5 fire severity classes are:
#  5 = Canopy burnt
#  4 = High canopy scorch (>80% canopy scorch)
#  3 = Medium canopy scorch (20% - 80% canopy scorch)
#  2 = Low canopy scorch (<20% canopy scorch)
#  1 = Unburnt
#  0 = No data

filter_by_fire_severity <- function(obs, taxon) {
  taxon_fire_rating <- taxon$fire_rating
  # For NAs in taxon fire rating we give the taxon the maximum value
  taxon_fire_rating[is.na(taxon_fire_rating)] <- 6

  # read in raster layer data:
  severity_rast <- terra::rast(FIRE_SEVERITY_RASTER_PATH)

  # Make subsets of geographic points only for both lists:
  locs <- dplyr::select(obs, "decimalLongitude", "decimalLatitude")
  # Note: order of Longitude and Latitude is reversed here.

  # Create SpatialPoints objects from lists of geographic points,
  # telling R to use lat-long coordinates projection (4326):
  sf_locs <- sf::st_as_sf(locs, coords = c("decimalLongitude",
              "decimalLatitude"), crs = 4326)
  # transform SpatialPoints objects to same projection as raster layer:
  points <- sf::st_transform(sf_locs, crs = terra::crs(severity_rast))

  # Extract raster layer values to new lists:
  fire_severity_points <- raster::extract(severity_rast, terra::vect(points))

  # Replace 'NA's with zeros:
  # Needed because extracted lists will have NAs for any locations
  # that fall outside of the raster layer boundaries)
  fire_severity_points[is.na(fire_severity_points)] <- 0

  # Combine extracted fire severity lists with initial data frames:
  severity_df <- add_column(obs,
            fire_severity = fire_severity_points$fire_severity)

  # Subset new data frames by likely survivors of a given fire intensity:
  survivors_df <- dplyr::filter(severity_df, fire_severity < taxon_fire_rating)
  return(survivors_df)
}
