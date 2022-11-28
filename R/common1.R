## Install and load required packages --------
packages <- c("galah","fpc","lubridate","sf","fs","tidyverse",
              "RcppTOML","terra","igraph","ini")
# Install packages not yet installed
## Note: "sf" is difficult and probably requires gdal to be installed
## outside of R for MacOS and Linux
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

library(galah)
library(fpc) # has dbscan
library(lubridate)
library(sf)
library(fs)
library(tidyverse)
library(RcppTOML)
library(terra)
library(igraph)
library(ini)

## this needs to be set each session:
galah_config(email = "peterkriesner@yahoo.com", download_reason_id = 0,
             verbose = TRUE)

## Load values --------

# Set primary path for all input data depending on computing environment
Ubuntu_home <- "/home/peter"
MacOS_home <- "/Users/bioinformatics"
Q_home <- "Q:"
PKWin11_home <- "C:/Users/peter"

HOME_dir <- fs::path_home()

set_GRI <- function(HOME_dir) {
  if (HOME_dir == Ubuntu_home) {
  GRI_set <- "/home/peter/GRI2022"
  } else if (HOME_dir == Q_home) {
    GRI_set <- "Q:/peterk/GRI2022"
  } else if (HOME_dir == MacOS_home) {
    GRI_set <- "/Users/bioinformatics/peterk/GRI2022"
  } else if (HOME_dir == PKWin11_home) {
    GRI_set <- "C:/Users/peter/Documents/R/GRI2022" # not set up yet
  } else {
    cat("Starting directory was wrong!")
  }
  return(GRI_set)
}
GRI_dir <- set_GRI(HOME_dir)
setwd(paste0(GRI_dir,"/R"))

set_DATA <- function(HOME_dir) {
  if (HOME_dir == Ubuntu_home) {
    DATA_set <- "/home/peter/data"
  } else if (HOME_dir == Q_home) {
    DATA_set <- "Q:/peterk/data"
  } else if (HOME_dir == MacOS_home) {
    DATA_set <- "/Users/bioinformatics/peterk/data"
  } else if (HOME_dir == PKWin11_home) {
    DATA_set <- "C:/Users/peter/data"
  } else {
    cat("Starting directory was wrong!")
  }
  return(DATA_set)
}
DATA_dir <- set_DATA(HOME_dir)

# set data path depending on computing environment
datapath <- file.path(DATA_dir)

## Paths to particular files and folders
taxapath <- file.path(datapath, "taxa")
groupingspath <- file.path(datapath, "groupings")
AoOpath <- file.path(datapath, "AoOfiles")
SMPpath <- file.path(datapath, "SMPfiles")
HIMpath <- file.path(datapath, "HIMfiles")

# metric reference system epsg code
METRIC_EPSG <- 3111
# lat/lon reference system epsg code
LATLON_EPSG <- 4326

BATCH_TAXA_CSV <- "batch_taxa.csv"
BATCH_TAXA_CSV_PATH <- file.path(datapath, BATCH_TAXA_CSV)
GALAH_MAXROWS <- 700

# Load main taxa dataframe from csv ------
taxa <- read.csv(BATCH_TAXA_CSV_PATH, header = TRUE)
# head(taxa[,c(1:3,38)])

# Access config.toml file and load 13 config variables --------
CONFIG_FILE <- "config.toml"
CONFIG_PATH <- file.path(datapath, CONFIG_FILE)
list2env(RcppTOML::parseTOML(file.path(datapath, CONFIG_FILE)), globalenv())
## expanded basisOfRecord (cannot put >1 option in TOML file?):
## other options are:
## "MACHINE_OBSERVATION","OBSERVATION","UNKNOWN"
BASIS3 <- c(BASIS, "MATERIAL_SAMPLE", "PRESERVED_SPECIMEN", "OBSERVATION")

HABITAT_RASTER <- "habitat.tif"
HABITAT_RASTER_PATH <- file.path(datapath, HABITAT_RASTER)
FIRE_SEVERITY_RASTER <- "fire_severity.tif"
FIRE_SEVERITY_RASTER_PATH <- file.path(datapath, FIRE_SEVERITY_RASTER)

RESISTANCE_RASTER <- "resistance.tif"
# Plot rasters to test them...
# HABITAT_RASTER_PATH |> terra::rast() |> plot()
# FIRE_SEVERITY_RASTER_PATH |> terra::rast() |> plot()

## an adjustment for buffering of observation record points
##  slightly < 0.5 (i.e. half of epsilon) to avoid ambiguity
EPSILON_SENSITIVITY_SCALAR <- 0.498
# Define timespan
TIMESPAN <- c(TIME_START:TIME_END)
# In case downloads run out of time
options(timeout=500)

mask_layer <- terra::rast(HABITAT_RASTER_PATH) < 0
# new code (Sept 2022):
### fix the projection for habitat.tif ###
mask_layer <- terra::project(mask_layer, paste0("EPSG:", METRIC_EPSG))
# cat(terra::crs(mask_layer))
# terra::plot(mask_layer) # errors first time this runs??


# functions --------

# Get the directory path for files relating to a specific taxon
taxon_path <- function(taxon, taxapath) {
  # We use underscores in the directory name
  underscored <- gsub(" ", "_", taxon$ala_search_term)[[1]]
  taxonpath <- file.path(taxapath, underscored)
  # Create the directory if it doesn't exist yet
  dir.create(taxonpath, recursive = TRUE)
  return(taxonpath)
}

# Download from a URL if the file doesn't exist already
maybe_download <- function(url, path) {
  if (!file.exists(path)) {
    download.file(url, path)
  }
}

# Download from s3 bucket if it is passed as a command line argument
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  bucket_url <- paste0(args[1], "/")
  batch_taxa_url <- paste0(bucket_url, BATCH_TAXA_CSV)
  config_url <- paste0(bucket_url, CONFIG_FILE)
  habitat_raster_url <- paste0(bucket_url, HABITAT_RASTER)
  fire_severity_raster_url <- paste0(bucket_url, FIRE_SEVERITY_RASTER)
  
  # Download
  maybe_download(fire_severity_raster_url, FIRE_SEVERITY_RASTER_PATH)
  maybe_download(habitat_raster_url, HABITAT_RASTER_PATH)
  # If we are on aws batch, always download updated taxa and config
  if (Sys.getenv("AWS_BATCH_CE_NAME") != "") {
    download.file(batch_taxa_url, BATCH_TAXA_CSV_PATH)
    download.file(config_url, CONFIG_PATH)
  } else {
    maybe_download(batch_taxa_url, BATCH_TAXA_CSV_PATH)
    maybe_download(config_url, CONFIG_PATH)
  }
}
