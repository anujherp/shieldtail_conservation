## Script 1: This script installs and loads the packages needed to run all subsequent scripts.
## Please change the paths to make them match your machine.
## Download all the files required from the following URL:

# Load and install required packages
needed_packages <- c("tidyverse","sf", "raster","exactextractr", "terra", "rgdal", "geosphere",
                     "rgeos", "ggrepel", "car", "ggextra","scales","lwgeom", "ggpattern", "rlang",
                     "purrr", "harmonicmeanp")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(list=ls()); gc()

# Download the datasets following instructions given on (ZENODO FILE)
# Make sure the parent folder is called 'shieldtail_conservation'.
wd <- list()
wd$data <- "C:/Users/Anuj Shinde/My Drive/Papers/1-work-in-progress/shieldtail distributions/Animal Conservation/shieldtail_conservation/dataframes/"
wd$shapefiles <- "C:/Users/Anuj Shinde/My Drive/Papers/1-work-in-progress/shieldtail distributions/Animal Conservation/shieldtail_conservation/shapefiles/"
wd$rasters <- "C:/Users/Anuj Shinde/My Drive/Papers/1-work-in-progress/shieldtail distributions/Animal Conservation/shieldtail_conservation/rasters/"
wd$output <- "C:/Users/Anuj Shinde/My Drive/Papers/1-work-in-progress/shieldtail distributions/Animal Conservation/shieldtail_conservation/output/"

## Turn of 'sf' package's spherical geometry feature
# This will assume planar geometry when running vector map functions.
sf::sf_use_s2(FALSE)

### Load Required Datafiles.
## Indian map with political boundaries
sf_india <- read_sf(paste0(wd$shapefiles, "India_State_Boundary.shp")) %>%
  st_transform(4326)

## Range polygon for the western ghats mountain range.
wg <- read_sf(paste0(wd$shapefiles, "wg_border_region11.shp")) %>%
  st_set_crs(4326)

## Global shoreline polygon (quite finescale)
shoreline <- read_sf(paste0(wd$shapefiles, "shores_penindia.shp")) %>%
  st_crop(c(xmin = 68.1167 - 3, xmax = 97.4167 + 3, ymin = 8.0667 - 3, ymax = 37.1 + 3), shoreline) %>% st_set_crs(4326)

## Geographic dataset of India's Protected Areas (as of 2016).
protected_areas <- read_sf(paste0(wd$shapefiles, "pa_india.shp")) %>%
  st_union() %>% st_set_crs(4326)

## Load raster layers (India DEM, built environment, Tree cover & loss).
elevation <- raster(paste0(wd$rasters, "Ind_raster.tiff"))
built <- terra::rast(paste0(wd$rasters, "Built_area_and_change.tif"))
tree_cover <- raster(paste0(wd$rasters, "pen_india_treecover.tif"))
tree_loss <- raster(paste0(wd$rasters, "treecover_loss_bin.tif"))