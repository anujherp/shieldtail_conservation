## Script 1: This script contains a list of required R packages and 
## paths that I will be using throughout the workflow.

# Load and install required packages
needed_packages <- c("tidyverse","sf", "raster","exactextractr", "terra", "rgdal", "geosphere",
                     "rgeos", "ggrepel", "car", "ggextra", "red","scales","lwgeom", "ggpattern", "rlang")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(list=ls()); gc()

# Download the datasets following instructions given on (ZENODO FILE)
# make sure your parent folder is called 'shieldtail_conservation'.
wd <- list()
wd$data <- ".../shieldtail_conservation/dataframes/"
wd$shapefiles <- ".../shieldtail_conservation/shapefiles/"
wd$rasters <- ".../shieldtail_conservation/rasters/"
wd$output <- ".../shieldtail_conservation/output/"

## Turning of 'sf' package's spherical geometry feature
# doing this will assume planar geometry when running vector map functions.
sf::sf_use_s2(FALSE)

### Loading Required Datafiles.
## Indian map with political boundaries
sf_india <- read_sf(paste0(wd$shapefiles, "India_State_Boundary.shp")) %>%
  st_transform(4326)

## Range polygon for the western ghats mountain range
wg <- read_sf(paste0(wd$shapefiles, "wg_border_region11.shp")) %>%
  st_set_crs(4326)

## Global shoreline polygon (quite finescale)
shoreline <- read_sf(paste0(wd$shapefiles, "shores_penindia.shp")) %>%
  st_crop(c(xmin = 68.1167 - 3, xmax = 97.4167 + 3, ymin = 8.0667 - 3, ymax = 37.1 + 3), shoreline) %>% st_set_crs(4326)

## Geographic dataset of India's Protected Areas (as of 2016)
protected_areas <- read_sf(paste0(wd$shapefiles, "pa_india.shp")) %>%
  st_union() %>% st_set_crs(4326)

## Load raster layers (India DEM, built environment, Tree cover & loss)
elevation <- raster(paste0(wd$rasters, "Ind_raster.tiff"))
built <- terra::rast(paste0(wd$rasters, "Built_area_and_change.tif"))
tree_cover <- raster(paste0(wd$rasters, "pen_india_treecover.tif"))
tree_loss <- raster(paste0(wd$rasters, "treecover_loss_bin.tif"))