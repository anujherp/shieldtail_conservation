## Script 4: This chunk uses range maps from IUCN, GARD, and Shieldtail Mapping Project and overlays them.

# List of shieldtail species of interest.
uro_list <- c("Melanophidium khairei","Platyplectrurus madurensis","Plectrurus perroteti",
              "Rhinophis karinthandani","Rhinophis melanoleucus","Teretrurus hewstoni",
              "Teretrurus sanguineus","Uropeltis ellioti","Uropeltis maculata",
              "Uropeltis myhendrae","Uropeltis pulneyensis","Uropeltis rubrolineata")

# Manually enter species name.
which_sp <- "Melanophidium khairei"

# Load the dataset
shieldtail_raw <- read_csv(paste0(wd$data,"shieldtail_records.csv"), show_col_types = FALSE)
shieldtail_long <- shieldtail_raw %>% 
  filter(scientific_name == which_sp) %>%
  drop_na(latitude) %>% distinct(latitude, .keep_all = T)

# Load the range maps (SMP, IUCN, and GARD respectively).
uro_rangemaps <- read_sf(paste0(wd$shapefiles,"uro_ranges.shp")) %>% filter(species == which_sp) %>% st_set_crs(4326)
uro_iucn <- read_sf(paste0(wd$shapefiles,"uro_iucn.shp")) %>% filter(species == which_sp)
uro_gard <- read_sf(paste0(wd$shapefiles,"uropeltid_gard.shp")) %>% filter(binomial == which_sp)

# Unioning all polygons per species to determine the overall geographic extent I will be dealing with.
poly <- st_union(uro_iucn, uro_gard) %>% st_convex_hull()
poly <- st_union(uro_rangemaps, poly, uro_gard) %>% st_convex_hull()

## For unasessed species (they don't have IUCN or GARD ranges)
# poly <- uro_rangemaps

# Calculate the bounding box with a buffer for aesthetic padding
bbox <- st_bbox(poly)

# Define a small padding factor to add some space around the polygons
padding_factor_lon <- 0.05  # Adjust this as needed for more or less space
padding_factor_lat <- 0.05 
# Adjust the bounding box with the padding factor
bbox[1] <- bbox[1] - (bbox[3] - bbox[1]) * padding_factor_lon  # xmin
bbox[3] <- bbox[3] + (bbox[3] - bbox[1]) * padding_factor_lon  # xmax
bbox[2] <- bbox[2] - (bbox[4] - bbox[2]) * padding_factor_lat  # ymin
bbox[4] <- bbox[4] + (bbox[4] - bbox[2]) * padding_factor_lat  # ymax

# Set the longitude and latitude extents
longs <- c(bbox[1], bbox[3])
long_breaks <- round(seq(bbox[1],bbox[3], length.out = 5),1)
lats <- c(bbox[2], bbox[4])
lats_breaks <- round(seq(bbox[2], bbox[4], length.out = 5),1)

# Plot occurrences and convex hull
map_fig <- ggplot() +
  # Background: Shoreline and India
  geom_sf(data = shoreline, colour = "#3b3b3b", fill = "#F5F5F5", alpha = 0.95, lwd = 1) +
  geom_sf(data = sf_india, colour = "#3b3b3b", fill = "#F5F5F5", alpha = 1, lwd = 2) +
  
  # Species range layers
  geom_sf(data = uro_gard, fill = "#FFB31A", colour = "black", alpha = 0.65, lwd = 0) +
  geom_sf(data = uro_iucn, fill = "#de1a24", colour = "black", alpha = 0.65, lwd = 0) +
  
  geom_sf(data = uro_rangemaps, fill = "#3f8f29", colour = "black", alpha = 0.65, lwd = 0) +
  
  # Species points
  geom_point(
    data = shieldtail_long,
    mapping = aes(x = longitude, y = latitude),
    fill = "#056517",
    size = 10,
    pch = 21
  ) +
  
  # Coordinate and scales
  coord_sf(
    xlim = c(longs[1] - 0.1, longs[2] + 0.1),
    ylim = c(lats[1] - 0.1, lats[2] + 0.1),
    expand = TRUE
  ) +
  scale_x_continuous(breaks = long_breaks) +
  scale_y_continuous(breaks = lats_breaks) +
  
  # Themes
  theme_bw(base_size = 70) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.text.x = element_text(size = 40),
    axis.ticks.length = unit(0.3, "cm"),
    panel.background = element_rect(fill = "lightblue"),
    plot.margin = margin(15, 15, 15, 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2) # <- updated size to linewidth
  )

## Save the plot as a JPEG file.
ggsave(map_fig,filename = paste0(wd$output,paste0("IUCN_GARD&SMP comparison/", 
       substr(shieldtail_long$scientific_name[1], 1, 1),"_", word(shieldtail_long$scientific_name[1],
       start = 2,sep=" ")),"_w_IUCN&GARD.png"), width = 800, height = 800, units = "mm",
       device='png', dpi=400)
