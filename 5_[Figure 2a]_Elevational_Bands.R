## Script 5: This chunk here calculates and plots the elevational ranges of species polygons.
# It also compares variance the elevations ranges of point and polygonal data using the Levene's Test.
# Figure 2a and Figure S1.

# Load the dataset
shieldtail_raw <- read_csv(paste0(wd$data,"shieldtail_records.csv"), show_col_types = FALSE)

# Load the rangemaps (SMP)
uro_rangemaps <- read_sf(paste0(wd$shapefiles,"uro_ranges.shp"))

# List of shieldtail species of interest.
uro_list <- c("Melanophidium khairei","Platyplectrurus madurensis","Plectrurus perroteti",
              "Rhinophis karinthandani","Rhinophis melanoleucus","Teretrurus hewstoni",
              "Teretrurus sanguineus","Uropeltis ellioti","Uropeltis maculata",
              "Uropeltis myhendrae","Uropeltis pulneyensis","Uropeltis rubrolineata")

#### POINT DATA
# Here I am isolating geographic point information from the dataframe
point_elev_df <- shieldtail_raw %>% drop_na(latitude) %>%
  filter(scientific_name %in% uro_list) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  dplyr::select(scientific_name)

# Extracting the elevation values of all point data
point_elev_df <- point_elev_df %>%
  mutate(point_elev = round(exactextractr::exact_extract(elevation, st_as_sf(point_elev_df),fun= c('mean')),3),
         species = paste0(substr(scientific_name, 1, 1),". ", word(scientific_name, 
         start = 2, sep=" "))) %>% st_drop_geometry() %>% dplyr::select(-scientific_name)

# This dataframe is for exporting the point values as a csv.
point_elev_csv <- point_elev_df %>%
  group_by(species) %>%
  summarise(min_alt = min(point_elev),
            max_alt = max(point_elev),
            med_alt = median(point_elev),
            avg_alt = mean(point_elev))

#### POLYGON DATA
# Here, I am extracting elevation variables from the polygon.

# Convert polygon to Spatial for raster::mask/crop
poly_sp <- as(uro_rangemaps, "Spatial")
uro_points <- st_coordinates()
# Extract elevation at points
pt_elev <- raster::extract(elevation, as(point_elev_df, "Spatial"))
pt_var  <- var(pt_elev, na.rm = TRUE)
pt_var


uro_sp <- as(uro_rangemaps, "Spatial") #convert the sf df to a spatpoly
rpts <- list(); rpts_sf <- list()
for (i in seq_along(uro_list)) {
  poly_i <- uro_sp[i,] 
  r_i <- raster::crop(elevation, poly_i)
  rpts <- rasterToPoints(mask(r_i, poly_i), spatial = TRUE)
  rpts_sf[[i]] <- st_as_sf(rpts)
} ## crop and mask the elevation layer to each of the species polygons

# Thin to reduce spatial autocorrelation
set.seed(123)
thin_sf <- list()
df <- data.frame(st_coordinates(rpts_sf[[i]][["geometry"]]), SPEC ="Poly")
for (i in 1:seq_along(uro_list)) {
  thin_out <- thin(
    loc.data = df,
    lat.col = "Y",
    long.col = "X",
    spec.col = "SPEC",
    thin.par = 2000,     # 5 km thinning radius
    reps = 1,          # one thinned set is enough
    write.files = FALSE,
    verbose = FALSE
  )
  
  sightings_thinned = spThin::thin(df, 
                                   lat.col = "Y", 
                                   long.col = "X", 
                                   spec.col = "SPEC", 
                                   thin.par = 10, 
                                   reps = 100, 
                                   locs.thinned.list.return = TRUE, 
                                   write.files = F, 
                                   write.log.file = FALSE)[[1]]#try thinning at different scales
  plot(sightings_thinned)
  plot(df$X,df$Y)
}

# Extract coordinates of thinned points
thin_indices <- thin_sf[[1]]$SampleID
thin_pts <- rpts_sf[thin_indices,]

# Get de-autocorrelated raster values
thin_vals <- thin_pts$layer




## Create a list with the elevation values and corresponding coordinates for each polygon
moran_results <- list()
for (i in seq_along(pol_elev)) {
  
  r <- pol_elev[[i]]
  # get values and coordinates
  vals <- getValues(r)
  coords <- xyFromCell(r, 1:ncell(r))
  
  # remove NA cells
  ok <- !is.na(vals)
  vals <- vals[ok]
  coords <- coords[ok, ]
  
  # store for downstream Moran’s I
  moran_results[[i]] <- list(values = vals, coords = coords)
}

## Build spatial weights using k-nearest neighbour, and compute the Moran's I

for (i in seq_along(moran_results[1:3])) {
  vals  <- moran_results[[i]]$values
  coords <- moran_results[[i]]$coords
  
  # Build 8-nearest neighbor graph
  knn <- spdep::knearneigh(coords, k = 8)
  nb <- spdep::knn2nb(knn)
  lw <- spdep::nb2listw(nb, style = "W")
  
  # Moran's I test
  m_test <- spdep::moran.test(vals, lw)
  
  moran_results[[i]]$moran <- m_test
}









pol_elev_df <- extract(elevation, uro_rangemaps)
names(pol_elev_df) <- uro_list

pol_elev_df <- pol_elev_df %>% 
  unlist(use.names = T) %>%
  as.data.frame()
pol_elev_df$species <- rownames(pol_elev_df)
rownames(pol_elev_df) <- NULL
colnames(pol_elev_df) <- c("pol_elev", "species")

pol_elev_df <-  pol_elev_df %>% mutate(species = gsub('[[:digit:]]+', '',species)) %>%
  dplyr::select(species, pol_elev) %>% mutate(species = paste0(substr(species, 1, 1),". ", word(species, start = 2, sep=" ")))

pol_elev_df <- pol_elev_df %>% mutate(species = fct_reorder(species, pol_elev, .fun = min, .desc = F))

# This dataframe is for exporting the polygonal values as a csv.
pol_elev_csv <- pol_elev_df %>%
  group_by(species) %>%
  summarise(min_alt = min(pol_elev),
            max_alt = max(pol_elev),
            med_alt = median(pol_elev),
            avg_alt = mean(pol_elev))

## A dataframe containing both point and polygonal elevation values.
elev_csv <- left_join(point_elev_csv, pol_elev_csv, by = c("species"))


## Plotting a violine plot (Figure 2a)

ggplot() +
  geom_violin(data = pol_elev_df, aes(x = species, y= pol_elev), fill = alpha("#3f8f29", 0.6), trim = T, scale = "width", adjust = 0.3, show.legend = F, lwd = 2) +
  geom_point(data = point_elev_df, aes(x = species, y = point_elev), fill = "#056517", size = 9, pch = 21, show.legend = F) +
  scale_y_continuous(breaks=seq(0,2400,300), limits = c(0, 2570)) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "Species", y = "Elevation [m]") +
  coord_flip() +
  theme_linedraw(base_size = 15) +
  theme(axis.text.y = element_text(face = "italic", size = 64))

# saving the plot as a jpeg file
ggsave(filename = "elevation_plots_n.jpg",path = paste0(wd$output), width = 900, 
       height = 1200, units = "mm", device='jpeg', dpi=300)

#### Levene's Test for Homogeneity of Variance

for_levene <- bind_rows(
  point_elev_df %>% mutate(elevation = point_elev, source = 1) %>% dplyr::select(species, elevation, source),
  pol_elev_df %>% mutate(elevation = pol_elev, source = 2) %>% dplyr::select(species, elevation, source)
) %>%
  mutate(source = dplyr::recode(source, "1" = "points", "2" = "polygon"))

levene_p <- list()
q <- 0
for (q in 1:length(unique(for_levene$species))) {
  variance = leveneTest(elevation ~ factor(source), subset(for_levene, species == paste0(unique(for_levene$species)[q])))
  levene_p[[q]] = list(species = unique(for_levene$species)[q],
                       F_val = round(variance$`F value`,3),
                       P_val = round(variance$`Pr(>F)`,3))
  variance <- NULL
}
levene_p <- bind_rows(levene_p) %>% drop_na(P_val)

## Plotting a box and whiskers plot (Figues S2)
ggplot() +
  geom_boxplot(data= for_levene, aes(x = species, y = elevation, fill = source),
               outlier.size = 7) +
  labs(x = "Species",
       y = "Elevation [m]") +
  theme_bw(base_size = 70) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_text(data = levene_p, 
            aes(x = species, y = max(for_levene$elevation) + 10, 
                label = paste("p =", round(P_val, 4))),
            color = "maroon", size = 15) +
  theme(axis.text.x = element_text(face = "italic"))

## Saving the plot as a jpeg file.
ggsave(filename = paste0(wd$output,"Elevation_Box.jpg"), width = 1200, 
       height = 800, units = "mm", device='jpg', dpi=300)