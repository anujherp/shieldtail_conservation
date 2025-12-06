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
# Extracting the elevation values of all point data
pts <- shieldtail_raw %>%
  drop_na(latitude) %>%
  filter(scientific_name %in% uro_list) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  mutate(point_elev = round(
      raster::extract(elevation,.,),3)) %>%
  st_drop_geometry()

# This dataframe is for exporting the point values as a csv.
point_elev_csv <- pts %>%
  group_by(scientific_name) %>%
  summarise(min_alt = min(point_elev),
            max_alt = max(point_elev),
            med_alt = median(point_elev),
            avg_alt = mean(point_elev))

#### POLYGON DATA
# Here, I am extracting elevation variables from the polygon.

# Convert polygon to Spatial for raster::mask/crop
uro_sp <- as(uro_rangemaps, "Spatial") #convert the sf df to a spatpoly

rpts_thinned <- setNames(vector("list", length(uro_list)), uro_list)
rpts_draw    <- setNames(vector("list", length(uro_list)), uro_list)

set.seed(130)

for (i in seq_along(uro_list)) {
  
  ## 1. Crop + mask raster to polygon
  poly_i <- uro_sp[i, ]
  r_i    <- raster::crop(elevation, poly_i)
  r_m    <- raster::mask(r_i, poly_i)
  
  ## 2. Raster → points → sf
  r_sf <- data.frame(st_as_sf(rasterToPoints(r_m, spatial = TRUE))) %>%
    rename(elevation = Ind_raster)
  
  rpts_thinned[[i]] <- r_sf
  
  ## 3. Random sampling (1000 samples, each size = species count)
  draw_n <- sum(pts$scientific_name == uro_list[i])
  
  rpts_draw[[i]] <- replicate(
    1000,
    r_sf[sample(nrow(r_sf), draw_n, replace = T), ],
    simplify = FALSE
  )
  print(i)
}

rpts_df <- rpts_draw %>% 
  map(~ bind_rows(.x, .id = "replicate")) %>%   # bind 1000 replicate data.frames per species
  bind_rows(.id = "species") %>%
  mutate(replicate = as.numeric(replicate)) %>%
  dplyr::select(-geometry)

pts <- pts %>% rename(elevation = point_elev, species = scientific_name) %>%
  mutate(replicate = 0) %>% dplyr::select(species, replicate, elevation)

elev_df <- bind_rows(rpts_df, pts) %>% mutate(source = case_when(replicate<1 ~ "point",
                                                                 replicate>0 ~ "polygon"))
## levene-test

levene_bootstrap_species <- function(df_species) {
  
  # separate point & polygon elevations
  obs <- df_species %>% filter(source == "point")
  poly <- df_species %>% filter(source == "polygon")
  
  # storage for 1000 tests
  lev_stats <- vector("list", 1000)
  
  for (i in 1:1000) {
    
    poly_i <- poly %>% filter(replicate == i)
    
    # combine obs + polygon replicate i 
    df_i <- bind_rows(obs, poly_i)
    
    # run Levene test
    lev <- leveneTest(elevation ~ source, data = df_i)
    
    lev_stats[[i]] <- tibble(
      replicate = i,
      F_value = lev$`F value`[1],
      df_num = lev$Df[1],
      df_den = lev$Df[2],
      p_value = lev$`Pr(>F)`[1]
    )
  }
  
  lev_df <- bind_rows(lev_stats)
  
  # summary across 1000 replicates
  summary_row <- lev_df %>%
    summarise(
      mean_F = round(mean(F_value),3),
      mean_p = round(mean(p_value),3),
      prop_p_less_0.05 = round(mean(p_value < 0.05),3)
    ) %>%
    mutate(
      species = unique(df_species$species))
  
  return(summary_row)
}

levene_bootstrap_results <- elev_df %>%
  group_by(species) %>%
  group_modify(~ levene_bootstrap_species(.x)) %>%
  ungroup()

## Plotting a box and whiskers plot (Figues S2)
elevation_box <- ggplot() +
  geom_boxplot(data= elev_df, aes(x = species, y = elevation, fill = source),
               outlier.size = 5, lwd = 2, coef = 2, show.legend = F) +
  labs(x = "Species",
       y = "Elevation [m]") +
  theme_bw(base_size = 60) +
  scale_x_discrete(labels = c("M. khairei","P. madurensis","P. perroteti",
                              "R. karinthandani","R. melanoleucus","T. hewstoni",
                              "T. sanguineus","U. ellioti","U. maculata",
                              "U. myhendrae","U. pulneyensis","U. rubrolineata"),
    guide = guide_axis(n.dodge = 2)) +
  geom_text(data = levene_bootstrap_results, 
            aes(x = species, y = max(elev_df$elevation) + 10, 
                label = paste("p =", round(mean_p, 3))),
            color = "maroon", size = 15) +
  theme(axis.text.x = element_text(face = "italic"))

## Saving the plot as a jpeg file.
ggsave(elevation_box, filename = paste0(wd$output,"Elevation_Box_n.png"), width = 1200, 
       height = 800, units = "mm", device='png', dpi=500)

p <- list()
for (sp in unique(elev_df$species)) {
  
  df_sp <- elev_df %>% filter(species == sp)
  
  p[[sp]] <- ggplot(df_sp, aes(x = elevation)) +
    geom_histogram(binwidth = 50, fill = "skyblue", color = "black") +
    geom_vline(aes(xintercept = mean(elevation[source == "point"])), color = "red", size = 1.2) +
    theme_minimal() +
    labs(
      x = "Elevation",
      y = "Count",
      title = paste0(sp),
      subtitle = "Red line = observed mean elevation",
      size = 2)
  }
combined_plot <- wrap_plots(p, ncol = 3)   # like facet_wrap


# This dataframe is for exporting the polygonal values as a csv.
pol_elev_csv <- pol_elev_df %>%
  group_by(species) %>%
  summarise(min_alt = min(pol_elev),
            max_alt = max(pol_elev),
            med_alt = median(pol_elev),
            avg_alt = mean(pol_elev))

## A dataframe containing both point and polygonal elevation values.
elev_csv <- left_join(point_elev_csv, pol_elev_csv, by = c("species"))

## Plotting a violin plot (Figure 2a)

rpts_elevation <- bind_rows(rpts_thinned, .id="species") %>% 
  dplyr::select(-geometry)

ggplot() +
  geom_violin(data = rpts_elevation, aes(x = species, y= elevation), fill = alpha("#3f8f29", 0.6), trim = T, scale = "width", adjust = 0.3, show.legend = F, lwd = 2) +
  geom_point(data = subset(elev_df, source == "point"), aes(x = species, y = elevation), fill = "#056517", size = 9, pch = 21, show.legend = F) +
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