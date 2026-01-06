## Script 5: Operations on point and polygonal elevation data.

## Load the data
# occurrence data
shieldtail_raw <- read_csv(paste0(wd$data,"shieldtail_records.csv"), show_col_types = FALSE)
# polygon shapefiles
uro_rangemaps <- read_sf(paste0(wd$shapefiles,"uro_ranges.shp")) %>% 
  mutate(species = paste0(substr(species, 1, 1),". ",
                   word(species, start = 2, sep=" ")))

# Extract elevation values from the cleaned occurrence point dataframe.
pts <- shieldtail_raw %>%
  drop_na(latitude) %>% # remove NA coordinates
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  mutate(elevation = round(
      raster::extract(elevation,.,),3)) %>%
  st_drop_geometry() %>% rename(species = scientific_name) %>%
  mutate(replicate = 0, 
         species = paste0(substr(species, 1, 1),". ",
                  word(species, start = 2, sep=" "))) %>%
  dplyr::select(species, replicate, elevation)

# Convert the polygon sf dataframe to a SpatialPolygonDataFrame
uro_sp <- as(uro_rangemaps, "Spatial")

# Extract elevation values from the polygonal data.
set.seed(130) #random seed
# empty lists for the loop.
rpts_thinned <- setNames(vector("list", length(uro_rangemaps$species)), uro_rangemaps$species)
rpts_draw    <- setNames(vector("list", length(uro_rangemaps$species)), uro_rangemaps$species)

# Loop that 1) crops and masks the elevation raster layer to range polygons.
# 2) Converts the masked raster grid to points.
# 3) Randomly draws a sample (based on occurrence data) for each species, 1000 times.
for (i in seq_along(uro_rangemaps$species)) {
  
  ## 1. Crop + mask raster to polygon
  poly_i <- uro_sp[i, ]
  r_i    <- raster::crop(elevation, poly_i)
  r_m    <- raster::mask(r_i, poly_i)
  
  ## 2. Raster → points → sf
  r_sf <- data.frame(st_as_sf(rasterToPoints(r_m, spatial = TRUE))) %>%
    rename(elevation = Ind_raster)
  
  rpts_thinned[[i]] <- r_sf # List with masked polygonal elevation points and geometry.
  
  ## 3. Random sampling (1000 samples, each size = species count)
  draw_n <- sum(pts$species == uro_rangemaps$species[i])
  
  rpts_draw[[i]] <- replicate( # List with sampled points and their geometry.
    1000,
    r_sf[sample(nrow(r_sf), draw_n, replace = T), ],
    simplify = FALSE
  )
  print(i)
}

# Convert the 1000 replicates of the random draw list to a dataframe with appropriate structure.
rpts_df <- rpts_draw %>% 
  map(~ bind_rows(.x, .id = "replicate")) %>%   # bind 1000 replicate data.frames per species
  bind_rows(.id = "species") %>%
  mutate(replicate = as.numeric(replicate)) %>%
  dplyr::select(-geometry)

## Merge the raster points and observed occurrence points dataframe, with the `source` column as an identifier.
elev_df <- bind_rows(rpts_df, pts) %>% 
  mutate(source = case_when(replicate<1 ~ "point",
                            replicate>0 ~ "polygon"))

# Function takes in a dataframe (eg., elev_df) and runs a Levene's test of homogeneity of variance between observed points and each of the 1000 replicate, for all 12 species.
# Returns a dataframe of mean Levene's test statistic for each species.
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
    
    lev_stats[[i]] <- tibble(replicate = i,
                            F_value = lev$`F value`[1],
                            p_value = lev$`Pr(>F)`[1])
  }
  
  lev_df <- bind_rows(lev_stats)
  
  # summary across 1000 replicates
  summary_row <- lev_df %>%
    summarise(mean_F = round(mean(F_value),3),
              mean_p = round(mean(p_value),3),
              prop_p_less_0.05 = round(mean(p_value < 0.05),3)) %>%
    mutate(species = unique(df_species$species))
  
  return(summary_row)
}

## Run the `levene_bootstrap_species()` function with the `elev_df` across all species.
levene_bootstrap_results <- elev_df %>%
  group_by(species) %>%
  group_modify(~ levene_bootstrap_species(.x)) %>%
  ungroup()

## Plot a box and whiskers plot, with mean Levene's test P values 
## (Figure S2)
elevation_box <- ggplot() +
  geom_boxplot(data= elev_df, aes(x = species, y = elevation, fill = source),
               outlier.size = 5, lwd = 2, coef = 2, show.legend = F) +
  labs(x = "Species",
       y = "Elevation [m]") +
  theme_bw(base_size = 70) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_text(data = levene_bootstrap_results, 
            aes(x = species, y = max(elev_df$elevation) + 10, 
                label = paste("p =", round(mean_p, 3))),
            color = "maroon", size = 15) +
  theme(axis.text.x = element_text(face = "italic"))

## Save the plot as a png file.
ggsave(elevation_box, filename = paste0(wd$output,"Elevation_Box_n.png"), width = 1200, 
       height = 800, units = "mm", device='png', dpi=500)

## Plot a violin plot (Figure 2a)
## Use the full raster points dataframe.
rpts_elevation <- bind_rows(rpts_thinned, .id = "species") %>%
  dplyr::select(-geometry) %>%
  mutate(species = fct_reorder(species, elevation,
                               .fun = min, .desc = F))

elev_df$species <- factor(elev_df$species, levels = levels(rpts_elevation$species))

# Plot a violin graph, with occurrence points overlayed.
violin <- ggplot() +
  geom_violin(data = rpts_elevation, aes(x = species, y= elevation), fill = alpha("#3f8f29", 0.5), trim = T, scale = "width", adjust = 0.3, show.legend = F, lwd = 2) +
  geom_point(data = subset(elev_df, source == "point"), aes(x = species, y = elevation), fill = "#056517", size = 11, pch = 21, show.legend = F) +
  scale_y_continuous(breaks=seq(0,2400,300), limits = c(0, 2530)) +
  labs(x = "Species", y = "Elevation [m]") +
  # coord_flip() +
  theme_bw(base_size = 70) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.75, vjust = 0.9, face = "italic")
  )
# save the plot as a jpeg file
ggsave(violin, filename = "elevation_plots_n.png",path = paste0(wd$output), width = 750, 
       height = 1000, units = "mm", device='png', dpi=400)

## creating Table (S5)
# This dataframe is for exporting the point values as a csv.
point_elev_csv <- pts %>%
  group_by(species) %>%
  summarise(min_elev = min(elevation),
            max_elev = max(elevation),
            med_elev = median(elevation),
            avg_elev = mean(elevation))

# This dataframe is for exporting the polygonal values as a csv.
pol_elev_csv <- rpts_elevation %>%
  group_by(species) %>%
  summarise(min_elev = min(elevation),
            max_elev = max(elevation),
            med_elev = median(elevation),
            avg_elev = mean(elevation))

## A dataframe containing both point and polygonal elevation values.
elev_csv <- left_join(point_elev_csv, pol_elev_csv, by = c("species"),
                      suffix = c("_point", "_polygon")) %>% mutate(across(2:9, ~round(.x,0)))
write_csv(elev_csv, paste0(wd$data, "elev_summary.csv"))
