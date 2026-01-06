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
  mutate(source = case_when(replicate==0 ~ "point",
                            replicate>0 ~ "polygon"))
## Group these data by species and replicates, and compute variances for each set.
elev_rangediff <- elev_df %>% 
  group_by(species, replicate, source) %>%
  summarise(elevational_range = max(elevation) - min(elevation),
              .groups = "drop")

elev_rangediff_ci <-  elev_rangediff %>%
  group_by(species) %>%
  summarise(
    obs_elev_range = elevational_range[replicate == 0],
    ci_low  = quantile(elevational_range[replicate > 0], 0.025),
    ci_high = quantile(elevational_range[replicate > 0], 0.975),
    within_95 = obs_elev_range >= ci_low & obs_elev_range <= ci_high,
      .groups = "drop")

### Test statistic using the t.test formula.
testing_elev_rangediff <- elev_rangediff %>%
  group_by(species) %>%
  summarise(
    ttest = list(
      t.test(
        x  = elevational_range[replicate > 0],  # random sample distribution
        mu = elevational_range[replicate == 0]  # observed range
      )
    ),
    t_stat  = round(ttest[[1]]$statistic, 5),
    p_value = round(ttest[[1]]$p.value, 10)
  ) %>% dplyr::select(-ttest)

## Run one-sample chi-squared variance test
chisq_stat <- elev_df_var %>%
  group_by(species) %>%
  summarise(
    chisq_test = list(
      EnvStats::varTest(
        x = sample_variances[replicate > 0],
        sigma.squared = sample_variances[replicate == 0]
      )
    ),
    chisq_val = chisq_test[[1]]$statistic[["Chi-Squared"]],
    p_val     = chisq_test[[1]]$p.value,
    .groups = "drop"
  ) %>% dplyr::select(-chisq_test)

## CI of the empirical distribution.
within_ci <- elev_df_var %>%
  group_by(species) %>%
  summarise(
    obs_var = sample_variances[replicate == 0],
    ci_low  = quantile(sample_variances[replicate > 0], 0.025),
    ci_high = quantile(sample_variances[replicate > 0], 0.975),
    within_95 = obs_var >= ci_low & obs_var <= ci_high, 
    .groups = "drop"
  )

## Plot histograms of the sample variance distribution along with the observed value and 95% CI.
ggplot(elev_rangediff, aes(x = elevational_range)) +
  geom_histogram(
    data = subset(elev_rangediff, replicate > 0),
    bins = 20,
    fill = "grey80",
    colour = "grey50"
  ) +
  # geom_rect(
  #   data = within_ci,
  #   inherit.aes = FALSE,
  #   aes(
  #     xmin = ci_low,
  #     xmax = ci_high,
  #     ymin = 0,
  #     ymax = Inf
  #   ),
  #   fill = "steelblue",
  #   alpha = 0.3
  # ) +
  geom_vline(
    data = subset(elev_rangediff, replicate == 0),
    aes(xintercept = elevational_range),
    colour = "red",
    linewidth = 0.8
  ) +
  facet_wrap(~ species, scales = "free_x") +
  theme_bw()
