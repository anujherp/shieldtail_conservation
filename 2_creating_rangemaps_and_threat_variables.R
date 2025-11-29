## Script 2: This chunk creates minimum convex polygons for the below listed 
# 12 shieldtail species and also calculates a number of variables.

# List of shieldtail species of interest.
uro_list <- c("Melanophidium khairei","Platyplectrurus madurensis","Plectrurus perroteti",
              "Rhinophis karinthandani","Rhinophis melanoleucus","Teretrurus hewstoni",
              "Teretrurus sanguineus", "Uropeltis ellioti", "Uropeltis maculata","Uropeltis myhendrae",
              "Uropeltis pulneyensis","Uropeltis rubrolineata")

# This loop runs twelve times, drawing a polygon and computing variables for each species.
uro_plots <- list()
uro_ranges <- list()
i <- 0
for (i in 10:length(uro_list)) {
  # Here, I fetch the occurrence data sheet stored in the dataframes folder.
  shieldtail_raw <- read_csv(paste0(wd$data,"shieldtail_records.csv")) %>% filter(scientific_name == uro_list[i])
  
  ## I create a new variable to store dataset without any NA values, and unique coordinates.
  shieldtail_long <- shieldtail_raw %>%
    drop_na(latitude) %>% distinct(latitude, .keep_all = T)
  
  ## To improve the speed of each loop iteration, I crop the raster layers to match the extents of my point data, thus working with smaller-sized files.
  box <- c(xmin = min(shieldtail_long$longitude) - 0.1, xmax = max(shieldtail_long$longitude) + 0.1, ymin = min(shieldtail_long$latitude) - 0.1, ymax = max(shieldtail_long$latitude) + 0.1)

  built_area <- crop(built, extent(box))
  tree_cover1 <- crop(tree_cover, extent(box))
  tree_loss1 <- crop(tree_loss, extent(box))
  
  ## The built environment raster layer has categorical values, where cell value of 0 means not built-up land, 1 corresponds to built-up land as of 2000, and 2 is the expanded built-up land since then.
  # I reclassifiy this raster layer so I have two different layers for pre-exisiting built env. and expanded built env.
  # raster reclassification (categorical to binary)
  rcl_exp <- matrix(c(0, 0.8, 0,   # Values < 1 → 0
                      0.9, 1.1, 1,   # Value = 1 → 1
                      1.8, 2.2, 0), # Values > 1 → 0
                    ncol = 3, byrow = TRUE)
  # b_exp is the expanded built-up area.
  b_exp <- raster(classify(built_area, rcl_exp))
  
  rcl_stable <- matrix(c(0, 0.8, 0,   # Values < 1 → 0
                         0.9, 1.1, 0,   # Value = 1 → 1
                         1.8, 2.2, 1),  # Values 2 to <3 → 2
                       ncol = 3, byrow = TRUE)
  # b_stb is the pre-exisitng or stable built-up area.
  b_stb <- raster(classify(built_area, rcl_stable, include.lowest = T))
  
  # I will convert raster resolution from degrees to meters
  resolution_meters <- round(res(tree_cover1)[1] * (111320 * cos(mean(min(shieldtail_long$latitude), max(shieldtail_long$latitude)) * (pi / 180))),3)
  resolution_meters1 <- round(res(built_area)[1] * (111320 * cos(mean(min(shieldtail_long$latitude), max(shieldtail_long$latitude)) * (pi / 180))),3)
  
  ## I will change raster values to represent area tree cover/ loss
  b_exp <- (resolution_meters1^2) * (b_exp)
  b_stb <- (resolution_meters1^2) * (b_stb)
  tree_cover1 <- (resolution_meters^2) * (tree_cover1/100)
  tree_loss1 <- (resolution_meters^2) * (tree_loss1)
  
  ## We have taken a decision to split Uropeltis ellioti's range based on a hclustering threshold, this loop does that.
  # Read the methods section for more information.
  if (shieldtail_long$scientific_name[1] == "Uropeltis ellioti") {
    uro_coord <- data.frame("name" = shieldtail_long$scientific_name, "longitude" = shieldtail_long$longitude, "latitude" = shieldtail_long$latitude)
    # Convert data to a SpatialPointsDataFrame object
    xy <- SpatialPointsDataFrame(matrix(c(uro_coord$longitude, uro_coord$latitude), ncol = 2), data.frame(ID = seq(1:length(uro_coord$longitude))), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
    # Generate a geodesic distance matrix in meters
    mdist <- distm(xy)
    # Cluster all points using a hierarchical clustering approach
    hc <- hclust(as.dist(mdist), method = "complete")
    # Define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
    # h is the height of the third split, see the hclust object
    uro_coord$cluster <- cutree(hc, k = 3)
    ## Creating a convex hull polygon
    uro_coord <- uro_coord %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
      st_buffer(dist = 0.0097) #Arc degrees, creates a buffer of ~3.5 sq km, buffer radius of 1km.
    
    ## To polygonise and join points
    pol1 <- uro_coord %>% filter(cluster == 1) %>%
      group_by(cluster) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_convex_hull() %>% st_intersection(shoreline) %>% st_combine()
    
    pol2 <- uro_coord %>% filter(cluster == 3) %>%
      group_by(cluster) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_convex_hull() %>% st_intersection(shoreline) %>% st_combine()
    
    pol3 <- uro_coord %>% filter(cluster == 2) %>%
      group_by(cluster) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_convex_hull() %>% st_intersection(shoreline) %>% st_combine()
    
    pol <- data.frame(st_union(c(pol1, pol2, pol3)))
    rm(pol1,pol2,pol3, uro_coord, mdist, hc, xy)
  } else {
    ## Here I create a range polygon with small, smoothening buffer.
    pol <- shieldtail_long %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
      st_buffer(dist = 0.0097) %>% #Arc degrees, creates a buffer of ~3.5 sq km, radius 1km
      summarise(geometry = st_combine(geometry)) %>%
      st_convex_hull() %>% st_intersection(shoreline) %>%
      summarise(geometry = st_combine(geometry))
  }
  
  ## The below lines compute the prevalence of threat parameters and other variables within the ranges
  pol <- pol %>% mutate(species = shieldtail_long$scientific_name[i],
                        area_km = round(as.numeric(st_area(pol$geometry)/1000000),3),
                        count = nrow(shieldtail_raw),
                        total_rm = sum(!is.na(shieldtail_raw$`dead/dor`)),
                        dor_count = sum(shieldtail_raw$`dead/dor` == "dor", na.rm = T),
                        dor_prop = round(dor_count / total_rm, 3),           
                        PA_cov = ifelse(st_intersects(st_as_sf(pol),
                                                      protected_areas, sparse =F)[1,1] == T,
                                        as.numeric(round(st_area(st_intersection(st_as_sf(pol),
                                                                                 protected_areas)) / 1000000, 3)), as.numeric(0)),
                        NPA_covp = round(100-(PA_cov/area_km)*100,3),
                        WG_cov = ifelse(st_intersects(st_as_sf(pol), wg, 
                                                      sparse =F)[1,1] == T,
                                        as.numeric(round(st_area(st_intersection(st_as_sf(pol),
                                                                                 wg)) / 1000000, 3)), as.numeric(0)),
                        WG_covp = round((WG_cov/area_km)*100,3),
                        b_stb_cov = round(exact_extract(b_stb, st_as_sf(pol), fun=                                                     c('sum'))/1000000,3),
                        b_exp_cov = round(exact_extract(b_exp, st_as_sf(pol), fun=                                                       c('sum'))/1000000,3),
                        b_exp_covp = round((b_exp_cov/b_stb_cov)*100,3),
                        tree_cov = round(exact_extract(tree_cover1, st_as_sf(pol), fun=                                                     c('sum'))/1000000,3),
                        tree_covp = round((tree_cov/area_km)*100,3),
                        lost_cov = round(exact_extract(tree_loss1, st_as_sf(pol), fun=                                                       c('sum'))/1000000,3),
                        lost_covp = round((lost_cov/tree_cov)*100,3),
  ) %>%
    dplyr::select(species, area_km, count, total_rm, dor_count, dor_prop, PA_cov, NPA_covp, WG_cov, WG_covp, b_stb_cov, b_exp_cov, b_exp_covp, tree_cov, tree_covp, lost_cov, lost_covp, geometry)
  
  uro_ranges[[i]] <- pol
  print(i)
}
# loop ends.

## Binding all list elements into a dataframe.
# Below dataframe has range geometry information and variable values
uro_ranges <- bind_rows(uro_ranges)
## New sf dataframe containing range geometry information.
uro_rangemaps <- uro_ranges %>% arrange(species) %>%
  dplyr::select(species, area_km, count)
