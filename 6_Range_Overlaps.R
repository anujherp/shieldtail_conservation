## Script 6: Calculating range size and position overlaps b/w IUCN, SMP, and GARD.

# List of shieldtail species of interest.
uro_list <- c("Melanophidium khairei","Platyplectrurus madurensis","Plectrurus perroteti",
              "Rhinophis karinthandani","Rhinophis melanoleucus","Teretrurus hewstoni",
              "Teretrurus sanguineus","Uropeltis ellioti","Uropeltis maculata",
              "Uropeltis myhendrae","Uropeltis pulneyensis","Uropeltis rubrolineata")

# Load the range maps (SMP, IUCN, and GARD).
uro_iucn <- read_sf(paste0(wd$shapefiles,"uro_iucn.shp")) %>% filter(species %in% uro_list) %>% 
  mutate(species = paste0(substr(species, 1, 1),". ",word(species, start = 2, sep=" ")))
uro_gard <- read_sf(paste0(wd$shapefiles,"uropeltid_gard.shp")) %>% filter(binomial %in% uro_list)
uro_rangemaps <- read_sf(paste0(wd$shapefiles,"uro_ranges.shp")) %>% filter(species %in% uro_iucn$species) %>%
  arrange(species)

## Calculating the EOO of IUCN and GARD ranges.
iucn_eoo <- uro_iucn %>% st_as_sf(st_coordinates()) %>% st_convex_hull() %>%
  mutate(area = round(as.numeric(st_area(.)/1000000),3))
gard_eoo <- uro_gard %>% st_as_sf(st_coordinates()) %>% st_convex_hull() %>%
  mutate(area = round(as.numeric(st_area(.)/1000000),3))

### Calculating overlap (Equation 1)
iucnsmp_overlap = c(); eoo_iucnsmp_overlap = c(); gardsmp_range_overlap = c(); 
eoo_gardsmp_overlap = c(); gardiucn_range_overlap = c(); eoo_gardiucn_overlap = c()
i <-0
for (i in 1:length(uro_iucn$species)) {
  iucnsmp_overlap[i] <- round(st_area(st_intersection(uro_rangemaps$geometry[i], uro_iucn$geometry[i]))/ st_area(st_union(uro_rangemaps$geometry[i],uro_iucn$geometry[i]))*100,3)
  eoo_iucnsmp_overlap[i] <- round(st_area(st_intersection(uro_rangemaps$geometry[i], iucn_eoo$geometry[i]))/ st_area(st_union(uro_rangemaps$geometry[i],iucn_eoo$geometry[i]))*100,3)
  gardsmp_range_overlap[i] <- round(st_area(st_intersection(uro_rangemaps$geometry[i], uro_gard$geometry[i]))/ st_area(st_union(uro_rangemaps$geometry[i],uro_gard$geometry[i]))*100,3)
  eoo_gardsmp_overlap[i] <- round(st_area(st_intersection(uro_rangemaps$geometry[i], gard_eoo$geometry[i]))/ st_area(st_union(uro_rangemaps$geometry[i],gard_eoo$geometry[i]))*100,3)
  gardiucn_range_overlap[i] <- round(st_area(st_intersection(uro_iucn$geometry[i], uro_gard$geometry[i]))/ st_area(st_union(uro_iucn$geometry[i],uro_gard$geometry[i]))*100,3)
  eoo_gardiucn_overlap[i] <- round(st_area(st_intersection(iucn_eoo$geometry[i], gard_eoo$geometry[i]))/ st_area(st_union(iucn_eoo$geometry[i],gard_eoo$geometry[i]))*100,3)
}

## variable containing overlap values
overlap <- data.frame(species = uro_iucn$species, 
                      iucnsmp_overlap = iucnsmp_overlap, 
                      eoo_iucnsmp_overlap = eoo_iucnsmp_overlap,
                      gardsmp_range_overlap = gardsmp_range_overlap,
                      eoo_gardsmp_overlap = eoo_gardsmp_overlap,
                      gardiucn_range_overlap = gardiucn_range_overlap,
                      eoo_gardiucn_overlap = eoo_gardiucn_overlap)

# Saving the overlap values as a csv file.
write_csv(overlap, paste0(wd$data, "range_eoo_overlaps.csv"))