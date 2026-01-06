## Script 6: Calculating range size and position overlaps b/w IUCN, SMP, and GARD.

## Load the data.
# List of shieldtail species of interest.
uro_list <- c("Melanophidium khairei","Platyplectrurus madurensis","Plectrurus perroteti",
              "Rhinophis karinthandani","Rhinophis melanoleucus","Teretrurus hewstoni",
              "Teretrurus sanguineus","Uropeltis ellioti","Uropeltis maculata",
              "Uropeltis myhendrae","Uropeltis pulneyensis","Uropeltis rubrolineata")

# Load the range maps (SMP, IUCN, and GARD).
uro_iucn <- read_sf(paste0(wd$shapefiles,"uro_iucn.shp")) %>% filter(species %in% uro_list) %>% 
  mutate(species = paste0(substr(species, 1, 1),". ",word(species, start = 2, sep=" ")))
uro_gard <- read_sf(paste0(wd$shapefiles,"uropeltid_gard.shp")) %>% filter(binomial %in% uro_list)%>% 
  mutate(species = paste0(substr(binomial, 1, 1),". ",word(binomial, start = 2, sep=" "))) %>% dplyr::select(-binomial)
uro_rangemaps <- read_sf(paste0(wd$shapefiles,"uro_ranges.shp")) %>%
  arrange(species) %>% mutate(species = paste0(substr(species, 1, 1),". ",word(species, start = 2, sep=" "))) %>% 
  filter(species %in% uro_iucn$species)

## Calculate the EOO of IUCN and GARD ranges.
iucn_eoo <- uro_iucn %>% st_as_sf(st_coordinates()) %>% st_convex_hull() %>%
  mutate(area = round(as.numeric(st_area(.)/1000000),3))
gard_eoo <- uro_gard %>% st_as_sf(st_coordinates()) %>% st_convex_hull() %>%
  mutate(area = round(as.numeric(st_area(.)/1000000),3))
smp_eoo <- uro_rangemaps %>% st_as_sf(st_coordinates()) %>% st_convex_hull() %>%
  mutate(area = round(as.numeric(st_area(.)/1000000),3))

## Calculate Range overlap (Equation 1)
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

# Save the overlap values as a csv file.
write_csv(overlap, paste0(wd$data, "range_eoo_overlaps.csv"))

## Range Size Differences
## Calculate the area (called extent here) of IUCN range and the overlap between SMP & IUCN. 
# % difference = (IUCN EOO - SMP EOO)/ mean(IUCN EOO,SMP EOO) X 100

smp_area <- uro_rangemaps %>% 
  mutate(smp_area = round(as.numeric(st_area(.)/1000000),3))
iucn_area <- uro_iucn %>% 
  mutate(iucn_area = round(as.numeric(st_area(.)/1000000),3))
gard_area <- uro_gard %>% 
  mutate(gard_area = round(as.numeric(st_area(.)/1000000),3))

areas <- data.frame(species = paste0(substr(uro_rangemaps$species, 1, 1),". ", 
                                     word(uro_rangemaps$species, start = 2, sep=" ")),
                    iucn= iucn_area$iucn_area, smp = smp_area$smp_area,
                    gard = gard_area$gard_area)

## Calculate difference in areas
areas <- areas %>% rowwise() %>% 
  mutate(iucn_diff = ((iucn - smp)/mean(c(iucn,smp))) *100,
         gard_diff = ((gard - smp)/mean(c(gard,smp))) *100) %>% 
  pivot_longer(cols = c(iucn_diff, gard_diff), names_to = "dataset", values_to = "difference") %>% 
  dplyr::select(-c(gard,iucn,smp))

## Plot a bar graph (Fgure 2b)
range_diff_plot <- ggplot(data = areas, aes(x = species, y = difference, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.51), color = "#57595B",
  width = 0.50, show.legend = F, alpha = 0.75, lwd = 3) +
  scale_fill_manual(values = c("#FFB31A","#de1a24")) +
  geom_hline(yintercept = 0, colour = "black", linetype = 1, lwd = 3) +
  xlab("Species") + ylab("Range Size Difference") + 
  scale_y_continuous(breaks = seq(-200, 200, 100), limits = c(-220, 270), labels = c("-200%", "-100%", "0%", "100%", "200%")) +
  theme_bw(base_size = 70) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.75, vjust = 0.9, face = "italic"))
## Save the plot as a jpeg file.
ggsave(range_diff_plot, filename = "range_size_differences.png",path = paste0(wd$output), width = 700, 
       height = 800, units = "mm", device='png', dpi=400)


## Exporting Table (S7).
## Calculate differences in EOO areas
eoo_areas <- data.frame(species = paste0(substr(uro_rangemaps$species, 1, 1),". ", 
                                         word(uro_rangemaps$species, start = 2, sep=" ")),
                        iucn_eoo= iucn_eoo$area, smp_eoo = smp_eoo$area,
                        gard_eoo = gard_eoo$area, 
                        iucn= iucn_area$iucn_area, smp = smp_area$smp_area,
                        gard = gard_area$gard_area)
eoo_areas <- eoo_areas %>% rowwise() %>% 
  mutate(
    iucn_diff      = ((iucn      - smp) / mean(c(iucn,      smp))) * 100,
    gard_diff      = ((gard      - smp) / mean(c(gard,      smp))) * 100,
    iucn_eoo_diff  = ((iucn_eoo  - smp_eoo) / mean(c(iucn_eoo,  smp_eoo))) * 100,
    gard_eoo_diff  = ((gard_eoo  - smp_eoo) / mean(c(gard_eoo,  smp_eoo))) * 100
  ) %>%
  mutate(
    across(2:11, ~ round(.x, 1))
  )

write.csv(eoo_areas, paste0(wd$data, "range size differences.csv"))
