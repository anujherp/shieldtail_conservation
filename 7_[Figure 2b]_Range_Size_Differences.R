## Script 7: This chunk computes Range size differences between rangemaps from different sources.
# It also plots Figure 2b

# List of shieldtail species of interest.
uro_list <- c("Melanophidium khairei","Platyplectrurus madurensis","Plectrurus perroteti",
              "Rhinophis karinthandani","Rhinophis melanoleucus","Teretrurus hewstoni",
              "Teretrurus sanguineus","Uropeltis ellioti","Uropeltis maculata",
              "Uropeltis myhendrae","Uropeltis pulneyensis","Uropeltis rubrolineata")

# Load the range maps (SMP, IUCN, and GARD).
uro_iucn <- read_sf(paste0(wd$shapefiles,"uro_iucn.shp")) %>% filter(species %in% uro_list)
uro_gard <- read_sf(paste0(wd$shapefiles,"uropeltid_gard.shp")) %>% filter(binomial %in% uro_list)
uro_rangemaps <- read_sf(paste0(wd$shapefiles,"uro_ranges.shp")) %>% filter(species %in% uro_iucn$species) %>%
  arrange(species)
## Calculating the area (called extent here) of IUCN range and the overlap between SMP & IUCN. 
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

## Calculating difference in areas
areas <- areas %>% rowwise() %>% 
  mutate(iucn_diff = ((iucn - smp)/mean(c(iucn,smp))) *100,
  gard_diff = ((gard - smp)/mean(c(gard,smp))) *100) %>% 
  pivot_longer(cols = c(iucn_diff, gard_diff), names_to = "dataset", values_to = "difference") %>% 
  dplyr::select(-c(gard,iucn,smp))

areas_csv <- pivot_wider(areas, names_from = dataset, values_from = difference)
write.csv(areas_csv, paste0(wd$data, "over and underestimation.csv"))

## Plotting a bar graph (Fgure 2b)
ggplot(data = areas, aes(x = species, y = difference, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, show.legend = F, alpha = 0.85) +
  scale_fill_manual(values = c("#FFB31A","#de1a24")) +
  geom_hline(yintercept = 0, colour = "black", linetype = 1, lwd = 3) +
  xlab("Species") + ylab("Range Size Difference") + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(breaks = seq(-200, 200, 100), limits = c(-300, 300), labels = c("-200%", "-100%", "0%", "100%", "200%")) +
  theme_light(base_size = 20) +
  theme(axis.text.x = element_text(face = "italic"))

## Saving the plot as a jpeg file.
ggsave(filename = "iucn_estimation.jpg",path = paste0(wd$output), width = 1000, 
       height = 500, units = "mm", device='jpeg', dpi=300)