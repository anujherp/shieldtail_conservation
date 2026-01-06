## Script 9: This chunk plots the contribution of threat parameters to the overall threat score, Figure 2B

# Load uro_threat, which I created in script 3.
uro_threat <- read_csv(paste0(wd$data, "uro_threat.csv"))

# threat cont will act as a variable storing individual threat contributions.
threat_cont <- uro_threat %>% rowwise() %>%
  mutate(built_cont = log_clust*built_clust,
         npa_cont = log_clust*NPA_covp_clust,
         treeloss_cont = log_clust*tl_clust,
         species = paste0(substr(species, 1, 1),". ", word(species, start = 2, sep=" "))) %>%
  dplyr::select(species, built_cont, npa_cont, treeloss_cont, threat_score, z_score, threat_cat) %>%
  ungroup()

## Threat contribution to export as a csv.
threat_cont_tb <- threat_cont %>% rowwise() %>%
  mutate(prop_built = round(built_cont/threat_score,2),
         prop_npa = round(npa_cont/threat_score,2),
         prop_treeloss = round(treeloss_cont/threat_score, 2),
         risk_cat = as.factor(threat_cat))

## Reordering the dataframe by descending order of total_threat, for plotting
threat_cont <- gather(threat_cont_tb, key = "threat", value = "val", 2:4)
threat_cont$species <- fct_reorder(threat_cont$species, threat_cont$threat_score, .desc = TRUE)

# I am trying to plot Z-score vertical axis in the same figure.
#so here I define transformation function for the secondary axis.
z_min <- min(threat_cont$z_score)
z_max <- max(threat_cont$z_score)
val_min <- 5.85  # Assuming y-axis starts from 0
val_max <- 27.90 # Same as limits in scale_y_continuous

sec_axis_transform <- function(y) {
  ((y - val_min) / (val_max - val_min)) * (z_max - z_min) + z_min
}
# Define z-score breaks manually (adjust as needed)
z_breaks <- seq(floor(z_min), ceiling(z_max), by = 0.5)  # Example: Breaks every 0.5


threat_plot <- ggplot(threat_cont, aes(x = species, y = val, fill = threat)) +
  geom_bar(
    stat = "identity",
    width = 0.75,
    color = "#57595B",
    lwd = 3,
    alpha = 0.95,
    show.legend = F
  ) +
  scale_fill_manual(values = c(
    "npa_cont"      = "#658C58",   # or any colour scheme you prefer
    "built_cont"    = "#F0E491",
    "treeloss_cont" = "#BBC863"
  )) +
  labs(
    x = "Species",
    y = "Threat Score",
    fill = "Threat Type"
  ) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(
    breaks = seq(0, 30, 3),
    limits = c(0, 30),
    sec.axis = sec_axis(~ sec_axis_transform(.), name = "Z-statistic")
  ) +
  theme_bw(base_size = 70) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.75, vjust = 0.9, face = "italic")
  )

# Saving the plot as a jpeg file.
ggsave(threat_plot, filename = "threat_cont.png",path = paste0(wd$output), width = 700, 
       height = 800, units = "mm", device='png', dpi=500)
