## Script 8: This chunk makes a bubble plot (from Figure 3a).
# It also plots Figure S2.

# Load the datafile
uro_ranges <- read_csv(paste0(wd$data,"uro_ranges.csv"), show_col_types = FALSE) %>%
  mutate(sl = 1:12)

# Scatter Bubble Plot (Figure 3a)
ggplot(uro_ranges, aes(y = b_exp_covp, x = lost_covp, size = NPA_covp)) +
  labs(y = "Built Environment Expansion", x = "Tree Cover Loss") +
  geom_text_repel(aes(label = paste0(substr(species, 1, 1),". ", word(species, start = 2, sep=" "))), size = 15, color = "black", 
                  box.padding = 10, fontface = "italic") +  # Italic species labels
  scale_x_continuous(breaks=seq(0,15,3), limits = c(0, 17), labels = c("0%", "3%", "6%", "9%","12%", "15%")) +
  scale_y_continuous(breaks=seq(0,100,20), limits = c(0,105), labels = c("0%", "20%", "40%", "60%","80%", "100%")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  geom_point(fill = alpha("#3f8f29", 0.3), colour = "#056517", shape = 21, stroke = 3.5, show.legend = F) +
  scale_size(range = c(1, 90)) +
  theme_bw(base_size = 70)

# Saving the plot as a jpeg file.
ggsave(filename = "pa_vs_be_lab.jpg",path = paste0(wd$output), width = 600, 
       height = 800, units = "mm", device='jpeg', dpi=300)

####
## Scatter plot from Figure S2
ggplot(uro_ranges, aes(y = b_exp_covp, x = NPA_covp)) +
  labs(y = "Built Environment Expansion [%]", x = "Outside Protected Areas [%]") +
  scale_x_continuous(breaks=seq(0,100,20), limits = c(0, 103)) +
  scale_y_continuous(breaks=seq(0,100,20), limits = c(0,105)) +
  coord_cartesian(expand = FALSE, clip = "off") +
  geom_point(fill = alpha("#3f8f29", 0.3), colour = "#056517", shape = 16, stroke = 3.5, show.legend = F, size = 16) +
  geom_smooth(method = "lm", se = F, size = 8) +
  scale_size(range = c(1, 90)) +
  theme_bw(base_size = 70)

# Saving the plot as a jpeg file.
ggsave(filename = "pa_vs_be_lm.jpg",path = paste0(wd$output), width = 650, 
       height = 750, units = "mm", device='jpeg', dpi=300)