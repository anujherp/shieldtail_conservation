## Script 3: This chunk contains the complete pathway to calculate the threat score (Equation 2), using threat parameter values and range sizes.

# I'll first load my dataset
uro_ranges <- read_csv(paste0(wd$data,"uro_ranges.csv")) %>% 
  mutate(log_range = round(log10(area_km),2))

# uro threat variable will have threat score related columns appended at the end of this list.
uro_threat <- uro_ranges %>%
  # Sort by area_km and apply clustering
  arrange(desc(log_range)) %>%
  mutate(log_clust = 6 - round(log_range, 2)) %>%
  arrange(b_exp_covp) %>%
  mutate(built_clust = cutree(hclust(dist(b_exp_covp, method = "euclidean"), method = "ward.D2"), k = 3)) %>%
  arrange(NPA_covp) %>%
  mutate(NPA_covp_clust = cutree(hclust(dist(NPA_covp, method = "euclidean"), method = "ward.D2"), k = 3)) %>%
  arrange(lost_covp) %>%
  mutate(tl_clust = cutree(hclust(dist(lost_covp, method = "euclidean"), method = "ward.D2"), k = 3)) %>%
  dplyr::select(species, area_km, log_range, log_clust, b_exp_covp, NPA_covp, lost_covp, built_clust, NPA_covp_clust, tl_clust)

## Threat score is called risk_score over here.
uro_threat_full <- uro_threat %>%
  rowwise() %>%
  mutate(threat_score = round(log_clust * built_clust + log_clust * tl_clust + log_clust * NPA_covp_clust, 3)) %>%
  ungroup() %>%  # Always ungroup after rowwise
  mutate(z_score = round((threat_score - mean(threat_score, na.rm = TRUE)) /
                           sd(threat_score, na.rm = TRUE), 3),
         threat_cat = case_when(z_score > 1 ~ "H",
                              z_score < -1 ~ "L",
                              z_score < 1 & z_score > -1 ~ "M")) %>%
  dplyr::select(species, area_km, log_clust, built_clust, NPA_covp_clust, tl_clust, threat_score, z_score, threat_cat)

# save this dataframe
write_csv(uro_threat_full,paste0(wd$data,"uro_threat.csv"))

## Sensitivity analyses: 1) Excl. Built Env. and 2) Excl. Non-Protected Area

## 1) Compute threat score, excluding the Built Up Area Expansion parameter.

uro_threat_wo_built <- uro_threat %>%
  rowwise() %>%
  mutate(threat_score = round(log_clust * tl_clust + log_clust * NPA_covp_clust, 3)) %>%
  ungroup() %>%  # Always ungroup after rowwise
  mutate(z_score = round((threat_score - mean(threat_score, na.rm = TRUE)) /
                           sd(threat_score, na.rm = TRUE), 3),
         threat_cat = case_when(z_score > 1 ~ "H",
                                z_score < -1 ~ "L",
                                z_score < 1 & z_score > -1 ~ "M")) %>%
  dplyr::select(species, area_km, log_clust, NPA_covp_clust, tl_clust, threat_score, z_score, threat_cat)

## 2) Compute threat score, excluding the Non-Protected Area Coverage parameter.

uro_threat_wo_NPA <- uro_threat %>%
  rowwise() %>%
  mutate(threat_score = round(log_clust * built_clust + log_clust * tl_clust, 3)) %>%
  ungroup() %>%  # Always ungroup after rowwise
  mutate(z_score = round((threat_score - mean(threat_score, na.rm = TRUE)) /
                           sd(threat_score, na.rm = TRUE), 3),
         threat_cat = case_when(z_score > 1 ~ "H",
                                z_score < -1 ~ "L",
                                z_score < 1 & z_score > -1 ~ "M")) %>%
  dplyr::select(species, area_km, log_clust, built_clust, tl_clust, threat_score, z_score, threat_cat)
