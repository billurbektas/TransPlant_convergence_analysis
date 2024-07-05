spools = pools %>%
  dplyr::select(comm_pool)%>%
  unnest(comm_pool)%>%
  dplyr::select(Region, originSiteID, destSiteID, SpeciesName, pool)%>%
  ungroup()%>%
  distinct()%>%
  filter(!is.na(pool))%>%
  filter(originSiteID !=destSiteID)

dfs = 
  dfx %>%
  mutate(experiment = paste0(originSiteID, " → ", destSiteID))%>%
  dplyr::select(Region, originSiteID, destSiteID, experiment, comm_wide)%>%
  mutate(comm_wide = map(comm_wide, ~{.} %>% pivot_longer(cols = !c(ODT, Year, destPlotID),  # Exclude the first three columns
                                                          names_to = "SpeciesName",         # Name of the new 'key' column
                                                          values_to = "Rel_Cover" )))%>%
  unnest(comm_wide)%>%
  left_join(spools)%>%
  rename(Gradient = Region)%>%
  left_join(exptab %>% dplyr::select(Gradient, destSiteID,`Transplantation year`))

sPRC = PRC %>%
  dplyr::select(Region, originSiteID, destSiteID, change, coef, sp.scores)%>%
  mutate(experiment = paste0(originSiteID, " → ", destSiteID))%>%
  dplyr::select(-originSiteID, -destSiteID)
  
ssp.pred = 
  sp.pred.w %>%
  mutate(pred = purrr::map(pred, ~{.} %>%
                             group_by(change, log_year_0)%>%
                             summarize_at(vars(fit:lwr), mean)%>%
                             ungroup()))%>%
  dplyr::select(axis, treatment, pool, pred)%>%
  unnest(pred)%>%
  mutate(experiment = "model")%>%
  mutate(Year_0 = exp(log_year_0))%>%
  mutate(change = ifelse(change == "Distance to destination controls"& treatment == "originControls", "Distance between controls",change))%>%
  mutate(change = ifelse(change == "Distance to origin controls"& treatment == "destControls", "remove",change))%>%
  filter(change !="remove")%>%
  mutate(change = ifelse(change == "Distance to origin controls", "origin controls", "destination controls"))%>%
  mutate(change = factor(change, levels = c("origin controls",
                                            "destination controls")))%>%
  mutate(treatment = recode(treatment, destControls = "Destination controls",
                            originControls = "Origin controls", 
                            warmed = "Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "PRC axis 1",
                       RDA2 = "PRC axis 2"))%>%
  mutate(pool = recode(pool, colonizing = "Colonizing",
                       overlapping = "Overlapping",
                       strictly_high_elevation = "Strictly high elevation"))%>%
  filter(treatment == "Warmed")

pvar = reg%>%
  mutate(pool = recode(pool, overall = "community",
                       strictly_high_elevation = "strictly high-elevation"))%>%
  
  mutate(pool = factor(pool, 
                       levels = c("community", "strictly high-elevation", 
                                  "overlapping", "colonizing")))%>%
  mutate(type = factor(type, levels = c("start", "end","rate")))%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  dplyr::select(-originSiteID, -destSiteID)%>%
  rename(`Destination precipitation` = destP,
         `Destination temperature` = destT, 
         `Experimental warming` = diffT, 
         `Destination CWM plant size` = destPS,
         `Origin CWM plant size` = oriPS,
         `Destination CWM resource-acquisition` = destRA,
         `Origin CWM resource-acquisition` = oriRA,
         `Experimental duration` = YearRange,
         `Plot size` = PlotSize
         )%>%
  mutate(pool = ifelse(pool != "community", paste0("Species pool: ", pool), "Community"))%>%
  pivot_longer(cols = `Destination precipitation`:`Plot size`, names_to = "variables")
  

saveRDS(sp.scores, file = "output/sp.scores.rds")
saveRDS(pred, file = "output/pred.rds")
saveRDS(slopes, file = "output/slopes.rds")
saveRDS(spools, file = "output/pools.rds")
saveRDS(dfs, file = "output/dfs.rds")
saveRDS(sPRC, file = "output/PRC.rds")
saveRDS(ssp.pred, file = "output/sp.pred.rds")
saveRDS(pvar, file = "output/pvar.rds")

