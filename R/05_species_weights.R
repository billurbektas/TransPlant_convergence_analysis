#### GET SPECIES WEIGHTS FROM PRCs #####
# Get species weights
sp.scores = 
  PRC %>%
  filter(type %in% c("reduced_treatment"))%>%
  dplyr::select(-Region, -originSiteID, -destSiteID, -ODT, -destPlotID, -Year, -comm_spp, -PRC, -prc.res, -type)%>%
  mutate(coef = purrr::map(coef, ~{.}%>%
                      filter(axis == "RDA1")%>%
                      summarize_at("contrast", mean)))%>%
  mutate(sp.scores = purrr::map(sp.scores, ~{.}%>%
                           filter(name == "RDA1")))%>%
  unnest(coef)%>%
  unnest(sp.scores)%>%
  dplyr::select(-name)%>%
  filter(pool != "control")%>%
  filter(Region != "US_Arizona")%>%
  mutate(contrast = as.numeric(sign(contrast)))%>%
  mutate(value = value*contrast, value)%>% # To have the warming effect all the time
  
  left_join(tr %>% dplyr::select(SpeciesName, plant_size, resource_acquisition), by = "SpeciesName")%>%
  mutate(plant_size = replace_na(plant_size, 0), #trees to get 0 weight
         resource_acquisition = replace_na(resource_acquisition,0))%>%
  nest(.by = c("Region", "originSiteID","destSiteID","pool"))%>%
  filter(pool %in% c("overlapping", "strictly_high_elevation", "colonizing"))%>%
  mutate(lm = purrr::map(data, ~{.} %>%
                    lm(value~change, data =.)))%>%
  mutate(em = purrr::map(lm, ~left_join(as.tibble(emmeans(.,c("change"))), as.tibble(emmeans::test(emmeans(.,c("change")))), by = c("change","emmean", "SE", "df"))))%>%
  mutate(em = purrr::map(em, ~{.} %>%
                    rename(p.value.em = p.value)%>%
                    dplyr::select(change, emmean, lower.CL, upper.CL, p.value.em)))%>%
  unnest(em)%>%
  dplyr::select(-lm)%>%
  mutate(pval.em = pval(p.value.em))

## Get species weights ----
sp.scores = 
  left_join(
    sp.scores %>%
      dplyr::select(Region, originSiteID, destSiteID, pool, data)%>%
      mutate(data = purrr::map(data, ~{.} %>%
                          group_by(change)%>%
                          summarize_at(vars(plant_size:resource_acquisition), mean)%>%
                          dplyr::select(change, plant_size, resource_acquisition))) %>%
      unnest(data) %>%
      distinct(),
    sp.scores %>%
      dplyr::select(Region, originSiteID, destSiteID, pool, change, emmean, lower.CL, upper.CL, p.value.em, pval.em),
    by = c("Region", "originSiteID","destSiteID","pool","change")
  )


p1=
  sp.scores %>%
  filter(pool %in% c("colonizing", "overlapping", "strictly_high_elevation"))%>%
  filter(!is.na(pval.em))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(pval.em = ifelse(pval.em == "non-significant", "non-significant", "significant"))%>%
  mutate(pool = recode(pool, colonizing = "Colonizing",
                       overlapping = "Overlapping",
                       strictly_high_elevation = "Strictly \nhigh elevation"))%>%
  group_by(change, pool)%>%
  mutate(mean = mean(emmean, na.rm = TRUE))%>%
  mutate(max = max(emmean, na.rm = TRUE))%>%
  ggplot(aes(pool, emmean, group = pool, color = pool))+
  geom_hline(yintercept = 0, color = "grey50")+
  TP_theme()+
  facet_grid(.~change)+
  geom_violin(position = position_dodge(width = 0.8))+
  geom_jitter(aes(fill = pval.em), shape = 21, size = 2, color = "white", position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1))+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, fill="red", position = position_dodge(width = 0.8))+
  geom_text(aes(label = round(mean, 3), y = max+0.02), position = position_dodge(width =0.3), show.legend = FALSE)+
  scale_fill_manual(values = c("grey50", "black"))+
  scale_color_manual(values = c( "#ff7f00", "#d25fff","#197af6"))+
  theme(legend.position = "bottom")+
  guides(color = "none")+
  labs(fill = "",color = "", x = "Species pools", y = "Average species weights")

pdf(here("plot", "species_weights.pdf"), height = 8, width = 10)
print(p1)
dev.off()
