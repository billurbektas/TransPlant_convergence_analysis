#### ASSESS EXPERIMENTAL DIFFERENCES ####
## Get climate ----
clim = 
  get.clim(sites = slopes %>% 
             dplyr::select(Region, originSiteID, destSiteID, Year) %>% 
             mutate(Year = as.numeric(Year)) %>%
             distinct(),
           climdata = climdata)%>%
  ungroup()%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))


## Get metadata ----
regmeta = 
  metadata %>%
  filter(Region %in% unique(df$Region))%>%
  mutate(YearRange = as.numeric(scale(YearRange, center = TRUE, scale = TRUE)))%>%
  mutate(PlotSize = as.numeric(scale(PlotSize, center = TRUE, scale = TRUE)))%>%
  dplyr::select(Region, destSiteID, YearRange, PlotSize)

# Assess experimental differences ----
reg = tax$exp[[2]] %>%
  mutate(pval = pval(p.value))%>%
  left_join(clim %>% mutate(across(destP:cumsumT, ~as.numeric(scale(., center = TRUE, scale = TRUE)))), by = c("experiment"))%>%
  left_join(cwm, by = c("Region", "originSiteID","destSiteID", "experiment"))%>%
  left_join(regmeta, by = c("Region","destSiteID"))

mf = as.formula("log_year_0.trend ~ 
                PlotSize + YearRange +
                diffT + diffT:YearRange +
                destP + destT +
                destPS + destRA + 
                oriPS + oriRA")

regmod = 
  reg %>%
  nest(.by = "change")%>%
  mutate(data = purrr::map(data, ~get.mod(data=., mf = mf)))%>%
  unnest(data)

pdf(here("plot","experimental_effects_slopes.pdf"), height = 13, width = 12)
pp = regmod %>%
  filter(rowname != "(Intercept)")%>%
  mutate(var = recode(rowname, PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      `YearRange:diffT` = "Experiment duration \n& warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM (plant size)",
                      oriRA = "Origin CWM (resource acquisition)",
                      destPS = "Destination CWM (plant size)",
                      destRA = "Destination CWM (resource acquisition)"))%>%
  mutate(type = case_when(var %in% c("Plot size", "Experiment duration", "Experimental warming", "Experiment duration \n& warming") ~ "Experimental effects",
                          var %in% c("Destination PET","Destination temperature")~"Climate effects",
                          var %in% c("Origin CWM (plant size)","Origin CWM (resource acquisition)",
                                     "Destination CWM (plant size)", "Destination CWM (resource acquisition)")~"Functional effects"))%>%
  mutate(type = factor(type, levels = c("Experimental effects","Climate effects","Functional effects")))%>%
  mutate(var = factor(var, levels = c("Destination CWM (plant size)", "Destination CWM (resource acquisition)",
                                      "Origin CWM (plant size)","Origin CWM (resource acquisition)",
                                      "Destination PET","Destination temperature","Experiment duration \n& warming","Experimental warming",
                                      "Plot size","Experiment duration")))%>%
  mutate(pval = factor(pval, levels = c("non-significant","*","**")))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  ggplot(aes(Estimate, var, alpha = pval))+
  geom_vline(xintercept = 0, color = "grey30")+
  TP_theme()+
  facet_grid(type~change, scales = "free_y")+
  geom_point(position = position_dodge(width = 0.5), size = 3)+
  geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`),
                 position = position_dodge(width = 0.5),
                 linewidth = 1)+
  scale_alpha_manual(values = c(0.2, 0.7, 1))+
  labs(y = "", x = "Effect sizes", alpha = "")
print(pp)
dev.off()

regmod =
  regmod %>%
  mutate(var = recode(rowname, 
                      `(Intercept)` = "Intercept",
                      PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      `YearRange:diffT` = "Experiment duration \n& warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM (plant size)",
                      oriRA = "Origin CWM (resource acquisition)",
                      destPS = "Destination CWM (plant size)",
                      destRA = "Destination CWM (resource acquisition)"))

regmod1 = regmod %>% dplyr::select(change, var, Estimate, `Std. Error`, `t value`,  
                            `2.5 %`, `97.5 %`, `Pr(>|t|)`)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

regmod2 = regmod %>% dplyr::select(change, R, Radj, F.stat, numdf, dendf, p.value)%>%distinct()%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

write.csv (regmod1, file = here("output", "regmod1.csv"))
write.csv (regmod2, file = here("output", "regmod2.csv"))


## Assess experimental differences in species weights ----
regsp = 
  sp.scores %>%
  left_join(regmeta, by = c("Region","destSiteID"))%>%
  left_join(clim %>% mutate(across(destP:cumsumT, ~as.numeric(scale(., center = TRUE, scale = TRUE)))), by = c("Region", "originSiteID","destSiteID"))

mf = as.formula("emmean~
                PlotSize + YearRange + 
                diffT + YearRange:diffT + 
                destP + destT + 
                plant_size + resource_acquisition")

regspmod = 
  regsp %>%
  nest(.by = c("pool","change"))%>%
  mutate(data = purrr::map(data, ~get.mod(data = ., mf = mf)))%>%
  unnest(data)

pdf(here("plot","experimental_effects_pool.pdf"), height = 11, width = 10)
pp1= regspmod %>%
  filter(rowname != "(Intercept)")%>%
  mutate(pool = recode(pool, overlapping = "Overlapping",
                       colonizing = "Colonizing",
                       strictly_high_elevation = "Strictly \nhigh elevation"))%>%
  mutate(var = recode(rowname, PlotSize = "Turf size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      `YearRange:diffT` = "Experiment duration \n& warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      plant_size = "Plant size traits",
                      resource_acquisition = "Resource acquisition traits"))%>%
  mutate(type = case_when(var %in% c("Turf size", "Experiment duration", "Experimental warming",
                                     "Experiment duration \n& warming") ~ "Experimental effects",
                          var %in% c("Destination PET","Destination temperature")~"Climate effects",
                          var %in% c("Plant size traits","Resource acquisition traits")~"Functional effects"))%>%
  mutate(type = factor(type, levels = c("Experimental effects","Climate effects","Functional effects")))%>%
  mutate(var = factor(var, levels = c("Plant size traits","Resource acquisition traits",
                                      "Destination PET","Destination temperature","Experiment duration \n& warming",
                                      "Experimental warming",
                                      "Turf size","Experiment duration")))%>%
  mutate(pval = factor(pval, levels = c("non-significant","*","**", "***")))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  ggplot(aes(Estimate, var, alpha = pval, color = pool))+
  geom_vline(xintercept = 0, color = "grey30")+
  TP_theme()+
  facet_grid(type~change, scales = "free_y")+
  geom_point(position = position_dodge(width = 0.5), size = 3)+
  geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`),
                 position = position_dodge(width = 0.5),
                 linewidth = 1)+
  scale_alpha_manual(values = c(0.2, 0.6,0.8, 1))+
  scale_color_manual(values = c( "#ff7f00", "#d25fff","#197af6"))+
  guides(color = guide_legend(order = 1, ncol = 1),
         alpha = guide_legend(order = 2, ncol = 1))+
  labs(y = "", x = "Effect sizes", alpha = "", color = "Species pools")+
  theme(strip.text = ggplot2::element_text(size  = 14,  hjust = 0))
print(pp1)
dev.off()

pdf(here("plot","experimental_effects_pools_scatter.pdf"), height = 10, width = 8)

pp2=regsp %>%
  dplyr::select(Region, originSiteID, destSiteID, pool, change, emmean,
         YearRange, PlotSize, diffT, destT, destP, 
         plant_size, resource_acquisition)%>%
  pivot_longer(cols = YearRange:resource_acquisition, names_to = "rowname")%>%
  left_join(regspmod, by = c("pool","change","rowname"))%>%
  filter(rowname != "(Intercept)")%>%
  mutate(pool = recode(pool, overlapping = "Overlapping",
                       colonizing = "Colonizing",
                       strictly_high_elevation = "Strictly \nhigh elevation"))%>%
  mutate(var = recode(rowname, PlotSize = "Turf size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      plant_size = "Plant size traits",
                      resource_acquisition = "Resource acquisition traits"))%>%
  mutate(type = case_when(var %in% c("Turf size", "Experiment duration", "Experimental warming") ~ "Experimental effects",
                          var %in% c("Destination PET","Destination temperature")~"Climate effects",
                          var %in% c("Plant size traits","Resource acquisition traits")~"Functional effects"))%>%
  mutate(type = factor(type, levels = c("Experimental effects","Climate effects","Functional effects")))%>%
  mutate(var = factor(var, levels = c("Plant size traits","Resource acquisition traits",
                                      "Destination PET","Destination temperature","Experimental warming",
                                      "Turf size","Experiment duration")))%>%
  mutate(pval = factor(pval, levels = c("non-significant","*","**", "***")))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  ggplot(aes(value, emmean, color = pool, alpha = pval, group = interaction(var, pool, pval)))+
  TP_theme()+
  geom_hline(yintercept = 0, color = "grey20")+
  geom_point(size = 0.8)+
  stat_smooth(fullrange = TRUE, method = "lm", geom = "line", se = FALSE, show.legend = TRUE, linewidth = 1.2)+
  facet_grid(type~change)+
  scale_alpha_manual(values = c(0.2, 0.7, 0.9, 1))+
  scale_color_manual(values = c( "#ff7f00", "#d25fff","#197af6"))+
  labs(y = "Average species weights", x = "Values of independent factors",
       color = "Species pools",
       alpha = "")+
  guides(color = guide_legend(order = 1, ncol = 1),
         alpha = guide_legend(order = 2, ncol = 1))+
  theme(strip.text = ggplot2::element_text(size  = 14,  hjust = 0))
print(pp2)
dev.off()  

pdf(here("plot", "experimental_effects_pools_merged.pdf"), height = 10, width = 16)
pp = ggarrange(pp1, pp2, widths = c(1.1,0.9), common.legend = TRUE, legend = "bottom")
print(pp)
dev.off()

regsptab =
  regspmod %>%
  mutate(var = recode(rowname, 
                      `(Intercept)` = "Intercept",
                      PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      cumsumT = "Experimental warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM (plant size)",
                      oriRA = "Origin CWM (resource acquisition)",
                      destPS = "Destination CWM (plant size)",
                      destRA = "Destination CWM (resource acquisition)"))
regspmod1 = regsptab %>% dplyr::select(change, pool, var, Estimate, `Std. Error`, `t value`,  
                                `2.5 %`, `97.5 %`, `Pr(>|t|)`)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

regspmod2 = regsptab %>% dplyr::select(change, pool, R, Radj, F.stat, numdf, dendf, p.value)%>%distinct()%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

write.csv (regspmod1, file = here("output", "regspmod1.csv"))
write.csv (regspmod2, file = here("output", "regspmod2.csv"))

