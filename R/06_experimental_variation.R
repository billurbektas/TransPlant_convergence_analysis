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

## Get species pool results -----
regsp = 
  sp.lm %>%
  filter(axis == "RDA1")%>%
  dplyr::select(pool, exp)%>%
  unnest(exp)%>%
  mutate(Region = str_extract(experiment, "^[A-Z]{2}_[^_]+"),
         originSiteID = str_extract(experiment, "(?<=_)[^_]+?(?=_[^_]+$)"),
         destSiteID = str_extract(experiment, "[^_]+$"))

regspend =
  sp.pred.w %>%
  mutate(pred = purrr::map(pred, ~{.} %>%
                             group_by(experiment, change, log_year_0)%>%
                             summarize_at(vars(fit:lwr), mean)%>%
                             ungroup()))%>%
  dplyr::select(axis, treatment, pool, pred)%>%
  unnest(pred)%>%
  mutate(Year_0 = exp(log_year_0))%>%
  filter(Year_0 %in% c(min(Year_0), max(Year_0)))%>%
  mutate(Year_0 = ifelse(Year_0 < 9, "start", "end"))%>%
  mutate(SE = (upr - lwr) / (2 * qnorm(0.975)))%>%
  mutate(Region = str_extract(experiment, "^[A-Z]{2}_[^_]+"),
         originSiteID = str_extract(experiment, "(?<=_)[^_]+?(?=_[^_]+$)"),
         destSiteID = str_extract(experiment, "[^_]+$"))%>%
  rename(type = Year_0, response = fit)%>%
  filter(treatment == "warmed" & axis == "RDA1")%>%
  dplyr::select(-treatment, -axis, -lwr, -upr, -log_year_0)
  

# Assess experimental differences ----
reg = tax %>%
  filter(treatment == "warmed" & axis == "RDA1")%>%
  dplyr::select(exp)%>%
  unnest(exp)%>%
  mutate(pool = "overall")%>%
  mutate(Region = str_extract(experiment, "^[A-Z]{2}_[^_]+"),
         originSiteID = str_extract(experiment, "(?<=_)[^_]+?(?=_[^_]+$)"),
         destSiteID = str_extract(experiment, "[^_]+$"))%>%
  bind_rows(regsp)%>%
  dplyr::select(-upper.CL, -lower.CL, -df, -t.ratio, -p.value)%>%
  mutate(type = "rate")%>%
  rename(response = log_year_0.trend)%>%
  bind_rows(regspend)%>%
  left_join(clim %>% mutate(across(destP:cumsumT, ~as.numeric(scale(., center = TRUE, scale = TRUE)))), by = c("Region","originSiteID","destSiteID","experiment"))%>%
  left_join(cwm, by = c("Region", "originSiteID","destSiteID", "experiment"))%>%
  left_join(regmeta, by = c("Region","destSiteID"))%>%
  dplyr::select(-oriP, -oriT, -diffP, -cumsumP, -cumsumT)

pdf(file = "plot/experimental_variation.pdf", width = 16, height = 10)
pvar = reg%>%
    mutate(pool = recode(pool, overall = "community",
                         strictly_high_elevation = "strictly high-elevation"))%>%
    
    mutate(pool = factor(pool, 
                         levels = c("community", "strictly high-elevation", 
                                    "overlapping", "colonizing")))%>%
    mutate(type = factor(type, levels = c("start", "end","rate")))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))
pvar1=
  pvar %>% filter(pool == "community")%>%
ggplot(aes(type, response, color = pool))+
  TP_theme()+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_violin(fill = "white")+
  geom_jitter(size = 0.9)+
  facet_grid(change~pool, scales = "free")+
  scale_color_manual(values = c("black"))+
  labs(y = "Changes in warmed communities in comparison to controls", x = "", color = "")+
  ylim(-0.35, 0.35)

pvar2 =
  pvar %>% filter(pool != "community")%>%
  ggplot(aes(type, response, color = pool))+
  TP_theme()+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_violin(fill = "white")+
  geom_jitter(size = 0.9)+
  facet_grid(change~pool, scales = "free")+
  scale_color_manual(values = c("#197af6", "#d25fff","#ff7f00"))+
  labs(y = "Changes in warmed communities in comparison to controls", x = "", color = "")
pvar = ggarrange(pvar1, pvar2, widths = c(0.5, 1.5))
print(pvar)
dev.off()

# Summary table
reg %>%
  dplyr::select(change, experiment, pool, response, type)%>%
  mutate(sign = sign(response))%>%
  group_by(change, pool, type, sign)%>%
  summarize(count = n()) %>% View


mf = as.formula("~PlotSize + YearRange +
                diffT + destP + destT +
                destPS + destRA + 
                oriPS + oriRA")
regmod = 
  reg %>%
  nest(.by = c("change","pool", "type"))%>%
  mutate(data = purrr::map(data, ~get.mod(data=., mf = mf)))%>%
  unnest(data)%>%
  mutate(var = recode(term, PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM \n(plant size)",
                      oriRA = "Origin CWM \n(resource acquisition)",
                      destPS = "Destination CWM \n(plant size)",
                      destRA = "Destination CWM \n(resource acquisition)"))%>%
  mutate(typex = case_when(var %in% c("Plot size", "Experiment duration", "Experimental warming", "Experiment duration \n& warming") ~ "Experimental effects",
                           var %in% c("Destination PET","Destination temperature")~"Climate effects",
                           var %in% c("Origin CWM \n(plant size)","Origin CWM \n(resource acquisition)",
                                      "Destination CWM \n(plant size)", "Destination CWM \n(resource acquisition)")~"Functional effects"))%>%
  mutate(typex = factor(typex, levels = c("Experimental effects","Climate effects","Functional effects")))%>%
  mutate(var = factor(var, levels = c("Destination CWM \n(plant size)", "Destination CWM \n(resource acquisition)",
                                      "Origin CWM \n(plant size)","Origin CWM \n(resource acquisition)",
                                      "Destination PET","Destination temperature","Experiment duration \n& warming","Experimental warming",
                                      "Plot size","Experiment duration")))%>%
  mutate(pval = factor(pval, levels = c("","*","**", "***")))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(type = factor(type, levels = c("start", "end","rate")))

regpred = 
  reg %>%
  filter(pool != "overall")%>%
  nest(.by = c("change","pool", "type"))%>%
  mutate(data = purrr::map(data, ~pred.mod(data=., mf = mf)))%>%
  unnest(data)%>%
  dplyr::select(-(se:pi.ub))%>%
  rename(term = moderator, response = pred, explanatory = value)%>%
  left_join(regmod, by = c("change","type", "pool","term"))%>%
  mutate(var = recode(term, PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM \n(plant size)",
                      oriRA = "Origin CWM \n(resource acquisition)",
                      destPS = "Destination CWM \n(plant size)",
                      destRA = "Destination CWM \n(resource acquisition)"))%>%
  mutate(typex = case_when(var %in% c("Plot size", "Experiment duration", "Experimental warming", "Experiment duration \n& warming") ~ "Experimental effects",
                           var %in% c("Destination PET","Destination temperature")~"Climate effects",
                           var %in% c("Origin CWM \n(plant size)","Origin CWM \n(resource acquisition)",
                                      "Destination CWM \n(plant size)", "Destination CWM \n(resource acquisition)")~"Functional effects"))%>%
  mutate(typex = factor(typex, levels = c("Experimental effects","Climate effects","Functional effects")))%>%
  mutate(var = factor(var, levels = c("Destination CWM \n(plant size)", "Destination CWM \n(resource acquisition)",
                                      "Origin CWM \n(plant size)","Origin CWM \n(resource acquisition)",
                                      "Destination PET","Destination temperature","Experiment duration \n& warming","Experimental warming",
                                      "Plot size","Experiment duration")))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(type = factor(type, levels = c("start", "end","rate")))



pdf(here("plot","experimental_effects_slopes.pdf"), height = 10, width = 9)
pp = regmod %>%
  filter(term != "intercept")%>%
  filter(pool == "overall")%>%
  ggplot(aes(estimate, var, alpha = pval, label = pval))+
  geom_vline(xintercept = 0, color = "grey30")+
  TP_theme()+
  facet_grid(typex~change, scales = "free_y")+
  geom_text(aes(x = conf.high +0.02), color = "#ff0000")+
  geom_point(size = 3, color = "#ff0000")+
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 linewidth = 1, color = "#ff0000")+
  guides(alpha = "none")+
  scale_alpha_manual(values = c(0.2, 1, 1, 1))+
  labs(y = "", x = "Effect sizes", alpha = "")
print(pp)
dev.off()

p.reg = 
  reg %>%
  pivot_longer(cols = c("PlotSize", "YearRange","diffT", "destP","destT","destPS","destRA","oriPS","oriRA") ,
               names_to = "term", values_to = "explanatory")%>%
  left_join(regmod, by = c("change","type", "pool","term"))%>%
  filter(pool != "overall")%>%
  mutate(var = recode(term, PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM \n(plant size)",
                      oriRA = "Origin CWM \n(resource acquisition)",
                      destPS = "Destination CWM \n(plant size)",
                      destRA = "Destination CWM \n(resource acquisition)"))%>%
  mutate(typex = case_when(var %in% c("Plot size", "Experiment duration", "Experimental warming", "Experiment duration \n& warming") ~ "Experimental effects",
                           var %in% c("Destination PET","Destination temperature")~"Climate effects",
                           var %in% c("Origin CWM \n(plant size)","Origin CWM \n(resource acquisition)",
                                      "Destination CWM \n(plant size)", "Destination CWM \n(resource acquisition)")~"Functional effects"))%>%
  mutate(typex = factor(typex, levels = c("Experimental effects","Climate effects","Functional effects")))%>%
  mutate(var = factor(var, levels = c("Destination CWM \n(plant size)", "Destination CWM \n(resource acquisition)",
                                      "Origin CWM \n(plant size)","Origin CWM \n(resource acquisition)",
                                      "Destination PET","Destination temperature","Experiment duration \n& warming","Experimental warming",
                                      "Plot size","Experiment duration")))%>%
  mutate(pval = factor(pval, levels = c("","*","**", "***")))%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(type = factor(type, levels = c("start", "end","rate")))

p1=
p.reg%>%
  filter(type %in% c("start", "end"))%>%
  filter(pval!= "")%>%
  mutate(change = ifelse(change == "Distance to origin controls", "In comparison to \norigin controls", "In comparison to \ndestination controls"))%>%
  mutate(change = factor(change, levels = c("In comparison to \norigin controls", "In comparison to \ndestination controls")))%>%
  ggplot(aes(explanatory, response, color = pool, group = interaction(var, pool)))+
  TP_theme()+
  geom_hline(yintercept = 0, color = "grey20")+
  geom_point(size = 0.8)+
  geom_line(data = regpred%>%
              filter(type %in% c("start", "end"))%>%
              mutate(change = ifelse(change == "Distance to origin controls", "In comparison to \norigin controls", "In comparison to \ndestination controls"))%>%
              mutate(change = factor(change, levels = c("In comparison to \norigin controls", "In comparison to \ndestination controls")))
              , aes(explanatory,response,  color = pool, alpha = pval, group = interaction(var, pool, pval)),
            linewidth = 1.2)+
  #stat_smooth(fullrange = TRUE, method = "lm", geom = "line", se = FALSE, show.legend = TRUE, linewidth = 1.2)+
  facet_nested(typex~change+type, scales= "free_x")+
  scale_alpha_manual(values = c(0.1, 1, 1, 1))+
  scale_color_manual(values = c( "#ff7f00", "#d25fff","#197af6"))+
  labs(y = "Proportional change in warmed communities", x = "Values of independent factors",
       color = "Species pools",
       alpha = "")+
  guides(color = "none",
         alpha = "none")+
  theme(strip.text = ggplot2::element_text(size  = 14,  hjust = 0))
p1

p2=
regmod %>%
  filter(term != "intercept")%>%
  filter(pool != "overall")%>%
  filter(type %in% c("start", "end"))%>%
  mutate(change = ifelse(change == "Distance to origin controls", "In comparison to \norigin controls", "In comparison to \ndestination controls"))%>%
  mutate(change = factor(change, levels = c("In comparison to \norigin controls", "In comparison to \ndestination controls")))%>%
  ggplot(aes(estimate, var, alpha = pval, label = pval, color = pool))+
  geom_vline(xintercept = 0, color = "grey30")+
  TP_theme()+
  geom_text(aes(x = conf.high +0.02), 
            position = position_dodge(width = 0.8),
            angle = 90, show.legend = FALSE)+
  geom_point(size = 1, position = position_dodge(width = 0.8))+
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 linewidth = 0.8, position = position_dodge(width = 0.8))+
  facet_nested(typex~change+type, scales = "free_y")+
  guides(alpha = "none", color = "none")+
  scale_alpha_manual(values = c(0.2, 1, 1, 1))+
  scale_color_manual(values = c("#ff7f00", "#d25fff","#197af6"))+
  scale_x_continuous(breaks = seq(-0.1,0.1, 0.1))+
  theme(axis.text.x = element_text(size =12))+theme(strip.text = ggplot2::element_text(size  = 14,  hjust = 0))+
  labs(y = "", x = "Effect sizes", alpha = "", color  = "Species pools")

pp = ggarrange(p2, p1)
pdf(here("plot","experimental_effects_species_1.pdf"), height = 8, width = 14)
print(pp)
dev.off()

p1=
  p.reg%>%
  filter(type %in% c("rate"))%>%
  mutate(change = ifelse(change == "Distance to origin controls", "In comparison to \norigin controls", "In comparison to \ndestination controls"))%>%
  mutate(change = factor(change, levels = c("In comparison to \norigin controls", "In comparison to \ndestination controls")))%>%
  ggplot(aes(explanatory, response, color = pool, alpha = pval, group = interaction(var, pool, pval, experiment)))+
  TP_theme()+
  geom_hline(yintercept = 0, color = "grey20")+
  geom_point(size = 0.8)+
  geom_line(data = regpred%>%
              filter(type %in% c("rate"))%>%
              mutate(change = ifelse(change == "Distance to origin controls", "In comparison to \norigin controls", "In comparison to \ndestination controls"))%>%
              mutate(change = factor(change, levels = c("In comparison to \norigin controls", "In comparison to \ndestination controls")))
            , aes(explanatory,response,  color = pool, alpha = pval, group = interaction(var, pool, pval)),
            linewidth = 1.2)+
  #stat_smooth(fullrange = TRUE, method = "lm", geom = "line", se = FALSE, show.legend = TRUE, linewidth = 1.2)+
  facet_nested(typex~change+type, scales= "free_x")+
  scale_alpha_manual(values = c(0.2, 1, 1, 1))+
  scale_color_manual(values = c( "#ff7f00", "#d25fff","#197af6"))+
  labs(y = "Change in proportional differences in warmed communities \nover experimental years (slopes)", x = "Values of independent factors",
       color = "Species pools",
       alpha = "")+
  guides(color = "none",
         alpha = "none")+
  theme(strip.text = ggplot2::element_text(size  = 14,  hjust = 0))
p1

p2=
  regmod %>%
  filter(term != "intercept")%>%
  filter(pool != "overall")%>%
  filter(type %in% c("rate"))%>%
  mutate(change = ifelse(change == "Distance to origin controls", "In comparison to \norigin controls", "In comparison to \ndestination controls"))%>%
  mutate(change = factor(change, levels = c("In comparison to \norigin controls", "In comparison to \ndestination controls")))%>%
  ggplot(aes(estimate, var, alpha = pval, label = pval, color = pool))+
  geom_vline(xintercept = 0, color = "grey30")+
  TP_theme()+
  geom_text(aes(x = conf.high +0.002), 
            position = position_dodge(width = 0.8),
            angle = 90, show.legend = FALSE)+
  geom_point(size = 1, position = position_dodge(width = 0.8))+
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 linewidth = 0.8, position = position_dodge(width = 0.8))+
  facet_nested(typex~change+type, scales = "free_y")+
  guides(alpha = "none", color = "none")+
  scale_alpha_manual(values = c(0.2, 1, 1, 1))+
  scale_color_manual(values = c("#ff7f00", "#d25fff","#197af6"))+
  #scale_x_continuous(breaks = seq(-0.1,0.1, 0.1))+
  theme(axis.text.x = element_text(size =12))+theme(strip.text = ggplot2::element_text(size  = 14,  hjust = 0))+
  labs(y = "", x = "Effect sizes", alpha = "", color  = "Species pools")
p2

pp = ggarrange(p2, p1)

pdf(here("plot","experimental_effects_species_2.pdf"), height = 8, width = 15)
print(pp)
dev.off()

regmod =
  regmod %>%
  mutate(var = recode(term, 
                      intercept = "Intercept",
                      PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      diffT = "Experimental warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM (plant size)",
                      oriRA = "Origin CWM (resource acquisition)",
                      destPS = "Destination CWM (plant size)",
                      destRA = "Destination CWM (resource acquisition)"))

regmod1 = regmod %>% dplyr::select(change, pool, type, var, estimate, std.error, statistic,  
                            conf.low, conf.high, p.value)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

regmod2 = regmod %>% dplyr::select(change, pool, type, i.squared, h.squared, tau.squared, tau.squared.se, 
                                   cochran.qe, p.value.cochran.qe, cochran.qm, p.value.cochran.qm,
                                   df.residual)%>%distinct()%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

write.csv (regmod1, file = here("output", "regmod1.csv"))
write.csv (regmod2, file = here("output", "regmod2.csv"))

# Additional appendix figures ----

pdf(file = "plot/cwm_differences.pdf", width = 10, height = 8)
cwm %>%
  mutate(plant_size = destPS-oriPS,
         resource_acquisition = destRA-oriRA)%>%
  dplyr::select(plant_size, resource_acquisition)%>%
  pivot_longer(cols = plant_size:resource_acquisition)%>%
  ggplot(aes(name, value))+
  geom_hline(yintercept = 0, color = "grey50")+
  TP_theme()+
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, color = "grey50")+
  labs(x = "", y= "Difference in CWMs between \ndestination and origin controls")
dev.off()  

# Relative abundance of the outsider species

pdf(file = "plot/relative_abundances_outsiders.pdf", width = 10, height = 8)
pools%>%
  dplyr::select(comm_pool)%>%
  unnest(comm_pool)%>%
  filter(ODT == "warmed")%>%
  filter(pool == "outsider")%>%
  group_by(Region, originSiteID, destSiteID, Year, UniqueID, pool)%>%
  summarize(tot = sum(Rel_Cover))%>%
  ggplot(aes(x = tot)) +
  TP_theme()+
  geom_histogram(fill = "white", color = "black")+
  geom_vline(aes(xintercept = mean(tot)), color = "pink", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median(tot)), color = "red", linetype = "dashed", size = 1) +
  geom_text(aes(x = mean(tot)+0.07, y = 170, label = sprintf("Mean: %.2f", mean(tot))), vjust = -1.5, color = "pink") +
  geom_text(aes(x = median(tot)-0.07, y = 170, label = sprintf("Median: %.2f", median(tot))), vjust = 1.5, color = "red") +
  
  
  labs(x = "Relative abundances", y = "Number of warmed communities across all experiments")
dev.off()

  
