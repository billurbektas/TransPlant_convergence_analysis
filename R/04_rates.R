#### BUILD LINEAR MODELS ####
## Distances along the years----
slopes = 
  PRC %>%
  filter(type == "all")%>%
  dplyr::select(Region, originSiteID, destSiteID, change, coef)%>%
  unnest(coef)%>%
  mutate(contrast = abs(contrast))%>%
  left_join(metadata %>% dplyr::select(Region, YearEstablished)%>% distinct(), by = "Region")%>%
  mutate(Year_0 = as.numeric(Year) - YearEstablished)%>%
  mutate(log_year_0 = log(Year_0))%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))

tax =
  slopes %>%
  rename(distance = contrast)%>%
  mutate(treatment = ifelse(treatment == "warmed", "warmed", "control"))%>%
  nest(.by = c("axis", "treatment"))%>%
  mutate(lmer.complete = purrr::map(data, ~lmer(distance~log_year_0*change*experiment + (1|Region),  data =.,control = lmerControl(optimizer = "bobyqa"))))%>%
  mutate(lmer.nonlog = purrr::map(data, ~lmer(distance~Year_0*change*experiment + (1|Region),  data =.,control = lmerControl(optimizer = "bobyqa"))))%>%
  mutate(R2 = purrr::map(lmer.complete, ~r.squaredGLMM(.)))%>%
  mutate(diff.AIC = map2(lmer.complete, lmer.nonlog, ~AIC(.x, .y)))%>%
  mutate(cont = purrr::map(lmer.complete, ~left_join(as.tibble(emtrends(., pairwise ~ change, var = "log_year_0")$emtrends), 
                                              as.tibble(emmeans::test(emtrends(., pairwise ~ change, var = "log_year_0")$emtrends)))))%>%
  mutate(exp = purrr::map(lmer.complete, ~left_join(as.tibble(emtrends(., pairwise ~ change*experiment, var = "log_year_0")$emtrends), 
                                             as.tibble(emmeans::test(emtrends(., pairwise ~ change*experiment, var = "log_year_0")$emtrends)))))%>%
  
  mutate(pred = purrr::map(data, ~{.} %>%
                      dplyr::select(Region, originSiteID, destSiteID, experiment, change) %>%
                      expand_grid(log_year_0 = c(seq(0,2.197225, 0.005), 2.197225))))%>%
  mutate(pred = map2(lmer.complete, pred, ~bind_cols(.y %>% dplyr::select(Region, originSiteID, destSiteID, experiment, change, log_year_0),
                                                     predictInterval(.x, newdata = as.data.frame(.y), n.sims = 1000, level = 0.95, type = "linear.prediction"))))%>%
  mutate(pred = purrr::map(pred, ~{.} %>%
                      group_by(change, log_year_0)%>%
                      summarize_at(vars(fit:lwr), mean)%>%
                      ungroup()))

# Predictions for the figure
pred =
  tax %>%
  dplyr::select(axis, treatment, pred)%>%
  unnest(pred)%>%
  mutate(treatment = case_when(treatment == "control"&change == "Distance to origin controls"~"destControls",
                               treatment == "control"&change == "Distance to destination controls"~"originControls",
                               treatment == "warmed"~"warmed"))%>%
  mutate(experiment = "model")%>%
  mutate(Year_0 = exp(log_year_0))%>%
  mutate(change = ifelse(change == "Distance to destination controls"& treatment == "originControls", "Distance between controls",change))%>%
  mutate(change = ifelse(change == "Distance to origin controls"& treatment == "destControls", "remove",change))%>%
  filter(change !="remove")%>%
  mutate(change = factor(change, levels = c("Distance between controls",
                                            "Distance to origin controls",
                                            "Distance to destination controls")))%>%
  mutate(treatment = recode(treatment, destControls = "Destination controls",
                            originControls = "Origin controls", 
                            warmed = "Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "Canonical coefficient (PRC axis 1)",
                       RDA2 = "Canonical coefficient (PRC axis 2)"))

## Get Figures----
pslope =
  slopes %>%
  mutate(change = ifelse(change == "Distance to destination controls"& treatment == "originControls", "Distance between controls",change))%>%
  mutate(change = ifelse(change == "Distance to origin controls"& treatment == "destControls", "remove",change))%>%
  filter(change !="remove")%>%
  mutate(change = factor(change, levels = c("Distance between controls",
                                            "Distance to origin controls",
                                            "Distance to destination controls")))%>%
  mutate(treatment = recode(treatment, destControls = "Destination controls",
                            originControls = "Origin controls", 
                            warmed = "Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "Canonical coefficient (PRC axis 1)",
                       RDA2 = "Canonical coefficient (PRC axis 2)"))%>%
  
  ggplot()+
  geom_hline(yintercept = 0, color = "black")+
  TP_theme()+
  geom_point(alpha = 0.3, aes(Year_0, contrast, group = interaction(experiment, change), color = change))+
  stat_smooth(aes(Year_0, contrast, group = interaction(experiment, change), color = change, linetype = change), 
              geom = "line", method = "lm", formula = y~log(x), se = FALSE, linewidth = 0.5, alpha = 0.3)+
  geom_line(data = pred, aes(Year_0, y= fit, group = interaction(experiment, change), color = change, linetype = change), linewidth = 1)+
  geom_ribbon(data = pred, aes(Year_0, ymin = lwr, ymax = upr, group = change, fill = change), alpha = 0.1, show.legend = FALSE)+
  facet_grid(axis~., scales = "free_y")+
  scale_color_manual(values = c("grey30","#ff0000","#ff0000"))+
  scale_fill_manual(values = c("grey30","#ff0000","#ff0000"))+
  scale_x_continuous(breaks = seq(0,9,1))+
  scale_y_continuous(breaks = seq(0,0.7,0.1))+
  scale_linetype_manual(values = c("solid","solid","dotted"))+
  labs(color = "", linetype = "", x= "Experimental years", y = "Taxonomic distance")
pslope

ptax=
  tax %>%
  dplyr::select(axis, treatment, cont)%>%
  unnest(cont)%>%
  mutate(treatment = case_when(treatment == "control"&change == "Distance to origin controls"~"destControls",
                               treatment == "control"&change == "Distance to destination controls"~"originControls",
                               treatment == "warmed"~"warmed"))%>%
  mutate(change = as.character(change))%>%
  mutate(change = ifelse(change == "Distance to destination controls"& treatment == "originControls", "Distance between controls",change))%>%
  mutate(change = ifelse(change == "Distance to origin controls"& treatment == "destControls", "remove",change))%>%
  filter(change !="remove")%>%
  mutate(change = factor(change, levels = c("Distance between controls",
                                            "Distance to origin controls",
                                            "Distance to destination controls")))%>%
  mutate(treatment = recode(treatment, destControls = "Destination controls",
                            originControls = "Origin controls", 
                            warmed = "Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "Canonical coefficient (PRC axis 1)",
                       RDA2 = "Canonical coefficient (PRC axis 2)"))%>%
  mutate(pval = factor(pval(p.value), levels = c("","*","***")))%>%
  
  
  ggplot(aes(change, log_year_0.trend, color = change, label = pval))+
  geom_hline(yintercept = 0, color = "grey50")+
  TP_theme()+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 1)+
  geom_text(aes(y = upper.CL+0.003))+
  facet_grid(axis~., scale = "free_x")+
  scale_color_manual(values = c("grey30","#ff0000","#ff0000"))+
  scale_alpha_manual(values = c(0.2,0.5,1))+
  guides(color = "none")+
  scale_x_discrete(labels = label_wrap(10))+
  theme(strip.text = ggplot2::element_text(size  = 16,  hjust = 0))+
  labs(x = "", linetype = "", y = "Change in taxonomic distances over experimental years (slopes)", alpha= "Significant slopes")
ptax

pdf(here("plot","taxonomic_distances.pdf"), height = 10, width = 13)
pp=ggarrange(pslope, ptax, widths = c(1.2, 0.9), common.legend = FALSE, legend = "bottom")
print(pp)
dev.off()

 ## Get results table
res.tax =
  tax %>%
  dplyr::select(axis, treatment, cont, R2)%>%
  unnest(R2)%>%
  unnest(cont)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 3)))

write.csv(res.tax, file = here ("output", "res.tax.csv")) 
# Get results per experiment
tax %>%
  filter(treatment == "warmed")%>%
  dplyr::select(treatment, axis, exp)%>%
  unnest(exp)%>%
  mutate(pval = ifelse(p.value<0.05, "s", "ns"))%>%
  mutate(sign = factor(as.character(sign(log_year_0.trend))))%>%
  group_by(axis, change, pval, sign)%>%
  summarize(n = n())

tax %>%
  filter(treatment == "warmed")%>%
  dplyr::select(treatment, axis, exp)%>%
  unnest(exp)%>%
  mutate(pval = ifelse(p.value<0.05, "s", "ns"))%>%
  group_by(axis, change)%>%
  summarize(max = max(log_year_0.trend), min = min(log_year_0.trend))
