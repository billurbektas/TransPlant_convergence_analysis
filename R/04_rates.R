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
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(treatment = case_when(treatment == "control"&change == "Distance to origin controls"~"Destination controls",
                               treatment == "control"&change == "Distance to destination controls"~"Origin controls",
                               treatment == "warmed"~"Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "Canonical coefficient (PRC axis 1)",
                       RDA2 = "Canonical coefficient (PRC axis 2)"))%>%
  mutate(experiment = "model")%>%
  mutate(Year_0 = exp(log_year_0))

## Get Figures----
pslope=
  slopes %>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(treatment = recode(treatment, destControls = "Destination controls",
                            originControls = "Origin controls", 
                            warmed = "Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "Canonical coefficient (PRC axis 1)",
                       RDA2 = "Canonical coefficient (PRC axis 2)"))%>%
  
  ggplot()+
  geom_hline(yintercept = 0, color = "black")+
  TP_theme()+
  geom_point(alpha = 0.3, aes(Year_0, contrast, group = interaction(experiment, treatment), color = treatment))+
  stat_smooth(aes(Year_0, contrast, group = interaction(experiment, treatment), color = treatment), geom = "line", method = "lm", formula = y~log(x), se = FALSE, linewidth = 0.5, alpha = 0.3)+
  geom_line(data = pred, aes(Year_0, y= fit, group = interaction(experiment, treatment), color = treatment), linewidth = 1)+
  geom_ribbon(data = pred, aes(Year_0, ymin = lwr, ymax = upr, group = interaction(experiment, treatment), fill = treatment), alpha = 0.1, show.legend = FALSE)+
  facet_grid(axis~change)+
  scale_color_manual(values = c("#ff7f00","#197af6","#ff0000"))+
  scale_fill_manual(values = c("#ff7f00","#197af6","#ff0000"))+
  scale_x_continuous(breaks = seq(0,9,1))+
  scale_y_continuous(breaks = seq(0,0.7,0.1))+
  labs(color = "", x= "Experimental years", y = "Taxonomic distance")

ptax=
  tax %>%
  dplyr::select(axis, treatment, cont)%>%
  unnest(cont)%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(treatment = case_when(treatment == "control"&change == "Distance to origin controls"~"Destination control",
                               treatment == "control"&change == "Distance to destination controls"~"Origin control",
                               treatment == "warmed"~"Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "Canonical coefficient (PRC axis 1)",
                       RDA2 = "Canonical coefficient (PRC axis 2)"))%>%
  mutate(pval = factor(pval(p.value), levels = c("non-significant","*","***")))%>%
  
  
  ggplot(aes(treatment, log_year_0.trend, color = treatment, alpha = pval))+
  geom_hline(yintercept = 0, color = "grey50")+
  TP_theme()+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 1)+
  facet_grid(axis~change, scale = "free_x")+
  scale_color_manual(values = c("#ff7f00","#197af6","#ff0000"))+
  scale_alpha_manual(values = c(0.3,0.6,1))+
  guides(color = "none")+
  scale_x_discrete(labels = label_wrap(10))+
  theme(strip.text = ggplot2::element_text(size  = 12,  hjust = 0))+
  labs(x = "", y = "Change in taxonomic distances over experimental years (slopes)", alpha= "Significant slopes")

pdf(here("plot","taxonomic_distances.pdf"), height = 9, width = 15)
pp=ggarrange(pslope, ptax, widths = c(1.2, 0.9), common.legend = TRUE, legend = "bottom")
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
