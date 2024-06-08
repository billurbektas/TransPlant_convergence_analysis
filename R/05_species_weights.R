#### GET SPECIES WEIGHTS FROM PRCs #####
# Get species weights
sp.coef = 
  PRC %>%
  filter(type == "all")%>%
  dplyr::select(Region, originSiteID, destSiteID, change, coef)%>%
  unnest(coef)

sp.wei = 
  PRC %>%
  filter(type == "all")%>%
  dplyr::select(change, sp.scores)%>%
  unnest(sp.scores)%>%
  rename(axis = name)
  
sp.scores = 
  left_join(sp.coef, sp.wei, by = c("Region", "originSiteID", "destSiteID", "change","axis"))%>%
  mutate(treatment = ifelse(treatment == "warmed", "warmed", "control"))%>%
  filter(pool %in% c("overlapping","strictly_high_elevation","colonizing"))%>%
  filter(Region != "US_Arizona")%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  mutate(sp.change = (contrast*value))%>% 
  group_by(Region, originSiteID, destSiteID, experiment, change, treatment, Year, axis)%>%
  mutate(contrib = sp.change/sum(abs(sp.change))*100)%>%
  ungroup()%>%
  left_join(metadata %>% dplyr::select(Region, YearEstablished) %>% distinct(), by = "Region")%>%
  mutate(Year_0 = as.numeric(Year) - YearEstablished)%>%
  mutate(log_year_0 = log(Year_0))%>%
  distinct()%>%
  nest(.by = c("axis", "treatment", "pool"))%>%
  filter(treatment == "warmed")
  

sp.lm = 
  sp.scores %>%
    mutate(lmer.complete = purrr::map(data, ~lmer(sp.change~log_year_0*change*experiment + (1|SpeciesName/Region),  
                                                data =.))) %>%
    mutate(cont = purrr::map(lmer.complete, ~left_join(as.tibble(emtrends(., pairwise ~ change, var = "log_year_0")$emtrends), 
                                                       as.tibble(emmeans::test(emtrends(., pairwise ~ change, var = "log_year_0")$emtrends)))))%>%
    mutate(exp = purrr::map(lmer.complete, ~left_join(as.tibble(emtrends(., pairwise ~ change*experiment, var = "log_year_0")$emtrends), 
                                                     as.tibble(emmeans::test(emtrends(., pairwise ~ change*experiment, var = "log_year_0")$emtrends)))))%>%
    mutate(R2 = purrr::map(lmer.complete, ~r.squaredGLMM(.)))

sp.pred.w = 
  sp.lm %>%
  filter(treatment =="warmed")%>%
  mutate(pred = purrr::map(data, ~{.} %>%
                             dplyr::select(Region, SpeciesName, experiment, change) %>%
                             expand_grid(log_year_0 = c(seq(0,2.197225, 0.02), 2.197225))))%>%
  mutate(pred = map2(lmer.complete, pred, ~bind_cols(.y %>% dplyr::select(Region, SpeciesName, experiment, change, log_year_0),
                                                     predictInterval(.x, newdata = as.data.frame(.y), n.sims = 1000, 
                                                                     level = 0.95, type = "linear.prediction",
                                                                     which = "fixed",
                                                                     include.resid.var = FALSE,
                                                                     .parallel = TRUE))))

  
 

# sp.pred.c = 
#   sp.lm %>%
#   filter(treatment =="control")%>%
#   mutate(pred = purrr::map(data, ~{.} %>%
#                              dplyr::select(Region, SpeciesName, experiment, change) %>%
#                              expand_grid(log_year_0 = c(seq(0,2.197225, 0.02), 2.197225))))%>%
#   mutate(pred = map2(lmer.complete, pred, ~bind_cols(.y %>% dplyr::select(Region, SpeciesName, experiment, change, log_year_0),
#                                                      predictInterval(.x, newdata = as.data.frame(.y), n.sims = 1000, 
#                                                                      level = 0.95, type = "linear.prediction",
#                                                                      which = "fixed",
#                                                                      include.resid.var = FALSE,
#                                                                      .parallel = TRUE))))%>%
#   mutate(pred = purrr::map(pred, ~{.} %>%
#                              group_by(change, log_year_0)%>%
#                              summarize_at(vars(fit:lwr), mean)%>%
#                              ungroup()))

for(axnum in c("PRC axis 1", "PRC axis 2")){
sp.pred =
#bind_rows(sp.pred.w, sp.pred.c)%>%
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
                       strictly_high_elevation = "Strictly \nhigh elevation"))%>%
  filter(treatment == "Warmed")%>%
  filter(axis == axnum)

sp.slope=
sp.scores %>%
  unnest(data)%>%
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
  group_by(Region, originSiteID, destSiteID, experiment, Year_0, change, axis, treatment, pool)%>%
  summarize(sp.change = mean(sp.change))%>%
  ungroup()%>%
  mutate(pool = recode(pool, colonizing = "Colonizing",
                       overlapping = "Overlapping",
                       strictly_high_elevation = "Strictly \nhigh elevation"))%>%
  filter(axis == axnum)%>%
  
ggplot()+
  geom_hline(yintercept = 0, color = "black")+
  TP_theme()+
  geom_point(alpha = 0.1, aes(Year_0, sp.change, group = interaction(experiment, change, pool), color = pool))+
  stat_smooth(aes(Year_0, sp.change, group = interaction(experiment, change, pool), color = pool, linetype = change), 
              geom = "line", method = "lm", formula = y~log(x), se = FALSE, linewidth = 0.5, alpha = 0.2)+
  geom_line(data = sp.pred, aes(Year_0, y= fit, group = interaction(experiment, change, pool), color = pool, linetype = change), linewidth = 1)+
  geom_ribbon(data = sp.pred, aes(Year_0, ymin = lwr, ymax = upr, group = interaction(change, pool), fill = pool), alpha = 0.1, show.legend = FALSE)+
  facet_grid(axis~., scales = "free_y")+
  scale_color_manual(values = c("#ff7f00", "#d25fff","#197af6"))+
  scale_fill_manual(values = c("#ff7f00", "#d25fff","#197af6"))+
  scale_linetype_manual(values = c("solid","dotted"))+
  scale_x_continuous(breaks = seq(0,9,1))+
  labs(color = "Species pools", linetype = "In comparison to", x= "Experimental years", y = "Proportional difference in warmed communities")

sp.rate=
  sp.lm %>%
  dplyr::select(axis, treatment, pool, cont)%>%
  unnest(cont)%>%
  mutate(treatment = case_when(treatment == "control"&change == "Distance to origin controls"~"destControls",
                               treatment == "control"&change == "Distance to destination controls"~"originControls",
                               treatment == "warmed"~"warmed"))%>%
  mutate(change = as.character(change))%>%
  mutate(change = ifelse(change == "Distance to destination controls"& treatment == "originControls", "Distance between controls",change))%>%
  mutate(change = ifelse(change == "Distance to origin controls"& treatment == "destControls", "remove",change))%>%
  filter(change !="remove")%>%
  mutate(change = ifelse(change == "Distance to origin controls", "In comparison to origin controls", 
                         "In comparison to destination controls"))%>%
  mutate(change = factor(change, levels = c("In comparison to origin controls",
                                          "In comparison to destination controls")))%>%
  mutate(treatment = recode(treatment, destControls = "Destination controls",
                            originControls = "Origin controls", 
                            warmed = "Warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "PRC axis 1",
                       RDA2 = "PRC axis 2"))%>%
  mutate(pval = factor(pval(p.value), levels = c("","*","***")))%>%
  filter(treatment == "Warmed")%>%
  filter(axis == axnum)%>%
  ggplot(aes(change, log_year_0.trend, color = pool, label = pval))+
  geom_hline(yintercept = 0, color = "grey50")+
  TP_theme()+
  geom_point(size = 3, position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.5, position = position_dodge(width = 0.9))+
  geom_text(position = position_dodge(width = 0.9), aes(y = upper.CL+0.001))+
  facet_grid(axis~., scale = "free_x")+
  scale_color_manual(values = c("#ff7f00", "#d25fff","#197af6"))+
  scale_alpha_manual(values = c(0.2,0.5,1))+
  guides(color = "none")+
  scale_x_discrete(labels = label_wrap(20))+
  theme(strip.text = ggplot2::element_text(size  = 16,  hjust = 0))+
  labs(x = "", y = "Change in proportional differences \nin warmed communities \nover experimental years (slopes)", alpha= "Significant slopes")
sp.rate

sp.end = 
  sp.pred %>% 
  filter(Year_0 %in% c(min(Year_0), max(Year_0)))%>%
  mutate(Year_0 = factor(ifelse(Year_0 < 9, "start", "end"), levels = c("start","end")))%>%
  mutate(change = ifelse(change == "origin controls", 
                         "In comparison to origin controls", 
                         "In comparison to destination controls"))%>%
  mutate(change = factor(change, levels = c("In comparison to origin controls",
                                            "In comparison to destination controls")))%>%
  filter(axis == axnum)%>%
  ggplot(aes(change, fit, color = interaction(Year_0,pool), label = Year_0))+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_text(position = position_dodge(width = 0.9),  aes(y = lwr - 0.01))+
  TP_theme()+
  geom_point(size = 3, position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = upr, ymax = lwr), linewidth = 0.5, position = position_dodge(width = 0.9))+
  facet_grid(axis~., scale = "free_x")+
  scale_color_manual(values = rep(c("#ff7f00", "#d25fff","#197af6"), each = 2))+
  guides(color = "none")+
  scale_x_discrete(labels = label_wrap(20))+
  theme(strip.text = ggplot2::element_text(size  = 16,  hjust = 0))+
  labs(x = "", y = "Proportional difference \nin warmed communities")
sp.end

pdf(here("plot",paste0("species_pool_changes_", axnum, ".pdf")), height = 10, width = 15)
pp=ggarrange(
  sp.slope, 
  ggarrange(sp.end, sp.rate, ncol = 1),
  ncol = 2, 
  common.legend = TRUE,
  legend = "bottom", widths = c(1.2, 0.9)
)
print(pp)
dev.off()
}
## Get results table
res.sp =
  sp.lm %>%
  dplyr::select(axis, treatment, pool, cont, R2)%>%
  unnest(R2)%>%
  unnest(cont)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 3)))

write.csv(res.sp, file = here ("output", "res.sp.csv")) 

# Get results per experiment
sp.lm %>%
  filter(treatment == "warmed")%>%
  dplyr::select(treatment, axis, pool, exp)%>%
  unnest(exp)%>%
  mutate(sign = sign(log_year_0.trend))%>%
  mutate(pval = ifelse(p.value<0.05, "s", "ns"))%>%
  group_by(axis, change, pool, pval, sign)%>%
  summarize(n = n()) %>% View
