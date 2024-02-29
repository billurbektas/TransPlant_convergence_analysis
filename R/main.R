library(patchwork)
library(tidyverse)
library(tidyr)
library(purrr) #map, map2, pmap
library(emmeans) #emmans, contrast
library(ggrepel)
library(ggpubr)
library(miceRanger) #miceRanger
library(lme4) #lmer
library(scales) #rescale
library(vegan) #prc, scores
library(MuMIn) #r.squaredGLMM

#detach("package:mapproj", unload = TRUE)
#detach("package:maps", unload = TRUE)

# Get the theme for the figures ----
source("R/theme_ggplot.R")
# Get the functions ----
source("R/functions.R")

# Load the clean data ----
## loads all the pre-processed data (vegetation surveys, traits, climate from CHELSA dataset)

load("data/TransPlantNetwork.RData")

# Organize the tibbles for the analysis ----
explist = list(CH_Calanda, CH_Calanda2, CH_Lavey, #Switzerland
               CN_Damxung, CN_Gongga, CN_Heibei, #China
               DE_Grainau, DE_Susalps, DE_TransAlps, #Germany
               FR_AlpeHuez, FR_Lautaret, #France
               IN_Kashmir, #India
               IT_MatschMazia1, IT_MatschMazia2, #Italy
               NO_Gudmedalen, NO_Lavisdalen, NO_Skjellingahaugen, NO_Ulvhaugen, #Norway
               SE_Abisko, #Sweden
               US_Arizona, US_Colorado, US_Montana #United States
) 
# Get the metadata per experiment
metadata = map_dfr(explist, ~ .x$meta)
metadata =
  metadata %>%
  mutate(YearEstablished = ifelse(Gradient == "FR_Lautaret",2016,YearEstablished))%>% #Correction to Lautaret experiment
  rename(Region = Gradient,
         PlotSize = PlotSize_m2)

# Get the gradients to be taken out from the analysis
takeoutregion = metadata %>%
  filter(YearRange < 3) %>%
  distinct(Region) %>%
  pull(Region)

# Get the undetermined taxa to be taken out from the analysis
takeoutsp = cleanedspecies %>%
  filter(is.na(matched_name)) %>%
  distinct(original_name) %>%
  pull(original_name)%>%
  c("Equ pal/pra")

# Re-calculate the relative cover only for vegetation
dd = dat %>%
  filter(!Region %in% takeoutregion)%>%
  select(-Elevation)%>%
  group_by(Region, Year, originSiteID, originBlockID, destSiteID, destBlockID, destPlotID, Treatment, turfID, UniqueID)%>%
  mutate(Total_Cover = sum(Cover, na.rm = TRUE), Rel_Cover = Cover/Total_Cover)%>%
  ungroup()%>%
  bind_rows(dat %>% filter(Region == "US_Arizona") %>% select(-Elevation))%>% # Only relative cover available for Arizona
  filter(!is.na(Rel_Cover)) %>% #4 NAs in DE_Susalps.
  filter(!is.na(SpeciesName))%>%
  filter(!(Region=="DE_Susalps" & Year==2016)) %>% # Not enough data
  filter(!(Region=="CH_Calanda2" & Year==2020))%>% # OriginControl not counted
  filter(!SpeciesName %in% takeoutsp) %>%
  left_join(metadata %>% select(Region, destSiteID, YearEstablished, YearMin),
            by = c("Region","destSiteID"))%>%
  mutate(Year_0 = Year-YearEstablished)%>%
  filter(Year_0 != 0)%>% # We do not take pre-transplantation vegetation surveys
  select(-(YearEstablished:Year_0))

# Standardize the species names 
cleanedspecies = cleanedspecies %>%
  filter(!(SpeciesName == "Ant alp" & matched_name == "Anthoxanthum odoratum subsp. nipponicum"))%>%
  filter(!(SpeciesName == "Equ pal/pra"))%>%
  filter(!(SpeciesName == "Rum ace subsp lap" & matched_name == "Rumex acetosa subsp. acetosa"))%>%
  filter(!(SpeciesName == "Sal gla" & matched_name == "Salix glauca glauca"))%>%
  filter(!(SpeciesName == "Ver alp" & matched_name == "Veronica alpina subsp. alpina"))

dd = dd %>%
  left_join(cleanedspecies %>% 
              select(Region, destSiteID, SpeciesName, matched_name) %>% 
              distinct(), by = c("SpeciesName", "Region", "destSiteID"))%>%
  mutate(matched_name = ifelse(is.na(matched_name), SpeciesName, matched_name))%>%
  select(-SpeciesName)%>%
  rename(SpeciesName = matched_name)%>%
  distinct()

# Unresolved botanical issues -> we assign genus names
# CN_Damxung: Stipa capillacea, Stipa capillata -> Stipa
#             Polygonaceae -> drives the pattern too much and all of a sudden disappears to be taken out. 
# DE_Susalps (GW-FE): Festuca pratensis (in destination control plots)? Elymus repens (in origin control plots)? -> drives the pattern too much and all of a sudden disappears to be taken out.
# DE_Susalps (GW-BT): Festuca pratensis (in destination control plots)? Elymus repens (in origin control plots)? -> drives the pattern too much and all of a sudden disappears to be taken out.
dd = 
  dd %>%
  mutate(SpeciesName = ifelse(Region == "CN_Damxung" & SpeciesName %in% c("Stipa capillacea","Stipa capillata"), "Stipa", SpeciesName))%>%
  mutate(Rel_Cover = ifelse(Region == "CN_Damxung" & SpeciesName == "Polygonaceae", 0,Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "FE" & destSiteID == "FE" &  SpeciesName == "Festuca pratensis", 0, Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "BT" & destSiteID == "BT" & SpeciesName == "Festuca pratensis", 0, Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "GW" & destSiteID == "GW" & SpeciesName == "Elymus repens", 0, Rel_Cover))%>%
  filter(Rel_Cover != 0)

# Nest the data per experiment
ddx = dd %>%
  select(Region, originSiteID, destSiteID, Treatment) %>% 
  distinct() %>% 
  filter(Treatment == "Warm") %>% 
  select(-Treatment) %>%
  mutate(comm = pmap(.l = list(R = Region, O = originSiteID, D = destSiteID), .f = function(R, O, D){
    bind_rows(
      originControls = dd %>% filter(Region == R & originSiteID == O & destSiteID == O & Treatment == "LocalControl"),
      destControls = dd %>% filter(Region == R &  originSiteID == D & destSiteID == D & Treatment == "LocalControl"),
      warmed =  dd %>% filter(Region == R & originSiteID == O & destSiteID == D & Treatment == "Warm"),
      .id = "ODT")
  })) 

# I take out species that occurred in communities only once over the whole experimental duration.
## This is because we cannot distinguish if there is an botanical identification error or true presence.
### However, I do not change the relative cover as the species contributed to the total cover. 
ddx = 
  ddx %>%
  mutate(comm = map(comm, ~{.} %>%
                      group_by(ODT, SpeciesName)%>%
                      mutate(n = n())%>%
                      filter(n>1)%>% 
                      ungroup()%>%
                      select(-n)))

# Get species pools ----
pools = ddx %>% 
  select(Region, originSiteID, destSiteID, comm) %>%
  mutate(low = map(comm, ~{.} %>% #species pool at low site controls across all years
                     select(ODT, SpeciesName) %>% 
                     filter(ODT == "destControls") %>%
                     distinct(.$SpeciesName) %>% 
                     flatten_chr(.)),
         high = map(comm, ~{.} %>% #species pool of high site controls across all years
                      select(ODT, SpeciesName) %>% 
                      filter(ODT =="originControls") %>%
                      distinct(.$SpeciesName) %>% 
                      flatten_chr(.)),
         warmed = map(comm, ~{.} %>% #species pool of transplanted turfs across all years
                        select(ODT, SpeciesName) %>% 
                        filter(ODT =="warmed") %>%
                        distinct(.$SpeciesName) %>% 
                        flatten_chr(.)),
         high_unique = map2(high, low, ~setdiff(.x, intersect(.x,.y))), #strictly high elevation species pool
         total = map2(high, low, ~unique(c(.x,.y))), #total
         overlap = map2(high, low, ~intersect(.x, .y)), # overlapping species pool
         resident = map2(high_unique, overlap, ~unique(c(.x,.y))), #high species pool and overlap for transplanted turfs
         non_resident = map2(warmed, resident, ~setdiff(.x, intersect(.x,.y))),
         colonizer = map2(non_resident, low, ~intersect(.x, .y)),
         outsider = map2(non_resident, low, ~setdiff(.x, intersect(.x,.y))),
         warmed_overlap = map2(overlap, warmed, ~intersect(.x, .y)), 
         warmed_high_unique = map2(high_unique, warmed, ~intersect(.x, .y)),
         warmed_colonizer = map2(colonizer, warmed, ~intersect(.x, .y)))%>% 
  
  mutate(comm_pool = pmap(list(comm, colonizer, overlap, high_unique, outsider), .f=function(comm, colonizer, overlap, high_unique, outsider){
    comm %>%
      mutate(pool = case_when(
        SpeciesName %in% unique(colonizer) ~ "colonizing",
        SpeciesName %in% unique(overlap) ~ "overlapping",
        SpeciesName %in% unique(high_unique) ~ "strictly_high_elevation",
        SpeciesName %in% unique(outsider) ~ "outsider"
      ))}))

# Get the percentage of each pool in warmed communities per experiment
sum.pools = pools %>% 
  select(-comm, -comm_pool)%>%
  group_by(Region, originSiteID, destSiteID)%>%
  mutate_at(vars(low:warmed_colonizer), ~ map(., length))%>%
  unnest(low:warmed_colonizer)%>%
  mutate_at(vars(high_unique:colonizer), ~(./total)*100)

sum.pools = sum.pools %>%
  select(Region, originSiteID, destSiteID, high_unique, overlap, colonizer)%>%
  pivot_longer(cols = high_unique:colonizer, names_to = "pools")%>%
  mutate(pools = recode(pools, colonizer = "Colonizing \nspecies pool",
                        high_unique = "Strictly high-elevation \nspecies pool",
                        overlap = "Overlapping \nspecies pool"))%>%
  mutate(pools = factor(pools, levels = c("Overlapping \nspecies pool", 
                                          "Strictly high-elevation \nspecies pool", 
                                          "Colonizing \nspecies pool")))%>%
  filter(Region != "US_Arizona")

# Anova
comp = contrast(emmeans(lm(value ~ pools, data = sum.pools), "pools"), "pairwise") %>%
  as_tibble() %>%
  mutate(xmin = map_chr(strsplit(contrast, " - "), ~ gsub("\\(|\\)", "", .x)[1]),
         xmax = map_chr(strsplit(contrast, " - "), ~ gsub("\\(|\\)", "", .x)[2]))

pdf(file = "plot/species_pool_percentages.pdf", height = 8, width = 8)
ggplot(sum.pools, aes(pools, value))+
  TP_theme()+
  geom_boxplot(aes(color = pools), show.legend = FALSE)+
  geom_jitter(aes(color = pools), show.legend = FALSE)+
  geom_bracket(data = comp, aes(xmin = xmin, xmax = xmax, y.position = c(120, 100, 90), label = paste0("p.value = ",round(p.value, 3))))+
  scale_color_manual( values = c("#197af6","#ff7f00","#d25fff"))+
  labs(x = "Species pools", y = "Percentage of species pool per experiment")

dev.off()

# Get the species pools per experiment per species
pools = pools %>%
  select(comm_pool)%>%
  unnest(comm_pool)%>%
  select(Region, originSiteID, destSiteID, ODT, SpeciesName, pool)%>%
  distinct()%>%
  mutate(pool=ifelse(is.na(pool), "low", pool))%>%
  select(-ODT)%>%
  distinct()

# Principal Response Curves----
# Prepare species sites matrices per experiment
ddx = ddx %>% 
  mutate(comm_wide = map(comm, ~{
    .x %>% dplyr::select(ODT, Year, SpeciesName, Rel_Cover, destPlotID) %>% 
      pivot_wider(names_from = SpeciesName, values_from = Rel_Cover, values_fill = list(Rel_Cover = 0), values_fn = list(Rel_Cover = sum)) %>%
      group_by(ODT)%>%
      complete(Year, destPlotID)%>% # There are cases where plots are missing. Create those plots in the dataset.
      group_by(ODT, Year)%>%
      mutate(across(where(is.numeric), ~replace_na(., mean(., na.rm=TRUE))))%>% #Fill the missing plots with the mean of ODT and Year.
      ungroup()
  })) 

# Get PRCs, species scores and coefficients
PRC = bind_rows(get.PRC(data = ddx, sel.ODT = c("warmed","originControls","destControls"), sel.control = "originControls") %>%
                  mutate(change = "Distance to origin controls",
                         type = "all",
                         sel.control = "originControls"),
                get.PRC(data = ddx, sel.ODT = c("warmed","originControls","destControls"), sel.control = "destControls") %>%
                  mutate(change = "Distance to destination controls",
                         type = "all",
                         sel.control = "destControls"),
                get.PRC(data = ddx, sel.ODT = c("warmed","originControls"), sel.control = "originControls") %>%
                  mutate(change = "Distance to origin controls",
                         type = "reduced_treatment",
                         sel.control = "originControls"),
                get.PRC(data = ddx, sel.ODT = c("warmed","destControls"), sel.control = "destControls") %>%
                  mutate(change = "Distance to destination controls",
                         type = "reduced_treatment",
                         sel.control = "destControls"),
                get.PRC(data = ddx, sel.ODT = c("originControls","destControls"), sel.control = "originControls") %>%
                  mutate(change = "Distance to origin controls",
                         type = "reduced_controls",
                         sel.control = "originControls"),
                get.PRC(data = ddx, sel.ODT = c("originControls","destControls"), sel.control = "destControls") %>%
                  mutate(change = "Distance to destination controls",
                         type = "reduced_controls",
                         sel.control = "destControls"))
perm.PRC = get.perm(PRC)
perm.df = 
  perm.PRC %>%
  mutate(perm = map(perm, ~rownames_to_column(.)))%>%
  unnest(perm)%>%
  rename(var = rowname)%>%
  filter(var %in% c("RDA1","RDA2"))%>%
  mutate(pval = pval(`Pr(>F)`))%>%
  select(-sel.control)%>%
  mutate(pval = ifelse(pval == "non-significant","", pval))
## Get PRC tables ----
res.df =
  PRC %>%
  filter(type != "reduced_controls")%>%
  select(Region, originSiteID, destSiteID, change, type, prc.res)%>%
  unnest(prc.res)%>%
  left_join(perm.df %>% select(-(Df:`Pr(>F)`)), by = c("Region","originSiteID","destSiteID","change","type","var"))%>%
  mutate(Variance = round(Variance*100,2))%>%
  mutate(R2 = round(R2, 2))%>%
  pivot_wider(names_from = "var",values_from = c("Variance", "pval"))%>%
  mutate(Variance_RDA1 = paste0(Variance_RDA1, pval_RDA1),
         Variance_RDA2 = paste0(Variance_RDA2, pval_RDA2))%>%
  select(-pval_RDA1, -pval_RDA2)

write.csv(res.df, file = "plot/res.df.csv")

## PRC figures ----
gg1 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "all", axisx = 1, sel.col = c("#ff7f00","#ff0000"))
gg2 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "all", axisx = 2, sel.col = c("#ff7f00","#ff0000"))
gg3 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "all", axisx = 1, sel.col = c("#197af6","#ff0000"))
gg4 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "all", axisx = 2, sel.col = c("#197af6","#ff0000"))

gg5 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "reduced_treatment", axisx = 1, sel.col = c("#ff0000"))
gg6 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "reduced_treatment", axisx = 2, sel.col = c("#ff0000"))
gg7 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "reduced_treatment", axisx = 1, sel.col = c("#ff0000"))
gg8 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "reduced_treatment", axisx = 2, sel.col = c("#ff0000"))

# Sp scores
ss1 = get.sp.scores.plot(data = PRC, changex = "Distance to origin controls", typex = "all", axisx = 1, sel.col = "grey50")
ss2 = get.sp.scores.plot(data = PRC, changex = "Distance to origin controls", typex = "all", axisx = 2, sel.col = "grey50")
ss3 = get.sp.scores.plot(data = PRC, changex = "Distance to destination controls", typex = "all", axisx = 1, sel.col = "grey50")
ss4 = get.sp.scores.plot(data = PRC, changex = "Distance to destination controls", typex = "all", axisx = 2, sel.col = "grey50")

ss5 = get.sp.scores.plot(data = PRC, changex = "Distance to origin controls", typex = "reduced_treatment", axisx = 1, sel.col = c("#ff0000"))
ss6 = get.sp.scores.plot(data = PRC, changex = "Distance to origin controls", typex = "reduced_treatment", axisx = 2, sel.col = c("#ff0000"))
ss7 = get.sp.scores.plot(data = PRC, changex = "Distance to destination controls", typex = "reduced_treatment", axisx = 1, sel.col = c("#ff0000"))
ss8 = get.sp.scores.plot(data = PRC, changex = "Distance to destination controls", typex = "reduced_treatment", axisx = 2, sel.col = c("#ff0000"))

pdf(file = "plot/experiment_PRCs.pdf", height = 12, width = 15)
for(i in 1:40){
  
  print((gg1[["plot.prc"]][[i]]+gg3[["plot.prc"]][[i]])/(gg2[["plot.prc"]][[i]]+gg4[["plot.prc"]][[i]]))
  print((ss1[["plot.sp"]][[i]]+ss3[["plot.sp"]][[i]])/(ss2[["plot.sp"]][[i]]+ss4[["plot.sp"]][[i]]))
  
  print((gg5[["plot.prc"]][[i]]+gg7[["plot.prc"]][[i]])/(gg6[["plot.prc"]][[i]]+gg8[["plot.prc"]][[i]]))
  print((ss5[["plot.sp"]][[i]]+ss7[["plot.sp"]][[i]])/(ss6[["plot.sp"]][[i]]+ss8[["plot.sp"]][[i]]))
}
dev.off()

## Distances along the years----
slopes = 
  PRC %>%
  filter(type == "all")%>%
  select(Region, originSiteID, destSiteID, change, coef)%>%
  unnest(coef)%>%
  mutate(contrast = abs(contrast))%>%
  left_join(metadata %>% select(Gradient, YearEstablished)%>% distinct()%>% rename(Region = Gradient), by = "Region")%>%
  mutate(Year_0 = as.numeric(Year) - YearEstablished)%>%
  mutate(log_year_0 = log(Year_0))%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))

library(merTools) #predictInterval
tax =
  slopes %>%
  rename(distance = contrast)%>%
  mutate(treatment = ifelse(treatment == "warmed", "warmed", "control"))%>%
  nest(.by = c("axis", "treatment"))%>%
  mutate(lmer.complete = map(data, ~lmer(distance~log_year_0*change*experiment + (1|Region),  data =.,control = lmerControl(optimizer = "bobyqa"))))%>%
  mutate(lmer.reduced = map(data, ~lmer(distance~log_year_0*change + (1|Region),  data =.,control = lmerControl(optimizer = "bobyqa"))))%>%
  mutate(R2 = map(lmer.complete, ~r.squaredGLMM(.)))%>%
  mutate(exp.chi = map2(lmer.complete, lmer.reduced, ~anova(.x, .y)))%>%
  mutate(cont = map(lmer.complete, ~left_join(as.tibble(emtrends(., pairwise ~ change, var = "log_year_0")$emtrends), 
                                              as.tibble(test(emtrends(., pairwise ~ change, var = "log_year_0")$emtrends)))))%>%
  mutate(exp = map(lmer.complete, ~left_join(as.tibble(emtrends(., pairwise ~ change*experiment, var = "log_year_0")$emtrends), 
                                             as.tibble(test(emtrends(., pairwise ~ change*experiment, var = "log_year_0")$emtrends)))))%>%
  
  mutate(pred = map(data, ~{.} %>%
                      dplyr::select(Region, originSiteID, destSiteID, experiment, change) %>%
                      expand_grid(log_year_0 = c(seq(0,2.197225, 0.005), 2.197225))))%>%
  mutate(pred = map2(lmer.complete, pred, ~bind_cols(.y %>% dplyr::select(Region, originSiteID, destSiteID, experiment, change, log_year_0),
                                                     predictInterval(.x, newdata = as.data.frame(.y), n.sims = 1000, level = 0.95, type = "linear.prediction"))))%>%
  mutate(pred = map(pred, ~{.} %>%
                      group_by(change, log_year_0)%>%
                      summarize_at(vars(fit:lwr), mean)%>%
                      ungroup()))
detach("package:merTools", unload = TRUE)
detach("package:arm", unload = TRUE)
detach("package:MASS", unload = TRUE)

pred =
  tax %>%
  select(axis, treatment, pred)%>%
  unnest(pred)%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(treatment = case_when(treatment == "control"&change == "Distance to origin controls"~"Destination controls",
                               treatment == "control"&change == "Distance to destination controls"~"Origin controls",
                               treatment == "warmed"~"warmed"))%>%
  mutate(axis = recode(axis, RDA1 = "Canonical coefficient (PRC axis 1)",
                       RDA2 = "Canonical coefficient (PRC axis 2)"))%>%
  mutate(experiment = "model")%>%
  mutate(Year_0 = exp(log_year_0))

## Get Figures----
pslope=
  slopes %>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(treatment = recode(treatment, destControls = "Destination controls",
                            originControls = "Origin controls"))%>%
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
  select(axis, treatment, cont)%>%
  unnest(cont)%>%
  mutate(change = factor(change, levels = c("Distance to origin controls","Distance to destination controls")))%>%
  mutate(treatment = case_when(treatment == "control"&change == "Distance to origin controls"~"Destination control",
                               treatment == "control"&change == "Distance to destination controls"~"Origin control",
                               treatment == "warmed"~"warmed"))%>%
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = ggplot2::element_text(size  = 12,  hjust = 0))+
  labs(x = "", y = "Change in taxonomic distances over experimental years (slopes)", alpha= "Significant slopes")

pdf("plot/taxonomic_distances.pdf", height = 9, width = 15)

ggarrange(pslope, ptax, widths = c(1.2, 0.9), common.legend = TRUE, legend = "bottom")

dev.off()

res.tax =
  tax %>%
  select(axis, treatment, cont, R2)%>%
  unnest(R2)%>%
  unnest(cont)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 3)))

write.csv(res.tax, file = "plot/res.tax.csv") 

## Functional traits ----
allsp = unique(dd$SpeciesName)
alltraits = alltraits %>% filter(!Region %in% takeoutregion)

ddt =
  alltraits %>%
  ungroup()%>%
  select(SpeciesName, TPlant_Veg_Height_cm, TSLA_cm2_g, TLeaf_Area_cm2, TN_percent, TC_percent, TP_percent, Seed_Mass)%>%
  rename(height = TPlant_Veg_Height_cm,
         SLA = TSLA_cm2_g,
         LA = TLeaf_Area_cm2,
         seed_mass = Seed_Mass,
         N_percent = TN_percent,
         C_percent = TC_percent,
         P_percent = TP_percent)%>%
  mutate(SLA = replace(SLA, SLA == Inf, NA))%>%
  mutate_at(.vars = c("height","SLA","LA","N_percent","C_percent","P_percent","seed_mass"), log) %>%
  mutate_at(.vars = c("height","SLA","LA","N_percent","C_percent","P_percent","seed_mass"), scale, center = TRUE, scale = TRUE)%>%
  mutate_at(.vars = c("height","SLA","LA","N_percent","C_percent","P_percent","seed_mass"), as.numeric)%>%
  distinct()

# Fill in the trait gaps with miceRanger algorithm
nb.NA = ddt %>% summarise_all(~ 100-round(sum(is.na(.))/length(.),2)*100)
print(nb.NA)
matLong = ddt[, c(-1)]
mrMeanMatch = miceRanger(matLong
                         , m = 10
                         , valueSelector = "meanMatch"
                         , returnModels = TRUE
                         , verbose = FALSE)
matLong_pred = as_tibble(completeData(mrMeanMatch)[[1]])

ddt = bind_cols(ddt %>% select(SpeciesName), matLong_pred)
nb.NA = ddt %>% summarise_all(~ 100-round(sum(is.na(.))/length(.),2)*100)
print(nb.NA)

trees = c("Larix decidua","Acer pseudoplatanus", "Pinus","Pinus sylvestris","Betula pubescens",
          "Sorbus aucuparia","Picea abies","Juniperus communis")
ddt = filter(ddt, !(SpeciesName %in% trees))

library(corrplot)
cor.mat = round(cor(ddt %>% dplyr::select(height:seed_mass)),2)
corrplot(cor.mat, type="upper", order="hclust", method = c("number"),
         tl.col="black", tl.srt=45)
library(paran)
ncomp=paran(ddt %>% dplyr::select(height, LA, SLA, seed_mass, N_percent), iterations = 5000, centile = 0, quietly = FALSE, 
            status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
            col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
            file = "", width = 640, height = 640, grdevice = "png", seed = 0)$Retained
detach("package:paran", unload = TRUE)
detach("package:MASS", unload = TRUE)

res.pca = PCA(ddt %>% dplyr::select(height, LA, SLA, seed_mass, N_percent), ncp = ncomp, scale.unit = TRUE)

pdf("plot/species_traits_pca.pdf", height = 6, width = 8)
fviz_pca_biplot(res.pca,
                geom ="point",
                col.var = "black",
                col.ind = "grey50")+
  labs(x = paste0("PC1 (",round(get_eigenvalue(res.pca)[1,2], 2), "%)"),
       y = paste0("PC2 (",round(get_eigenvalue(res.pca)[2,2], 2), "%)"),
       title = "")+
  xlim(-5,5)+
  ylim(-5, 5)+
  coord_equal()
dev.off()

print(res.pca$var$contrib)

ddt =
  ddt %>%
  mutate(plant_size = res.pca$ind$coord[,1],
         resource_acquisition = res.pca$ind$coord[,2])

cwm = 
  ddx %>%
  select(Region, originSiteID, destSiteID, comm)%>%
  mutate(comm_wide = map(comm, ~{.}%>%
                           left_join(ddt %>% select(SpeciesName, plant_size, resource_acquisition) , by = "SpeciesName", keep = FALSE)%>%
                           mutate(plant_size = replace_na(plant_size,0),
                                  resource_acquisition = replace_na(resource_acquisition, 0))%>%
                           mutate(across(plant_size:resource_acquisition, ~ .*Rel_Cover), na.rm=TRUE)%>%
                           group_by(Region, destSiteID, originSiteID, ODT, Year, destPlotID)%>%
                           summarize_at(vars(plant_size:resource_acquisition), sum, na.rm = TRUE) %>%
                           drop_na(plant_size:resource_acquisition)%>%
                           ungroup()%>%
                           select(ODT, Year, destPlotID, plant_size, resource_acquisition)))%>%
  mutate(comm_wide = map(comm_wide, ~{.}%>%
                           filter(ODT %in% c("originControls", "destControls")) %>%
                           group_by(ODT)%>%
                           summarize_at(vars(plant_size:resource_acquisition), mean)))%>%
  select(-comm)%>%
  unnest(comm_wide)%>%
  pivot_wider(names_from = "ODT", values_from = c("plant_size", "resource_acquisition"))%>%
  rename(destPS = plant_size_destControls,
         destRA = resource_acquisition_destControls,
         oriPS = plant_size_originControls,
         oriRA = resource_acquisition_originControls)%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  ungroup()

## Get climate ----
clim = 
  get.clim(sites = slopes %>% 
             select(Region, originSiteID, destSiteID, Year) %>% 
             mutate(Year = as.numeric(Year)) %>%
             distinct(),
           climdata = climdata)%>%
  ungroup()%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))


## Get metadata ----
regmeta = 
  metadata %>%
  filter(Region %in% unique(dd$Region))%>%
  mutate(YearRange = as.numeric(scale(YearRange, center = TRUE, scale = TRUE)))%>%
  mutate(PlotSize = as.numeric(scale(PlotSize, center = TRUE, scale = TRUE)))%>%
  select(Region, destSiteID, YearRange, PlotSize)

# Assess experimental differences ----
reg = tax$exp[[2]] %>%
  mutate(pval = pval(p.value))%>%
  left_join(clim %>% mutate(across(destP:cumsumT, ~as.numeric(scale(., center = TRUE, scale = TRUE)))), by = c("experiment"))%>%
  left_join(cwm, by = c("Region", "originSiteID","destSiteID", "experiment"))%>%
  left_join(regmeta, by = c("Region","destSiteID"))

mf = as.formula("log_year_0.trend ~ 
                PlotSize + YearRange +
                cumsumT +
                destP + destT +
                destPS + destRA + 
                oriPS + oriRA")

regmod = 
  reg %>%
  nest(.by = "change")%>%
  mutate(data = map(data, ~get.mod(data=., mf = mf)))%>%
  unnest(data)

pdf("plot/experimental_effects_slopes.pdf", height = 13, width = 12)
regmod %>%
  filter(rowname != "(Intercept)")%>%
  mutate(var = recode(rowname, PlotSize = "Plot size",
                      YearRange = "Experiment duration",
                      cumsumT = "Experimental warming",
                      destP = "Destination PET",
                      destT = "Destination temperature",
                      oriPS = "Origin CWM (plant size)",
                      oriRA = "Origin CWM (resource acquisition)",
                      destPS = "Destination CWM (plant size)",
                      destRA = "Destination CWM (resource acquisition)"))%>%
  mutate(type = case_when(var %in% c("Plot size", "Experiment duration", "Experimental warming") ~ "Experimental effects",
                          var %in% c("Destination PET","Destination temperature")~"Climate effects",
                          var %in% c("Origin CWM (plant size)","Origin CWM (resource acquisition)",
                                     "Destination CWM (plant size)", "Destination CWM (resource acquisition)")~"Functional effects"))%>%
  mutate(type = factor(type, levels = c("Experimental effects","Climate effects","Functional effects")))%>%
  mutate(var = factor(var, levels = c("Destination CWM (plant size)", "Destination CWM (resource acquisition)",
                                      "Origin CWM (plant size)","Origin CWM (resource acquisition)",
                                      "Destination PET","Destination temperature","Experimental warming",
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
dev.off()

regmod =
  regmod %>%
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

regmod1 = regmod %>% select(change, var, Estimate, `Std. Error`, `t value`,  
                            `2.5 %`, `97.5 %`, `Pr(>|t|)`)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

regmod2 = regmod %>% select(change, R, Radj, F.stat, numdf, dendf, p.value)%>%distinct()%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

write.csv (regmod1, file = "plot/regmod1.csv")
write.csv (regmod2, file = "plot/regmod2.csv")

sp.scores = 
  PRC %>%
  filter(type %in% c("reduced_treatment"))%>%
  select(-Region, -originSiteID, -destSiteID, -ODT, -destPlotID, -Year, -comm_spp, -PRC, -prc.res, -type)%>%
  mutate(coef = map(coef, ~{.}%>%
                      filter(axis == "RDA1")%>%
                      summarize_at("contrast", mean)))%>%
  mutate(sp.scores = map(sp.scores, ~{.}%>%
                           filter(name == "RDA1")))%>%
  unnest(coef)%>%
  unnest(sp.scores)%>%
  select(-name)%>%
  filter(pool != "control")%>%
  filter(Region != "US_Arizona")%>%
  mutate(contrast = as.numeric(sign(contrast)))%>%
  mutate(value = value*contrast, value)%>% # To have the warming effect all the time
  
  left_join(ddt %>% select(SpeciesName, plant_size, resource_acquisition), by = "SpeciesName")%>%
  mutate(plant_size = replace_na(plant_size, 0),
         resource_acquisition = replace_na(resource_acquisition,0))%>%
  nest(.by = c("Region", "originSiteID","destSiteID","pool"))%>%
  mutate(lm = map(data, ~{.} %>%
                    lm(value~change, data =.)))%>%
  mutate(em = map(lm, ~left_join(as.tibble(emmeans(.,c("change"))), as.tibble(test(emmeans(.,c("change")))), by = c("change","emmean", "SE", "df"))))%>%
  mutate(em = map(em, ~{.} %>%
                    rename(p.value.em = p.value)%>%
                    select(change, emmean, lower.CL, upper.CL, p.value.em)))%>%
  unnest(em)%>%
  select(-lm)%>%
  mutate(pval.em = pval(p.value.em))

## Get species weights ----
sp.scores = 
  left_join(
    sp.scores %>%
      select(Region, originSiteID, destSiteID, pool, data)%>%
      mutate(data = map(data, ~{.} %>%
                          group_by(change)%>%
                          summarize_at(vars(plant_size:resource_acquisition), mean)%>%
                          select(change, plant_size, resource_acquisition))) %>%
      unnest(data) %>%
      distinct(),
    sp.scores %>%
      select(Region, originSiteID, destSiteID, pool, change, emmean, lower.CL, upper.CL, p.value.em, pval.em),
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

pdf("plot/species_pool_weights.pdf", height = 8, width = 10)
p1
dev.off()

## Assess experimental differences in species weights ----
regsp = 
  sp.scores %>%
  filter(pool %in% c("colonizing", "overlapping", "strictly_high_elevation"))%>%
  left_join(regmeta, by = c("Region","destSiteID"))%>%
  left_join(clim %>% mutate(across(destP:cumsumT, ~as.numeric(scale(., center = TRUE, scale = TRUE)))), by = c("Region", "originSiteID","destSiteID"))

mf = as.formula("emmean~
                PlotSize + YearRange + 
                cumsumT + destP + destT + 
                plant_size + resource_acquisition")

regspmod = 
  regsp %>%
  nest(.by = c("pool","change"))%>%
  mutate(data = map(data, ~get.mod(data = ., mf = mf)))%>%
  unnest(data)

pdf("plot/experimental_effects_pools.pdf", height = 11, width = 10)
pp1= regspmod %>%
  filter(rowname != "(Intercept)")%>%
  mutate(pool = recode(pool, overlapping = "Overlapping",
                       colonizing = "Colonizing",
                       strictly_high_elevation = "Strictly \nhigh elevation"))%>%
  mutate(var = recode(rowname, PlotSize = "Turf size",
                      YearRange = "Experiment duration",
                      cumsumT = "Experimental warming",
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
pp1
dev.off()

pdf("plot/experimental_effects_pools_scatter.pdf", height = 10, width = 8)

pp2=regsp %>%
  select(Region, originSiteID, destSiteID, pool, change, emmean,
         YearRange, PlotSize, cumsumT, destT, destP, 
         plant_size, resource_acquisition)%>%
  pivot_longer(cols = YearRange:resource_acquisition, names_to = "rowname")%>%
  left_join(regspmod, by = c("pool","change","rowname"))%>%
  filter(rowname != "(Intercept)")%>%
  mutate(pool = recode(pool, overlapping = "Overlapping",
                       colonizing = "Colonizing",
                       strictly_high_elevation = "Strictly \nhigh elevation"))%>%
  mutate(var = recode(rowname, PlotSize = "Turf size",
                      YearRange = "Experiment duration",
                      cumsumT = "Experimental warming",
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
  geom_point(size = 0.8)+
  stat_smooth(fullrange = TRUE, method = "lm", geom = "line", se = FALSE, show.legend = FALSE, linewidth = 1.2)+
  facet_grid(type~change)+
  scale_alpha_manual(values = c(0.2, 0.7, 1))+
  scale_color_manual(values = c( "#ff7f00", "#d25fff","#197af6"))+
  labs(y = "Average species weights", x = "Values of independent factors",
       color = "Species pools",
       alpha = "")+
  guides(color = guide_legend(order = 1, ncol = 1),
         alpha = guide_legend(order = 2, ncol = 1))+
  theme(strip.text = ggplot2::element_text(size  = 14,  hjust = 0))
pp2
dev.off()  

pdf("plot/experimental_effects_pools_merged.pdf", height = 10, width = 16)
ggarrange(pp1, pp2, widths = c(1.1, 0.9), common.legend = TRUE, legend = "bottom")
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
regspmod1 = regsptab %>% select(change, pool, var, Estimate, `Std. Error`, `t value`,  
                                `2.5 %`, `97.5 %`, `Pr(>|t|)`)%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

regspmod2 = regsptab %>% select(change, pool, R, Radj, F.stat, numdf, dendf, p.value)%>%distinct()%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))

write.csv (regspmod1, file = "plot/regspmod1.csv")
write.csv (regspmod2, file = "plot/regspmod2.csv")


# PRC example experiment ----
i = 1
gg1 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "all", axisx = 1:2, sel.col = c("#ff7f00","#ff0000"))
gg4 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "all", axisx = 1:2, sel.col = c("#197af6","#ff0000"))
gg5 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "reduced_treatment", axisx = 1:2, sel.col = c("#ff0000"))
gg8 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "reduced_treatment", axisx = 1:2, sel.col = c("#ff0000"))

sp.ex = PRC %>%
  filter(Region == "NO_Ulvhaugen" & originSiteID == "Ulvhaugen" & destSiteID == "Alrust" & type %in% c("all","reduced_treatment"))%>%
  select(Region, originSiteID, destSiteID, change, type, sp.scores)%>%
  mutate(p = map(sp.scores, ~{.} %>%
                   filter(value > 0.05 | value < -0.05)%>%
                   mutate(pool = ifelse(pool %in% c("control"), "destination control", pool))%>%
                   ggplot(aes(x = name, y = value, color = pool))+
                   geom_hline(yintercept = 0, color = "grey30")+
                   geom_vline(xintercept = c("RDA1", "RDA2"), color = "grey30")+
                   TP_theme()+
                   geom_point()+
                   geom_text_repel(aes(label = SpeciesName), hjust = 1, direction = "y", nudge_x = 0.5, show.legend = FALSE)+
                   scale_color_manual(values = c( "#ff7f00", "grey20","grey50","#d25fff","#197af6"))+
                   labs(y = "Species weights", x = "", color = "")
  ))


pdf("plot/ex_experiment.pdf", height = 15, width = 17)
(gg1$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[1]])/(gg4$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[2]])

(gg5$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[3]])/(gg8$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[4]])
dev.off()

# Paper Figure 2----
library(maps)
library(mapproj)
mapdata = left_join(clim, metadata %>% select(Region, destSiteID, Longitude, Latitude), by = c("Region","destSiteID"))
world_map = map_data("world")

#Create a base plot with gpplot2
p = ggplot() + coord_fixed() +
  xlab("") + ylab("")

#Add map to base plot
base_world_messy = p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                    colour="white", fill="#dddddd") +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  coord_map("ortho",ylim=c(25,180), orientation=c(61, 0, 0))

cleanup = 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = ggplot2::element_blank())

base_world = base_world_messy + cleanup

map_data = 
  base_world +
  geom_point(data=mapdata, 
             aes(x=as.numeric(Longitude), y=as.numeric(Latitude), fill = cumsumT), 
             pch=21, alpha=I(0.9), size = 4) + 
  scale_fill_gradient(low = "grey90", high = "#CD0000")+
  theme(legend.position="bottom", 
        legend.text.align = 0) +
  labs(fill= "Cumulative experimental warming (°C)")

map_data

zoom_to <- c(8.8, 46.3)
zoom_level <- 5
lon_span <- 360 / 2^zoom_level
lat_span <- 180 / 2^zoom_level

lon_bounds <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
lat_bounds <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

p2=
  base_world +
  geom_jitter(data=mapdata, 
              aes(x=as.numeric(Longitude), y=as.numeric(Latitude), fill = cumsumT), 
              pch=21, 
              alpha=I(0.9),
              stroke = 1,
              size = 6) + 
  guides(color = "none", fill = "none", size = guide_legend(title.position = "top"))+
  theme(legend.position="bottom", 
        legend.text.align = 0) + # omit plot title saying 'color'
  scale_fill_gradient(low = "grey90", high = "#CD0000")+
  coord_sf(xlim = lon_bounds, ylim = lat_bounds)+
  scale_size_manual(values = c(1:9)*1.5)
p2
zoom_to <- c(6.2, 60.8)
zoom_level <- 5
lon_span <- 360 / 2^zoom_level
lat_span <- 180 / 2^zoom_level

lon_bounds <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
lat_bounds <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

p3 =
  base_world +
  geom_jitter(data=mapdata, 
              aes(x=as.numeric(Longitude), y=as.numeric(Latitude), fill = cumsumT), 
              pch=21, 
              alpha=I(0.9),
              stroke = 1,
              size = 6) + 
  guides(color = "none", fill = "none", size = guide_legend(title.position = "top"))+
  theme(legend.position="bottom", 
        legend.text.align = 0) + # omit plot title saying 'color'
  scale_fill_gradient(low = "grey90", high = "#CD0000")+
  coord_sf(xlim = lon_bounds, ylim = lat_bounds)+
  scale_size_manual(values = c(1:9)*1.5)
p3

pdf(file = paste0("plot/experiments_map.pdf"), width = 6, height = 6)
map_data
p2
p3
dev.off()

detach("package:mapproj", unload = TRUE)
detach("package:maps", unload = TRUE)

ggclim =
  clim%>%
  mutate(Region = recode(Region, CH_Calanda = "CH_Calanda1"))%>%
  group_by(Region)%>%
  mutate(num = n())%>%
  mutate(Region_txt = paste0(Region, " (",num, ")"))
ppclim = 
  ggclim %>%
  select(Region, Region_txt, originSiteID, destSiteID, destP, destT, oriP, oriT, experiment)%>%
  pivot_longer(cols = destP:oriT)%>%
  mutate(var = str_extract(name, "[TP]"),
         name = str_extract(name, "(ori|dest)"))%>%
  mutate(name = recode(name, dest = "Destination site",
                       ori = "Origin site"))%>%
  pivot_wider(names_from = var, values_from = value)

pdf("plot/experiments.pdf", height = 7, width = 8)
ppclim %>%
  ggplot(aes(`T`, P, color = Region_txt))+
  TP_theme()+
  geom_segment(data = ggclim, aes(x=oriT, xend = destT,
                                  y = oriP, yend = destP), alpha = 0.7, show.legend = TRUE)+
  geom_segment(data = ggclim, aes(x=oriT, xend = destT,
                                  y = oriP, yend = destP),
               arrow = arrow(length=unit(0.30,"cm"), type = "closed"), alpha = 0.7, show.legend = FALSE)+
  geom_point(aes(shape = name),size = 5, alpha = 0.5)+
  xlab("Summer temperature (°C)")+
  ylab("Annual PET (mm)")+
  theme(legend.position = "right")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values = c("#332288","#332288","#52497D",#CH
                                "#117733","#209648","#52B372", #CN
                                "#88CCEE","#5BB1DC","#1E5F7F", #DE
                                "#CC6677","#EC4D68",#FR
                                "#AA4499","#AA4499", #IT
                                "#44AA99","#5D9087","#99DCD1","#4E6964",#NO,
                                "#DDCC77",#SE,
                                "#882255","#652444"#US
  ))+
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))+
  labs(color = "Gradient", shape = "Site", linetype = "Treatment")
dev.off()
# Experiments table ----
exptab =
  left_join(clim, metadata, by = c("Region", "destSiteID"))%>%
  rename(destLon = Longitude,
         destLat = Latitude,
         destElev = Elevation)%>%
  left_join(metadata %>% 
              select(Region, destSiteID, Longitude, Latitude, Elevation) %>%
              rename(originSiteID = destSiteID), by = c("Region","originSiteID"))%>%
  
  rename(oriLon = Longitude,
         oriLat = Latitude,
         oriElev = Elevation)%>%
  left_join(dd %>% 
              group_by(Region, originSiteID, destSiteID, Treatment, Year)%>% 
              summarize(Replicate = length(unique(destPlotID))) %>%
              ungroup()%>%
              select(Region, destSiteID, Replicate) %>%
              group_by(Region, destSiteID)%>% 
              summarize(Replicate = max(Replicate)), by = c("Region","destSiteID"))%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))%>%
  mutate(origin = paste0("(",oriLat, ", ", oriLon,")"))%>%
  mutate(dest = paste0("(",destLat, ",", destLon,")"))%>%
  select(Region, origin, dest, oriElev, destElev, YearEstablished, YearMin, YearRange, PlotSize, Replicate, oriT,destT, oriP, destP)%>%
  rename(Gradient = Region, 
         `Origin coordinates` = origin,
         `Destination coordinates`= dest,
         `Origin elevation` = oriElev,
         `Destination elevation` = destElev,
         `Transplantation year` = YearEstablished,
         `First year of vegetation survey` = YearMin, 
         `Experimental duration` = YearRange,
         `Turf size (m)` = PlotSize,
         `Number of replicates` = Replicate,
         `Origin temperature (°C)` = oriT,
         `Destination temperature(°C)` = destT, 
         `Origin PET(mm)` = oriP,
         `Destination PET(mm)` = destP)
write.csv(exptab, file= "plot/exptab.csv")
