#### BUILD PRINCIPAL RESPONSE CURVES ####
# Principal Response Curves----
## Prepare species sites matrices per experiment
dfx = dfx %>% 
  mutate(comm_wide = purrr::map(comm, ~{
    .x %>% dplyr::select(ODT, Year, SpeciesName, Rel_Cover, destPlotID) %>% 
      pivot_wider(names_from = SpeciesName, values_from = Rel_Cover, values_fill = list(Rel_Cover = 0), values_fn = list(Rel_Cover = sum)) %>%
      group_by(ODT)%>%
      complete(Year, destPlotID)%>% # There are cases where plots are missing. Create those plots in the dataset.
      group_by(ODT, Year)%>%
      mutate(across(where(is.numeric), ~replace_na(., mean(., na.rm=TRUE))))%>% #Fill the missing plots with the mean of ODT and Year.
      ungroup()
  })) 

## Get PRCs, species scores and coefficients ----
PRC = bind_rows(get.PRC(data = dfx, sel.ODT = c("warmed","originControls","destControls"), sel.control = "originControls") %>%
                  mutate(change = "Distance to origin controls",
                         type = "all",
                         sel.control = "originControls"),
                get.PRC(data = dfx, sel.ODT = c("warmed","originControls","destControls"), sel.control = "destControls") %>%
                  mutate(change = "Distance to destination controls",
                         type = "all",
                         sel.control = "destControls"),
                get.PRC(data = dfx, sel.ODT = c("warmed","originControls"), sel.control = "originControls") %>%
                  mutate(change = "Distance to origin controls",
                         type = "reduced_treatment",
                         sel.control = "originControls"),
                get.PRC(data = dfx, sel.ODT = c("warmed","destControls"), sel.control = "destControls") %>%
                  mutate(change = "Distance to destination controls",
                         type = "reduced_treatment",
                         sel.control = "destControls"),
                get.PRC(data = dfx, sel.ODT = c("originControls","destControls"), sel.control = "originControls") %>%
                  mutate(change = "Distance to origin controls",
                         type = "reduced_controls",
                         sel.control = "originControls"),
                get.PRC(data = dfx, sel.ODT = c("originControls","destControls"), sel.control = "destControls") %>%
                  mutate(change = "Distance to destination controls",
                         type = "reduced_controls",
                         sel.control = "destControls"))
## Get PRC figures per experiment ----
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

pdf(file = here("plot", "experiment_PRCs.pdf"), height = 12, width = 15)
for(i in 1:40){
  
  print((gg1[["plot.prc"]][[i]]+gg3[["plot.prc"]][[i]])/(gg2[["plot.prc"]][[i]]+gg4[["plot.prc"]][[i]]))
  print((ss1[["plot.sp"]][[i]]+ss3[["plot.sp"]][[i]])/(ss2[["plot.sp"]][[i]]+ss4[["plot.sp"]][[i]]))
  
  print((gg5[["plot.prc"]][[i]]+gg7[["plot.prc"]][[i]])/(gg6[["plot.prc"]][[i]]+gg8[["plot.prc"]][[i]]))
  print((ss5[["plot.sp"]][[i]]+ss7[["plot.sp"]][[i]])/(ss6[["plot.sp"]][[i]]+ss8[["plot.sp"]][[i]]))
}
dev.off()

## RUN permutation tests ----
perm.PRC = get.perm(PRC, set.nperm = set.nperm)
perm.df = 
  perm.PRC %>%
  mutate(perm = purrr::map(perm, ~rownames_to_column(.)))%>%
  unnest(perm)%>%
  rename(var = rowname)%>%
  filter(var %in% c("RDA1","RDA2"))%>%
  mutate(pval = pval(`Pr(>F)`))%>%
  dplyr::select(-sel.control)%>%
  mutate(pval = ifelse(pval == "non-significant","", pval))

## Get PRC tables ----
res.df =
  PRC %>%
  filter(type != "reduced_controls")%>%
  dplyr::select(Region, originSiteID, destSiteID, change, type, prc.res)%>%
  unnest(prc.res)%>%
  left_join(perm.df %>% dplyr::select(-(Df:`Pr(>F)`)), by = c("Region","originSiteID","destSiteID","change","type","var"))%>%
  mutate(Variance = round(Variance*100,2))%>%
  mutate(R2 = round(R2, 2))%>%
  pivot_wider(names_from = "var",values_from = c("Variance", "pval"))%>%
  mutate(Variance_RDA1 = paste0(Variance_RDA1, pval_RDA1),
         Variance_RDA2 = paste0(Variance_RDA2, pval_RDA2))%>%
  dplyr::select(-pval_RDA1, -pval_RDA2)

write.csv(res.df, file = here("output", "res.df.csv"))

