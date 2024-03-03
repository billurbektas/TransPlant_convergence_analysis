# PRC ANALYSIS ----
get.PRC = function(data, sel.ODT, sel.control){
data =
  data %>%
  mutate(comm_meta = purrr::map(comm_wide, ~{.} %>% 
                               filter(ODT %in% sel.ODT) %>%
                               dplyr::select(ODT, destPlotID, Year)))%>%
  mutate(comm_spp = purrr::map(comm_wide, ~{.} %>% 
                              filter(ODT %in% sel.ODT) %>%
                              dplyr::select(-ODT, -destPlotID, -Year)))%>%
  dplyr::select(-comm_wide, -comm)%>%
  mutate(comm_spp = purrr::map(comm_spp, ~ sqrt(.x)))%>%
  dplyr::select(Region, originSiteID,destSiteID, comm_meta, comm_spp)%>%
  unnest_wider(comm_meta)%>%
  # Get the principal response curves to determine the change from the selected controls
   mutate(PRC = pmap(.l = list(ODT=ODT, Year=Year, destPlotID = destPlotID, comm_spp=comm_spp), 
                      .f = function(ODT, Year, destPlotID, comm_spp){
                            ODT = relevel(as.factor(ODT), ref = sel.control)
                            Year = as.factor(Year)
                            prc(response = comm_spp, treatment = ODT, time = Year)}
                      ))%>%
   mutate(prc.res = pmap(.l = list(ODT=ODT, Year=Year, destPlotID = destPlotID, comm_spp = comm_spp, PRC = PRC), 
                         .f = function(ODT, Year, destPlotID, comm_spp, PRC){
                              res = tibble(var =c("Conditional","Constrained","Unconstrained", "RDA1","RDA2"),
                                           Variance = c(PRC$pCCA$tot.chi/PRC$tot.chi, #Conditional
                                                        PRC$CCA$tot.chi/PRC$tot.chi, #Constrained
                                                        PRC$CA$tot.chi/PRC$tot.chi, #Unconstrained
                                                        PRC$CCA$eig[1]/PRC$CCA$tot.chi, #RDA1
                                                        PRC$CCA$eig[2]/PRC$CCA$tot.chi #RDA2
                                                        ))%>%
                                mutate(R2 = RsquareAdj(PRC)$adj.r.squared)
                              return(res)}
                        ))%>%
  mutate(coef = pmap(.l = list(ODT=ODT, Year=Year, destPlotID = destPlotID, comm_spp = comm_spp, PRC = PRC),
                         .f = function(ODT, Year, destPlotID, comm_spp, PRC){
                               ODT = relevel(as.factor(ODT), ref = sel.control)
                               Year = as.factor(Year)
                               # anova can only read the object from the global environment
                               assign("ODT", ODT, envir = .GlobalEnv)
                               assign("Year", Year, envir = .GlobalEnv)
                               assign("destPlotID", destPlotID, envir = .GlobalEnv)
                               assign("comm_spp", comm_spp, envir = .GlobalEnv)

                               res = as.data.frame(summary(PRC,  axis = 1, scaling = "symmetric")$coefficients) %>%
                                 rownames_to_column(var = "treatment")%>%
                                 pivot_longer(cols = where(is.numeric), names_to = "Year", values_to = "contrast")%>%
                                 mutate(axis = "RDA1")%>%
                                 bind_rows(as.data.frame(summary(PRC, axis = 2)$coefficients) %>%
                                             rownames_to_column(var = "treatment")%>%
                                             pivot_longer(cols = where(is.numeric), names_to = "Year", values_to = "contrast")%>%
                                             mutate(axis = "RDA2"))
                               return(res)}
                         ))%>%
  mutate(sp.scores  = pmap(.l = list(Region = Region, originSiteID = originSiteID, destSiteID = destSiteID,
                                            ODT=ODT, Year=Year, destPlotID = destPlotID, comm_spp = comm_spp, PRC = PRC),
                              .f = function(Region, originSiteID, destSiteID, ODT, Year, destPlotID, comm_spp, PRC){
                                    ODT = relevel(as.factor(ODT), ref = sel.control)
                                    Year = as.factor(Year)
                                    # anova can only read the object from the global environment
                                    assign("ODT", ODT, envir = .GlobalEnv)
                                    assign("Year", Year, envir = .GlobalEnv)
                                    assign("destPlotID", destPlotID, envir = .GlobalEnv)
                                    assign("comm_spp", comm_spp, envir = .GlobalEnv)

                                    res =
                                      as.data.frame(scores(PRC, scaling = "symmetric", display = "species"))%>%
                                      rownames_to_column(var = "SpeciesName")%>%
                                      mutate(Region = Region, originSiteID = originSiteID, destSiteID = destSiteID)%>%
                                      left_join(pools, by = c("Region","originSiteID","destSiteID", "SpeciesName"))%>%
                                      mutate(pool = ifelse(is.na(pool), "control", pool))%>%
                                      pivot_longer(cols = RDA1:RDA2)
                                    return(res)}
  ))
return(data)
}

get.perm = function(data, set.nperm){
  data =
  data %>%
  mutate(perm = pmap(.l = list(ODT=ODT, Year=Year, destPlotID = destPlotID, comm_spp = comm_spp, PRC = PRC, sel.control = sel.control),
                     .f = function(ODT, Year, destPlotID, comm_spp, PRC, sel.control){
                       ODT = relevel(as.factor(ODT), ref = sel.control)
                       Year = as.factor(Year)
                       # anova can only read the object from the global environment
                       assign("ODT", ODT, envir = .GlobalEnv)
                       assign("Year", Year, envir = .GlobalEnv)
                       assign("destPlotID", destPlotID, envir = .GlobalEnv)
                       assign("comm_spp", comm_spp, envir = .GlobalEnv)
                       ctrl = how(plots = Plots(strata = destPlotID, type = "free"),
                                  within = Within(type = "none"), nperm = set.nperm)
                       res = anova(PRC, permutations = ctrl, by = "axis")
                       return(res)}
  ))
  data = data %>% dplyr::select(Region, originSiteID, destSiteID, change, type, sel.control, perm)
  return(data)
}

# PLOT PRC ----
get.PRC.plot = function(data, changex, typex, axisx, sel.col){
  return(
  data %>%
    filter(change == changex & type == typex)%>%
    mutate(plot.prc = pmap(.l = list(Region = Region, originSiteID = originSiteID, destSiteID = destSiteID, coef = coef),
                         .f = function(Region, originSiteID, destSiteID, coef){
                           coef %>%
                             filter(axis %in% paste0("RDA",axisx))%>%
                             mutate(treatment = recode(treatment, warmed = "warmed",
                                                       destControls = "Destination controls",
                                                       originControls = "Origin controls"))%>%
                             ggplot(aes(Year, contrast, color = treatment, linetype = axis))+
                             TP_theme()+
                             geom_hline(yintercept = 0, color = "grey50")+
                             geom_point()+
                             geom_line(aes(group = interaction(treatment, axis)))+
                             scale_color_manual(values = sel.col)+
                             labs(x="Experimental years",y=paste0("Canonical coefficient (PRC axis ",axisx,")"), color ="",
                                  title = paste0(Region, " (Origin: ", originSiteID, " Destination: ", destSiteID,")" ),
                                  subtitle = changex)+
                             ylim(-1,1)}))%>%
  dplyr::select(Region, originSiteID, destSiteID, plot.prc)
)
}

get.sp.scores.plot = function(data, changex, typex, axisx, sel.col){
  return(
    data %>%
      filter(change == changex & type == typex)%>%
      mutate(plot.sp = pmap(.l = list(Region = Region, originSiteID = originSiteID, destSiteID = destSiteID, sp.scores = sp.scores),
                             .f = function(Region, originSiteID, destSiteID, sp.scores){
                               sp.scores %>%
                                 rename(axis = name)%>%
                                 filter(axis == paste0("RDA",axisx))%>%
                                 ggplot(aes(pool, value))+
                                 TP_theme()+
                                 geom_hline(yintercept = 0, color = "grey50")+
                                 geom_boxplot(color = sel.col)+
                                 labs(x="Species pools",y=paste0("Species weights (PRC axis ",axisx,")"), color ="",
                                      title = paste0(Region, " (Origin: ", originSiteID, " Destination: ", destSiteID,")" ),
                                      subtitle = changex)+
                                 ylim(-1,1)}))%>%
      dplyr::select(Region, originSiteID, destSiteID, plot.sp)
  )
}

# LM----
get.mod = function(data, mf){
  m1 = lm(mf, data)
  m1sum = as.data.frame(cbind(summary(m1)$coefficients, confint(m1)))%>%
    rownames_to_column()%>%
    as.tibble()%>%
    mutate(pval = pval(`Pr(>|t|)`))%>%
    mutate(R = summary(m1)$r.squared, 
           Radj = summary(m1)$adj.r.squared,
           F.stat = summary(m1)$fstatistic[1],
           numdf = summary(m1)$fstatistic[2],
           dendf = summary(m1)$fstatistic[3],
           p.value = pf(F.stat, numdf, dendf)
           
    )
  return(m1sum)
}

# ORGANIZE climate data -----
get.clim = function(sites, climdata){
  clim = 
    climdata%>%
    ungroup()%>%
    filter(var != "precipitation")%>%
    dplyr::select(gradient, destSiteID, year, month, var, value)%>%
    filter(!(var == "temperature" & !month %in% c(6:8)))%>%
    pivot_wider(values_from = value, names_from = var)%>%
    group_by(gradient, destSiteID, year)%>%
    summarize(temperature = mean(temperature, na.rm = TRUE),
              pet = sum(pet, na.rm = TRUE))%>%
    rename(Year = year,
           Region = gradient,
           site = destSiteID)%>%
    ungroup()
  
  clim = 
    bind_rows(clim %>% right_join(sites %>% 
                                    ungroup() %>% 
                                    dplyr::select(Region, destSiteID, Year) %>% 
                                    rename(site = destSiteID)%>%
                                    distinct, 
                                  by = c("Region","site", "Year")),
              clim %>% right_join(sites %>% 
                                    ungroup() %>% 
                                    dplyr::select(Region, originSiteID, Year) %>% 
                                    rename(site = originSiteID)%>%
                                    distinct, 
                                  by = c("Region","site", "Year")))%>%
    distinct()%>%
    na.omit()%>% #We do not have climate data for 2020-2021
    ungroup()
  
  destclim =
    clim %>%
    rename(destSiteID = site)%>%
    right_join(sites %>% 
                 ungroup() %>% 
                 dplyr::select(Region, destSiteID, Year) %>% 
                 distinct, 
               by = c("Region","destSiteID", "Year"))%>%
    rename(dest_pet = pet,
           dest_temperature = temperature)
  
  oriclim =
    clim %>%
    rename(originSiteID = site)%>%
    right_join(sites %>% 
                 ungroup() %>% 
                 dplyr::select(Region, originSiteID, Year) %>% 
                 distinct, 
               by = c("Region","originSiteID", "Year"))%>%
    rename(ori_pet = pet,
           ori_temperature = temperature)
  
  distclim =
    sites %>% 
    ungroup() %>% 
    dplyr::select(Region, destSiteID, originSiteID, Year) %>% 
    distinct()%>%
    left_join(destclim, by = c("Region", "destSiteID", "Year"))
  
  distclim =
    distclim %>%
    left_join(oriclim %>% dplyr::select(Region, originSiteID, Year, 
                                 ori_pet, ori_temperature), by = c("Region", "Year","originSiteID"))
  
  distclim = 
    distclim %>%
    mutate(diff_temperature = dest_temperature-ori_temperature,
           diff_pet = dest_pet-ori_pet)
  
  distclim =
    distclim %>%
    group_by(Region, originSiteID, destSiteID)%>%
    mutate(cumsum_temp = cumsum(diff_temperature),
           cumsum_precip = cumsum(diff_pet))
  
  distclim =
    distclim%>%
    na.omit()%>%
    group_by(Region, originSiteID, destSiteID)%>%
    summarize(destP = mean(dest_pet),
              destT = mean(dest_temperature),
              oriP = mean(ori_pet),
              oriT = mean(ori_temperature),
              diffP = mean(diff_pet),
              diffT = mean(diff_temperature),
              cumsumP = max(cumsum_precip),
              cumsumT = max(cumsum_temp))
  return(distclim)
}

# UTILS ----
setRowNames = function(x){
  rownames(x)= x$pool
  return(x)
}
namedist = function(x){
  nameDist = c()
  x = as.vector(x)
  for(i in 1:(length(x)-1)){
    xc = paste0(x[i], sep = "-",x[(i+1):length(x)])
    nameDist = c(nameDist, xc)
  }
  return(nameDist)
}
pval = function(x){
  return(ifelse(x<0.001, "***", ifelse(x<0.01, "**", ifelse(x<0.05, "*", "non-significant"))))
}

