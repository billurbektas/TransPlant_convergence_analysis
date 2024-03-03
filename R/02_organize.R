#### ORGANIZE THE DATA ####
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
  mutate(YearEstablished = ifelse(Gradient == "FR_Lautaret",2016,YearEstablished))%>% # Correction to Lautaret experiment
  mutate(Elevation = ifelse(Gradient == "CH_Calanda" & Elevation == 2800, 2600, Elevation))%>% # Correction to Calanda experiment
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
df = dat %>%
  filter(!Region %in% takeoutregion)%>%
  dplyr::select(-Elevation)%>%
  group_by(Region, Year, originSiteID, originBlockID, destSiteID, destBlockID, destPlotID, Treatment, turfID, UniqueID)%>%
  mutate(Total_Cover = sum(Cover, na.rm = TRUE), Rel_Cover = Cover/Total_Cover)%>%
  ungroup()%>%
  bind_rows(dat %>% filter(Region == "US_Arizona") %>% dplyr::select(-Elevation))%>% # Only relative cover available for Arizona
  filter(!is.na(Rel_Cover)) %>% #4 NAs in DE_Susalps.
  filter(!is.na(SpeciesName))%>%
  filter(!(Region=="DE_Susalps" & Year==2016)) %>% # Not enough data
  filter(!(Region=="CH_Calanda2" & Year==2020))%>% # OriginControl not counted
  filter(!SpeciesName %in% takeoutsp) %>%
  left_join(metadata %>% dplyr::select(Region, destSiteID, YearEstablished, YearMin),
            by = c("Region","destSiteID"))%>%
  mutate(Year_0 = Year-YearEstablished)%>%
  filter(Year_0 != 0)%>% # We do not take pre-transplantation vegetation surveys
  dplyr::select(-(YearEstablished:Year_0))

# Standardize the species names 
cleanedspecies = cleanedspecies %>%
  filter(!(SpeciesName == "Ant alp" & matched_name == "Anthoxanthum odoratum subsp. nipponicum"))%>%
  filter(!(SpeciesName == "Equ pal/pra"))%>%
  filter(!(SpeciesName == "Rum ace subsp lap" & matched_name == "Rumex acetosa subsp. acetosa"))%>%
  filter(!(SpeciesName == "Sal gla" & matched_name == "Salix glauca glauca"))%>%
  filter(!(SpeciesName == "Ver alp" & matched_name == "Veronica alpina subsp. alpina"))

df = df %>%
  left_join(cleanedspecies %>% 
              dplyr::select(Region, destSiteID, SpeciesName, matched_name) %>% 
              distinct(), by = c("SpeciesName", "Region", "destSiteID"))%>%
  mutate(matched_name = ifelse(is.na(matched_name), SpeciesName, matched_name))%>%
  dplyr::select(-SpeciesName)%>%
  rename(SpeciesName = matched_name)%>%
  distinct()

# Unresolved botanical issues -> we assign genus names
# CN_Damxung: Stipa capillacea, Stipa capillata -> Stipa
#             Polygonaceae -> drives the pattern too much and all of a sudden disappears to be taken out. 
# DE_Susalps (GW-FE): Festuca pratensis (in destination control plots)? Elymus repens (in origin control plots)? -> drives the pattern too much and all of a sudden disappears to be taken out.
# DE_Susalps (GW-BT): Festuca pratensis (in destination control plots)? Elymus repens (in origin control plots)? -> drives the pattern too much and all of a sudden disappears to be taken out.
df = 
  df %>%
  mutate(SpeciesName = ifelse(Region == "CN_Damxung" & SpeciesName %in% c("Stipa capillacea","Stipa capillata"), "Stipa", SpeciesName))%>%
  mutate(Rel_Cover = ifelse(Region == "CN_Damxung" & SpeciesName == "Polygonaceae", 0,Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "FE" & destSiteID == "FE" &  SpeciesName == "Festuca pratensis", 0, Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "BT" & destSiteID == "BT" & SpeciesName == "Festuca pratensis", 0, Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "GW" & destSiteID == "GW" & SpeciesName == "Elymus repens", 0, Rel_Cover))%>%
  filter(Rel_Cover != 0)

# Nest the data per experiment
dfx = df %>%
  dplyr::select(Region, originSiteID, destSiteID, Treatment) %>% 
  distinct() %>% 
  filter(Treatment == "Warm") %>% 
  dplyr::select(-Treatment) %>%
  mutate(comm = pmap(.l = list(R = Region, O = originSiteID, D = destSiteID), .f = function(R, O, D){
    bind_rows(
      originControls = df %>% filter(Region == R & originSiteID == O & destSiteID == O & Treatment == "LocalControl"),
      destControls = df %>% filter(Region == R &  originSiteID == D & destSiteID == D & Treatment == "LocalControl"),
      warmed =  df %>% filter(Region == R & originSiteID == O & destSiteID == D & Treatment == "Warm"),
      .id = "ODT")
  })) 

# I take out species that occurred in communities only once over the whole experimental duration.
## This is because we cannot distinguish if there is an botanical identification error or true presence.
### However, I do not change the relative cover as the species contributed to the total cover. 
dfx = 
  dfx %>%
  mutate(comm = purrr::map(comm, ~{.} %>%
                      group_by(ODT, SpeciesName)%>%
                      mutate(n = n())%>%
                      filter(n>1)%>% 
                      ungroup()%>%
                      dplyr::select(-n)))

# Get species pools ----
pools = dfx %>% 
  dplyr::select(Region, originSiteID, destSiteID, comm) %>%
  mutate(low = purrr::map(comm, ~{.} %>% #species pool at low site controls across all years
                     dplyr::select(ODT, SpeciesName) %>% 
                     filter(ODT == "destControls") %>%
                     distinct(.$SpeciesName) %>% 
                     flatten_chr(.)),
         high = purrr::map(comm, ~{.} %>% #species pool of high site controls across all years
                      dplyr::select(ODT, SpeciesName) %>% 
                      filter(ODT =="originControls") %>%
                      distinct(.$SpeciesName) %>% 
                      flatten_chr(.)),
         warmed = purrr::map(comm, ~{.} %>% #species pool of transplanted turfs across all years
                        dplyr::select(ODT, SpeciesName) %>% 
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
  dplyr::select(-comm, -comm_pool)%>%
  group_by(Region, originSiteID, destSiteID)%>%
  mutate_at(vars(low:warmed_colonizer), ~ purrr::map(., length))%>%
  unnest(low:warmed_colonizer)%>%
  mutate_at(vars(high_unique:colonizer), ~(./total)*100)

sum.pools = sum.pools %>%
  dplyr::select(Region, originSiteID, destSiteID, high_unique, overlap, colonizer)%>%
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

pdf(file = here("plot","species_pool_percentages.pdf"), height = 8, width = 8)
pp= ggplot(sum.pools, aes(pools, value))+
  TP_theme()+
  geom_boxplot(aes(color = pools), show.legend = FALSE)+
  geom_jitter(aes(color = pools), show.legend = FALSE)+
  geom_bracket(data = comp, aes(xmin = xmin, xmax = xmax, y.position = c(120, 100, 90), label = paste0("p.value = ",round(p.value, 3))))+
  scale_color_manual( values = c("#d25fff","#197af6","#ff7f00"))+
  labs(x = "Species pools", y = "Percentage of species pool per experiment")
print(pp)
dev.off()

# Get the species pools per experiment per species
pools = pools %>%
  dplyr::select(comm_pool)%>%
  unnest(comm_pool)%>%
  dplyr::select(Region, originSiteID, destSiteID, ODT, SpeciesName, pool)%>%
  distinct()%>%
  mutate(pool=ifelse(is.na(pool), "low", pool))%>%
  dplyr::select(-ODT)%>%
  distinct()

# Get functional traits
## Functional traits ----
allsp = unique(df$SpeciesName)
alltraits = alltraits %>% filter(!Region %in% takeoutregion)

tr =
  alltraits %>%
  ungroup()%>%
  dplyr::select(SpeciesName, TPlant_Veg_Height_cm, TSLA_cm2_g, TLeaf_Area_cm2, TN_percent, TC_percent, TP_percent, Seed_Mass)%>%
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
nb.NA = tr %>% summarise_all(~ 100-round(sum(is.na(.))/length(.),2)*100)
print(nb.NA)
matLong = tr[, c(-1)]
mrMeanMatch = miceRanger(matLong
                         , m = 10
                         , valueSelector = "meanMatch"
                         , returnModels = TRUE
                         , verbose = FALSE)
matLong_pred = as_tibble(completeData(mrMeanMatch)[[1]])

tr = bind_cols(tr %>% dplyr::select(SpeciesName), matLong_pred)
nb.NA = tr %>% summarise_all(~ 100-round(sum(is.na(.))/length(.),2)*100)
print(nb.NA)

trees = c("Larix decidua","Acer pseudoplatanus", "Pinus","Pinus sylvestris","Betula pubescens",
          "Sorbus aucuparia","Picea abies","Juniperus communis")
tr = filter(tr, !(SpeciesName %in% trees))
ncomp=paran(tr %>% dplyr::select(height, LA, SLA, seed_mass, N_percent), iterations = 5000, centile = 0, quietly = FALSE, 
            status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
            col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
            file = "", width = 640, height = 640, grdevice = "png", seed = 0)$Retained

res.pca = PCA(tr %>% dplyr::select(height, LA, SLA, seed_mass, N_percent), ncp = ncomp, scale.unit = TRUE)

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

tr =
  tr %>%
  mutate(plant_size = res.pca$ind$coord[,1],
         resource_acquisition = res.pca$ind$coord[,2])

cwm = 
  dfx %>%
  dplyr::select(Region, originSiteID, destSiteID, comm)%>%
  mutate(comm_wide = purrr::map(comm, ~{.}%>%
                           left_join(tr %>% dplyr::select(SpeciesName, plant_size, resource_acquisition) , by = "SpeciesName", keep = FALSE)%>%
                           mutate(plant_size = replace_na(plant_size,0),
                                  resource_acquisition = replace_na(resource_acquisition, 0))%>%
                           mutate(across(plant_size:resource_acquisition, ~ .*Rel_Cover), na.rm=TRUE)%>%
                           group_by(Region, destSiteID, originSiteID, ODT, Year, destPlotID)%>%
                           summarize_at(vars(plant_size:resource_acquisition), sum, na.rm = TRUE) %>%
                           drop_na(plant_size:resource_acquisition)%>%
                           ungroup()%>%
                           dplyr::select(ODT, Year, destPlotID, plant_size, resource_acquisition)))%>%
  mutate(comm_wide = purrr::map(comm_wide, ~{.}%>%
                           filter(ODT %in% c("originControls", "destControls")) %>%
                           group_by(ODT)%>%
                           summarize_at(vars(plant_size:resource_acquisition), mean)))%>%
  dplyr::select(-comm)%>%
  unnest(comm_wide)%>%
  pivot_wider(names_from = "ODT", values_from = c("plant_size", "resource_acquisition"))%>%
  rename(destPS = plant_size_destControls,
         destRA = resource_acquisition_destControls,
         oriPS = plant_size_originControls,
         oriRA = resource_acquisition_originControls)%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  ungroup()
