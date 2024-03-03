# PRC example experiment ----
i = 1
gg1 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "all", axisx = 1:2, sel.col = c("#ff7f00","#ff0000"))
gg4 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "all", axisx = 1:2, sel.col = c("#197af6","#ff0000"))
gg5 = get.PRC.plot(data = PRC, changex = "Distance to origin controls", typex = "reduced_treatment", axisx = 1:2, sel.col = c("#ff0000"))
gg8 = get.PRC.plot(data = PRC, changex = "Distance to destination controls", typex = "reduced_treatment", axisx = 1:2, sel.col = c("#ff0000"))

sp.ex = PRC %>%
  filter(Region == "NO_Ulvhaugen" & originSiteID == "Ulvhaugen" & destSiteID == "Alrust" & type %in% c("all","reduced_treatment"))%>%
  dplyr::select(Region, originSiteID, destSiteID, change, type, sp.scores)%>%
  mutate(p = purrr::map(sp.scores, ~{.} %>%
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


pdf(here("plot","ex_experiment.pdf"), height = 15, width = 17)
pp = (gg1$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[1]])/(gg4$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[2]])
print(pp)
pp = (gg5$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[3]])/(gg8$plot.prc[[i]]+ylab("Canonical coefficients") + sp.ex$p[[4]])
print(pp)
dev.off()

# Paper Figure 2----

mapdata = left_join(clim, metadata %>% dplyr::select(Region, destSiteID, Longitude, Latitude), by = c("Region","destSiteID"))
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
             aes(x=as.numeric(Longitude), y=as.numeric(Latitude), fill = diffT), 
             pch=21, alpha=I(0.9), size = 4) + 
  scale_fill_gradient(low = "grey90", high = "#CD0000")+
  theme(legend.position="bottom", 
        legend.text.align = 0) +
  labs(fill= "Experimental warming (°C)")

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
              aes(x=as.numeric(Longitude), y=as.numeric(Latitude), fill = diffT), 
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
              aes(x=as.numeric(Longitude), y=as.numeric(Latitude), fill = diffT), 
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

pdf(file = paste0(here("plot","experiment_maps.pdf")), width = 6, height = 6)
print(map_data)
print(p2)
print(p3)
dev.off()

ggclim =
  clim%>%
  mutate(Region = recode(Region, CH_Calanda = "CH_Calanda1"))%>%
  group_by(Region)%>%
  mutate(num = n())%>%
  mutate(Region_txt = paste0(Region, " (",num, ")"))
ppclim = 
  ggclim %>%
  dplyr::select(Region, Region_txt, originSiteID, destSiteID, destP, destT, oriP, oriT, experiment)%>%
  pivot_longer(cols = destP:oriT)%>%
  mutate(var = str_extract(name, "[TP]"),
         name = str_extract(name, "(ori|dest)"))%>%
  mutate(name = recode(name, dest = "Destination site",
                       ori = "Origin site"))%>%
  pivot_wider(names_from = var, values_from = value)

pdf(here("plot","experiments_clim.pdf"), height = 7, width = 8)
pp = ppclim %>%
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
print(pp)
dev.off()
# Experiments table ----
exptab =
  left_join(clim, metadata, by = c("Region", "destSiteID"))%>%
  rename(destLon = Longitude,
         destLat = Latitude,
         destElev = Elevation)%>%
  left_join(metadata %>% 
              dplyr::select(Region, destSiteID, Longitude, Latitude, Elevation) %>%
              rename(originSiteID = destSiteID), by = c("Region","originSiteID"))%>%
  
  rename(oriLon = Longitude,
         oriLat = Latitude,
         oriElev = Elevation)%>%
  left_join(df %>% 
              group_by(Region, originSiteID, destSiteID, Treatment, Year)%>% 
              summarize(Replicate = length(unique(destPlotID))) %>%
              ungroup()%>%
              dplyr::select(Region, destSiteID, Replicate) %>%
              group_by(Region, destSiteID)%>% 
              summarize(Replicate = max(Replicate)), by = c("Region","destSiteID"))%>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), NA, round(., 2)))%>%
  mutate(origin = paste0("(",oriLat, ", ", oriLon,")"))%>%
  mutate(dest = paste0("(",destLat, ",", destLon,")"))%>%
  dplyr::select(Region, destSiteID, originSiteID, origin, dest, oriElev, destElev, YearEstablished, YearMin, YearRange, PlotSize, Replicate, oriT,destT, oriP, destP)%>%
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
write.csv(exptab, file= here("output", "exptab.csv"))
