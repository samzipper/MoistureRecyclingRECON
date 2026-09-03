## Figure_Map-IrrStatus+State.R

source(file.path("code", "paths+packages.R"))

sf_grid_lanid <- st_read(file.path("data", "ERA5-HPA_LANID_AverageIrrigation.shp"))
sf_conus <- st_read(file.path("data", "boundaries", "CONUS-States_TIGRIS.shp"))

# location with CONUS
p_state <-
  ggplot() +
  geom_sf(data = sf_conus) +
  geom_sf(data = sf_grid_lanid, aes(fill = state)) +
  labs(title = "HPA study domain",
       subtitle = "State assigned by largest area within grid cell",
       fill = "State") +
  theme(legend.position = "bottom")

p_irr <- 
  ggplot(sf_grid_lanid, aes(fill = irr_avgPrc)) +
  geom_sf() +
  scale_fill_viridis_c(direction = -1,
                       limits = c(0, max(sf_grid_lanid$irr_avgPrc)),
                       labels = scales::percent) +
    labs(title = "Percent of grid cell irrigated",
         subtitle = "LANID, 2008-2017 average",
         fill = NULL) +
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

(p_state + p_irr +
  plot_layout(nrow = 1, widths = c(1, 0.4))) |> 
  ggsave(file.path("figures+tables", "Figure_Map-IrrStatus+State.png"),
       plot = _, width = 210, height = 120, units = "mm")
