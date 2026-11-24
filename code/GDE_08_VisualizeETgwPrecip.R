## GDE_08_VisualizeETgwPrecip.R


source(file.path("code", "paths+packages.R"))

polygon_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_Global-Land-NoPoles.shp")
path_trackETgw <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_07_TrackETgw")

# load data used in all months
sf_grid <- st_read(polygon_grid_file)
sf_grid$area_m2 <- as.numeric(st_area(sf_grid))

# loop through months
for (m in seq(1,12)){
  
  # testing - july
  #m <- 7
  mo <- sprintf("%02d", m)
  df_mo_sink <- 
    read_csv(file.path(path_trackETgw, paste0("GDE_07_TrackETgwToSink_", mo, ".csv"))) |> 
    dplyr::select(cellid = cellid_sink, precipAtSinkFromETgw_m3, precipitation_m3_ERA, month)
  
  if (m == 1){
    df_sink <- df_mo_sink
  } else {
    df_sink <- bind_rows(df_sink, df_mo_sink)
  }
}

# clean up df_sink weird values
sf_sink <- 
  df_sink |> 
  right_join(sf_grid, by = "cellid") |> 
  st_as_sf() |> 
  mutate(precipAtSinkFromETgw_mm = 1000*precipAtSinkFromETgw_m3 / area_m2,
         precipitation_mm_ERA = 1000*precipitation_m3_ERA / area_m2)

# minimum depth to consider in mm; everything below this, set to 0
depth_tol_mm <- 1
sf_sink$precipAtSinkFromETgw_mm[sf_sink$precipAtSinkFromETgw_mm < depth_tol_mm] <- 0
sf_sink$precipitation_mm_ERA[sf_sink$precipitation_mm_ERA < depth_tol_mm] <- 0
sf_sink$precipAtSinkFromETgw_prc <- sf_sink$precipAtSinkFromETgw_mm / sf_sink$precipitation_mm_ERA
sf_sink$precipAtSinkFromETgw_prc[sf_sink$precipitation_mm_ERA < depth_tol_mm] <- NA

hist(sf_sink$precipAtSinkFromETgw_prc)

# add spatial information
mo_test <- seq(1, 12)
sf_plot <- 
  sf_sink |> 
  subset(month %in% mo_test)

# plot - percent of precipitation at cell from dryland GDE ETgw region
ggplot() +
  geom_sf(data = sf_plot, 
          aes(fill =  cut(precipAtSinkFromETgw_prc, c(0, 0.01, 0.05, 0.10, 0.25, 0.5), include.lowest = T)), 
          color = NA) +
  #geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_d(name = "Precip from ETgw [%]", direction = -1) +
  facet_wrap(~month) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggsave(file.path("figures+tables", "GDE_07_MapPrecipFromETgw-prc.png"),
       width = 260, height = 120, units = "mm")

ggplot() +
  geom_sf(data = sf_plot, 
          aes(fill =  precipAtSinkFromETgw_mm), 
          color = NA) +
  #geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_c(name = "Precip from ETgw [mm]", direction = -1) +
  facet_wrap(~month) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggsave(file.path("figures+tables", "GDE_07_MapPrecipFromETgw-mm.png"),
       width = 260, height = 120, units = "mm")

ggplot() +
  geom_sf(data = sf_plot, 
          aes(fill =  precipAtSinkFromETgw_prc), 
          color = NA) +
  #geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_b(name = "Precip from ETgw [%]", direction = -1) +
  facet_wrap(~month, nrow = 2) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggplot() +
  geom_sf(data = sf_plot, aes(fill =  precipitation_m3_ERA), color = NA) +
  #geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_c(name = "Precip from ERA [m^3]") +
  facet_wrap(~month, nrow = 2) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
