## Figure_Map-AnnualPrecip.R

source(file.path("code", "paths+packages.R"))

# load data
sf_grid_CONUS <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp"))
sf_states <- st_read(file.path("data", "boundaries", "CONUS-States_TIGRIS.shp"))
df_precip <- read_csv(file.path("data", "RECON_01_CONUStotalPrecipitationVolume.csv"))

# join, calculate area, convert to precip_mm
sf_grid_CONUS_precip <-
  left_join(sf_grid_CONUS, df_precip, by = "cellid")
sf_grid_CONUS_precip$area_m2 <- as.numeric(st_area(sf_grid_CONUS_precip))
sf_grid_CONUS_precip$precip_mm <- 
  1000 * (sf_grid_CONUS_precip$precipitation_m3 / sf_grid_CONUS_precip$area_m2)
sf_grid_CONUS_precip$precip_mm_cut <- cut(sf_grid_CONUS_precip$precip_mm, 
                                          breaks = c(floor(min(sf_grid_CONUS_precip$precip_mm)),
                                                     seq(200, 2000, 200), 
                                                     ceiling(max(sf_grid_CONUS_precip$precip_mm))),
                                          dig.lab = 4)

# plot 
ggplot() + 
  geom_sf(data = sf_grid_CONUS_precip, aes(fill = precip_mm_cut), color = NA) +
  geom_sf(data = sf_states, fill = NA) +
  scale_fill_viridis_d(name = "Precipitation [mm]", 
                       option = "G", direction = -1) +
  labs(title = "Mean Annual Precipitation, ERA5")
ggsave(file.path("figures+tables", "Figure_Map-AnnualPrecip.png"),
       width = 190, height = 90, units = "mm")
