source(file.path("code", "paths+packages.R"))

gee_out <- 
  read_csv("G:/My Drive/GEEexports/ERA5-CONUS_LANID_PolygonCounts_Annual.csv") |> 
  dplyr::select(-`system:index`, -`.geo`)
sf_era <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp"))

sf_era <- 
  left_join(sf_era, subset(gee_out, year == 2008))

ggplot(sf_era, aes(fill = pct_irrigated)) +
  geom_sf()

(as.numeric(st_area(sf_era))) |> hist()
hist(gee_out$land_area_m2)
max(gee_out$n_land)

mean(st_area(sf_era))
mean(gee_out$land_area_m2)
