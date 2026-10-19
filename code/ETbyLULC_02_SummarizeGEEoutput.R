## ETbyLULC_02_SummarizeGEEoutput.R
# output from GEE script `ETbyLULC_01_ETbyCropIrrigation` must be copied over to 
#  the `data/ETbyLULC` folder in the working directory

# set up workspace
source(file.path("code", "paths+packages.R"))

# load data
files_yearly <- list.files(file.path("data", "ETbyLULC"))
sf_cells <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_HPA.shp"))

f <- 1

y <- as.numeric(str_sub(files_yearly[f], -8, -5))

# load GEE output
df_y <- read_csv(file.path("data", "ETbyLULC", files_yearly[f]))

# calculate cubic meters of ET for each row
df_y$et_m3 <- (df_y$et_mm_mean/1000) * df_y$area_ha * 10000

# for each cell/month, calculate actual ET in m3
df_ym_ETa <-
  df_y |> 
  mutate(et_m3 = (et_mm_mean/1000)*area_ha*10000) |> 
  group_by(cellid, year, month) |> 
  summarize(ET_m3_total = sum(et_m3))

# alternate reality ET: calculate ET in m3 if all irrigated cropland was rainfed of same crop type
df_y_cropETdry <-
  df_y |> 
  filter(irrigated == 0)

df_y_cropArea <- 
  df_y |> 
  group_by(cellid, year, month, cdl) |> 
  summarize(area_ha_total = sum(area_ha))

df_ym_ETalt <- 
  left_join(select(df_y_cropArea, cellid, year, month, cdl, area_ha_total), df_y_cropETdry, by = c("cellid", "year", "month", "cdl")) |> 
  replace_na(list(et_mm_mean = 0)) |> # a couple land cover classes with basically no area don't have any ET, just set to 0  
  mutate(et_m3 = (et_mm_mean/1000)*area_ha_total*10000) |> 
  group_by(cellid, year, month) |> 
  summarize(ETalt_m3_total = sum(et_m3))

# join to cells
sf_cellET <-
  left_join(df_ym_ETa, df_ym_ETalt, by = c("cellid", "year", "month")) |> 
  right_join(sf_cells, by = "cellid") |> 
  st_as_sf()

ggplot(sf_cellET, aes(fill = ET_m3_total)) +
  geom_sf() +
  facet_wrap(~month) +
  scale_fill_viridis_c(name = "ET [m^3]", direction = -1)

ggplot(sf_cellET, aes(fill = (ETalt_m3_total - ET_m3_total)/ET_m3_total)) +
  geom_sf() +
  facet_wrap(~month) +
  scale_fill_gradient2(name = "ET Loss [%]", labels = scales::percent)
