## ETbyLULC_02_SummarizeGEEoutput.R
# output from GEE script `ETbyLULC_01_ETbyCropIrrigation` must be copied over to 
#  the `data/ETbyLULC` folder in the working directory

# set up workspace
source(file.path("code", "paths+packages.R"))

# load data
files_yearly <- list.files(file.path("data", "ETbyLULC"))
sf_cells <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_HPA.shp"))

for (f in 1:length(files_yearly)){
  #f <- 1
  
  y <- as.numeric(str_sub(files_yearly[f], -8, -5))
  
  # load GEE output
  df_y <- read_csv(file.path("data", "ETbyLULC", files_yearly[f]))
  
  # calculate cubic meters of ET for each row (et_mm_mean is in mm)
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
  
  df_ym_ETboth <- 
    left_join(df_ym_ETa, df_ym_ETalt, by = c("cellid", "year", "month"))
  
  if (f == 1){
    df_cellETbyYear <- df_ym_ETboth
  } else {
    df_cellETbyYear <- bind_rows(df_cellETbyYear, df_ym_ETboth)
  }
  
}

# summarize to mean monthly and annual
df_cellETmonthlyMean <-
  df_cellETbyYear |> 
  subset(year == 2008) |>  #  issue with NA values - GEE problem? just use 2008 for now
  group_by(cellid, month) |> 
  summarize(ETa_m3_mean = mean(ET_m3_total),
            ETalt_m3_mean = mean(ETalt_m3_total),
            ETirr_m3_mean = ETa_m3_mean - ETalt_m3_mean,
            ETirr_m3_prc = ETirr_m3_mean / ETa_m3_mean)

df_cellETannualMean <-
  df_cellETbyYear |> 
  subset(year == 2008) |>  #  issue with NA values - GEE problem? just use 2008 for now
  group_by(cellid) |> 
  summarize(ETa_m3_sum = sum(ET_m3_total),
            ETalt_m3_sum = sum(ETalt_m3_total),
            ETirr_m3_sum = ETa_m3_sum - ETalt_m3_sum,
            ETirr_m3_prc = ETirr_m3_sum / ETa_m3_sum)


# join to cells
sf_cellETmonthly <-
  right_join(df_cellETmean, sf_cells, by = "cellid") |> 
  st_as_sf()

sf_cellETannual <- right_join(df_cellETannualMean, sf_cells, by = "cellid") |> 
  st_as_sf()

# annual plots
ggplot(sf_cellETannual, aes(fill = ETa_m3_sum)) +
  geom_sf() +
  scale_fill_viridis_c(name = "ET [m^3]", direction = -1) +
  theme(legend.position = "bottom")

ggplot(sf_cellETannual, aes(fill = ETirr_m3_prc)) +
  geom_sf() +
  scale_fill_viridis_c(name = "ET from Irrigation [%]", direction = -1, labels = scales::percent) +
  theme(legend.position = "bottom") +
  labs(title = "Annual Average")
ggsave(file.path("figures+tables", "ETbyLULC_02_AnnualETfromIrr-Map.png"),
       width = 4, height = 6, units = "in")

ggplot(sf_cellETannual, aes(x = ETirr_m3_prc)) +
  geom_histogram(breaks = seq(-0.1, 0.3, 0.01)) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Annual Average",
       y = "# of Grid Cells",
       x = "ET from Irrigation [%]")
ggsave(file.path("figures+tables", "ETbyLULC_02_AnnualETfromIrr-Hist.png"),
       width = 4, height = 3, units = "in")

sf_cellETannual |> 
  summary()

# monthly plots
ggplot(sf_cellETmonthly, aes(fill = ETa_m3_mean)) +
  geom_sf() +
  facet_wrap(~month) +
  scale_fill_viridis_c(name = "ET [m^3]", direction = -1)

ggplot(sf_cellETmonthly, aes(fill = ETirr_m3_prc)) +
  geom_sf() +
  facet_wrap(~month) +
  scale_fill_gradient2(name = "ET from Irrigation [%]", labels = scales::percent)


sf_cellETmonthly |> 
  subset(month %in% seq(6,9)) |> 
  ggplot(aes(fill = ETirr_m3_prc)) +
  geom_sf() +
  facet_wrap(~month) +
  scale_fill_viridis_c(name = "ET from Irrigation [%]", 
                       breaks = seq(0, 0.2, 0.1),
                       labels = scales::percent, 
                       direction = -1) +
  theme(legend.position = "bottom") +
  labs(title = "Monthly Average, Peak Irrigation Season")

ggsave(file.path("figures+tables", "ETbyLULC_02_MonthlyETfromIrr-GrowingSeason-Map.png"),
       width = 6, height = 9, units = "in")


sf_cellETmonthly |> 
  subset(month %in% seq(6,9)) |> 
  ggplot(aes(x = ETirr_m3_prc)) +
  geom_histogram(breaks = seq(-0.1, 0.3, 0.01)) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Monthly Average, Peak Irrigation Season",
       y = "# of Grid Cells",
       x = "ET from Irrigation [%]")
ggsave(file.path("figures+tables", "ETbyLULC_02_MonthlyETfromIrr-GrowingSeason-Hist.png"),
       width = 4, height = 3, units = "in")

sf_cellETmonthly |> 
  subset(month %in% seq(6,9)) |> 
  summary()
