## ETbyLULC_02_SummarizeGEEoutput.R
# output from GEE script `ETbyLULC_01_ETbyCropIrrigation` must be copied over to 
#  the `data/ETbyLULC` folder in the working directory
# script UTrack-Polygon_02_MonthlyERASummary.R must be run already for comparison.

# set up workspace
source(file.path("code", "paths+packages.R"))

# load data
files_yearly <- list.files(file.path("data", "ETbyLULC"))
sf_cells <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_HPA.shp"))
sf_cells$area_m2 <- as.numeric(st_area(sf_cells))

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
            ETirr_prc = ETirr_m3_mean / ETa_m3_mean)

df_cellETannualMean <-
  df_cellETbyYear |> 
  subset(year == 2008) |>  #  issue with NA values - GEE problem? just use 2008 for now
  group_by(cellid) |> 
  summarize(ETa_m3_sum = sum(ET_m3_total),
            ETalt_m3_sum = sum(ETalt_m3_total),
            ETirr_m3_sum = ETa_m3_sum - ETalt_m3_sum,
            ETirr_prc = ETirr_m3_sum / ETa_m3_sum)

# save output
write_csv(df_cellETmonthlyMean, file.path("data", "ETbyLULC", "ETbyLULC_02_SummarizeMonthly.csv"))
write_csv(df_cellETannualMean, file.path("data", "ETbyLULC", "ETbyLULC_02_SummarizeAnnual.csv"))

# load ERA ET
df_ERA <- read_csv(file.path("data", "UTrack-Polygon_02_MonthlyERASummary.csv"))

df_ERA.OpenET <- left_join(df_cellETmonthlyMean, df_ERA, by = c("cellid", "month"))

# join to cells
sf_cellETmonthly <-
  right_join(df_ERA.OpenET, sf_cells, by = "cellid") |> 
  st_as_sf()

sf_cellETannual <- right_join(df_cellETannualMean, sf_cells, by = "cellid") |> 
  st_as_sf()

# ERA comparison
lm_ERA.OpenET <- lm(ETa_m3_mean ~ evaporation_m3, data = df_ERA.OpenET)
ggplot(sf_cellETmonthly, aes(x = 1000*evaporation_m3 / area_m2, 
                            y = 1000*ETa_m3_mean / area_m2, 
                            color = factor(month))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal() +
  scale_x_continuous(name = "RECON-corrected Monthly ERA Evaporation [mm]") +
  scale_y_continuous(name = "OpenET Ensemble Monthly ET [mm]") +
  scale_color_viridis_d(name = "Month", option = "C") +
  labs(title = "Comparison of Evaporation Source by Grid Cell",
       subtitle = paste0("R\u00b2 = ", round(summary(lm_ERA.OpenET)$r.squared, 3)))
ggsave(file.path("figures+tables", "ETbyLULC_02_CompareOpenET.ERA-Scatter.png"),
       width = 110, height = 95, units = "mm")

ggplot(sf_cellETmonthly, aes(fill = 1000*((ETa_m3_mean - evaporation_m3) / area_m2))) +
  facet_wrap(~month, nrow = 2) +
  geom_sf() +
  scale_fill_gradient2(name = "Residual [mm]") +
  labs(title = "Comparison of Evaporation by Grid Cell",
       subtitle = "Residual = OpenET - ERA") +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
ggsave(file.path("figures+tables", "ETbyLULC_02_CompareOpenET.ERA-Map.png"),
       width = 190, height = 130, units = "mm")

ggplot(sf_cellETmonthly, aes(fill = ((ETa_m3_mean - evaporation_m3) / evaporation_m3))) +
  facet_wrap(~month, nrow = 2) +
  geom_sf() +
  scale_fill_gradient2(name = "Residual [%]", labels = scales::percent) +
  labs(title = "Comparison of Evaporation by Grid Cell",
       subtitle = "Residual [%] = (OpenET - ERA) / ERA") +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
ggsave(file.path("figures+tables", "ETbyLULC_02_CompareOpenET.ERA-MapPrc.png"),
       width = 190, height = 130, units = "mm")

(p_ERA.OpenET_scatter + p_ERA.OpenET_map) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom")

# annual plots
ggplot(sf_cellETannual, aes(fill = ETa_m3_sum)) +
  geom_sf() +
  scale_fill_viridis_c(name = "ET [m^3]", direction = -1) +
  theme(legend.position = "bottom")

ggplot(sf_cellETannual, aes(fill = ETirr_prc)) +
  geom_sf() +
  scale_fill_viridis_c(name = "ET from Irrigation [%]", direction = -1, labels = scales::percent) +
  theme(legend.position = "bottom") +
  labs(title = "Annual Average")
ggsave(file.path("figures+tables", "ETbyLULC_02_AnnualETfromIrr-Map.png"),
       width = 4, height = 6, units = "in")

ggplot(sf_cellETannual, aes(x = ETirr_prc)) +
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

ggplot(sf_cellETmonthly, aes(fill = ETirr_prc)) +
  geom_sf() +
  facet_wrap(~month) +
  scale_fill_gradient2(name = "ET from Irrigation [%]", labels = scales::percent)


sf_cellETmonthly |> 
  subset(month %in% seq(6,9)) |> 
  ggplot(aes(fill = ETirr_prc)) +
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
  ggplot(aes(x = ETirr_prc)) +
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
