## GDE_06_TrackETgwToSink.R
# Track the contribution of ETgw from GDEs to precipitation at each grid cell.

source(file.path("code", "paths+packages.R"))

polygon_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_Global-Land-NoPoles.shp")
etgw_file <- file.path(
  path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_04_ETgwByERA",
  "GDE_04_ETgwByERA_AllCells.csv"
)
cell2cell_directory <- file.path(
  path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_05_UTrack-CellToCell"
)

# load data used in all months
sf_grid <- st_read(polygon_grid_file)
df_etgw <- read_csv(etgw_file)
df_ERA <- read_csv(file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_06_MonthlyERASummary.csv"))

# minimum fraction to retain a source-sink relationship
min_frac <- 1e-5


# loop through months
for (m in seq(1,12)){
  
  # testing - july
  #m <- 7
  mo <- sprintf("%02d", m)
  
  # get ERA data
  df_ERA_mo <- subset(df_ERA, month == m)
  
  # subset ETgw to month and calculate ETgw percent
  df_etgw$ETgw_mo_prc <- df_etgw[[paste0("ETgw_m3_", mo)]] / df_etgw[[paste0("ET_m3_", mo)]]
  sources_month <- unique(subset(df_etgw, ETgw_mo_prc >= min_frac)$cellid)
  
  # # map to inspect
  # left_join(sf_grid, dplyr::select(df_etgw, cellid, ETgw_mo_prc), by = "cellid") |> 
  #   ggplot(aes(fill = ETgw_mo_prc)) +
  #   geom_sf(color = NA) +
  #   scale_fill_viridis_c(direction = -1, labels = scales::percent) +
  #   coord_sf(xlim = c(-130, -80),
  #            ylim = c(20, 50),
  #            expand = FALSE)
  # 
  # left_join(sf_grid, dplyr::select(df_etgw, cellid, ETgw_m3_07), by = "cellid") |> 
  #   ggplot(aes(fill = ETgw_m3_07)) +
  #   geom_sf(color = NA) +
  #   scale_fill_viridis_c(direction = -1) +
  #   coord_sf(xlim = c(-130, -80),
  #            ylim = c(20, 50),
  #            expand = FALSE)
  # 
  # left_join(sf_grid, dplyr::select(df_etgw, cellid, ET_m3_07), by = "cellid") |> 
  #   ggplot(aes(fill = ET_m3_07)) +
  #   geom_sf(color = NA) +
  #   scale_fill_viridis_c(direction = -1) +
  #   coord_sf(xlim = c(-130, -80),
  #            ylim = c(20, 50),
  #            expand = FALSE)
  
  
  # load source-to-sink and join
  df_cell2cell_mo <- 
    file.path(cell2cell_directory,
              sprintf("GDE_05_UTrack-CellToCell_%02d.csv", m)) |> 
    read_csv() |> 
    subset(cellid_source %in% sources_month) |> 
    subset(cellid_sink != -9999)
  
  # create all source-sink pairs for this month
  df_allPairs <- expand.grid(cellid_source = sources_month, cellid_sink = unique(df_cell2cell_mo$cellid_sink))
  
  # summarize evaporation from source cells
  df_mo_evap <- 
    df_cell2cell_mo |>
    # add all source/sink combos
    right_join(df_allPairs, by = c("cellid_source", "cellid_sink")) |> 
    # NA means no source-sink relationship, so set to 0
    replace_na(list(fraction_sourceEvapToSink = 0)) |> 
    # add evaporation data
    left_join(dplyr::select(df_ERA_mo, cellid, evaporation_m3), by = c("cellid_source" = "cellid")) |> 
    rename(evaporation_m3_source = evaporation_m3) |> 
    left_join(dplyr::select(df_etgw, cellid, ETgw_mo_prc), by = c("cellid_source" = "cellid")) |> 
    rename(ETgw_prc_source = ETgw_mo_prc) |> 
    mutate(precipAtSinkFromSource_m3 = fraction_sourceEvapToSink*evaporation_m3_source,
           precipAtSinkFromSourceETgw_m3 = precipAtSinkFromSource_m3*ETgw_prc_source)
  
  # summarize precipitation at each sink cell
  df_mo_sink <- 
    df_mo_evap |> 
    group_by(cellid_sink) |> 
    summarize(precipAtSinkFromETgw_m3 = sum(precipAtSinkFromSourceETgw_m3, na.rm = F)) |> 
    right_join(dplyr::select(df_ERA_mo, cellid, precipitation_m3), by = c("cellid_sink" = "cellid")) |> 
    rename(precipitation_m3_ERA = precipitation_m3) |> 
    mutate(month = m) |> 
    # NA means no source-sink relationship, so set to 0
    replace_na(list(precipAtSinkFromETgw_m3 = 0))
  
  if (m == 1){
    df_sink <- df_mo_sink
  } else {
    df_sink <- bind_rows(df_sink, df_mo_sink)
  }
}


# add spatial information
mo_test <- 7
sf_sink <- 
  df_sink |> 
  subset(month == mo_test) |> 
  right_join(sf_grid, by = c("cellid_sink" = "cellid")) |> 
  st_as_sf()

# plot - percent of precipitation at cell from HPA region
ggplot() +
  geom_sf(data = sf_sink, aes(fill =  precipAtSinkFromETgw_m3 / precipitation_m3_ERA), color = NA) +
  #geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_c(name = "Precip from ETgw [%]", 
                       labels = scales::percent) +
  facet_wrap(~month, nrow = 2) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggplot() +
  geom_sf(data = sf_sink, aes(fill =  precipitation_m3_ERA), color = NA) +
  #geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_c(name = "Precip from ERA [m^3]",
                       trans = "log10") +
  facet_wrap(~month, nrow = 2) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
