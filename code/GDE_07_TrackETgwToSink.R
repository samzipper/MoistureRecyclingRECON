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

# testing - july
#m <- 7
mo <- sprintf("%02d", m)

# subset ETgw to month and calculate ETgw percent
df_etgw$ETgw_mo_prc <- df_etgw[[paste0("ETgw_m3_", mo)]] / df_etgw[[paste0("ET_m3_", mo)]]

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
df_mo <- 
  file.path(cell2cell_directory,
            sprintf("GDE_05_UTrack-CellToCell_%02d.csv", m)) |> 
  read_csv() 
