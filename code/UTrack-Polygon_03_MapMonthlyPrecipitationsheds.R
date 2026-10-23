## Utrack-Polygon_03_MapMonthlyPrecipitationsheds.R
# requires output from UTrack-Polygon_01 and _02

source(file.path("code", "paths+packages.R"))

## inputs and outputs
# Source IDs come from the HPA grid; sink IDs use the CONUS grid or -9999 outside CONUS.
polygon_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_HPA.shp")
conus_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp")
input_directory <- file.path(path_data, "Processed", "UTrack-Polygon_01_MonthlySourceToSink")

# load grids and state boundary
grid_hpa <- st_read(polygon_grid_file)
grid_conus <- st_read(conus_grid_file)
grid_conus$area_m2 <- as.numeric(st_area(grid_conus))
state_boundaries <- st_read(file.path("data", "boundaries", "CONUS-States_TIGRIS.shp"))

# create table with all source-sink pairs
df_allPairs <- expand.grid(cellid_source = grid_hpa$cellid, cellid_sink = grid_conus$cellid)

# load OpenET and ERA data
df_OpenET <- read_csv(file.path("data", "ETbyLULC", "ETbyLULC_02_SummarizeMonthly.csv"))
df_ERA <- read_csv(file.path("data", "UTrack-Polygon_02_MonthlyERASummary.csv"))

# loop through months
for (m in seq(1,12)){
  # subset ERA and OpenET to month
  df_ERA_mo <- subset(df_ERA, month == m)
  df_OpenET_mo <- subset(df_OpenET, month == m)
  
  # load source-to-sink and join
  df_mo <- 
    file.path(input_directory,
              sprintf("UTrack-Polygon_01_MonthlySourceToSink_%02d.csv", m)) |> 
    read_csv() 
  
  df_mo_evap <- 
    df_mo |>
    # add all source/sink combos
    right_join(df_allPairs, by = c("cellid_source", "cellid_sink")) |> 
    # NA means no source-sink relationship, so set to 0
    replace_na(list(fraction_sourceEvapToSink = 0)) |> 
    # add evaporation data
    left_join(dplyr::select(df_ERA_mo, cellid, evaporation_m3), by = c("cellid_source" = "cellid")) |> 
    rename(evaporation_m3_source = evaporation_m3) |> 
    left_join(dplyr::select(df_OpenET_mo, cellid, ETirr_prc), by = c("cellid_source" = "cellid")) |> 
    rename(ETirr_prc_source = ETirr_prc) |> 
    mutate(precipAtSinkFromSource_m3 = fraction_sourceEvapToSink*evaporation_m3_source,
           precipAtSinkFromSource_m3_noIrr = fraction_sourceEvapToSink*evaporation_m3_source*(1-ETirr_prc_source))
  
  # summarize precipitation at each sink cell
  df_mo_sink <- 
    df_mo_evap |> 
    subset(cellid_sink != -9999) |> 
    group_by(cellid_sink) |> 
    # FOR NOW- remove NAs. but need to fix so there aren't any NAs from GEE
    summarize(precipAtSinkFromHPA_m3 = sum(precipAtSinkFromSource_m3, na.rm = T),
              precipAtSinkFromHPA_m3_noIrr = sum(precipAtSinkFromSource_m3_noIrr, na.rm = T)) |> 
    mutate(month = m) |> 
    left_join(dplyr::select(df_ERA_mo, cellid, precipitation_m3), by = c("cellid_sink" = "cellid")) |> 
    rename(precipitation_m3_ERA = precipitation_m3)
  
  if (m == 1){
    df_sink <- df_mo_sink
  } else {
    df_sink <- bind_rows(df_sink, df_mo_sink)
  }
  
}

# add spatial information
sf_sink <- 
  right_join(df_sink, grid_conus, by = c("cellid_sink" = "cellid")) |> 
  st_as_sf() |> 
  # FOR TESTING
  subset(month %in% seq(5, 10))

# plot - percent of precipitation at cell from HPA region
ggplot() +
  geom_sf(data = sf_sink, aes(fill =  precipAtSinkFromHPA_m3 / precipitation_m3_ERA), color = NA) +
  geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_b(name = "Precip from HPA cells [%]", 
                       labels = scales::percent) +
  facet_wrap(~month, nrow = 2) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())


# plot - precipitation loss as percentage if no irrigation
ggplot() +
  geom_sf(data = sf_sink, aes(fill = (precipAtSinkFromHPA_m3 - precipAtSinkFromHPA_m3_noIrr) / precipitation_m3_ERA), 
          color = NA) +
  geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_c(name = "Precip Loss if No HPA Irrigation [%]", 
                       labels = scales::percent) +
  facet_wrap(~month, nrow = 2) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

# plot - precipitation loss as mm if no irrigation
ggplot() +
  geom_sf(data = sf_sink, aes(fill = 1000*(precipAtSinkFromHPA_m3 - precipAtSinkFromHPA_m3_noIrr) / area_m2), 
          color = NA) +
  geom_sf(data = state_boundaries, fill = NA, color = "white") +
  scale_fill_viridis_c(name = "Precip Loss if No HPA Irrigation [mm]") +
  facet_wrap(~month, nrow = 2) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
