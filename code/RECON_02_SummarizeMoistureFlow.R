## RECON_02_SummarizeMoistureFlow.R

source(file.path("code", "paths+packages.R"))

# load files
sf_hpa <- st_read(file.path("data", "ERA5-HPA_LANID_AverageIrrigation.shp"))
sf_conus <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp"))
sf_states <- st_read(file.path("data", "boundaries", "CONUS-States_TIGRIS.shp"))
df_flows <- read_csv(file.path("data", "RECON_01_HPAirrigationMoistureFlow.csv"))

# # inspect output
# ggplot(sf_hpa, aes(fill = irr_avgPrc > 0.1, color = cellid %in% df_flows$cellid_source)) +
#   geom_sf() +
#   scale_fill_manual(values = c("FALSE" = col.cat.org, "TRUE" = col.cat.grn)) +
#   scale_color_manual(values = c("TRUE" = col.cat.blu, "FALSE" = col.cat.yel))
# 
# length(unique(df_flows$cellid_source))
# length(unique(df_flows$cellid_sink))
# sum(sf_hpa$irr_avgPrc > 0.1)

# join state to output
df_flows_states <-
  df_flows |> 
  left_join(dplyr::select(st_drop_geometry(sf_conus), cellid, state), by = c("cellid_source" = "cellid")) |> 
  rename(state_source = state) |> 
  left_join(dplyr::select(st_drop_geometry(sf_conus), cellid, state), by = c("cellid_sink" = "cellid")) |> 
  rename(state_sink = state)

# inspect output
table(df_flows_states$state_source)
table(df_flows_states$state_sink)
length(unique(df_flows_states$state_sink))

# moisture flows: summarize state to state sums (by state_source and state_sink)
df_flows_s2s <-
  df_flows_states |> 
  group_by(state_source, state_sink) |> 
  summarize(moistureflow_s2s = sum(moistureflow_volume))

# plot for each state
states_hpa <- unique(df_flows_s2s$state_source)
for (state_plot in states_hpa){
  
  sf_flows_s2s_1state <- 
    sf_states |> 
    left_join(subset(df_flows_s2s, state_source == state_plot), by = c("STUSPS" = "state_sink")) |> 
    mutate(moistureflow_s2s = replace_na(moistureflow_s2s, 0))
  
  p_map <-
    ggplot(sf_flows_s2s_1state, aes(fill = moistureflow_s2s)) +
    geom_sf() +
    scale_fill_viridis_c(name = "Moisture Flow [m^3]", direction = -1, guide = "none") +
    labs(title = paste0("Fate of Evaporation from Irrigated HPA Pixels in ", state_plot))
  
  p_bar <- 
    ggplot(sf_flows_s2s_1state, aes(x = reorder(STUSPS, moistureflow_s2s), y = moistureflow_s2s, fill = moistureflow_s2s)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(name = "Moisture Flow [m^3]", expand = expansion(mult = c(0, 0.02))) +
    scale_x_discrete(name = "Receiving State") +
    scale_fill_viridis_c(name = "Moisture Flow [m^3]", direction = -1) +
    labs(title = paste0("Fate of Evaporation from Irrigated HPA Pixels in ", state_plot)) +
    theme(legend.position = c(0.95, 0.05),
          legend.justification = c(1, 0))
  
  ((p_map + p_bar) +
      plot_layout(ncol = 1, 
                  heights = c(1, 3))) |> 
    ggsave(file.path("figures+tables", paste0("RECON_02_SummarizeMoistureFlow-", state_plot, ".png")),
           plot = _, width = 150, height = 230, units = "mm")
}

## chord diagram
# for simplicity: get rid of states that are very small sinks
df_flows_sinkTotals <- 
  df_flows_s2s |> 
  group_by(state_sink) |> 
  summarize(moistureflow_sinkTotal = sum(moistureflow_s2s))
sink_max <- max(df_flows_sinkTotals$moistureflow_sinkTotal)
states_keep <- df_flows_sinkTotals$state_sink[df_flows_sinkTotals$moistureflow_sinkTotal > sink_max*0.1]

df_flows_s2s |> 
  subset(state_sink %in% states_keep) |> 
  rename(from = state_source, to = state_sink, value = moistureflow_s2s) |> 
  chordDiagram()
