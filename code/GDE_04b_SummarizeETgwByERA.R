## GDE_04b_SummarizeETgwByERA.R
# Run after GDE_04a - take individual tile outputs and combine to one entry per ERA cellid.
#  (this is necessary due to some ERA cells being split between multiple tiles)
#
# Author: Sam Zipper
# Date: 11/16/2026

# set up workspace
source(file.path("code", "paths+packages.R"))

# source for ETgw data
path_ETgwByERA <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_04_ETgwByERA")

# get list of tiles
tiles_all <- list.files(path_ETgwByERA, pattern = "GDE_04_ETgwByERA_")
tiles <- str_sub(sapply(str_split(tiles_all, pattern = "_"), `[`, 4), end = -5)

# load by tile and combine
for (tile in tiles){
  #tile <- "n35w110"
  df_t <- read_csv(file.path(path_ETgwByERA, paste0("GDE_04_ETgwByERA_", tile, ".csv")))
  
  # combine
  if (tile == tiles[1]){
    df_all <- df_t
  } else {
    df_all <- bind_rows(df_all, df_t)
  }
}

# summarize by cellid
df_cells <-
  df_all |> 
  group_by(cellid) |> 
  summarize(across(-tile, \(x) sum(x, na.rm = FALSE)))

# save
write_csv(df_cells, file.path(path_ETgwByERA, "GDE_04_ETgwByERA-AllCells.csv"))

# pivot longer
df_cells_long <-
  df_cells |> 
  pivot_longer(
    cols = -c(cellid, area_m2),
    names_to = c(".value", "month"),
    names_pattern = "(ETgw_m3|ET_m3)_(\\d{2})"
  ) %>%
  arrange(cellid, month) |> 
  mutate(ETgw_prc = ETgw_m3 / ET_m3)

# take a look
sf_era_native <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_Global-Land-NoPoles.shp"))

sf_era_tile_withETgw <- 
  right_join(df_cells_long, sf_era_native, by = "cellid") |> 
  st_as_sf() |> 
  subset(is.finite(ETgw_prc))

ggplot(sf_era_tile_withETgw, aes(fill = ETgw_prc)) +
  geom_sf(data = sf_era_native, color = NA, fill = "grey85") +
  geom_sf(color = NA) +
  scale_fill_viridis_c(name = "ETgw [%]", direction = -1) +
  facet_wrap(~month) +
  labs(title = "Dryland GDE ETgw as % of monthly ET (grey = not analyzed)") +
  theme(axis.text = element_blank(),
        legend.position = "bottom")
ggsave(file.path("figures+tables", "GDE_04b_MapETgwByPixel.png"),
       width = 260, height = 120, units = "mm")

