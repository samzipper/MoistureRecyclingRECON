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
tiles_all <- list.files(path_ETgwByERA, pattern = "GDE_04_ETgwByERA_n")
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
write_csv(df_cells, file.path(path_ETgwByERA, "GDE_04_ETgwByERA_AllCells.csv"))
