## GDE_04a_SummarizeETgwByERAbyTile.R
# This script is intended to summarize the ETgw estimated by GDE_03_CalculateETgw.R
# for each ERA grid cell to align with UTrack, with a separate output for each tile.
#
# Author: Sam Zipper
# Date: 11/10/2026

# set up workspace
source(file.path("code", "paths+packages.R"))

# source for ETgw data
path_ETgw <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_03_CalculateETgw")
path_ET.P <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GEEexport_ET+Pdata")
path_ETgwByERA <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_04_ETgwByERA")

# load ERA polygons
v_era_native <- terra::vect(file.path("data", "ERA5_grid", "ERA5_grid_Global-Land-NoPoles.shp"))

# load a random tile to reproject sf_era
r_tile <- rast(file.path(path_ETgw, "GDE_03_ETgw_prc_n35w110.tif"))
v_era <- project(v_era_native, crs(r_tile))

# get list of tiles
tiles_all <- list.files(path_ETgw, pattern = "GDE_03_ETgw_prc_*")
tiles <- str_sub(sapply(str_split(tiles_all, pattern = "_"), `[`, 5), end = -5)

for (tile in tiles){
  #tile <- "n35w110"
  
  # load ET and ETgw
  r_ET_mm <- rast(file.path(path_ET.P, paste0("GDE_ET_monthly_climatology_2008_2017_", tile, ".tif")))
  r_ETgw_prc <- rast(file.path(path_ETgw, paste0("GDE_03_ETgw_prc_", tile, ".tif")))
  
  # get area of each cell and convert to m^3
  r_area <- cellSize(r_ET_mm)
  r_ET_m3 <- (r_ET_mm/1000)*r_area
  r_ETgw_m3 <- r_ETgw_prc*r_ET_m3
  
  # for each ERA grid cell calculate:
  #  area
  #  total ET
  #  ET from GW
  v_era_tile <- crop(v_era, r_ETgw_prc)
  
  # get area of each ERA tile
  cell_area <- expanse(v_era_tile)
  
  # Sum raster values within each polygon
  cell_ET_m3 <- terra::extract(
    r_ET_m3,
    v_era_tile,
    fun = sum,
    na.rm = TRUE
  )
  names(cell_ET_m3) <- str_replace(names(cell_ET_m3), "ET_", "ET_m3_")
  
  cell_ETgw_m3 <- terra::extract(
    r_ETgw_m3,
    v_era_tile,
    fun = sum,
    na.rm = TRUE
  )
  
  names(cell_ETgw_m3)[2:13] <- paste0("ETgw_m3_", names(cell_ETgw_m3)[2:13])
  
  # compile output for each cell
  cell_out <- tibble(
    tile = tile,
    cellid = v_era_tile$cellid,
    area_m2 = round(cell_area, 1)
  ) |> 
    bind_cols(dplyr::select(cell_ET_m3, -ID),
              dplyr::select(cell_ETgw_m3, -ID))
  
  write_csv(cell_out, file.path(path_ETgwByERA, paste0("GDE_04_ETgwByERA_", tile, ".csv")))
  
  print(paste0(tile, " complete, ", which(tiles==tile)))
}
     

# # get ETgw as % of total ET
# cell_ETgw_prc <- cell_ETgw_mm / cell_ET_mm
# cell_ETgw_prc$cellid <- v_era_tile$cellid
# cell_ETgw_prc$ID <- NULL
# names(cell_ETgw_prc) <- str_replace(names(cell_ETgw_prc), "ET_", "ETgw_prc_")
# 
# ## visualize
# v_era_tile_withETgw <- merge(v_era_tile, cell_ETgw_prc, by.x = "cellid", by.y = "cellid")
# 
# sf_era_tile_withETgw <- st_as_sf(v_era_tile_withETgw)
# ggplot(sf_era_tile_withETgw, aes(fill = ETgw_prc_07)) +
#   geom_sf()
