## GDE_04_SummarizeETgwByERA.R
# This script is intended to summarize the ETgw estimated by GDE_03_CalculateETgw.R
# for each ERA grid cell to align with UTrack.
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

ggplot(sf_era, aes(fill = landPrc)) +
  geom_sf(color = "transparent")

# load a random tile to reproject sf_era
r_tile <- rast(file.path(path_ETgw, "GDE_03_ETgw_prc_n35w110.tif"))
v_era <- project(v_era_native, crs(r_tile))

# load tile
tile <- "n35w110"

# load ET and ETgw
r_ET_mm <- rast(file.path(path_ET.P, paste0("GDE_ET_monthly_climatology_2008_2017_", tile, ".tif")))
r_ETgw_prc <- rast(file.path(path_ETgw, paste0("GDE_03_ETgw_prc_", tile, ".tif")))
r_ETgw_mm <- r_ET_mm*r_ETgw_prc


# for each ERA grid cell calculate:
#  total ET
#  ET from GW
v_era_tile <- crop(v_era, r_ETgw)

# Sum raster values within each polygon
cell_ET_mm <- terra::extract(
  r_ET_mm,
  v_era_tile,
  fun = sum,
  na.rm = TRUE
)

cell_ETgw_mm <- terra::extract(
  r_ETgw_mm,
  v_era_tile,
  fun = sum,
  na.rm = TRUE
)

# get ETgw as % of total ET
cell_ETgw_prc <- cell_ETgw_mm / cell_ET_mm
cell_ETgw_prc$cellid <- v_era_tile$cellid
cell_ETgw_prc$ID <- NULL
names(cell_ETgw_prc) <- str_replace(names(cell_ETgw_prc), "ET_", "ETgw_prc_")

## NEED TO FIGURE OUT HOW TO DEAL WITH EDGES OF TILES

## visualize
v_era_tile_withETgw <- merge(v_era_tile, cell_ETgw_prc, by.x = "cellid", by.y = "cellid")

sf_era_tile_withETgw <- st_as_sf(v_era_tile_withETgw)
ggplot(sf_era_tile_withETgw, aes(fill = ETgw_prc_07)) +
  geom_sf()
