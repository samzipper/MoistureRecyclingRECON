## GDE_00_CreatePolygonsERA5.R
# This script is intended to create a polygon shapefile from the RECON ERA5 NetCDFs.

source(file.path("code", "paths+packages.R"))

## load data
# Open one of the RECON or ERA5 NetCDF files
r <- rast(file.path(path_data, "RECON", "RECON_ERA5_avgYear_0.5_volumes.nc"), subds = 1)
r_rotate <- rotate(r)

# load ERA5 land-sea mask (downloaded from copernicus climate data store)
era_lsm <- rast(file.path("data", "ERA5_grid", "ERA5-LSM.grib"))

# trim polar regions (> 60 lat)
era_lsm_landSea <-
  era_lsm |> 
  crop(ext(xmin(era_lsm), xmax(era_lsm), -60, 60))

plot(era_lsm_landSea)

# Convert to polygons
grid_lsm <- as.polygons(
  era_lsm_landSea,
  aggregate = FALSE,
  values = TRUE
)

grid_lsm$landPrc <- round(grid_lsm$`LSM_0-SFC`, 3)
grid_lsm$`LSM_0-SFC` <- NULL

grid_lsm$cellid <- 1:nrow(grid_lsm)

# subset to land only
grid_lsm_land <- grid_lsm[grid_lsm$landPrc > 0.5]

# Write output
writeVector(grid_lsm, 
            file.path("data", "ERA5_grid", "ERA5_grid_Global-NoPoles.shp"),
            overwrite = TRUE)

writeVector(grid_lsm_land, 
            file.path("data", "ERA5_grid", "ERA5_grid_Global-Land-NoPoles.shp"),
            overwrite = TRUE)
