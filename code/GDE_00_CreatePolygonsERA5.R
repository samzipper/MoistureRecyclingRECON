## GDE_00_CreatePolygonsERA5.R
# This script is intended to create a polygon shapefile from the RECON ERA5 NetCDFs.

source(file.path("code", "paths+packages.R"))

## load data
# Open one of the RECON or ERA5 NetCDF files
r <- rast(file.path(path_data, "RECON", "RECON_ERA5_avgYear_0.5_volumes.nc"), subds = 1)
r_rotate <- rotate(r)

# Set explicit CRS to WGS84 (EPSG:4326) - ERA5 data uses geographic coordinates
crs(r_rotate) <- "EPSG:4326"

# load ERA5 land-sea mask (downloaded from copernicus climate data store)
era_lsm_in <- rast(file.path("data", "ERA5_grid", "ERA5-LSM.grib"))

# Shift longitude coordinates by 180 degrees to correct hemisphere alignment
# (converts -180 to 0 to 0 to 180 and 0 to 180 to -180 to 0)
era_lsm_shifted <- 
  era_lsm_in |> 
  rotate() |> 
  shift(dx = -180)

# Set explicit CRS to WGS84 for the shifted LSM
crs(era_lsm_shifted) <- "EPSG:4326"

# Project shifted LSM to match the rotated RECON raster
era_lsm <- project(era_lsm_shifted, r_rotate, method = "average")

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

# Ensure explicit CRS on output shapefiles
crs(grid_lsm) <- "EPSG:4326"
crs(grid_lsm_land) <- "EPSG:4326"

# Write output with explicit CRS
writeVector(grid_lsm, 
            file.path("data", "ERA5_grid", "ERA5_grid_Global-NoPoles.shp"),
            overwrite = TRUE)

writeVector(grid_lsm_land, 
            file.path("data", "ERA5_grid", "ERA5_grid_Global-Land-NoPoles.shp"),
            overwrite = TRUE)
