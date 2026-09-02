## DataPrep_00_CreatePolygonsERA5.R
# This script is intended to create a polygon shapefile from the RECON ERA5 NetCDFs
# which can then be used in GEE to extract other datasets.

source(file.path("code", "paths+packages.R"))

## load data
# Open one of the RECON or ERA5 NetCDF files
r <- rast(file.path(path_data_recon, "RECON_ERA5_avgYear_0.5_volumes.nc"), subds = 1)
r_rotate <- rotate(r)

# CONUS states
conus <- states(cb = TRUE, year = 2024)
conus <- conus[!conus$STUSPS %in%
                 c("AK","HI","PR","GU","VI","MP","AS"), ]
conus_transform <- 
  st_transform(conus, crs = st_crs(r_rotate)) |> 
  st_make_valid()

conus_union <- 
  conus_transform |> 
  st_union() |> 
  st_make_valid()

# write CONUS boundary
st_write(conus_transform, file.path("data", "boundaries", "CONUS-States_TIGRIS.shp"),
         delete_dsn = T)
st_write(conus_union, file.path("data", "boundaries", "CONUS-Outline_TIGRIS.shp"),
         delete_dsn = T)

## trim raster to north america and to CONUS
# North America bounding box
na_ext <- ext(
  -171.7911,  # xmin
  -66.9647,   # xmax
  18.9162,   # ymin
  71.3578    # ymax
)

# Crop raster
r_na <- crop(r_rotate, na_ext)

# Convert to polygons
grid_na <- as.polygons(
    r_na,
    aggregate = FALSE,
    values = FALSE
  )

grid_na$cellid <- 1:nrow(grid_na)

# Write output
writeVector(
  grid_na,
  file.path("data", "ERA5_grid", "ERA5_grid_NorthAmerica.shp"),
  overwrite = TRUE
)

## crop to just cells wholly contained within CONUS
# Returns a logical vector, one value per grid cell
conus_terra <- vect(conus_union)
inside <- relate(grid_na, conus_terra, relation = "within")

# Subset to cells wholly contained within CONUS
grid_conus <- grid_na[inside, ]

# check
plot(conus_terra)
plot(grid_conus, add = TRUE)

writeVector(
  grid_conus,
  file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp"),
  overwrite = TRUE
)
