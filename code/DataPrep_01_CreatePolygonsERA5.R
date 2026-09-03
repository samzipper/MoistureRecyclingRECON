## DataPrep_01_CreatePolygonsERA5.R
# This script is intended to create a polygon shapefile from the RECON ERA5 NetCDFs
# which can then be used in GEE to extract other datasets.

source(file.path("code", "paths+packages.R"))

## load data
# Open one of the RECON or ERA5 NetCDF files
r <- rast(file.path(path_data, "RECON", "RECON_ERA5_avgYear_0.5_volumes.nc"), subds = 1)
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

# inspect
plot(r_rotate)
plot(grid_na, add = TRUE)

# Write output
writeVector(
  grid_na,
  file.path("data", "ERA5_grid", "ERA5_grid_NorthAmerica.shp"),
  overwrite = TRUE
)

## subset to CONUS grid cells
conus_terra <- vect(conus_union)

# create a raster with same geometry as ERA5 grid
r_template <- r_na[[1]]

# fraction of each cell covered by CONUS
frac <- rasterize(
  conus_terra,
  r_template,
  field = 1,
  cover = TRUE
)

grid_na$pct_conus <- values(frac)

grid_conus <- grid_na[
  grid_na$pct_conus >= 0.75,
]

# check
plot(conus_terra)
plot(grid_conus, add = TRUE)

# determine state
grid_sf <- st_as_sf(grid_conus)

# Make sure both layers have valid geometries
grid_sf <- st_make_valid(grid_sf)
conus_transform <- st_make_valid(conus_transform)

# Intersect cells and states
cell_state_int <- st_intersection(
  grid_sf |> select(cellid),
  conus_transform |> select(STUSPS)
)

# Calculate overlap area
cell_state_int$area <- st_area(cell_state_int)

# For each cell, keep the state with the largest overlap
cell_state_max <- cell_state_int |>
  st_drop_geometry() |>
  group_by(cellid) |>
  slice_max(area, n = 1, with_ties = FALSE) |>
  ungroup()

# Join back to grid
grid_conus_state <- grid_sf |>
  left_join(cell_state_max, by = "cellid") |> 
  dplyr::select(cellid, pct_conus, state = STUSPS)

# check
ggplot(grid_conus_state, aes(fill = state)) + geom_sf()

st_write(
  grid_conus_state,
  file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp"),
  delete_dsn = TRUE
)

## subset to HPA grid cells
hpa_vect <- vect(file.path("data", "boundaries", "HPA_extent.shp")) |> 
  project(crs(r_template))

# fraction of each cell covered by HPA
hpa_frac <- rasterize(
  hpa_vect,
  r_template,
  field = 1,
  cover = TRUE
)

grid_na$pct_hpa <- values(hpa_frac)

grid_hpa <- grid_na[
  grid_na$pct_hpa >= 0.25,
]
# check
#plot(conus_terra)
plot(hpa_vect)
plot(grid_hpa, add = TRUE)

# determine state
grid_hpa_sf <- st_as_sf(grid_hpa)

# Make sure both layers have valid geometries
grid_hpa_sf <- st_make_valid(grid_hpa_sf)

# Intersect cells and states
cell_hpa_state_int <- st_intersection(
  grid_hpa_sf |> select(cellid),
  conus_transform |> select(STUSPS)
)

# Calculate overlap area
cell_hpa_state_int$area <- st_area(cell_hpa_state_int)

# For each cell, keep the state with the largest overlap
cell_hpa_state_max <- cell_hpa_state_int |>
  st_drop_geometry() |>
  group_by(cellid) |>
  slice_max(area, n = 1, with_ties = FALSE) |>
  ungroup()

# Join back to grid
grid_hpa_state <- grid_hpa_sf |>
  left_join(cell_hpa_state_max, by = "cellid") |> 
  dplyr::select(cellid, pct_conus, state = STUSPS)

# check
ggplot(grid_hpa_state, aes(fill = state)) + geom_sf()

st_write(
  grid_hpa_state,
  file.path("data", "ERA5_grid", "ERA5_grid_HPA.shp"),
  delete_dsn = TRUE
)
