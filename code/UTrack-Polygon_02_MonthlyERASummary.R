## Monthly ERA5 evaporation and precipitation volumes for CONUS polygon-grid cells

source(file.path("code", "paths+packages.R"))

## inputs and outputs
conus_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp")
era5_directory <- file.path(path_data, "RECON")
output_file <- file.path("data", "UTrack-Polygon_02_MonthlyERASummary.csv")

## helpers
# ERA5 uses 0-360 longitudes, while the shapefile uses -180 to 180.
coordinate_key <- function(latitude, longitude) {
  paste(
    formatC(latitude, format = "f", digits = 6),
    formatC((longitude + 180) %% 360 - 180, format = "f", digits = 6),
    sep = "|"
  )
}

read_conus_cells <- function(grid_file) {
  if (!file.exists(grid_file)) {
    stop("Grid shapefile not found: ", grid_file)
  }

  grid <- st_read(grid_file, quiet = TRUE) |> st_transform(4326)
  if (!("cellid" %in% names(grid))) {
    stop(grid_file, " must contain a cellid field.")
  }
  cellid <- grid$cellid
  if (nrow(grid) == 0 || anyNA(cellid) ||
      any(cellid != as.integer(cellid)) || anyDuplicated(cellid)) {
    stop(grid_file, " must contain unique, nonmissing integer cellid values.")
  }

  # Use the bounding-box midpoint for each regular grid cell, as in script 01.
  centres <- t(vapply(st_geometry(grid), function(cell) {
    bbox <- st_bbox(cell)
    c(
      latitude = (bbox[["ymin"]] + bbox[["ymax"]]) / 2,
      longitude = (bbox[["xmin"]] + bbox[["xmax"]]) / 2
    )
  }, numeric(2)))
  if (any(!is.finite(centres))) {
    stop("CONUS shapefile contains nonfinite grid-cell centres.")
  }

  cells <- tibble(
    cellid = as.integer(cellid),
    coordinate_key = coordinate_key(
      centres[, "latitude"], centres[, "longitude"]
    )
  )
  if (anyDuplicated(cells$coordinate_key)) {
    stop("CONUS shapefile contains duplicate grid-cell centres.")
  }
  cells
}

read_month <- function(month, conus_cells) {
  era5_file <- file.path(
    era5_directory,
    sprintf("ERA5_%02d_0.5_volumes_corrected.nc", month)
  )
  if (!file.exists(era5_file)) {
    stop("ERA5 input file not found: ", era5_file)
  }

  era5 <- tidync::tidync(era5_file)
  required_variables <- c("ERA5_ET", "ERA5_TP")
  if (!all(required_variables %in% era5$variable$name)) {
    stop(era5_file, " must contain ERA5_ET and ERA5_TP.")
  }
  for (variable in required_variables) {
    units <- era5$attribute |>
      filter(.data$variable == .env$variable, name == "units")
    if (nrow(units) != 1 || !identical(as.character(units$value[[1]]), "m^3")) {
      stop(era5_file, ": ", variable, " must have units m^3.")
    }
  }

  volumes <- era5 |> tidync::hyper_tibble()
  if (!all(c("lat", "lon", required_variables) %in% names(volumes))) {
    stop(era5_file, " must provide ERA5_ET and ERA5_TP on a lat/lon grid.")
  }
  keys <- coordinate_key(volumes$lat, volumes$lon)
  if (anyDuplicated(keys)) {
    stop(era5_file, " contains duplicate grid-cell coordinates.")
  }
  indices <- match(conus_cells$coordinate_key, keys)
  if (anyNA(indices)) {
    stop(
      era5_file, ": ", sum(is.na(indices)),
      " CONUS grid cell centres do not match the ERA5 grid (within 1e-6 degrees)."
    )
  }

  # These are already monthly full-cell volumes in m^3, not depths or daily rates.
  # Retain the full volume; do not multiply by pct_conus or cell area.
  result <- tibble(
    cellid = conus_cells$cellid,
    month = as.integer(month),
    evaporation_m3 = volumes$ERA5_ET[indices],
    precipitation_m3 = volumes$ERA5_TP[indices]
  )
  if (any(!is.finite(result$evaporation_m3)) ||
      any(!is.finite(result$precipitation_m3))) {
    stop(era5_file, " contains missing or nonfinite volumes for CONUS cells.")
  }

  message("Month ", sprintf("%02d", month), ": indexed ", nrow(result), " CONUS cells.")
  result
}

if (sys.nframe() == 0) {
  conus_cells <- read_conus_cells(conus_grid_file)

  # Run all climatological months; each cellid has one row per month.
  months_run <- seq(1, 12)
  monthly_summary <- bind_rows(lapply(months_run, read_month, conus_cells = conus_cells)) |>
    arrange(cellid, month)

  # Write only after all requested months have been read and validated.
  write_csv(monthly_summary, output_file)
  message("Wrote ", nrow(monthly_summary), " cell-month records to ", output_file)
}

## plot and inspect
head(monthly_summary)
sf_grid <- st_read(conus_grid_file)

test_month <- 1
monthly_test <- monthly_summary |> filter(month == test_month)
sf_test <- left_join(sf_grid, monthly_test, by = "cellid")
ggplot(sf_test, aes(fill = evaporation_m3)) + 
  geom_sf() +
  scale_fill_viridis_c(trans = "log10", direction = -1)
ggplot(sf_test, aes(fill = precipitation_m3)) + geom_sf() +
  scale_fill_viridis_c(trans = "log10", direction = -1)

ggplot(sf_test, aes(fill = precipitation_m3/evaporation_m3)) + geom_sf() +
  scale_fill_viridis_c(direction = -1)
