## UTrack source-to-sink fractions for polygon-grid source cells

source(file.path("code", "paths+packages.R"))

## inputs and outputs
# Source IDs come from the HPA grid; sink IDs use the CONUS grid or -9999 outside CONUS.
polygon_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_HPA.shp")
conus_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_CONUS.shp")
utrack_directory <- file.path(path_data, "UTrack")
output_directory <- file.path(path_data, "Processed", "UTrack-Polygon_01_MonthlySourceToSink")
flow_encoding <- "encoded-log" # use "fraction" if moisture_flow is already a fraction
source_lat_chunk_size <- 1

## helpers
# UTrack uses 0-360 longitudes, while the shapefiles use -180 to 180.
normalize_longitude <- function(longitude) {
  (longitude + 180) %% 360 - 180
}

# Format coordinates consistently so polygon cells can be matched to UTrack cells.
coordinate_key <- function(latitude, longitude) {
  paste(
    formatC(latitude, format = "f", digits = 6),
    formatC(normalize_longitude(longitude), format = "f", digits = 6),
    sep = "|"
  )
}

read_grid_centres <- function(grid_file) {
  if (!file.exists(grid_file)) {
    stop("Grid shapefile not found: ", grid_file)
  }

  # Use WGS84 coordinates and the bounding-box midpoint for each regular grid cell.
  grid <- st_read(grid_file, quiet = TRUE) |> st_transform(4326)
  if (!("cellid" %in% names(grid))) {
    stop(grid_file, " must contain a cellid field.")
  }

  cellid <- grid$cellid
  if (anyNA(cellid) || any(cellid != as.integer(cellid)) || anyDuplicated(cellid)) {
    stop(grid_file, " must contain unique, nonmissing integer cellid values.")
  }

  # set centre based on bounding box of cell
  centres <- t(vapply(st_geometry(grid), function(cell) {
    bbox <- st_bbox(cell)
    c(
      latitude = (bbox[["ymin"]] + bbox[["ymax"]]) / 2,
      longitude = (bbox[["xmin"]] + bbox[["xmax"]]) / 2
    )
  }, numeric(2)))
  
  
  tibble(
    cellid = as.integer(cellid),
    latitude = centres[, "latitude"],
    longitude = centres[, "longitude"],
    coordinate_key = coordinate_key(
      centres[, "latitude"],
      centres[, "longitude"])
    )
}

decode_flow <- function(moisture_flow, encoding) {
  # Encoded UTrack values are -10 times the log of the flow fraction; zero means no flow.
  if (encoding == "fraction") {
    return(ifelse(moisture_flow > 0, moisture_flow, 0))
  }
  if (encoding == "encoded-log") {
    return(ifelse(moisture_flow == 0, 0, exp(-0.1 * moisture_flow)))
  }
  stop("flow_encoding must be 'encoded-log' or 'fraction'.")
}

validate_grid_keys <- function(grid, available_keys, grid_name) {
  # Ensure every CONUS sink center occurs in the target coordinates for this month.
  missing <- !grid$coordinate_key %in% available_keys
  if (any(missing)) {
    stop(
      sum(missing), " ", grid_name,
      " grid cell centres do not match the UTrack grid (within 1e-6 degrees)."
    )
  }
}

run_month <- function(month, source_cells, conus_cells, min_frac) {
  utrack_file <- file.path(
    utrack_directory,
    sprintf("utrack_climatology_0.5_%02d.nc", month)
  )
  if (!file.exists(utrack_file)) {
    stop("UTrack input file not found: ", utrack_file)
  }
  if (!requireNamespace("tidync", quietly = TRUE)) {
    stop("Package 'tidync' is required. Install it with install.packages('tidync').")
  }

  # Read only the selected source rows; each row retains its global target grid.
  ut <- tidync::tidync(utrack_file)
  source_lats <- sort(unique(source_cells$latitude))
  latitude_chunks <- split(
    source_lats,
    ceiling(seq_along(source_lats) / source_lat_chunk_size)
  )
  source_longitudes_utrack <- source_cells$longitude %% 360
  source_lon_min <- min(source_longitudes_utrack)
  source_lon_max <- max(source_longitudes_utrack)
  source_ids_found <- integer()
  target_keys <- character()
  
  # create empty output file that will be added to
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_directory,
    sprintf("UTrack-Polygon_01_MonthlySourceToSink_%02d.csv", month)
  )
  write_csv(
    tibble(
      cellid_source = integer(),
      cellid_sink = integer(),
      fraction_sourceEvapToSink = double()
    ),
    output_file
  )

  # Process source latitude bands in chunks to limit memory use.
  for (chunk_i in seq_along(latitude_chunks)) {
    chunk_lats <- latitude_chunks[[chunk_i]]
    chunk_source_cells <- source_cells |>
      filter(latitude %in% chunk_lats)
    flows <- ut |>
      tidync::hyper_filter(
        sourcelat = sourcelat >= min(chunk_lats) - 0.26 &
          sourcelat <= max(chunk_lats) + 0.26,
        sourcelon = sourcelon >= source_lon_min - 0.26 &
          sourcelon <= source_lon_max + 0.26
      ) |>
      tidync::hyper_tibble() |>
      mutate(
        source_key = coordinate_key(sourcelat, sourcelon),
        sink_key = coordinate_key(targetlat, targetlon)
      ) |>
      # Keep only source cells represented by the polygon grid.
      filter(source_key %in% chunk_source_cells$coordinate_key)

    if (nrow(flows) == 0) {
      next
    }

    source_ids_found <- union(
      source_ids_found,
      chunk_source_cells$cellid[match(unique(flows$source_key), chunk_source_cells$coordinate_key)]
    )
    target_keys <- union(target_keys, unique(flows$sink_key))

    chunk_result <- flows |>
      mutate(
        cellid_source = chunk_source_cells$cellid[
          match(source_key, chunk_source_cells$coordinate_key)
        ],
        cellid_sink = conus_cells$cellid[
          match(sink_key, conus_cells$coordinate_key)
        ],
        # Any sink not found in the CONUS grid is grouped under the outside-CONUS ID.
        cellid_sink = replace_na(cellid_sink, -9999L),
        fraction_sourceEvapToSink = decode_flow(
          moisture_flow, flow_encoding
        )
      ) |>
      filter(
        !is.na(fraction_sourceEvapToSink),
        fraction_sourceEvapToSink > 0
      ) |>
      group_by(cellid_source, cellid_sink) |>
      summarize(
        fraction_sourceEvapToSink = sum(
          fraction_sourceEvapToSink
        ),
        .groups = "drop"
      ) |>
      # Normalize after combining sinks so each source's fractions sum to one.
      group_by(cellid_source) |>
      mutate(
        source_total = sum(fraction_sourceEvapToSink),
        fraction_sourceEvapToSink =
          fraction_sourceEvapToSink / source_total
      ) |>
      ungroup() |>
      select(-source_total)
    if (nrow(chunk_result) > 0) {
      # Omit negligible links and append each chunk to the monthly CSV.
      chunk_result |> 
        subset(fraction_sourceEvapToSink > min_frac) |> 
        write_csv(output_file, append = TRUE, col_names = FALSE)
    }
  }

  if (!setequal(source_ids_found, source_cells$cellid)) {
    missing_ids <- setdiff(source_cells$cellid, source_ids_found)
    stop(
      "Some polygon source cells were not found on the UTrack source grid: ",
      paste(head(missing_ids, 10), collapse = ", ")
    )
  }
  validate_grid_keys(conus_cells, target_keys, "CONUS")

  message("Month ", sprintf("%02d", month), ": wrote ", output_file)
  invisible(output_file)
}

if (sys.nframe() == 0) {
  if (!flow_encoding %in% c("encoded-log", "fraction")) {
    stop("flow_encoding must be 'encoded-log' or 'fraction'.")
  }
  if (source_lat_chunk_size < 1) {
    stop("source_lat_chunk_size must be positive.")
  }

  # Match polygon cell IDs to UTrack source coordinates and CONUS sink IDs.
  source_cells <- read_grid_centres(polygon_grid_file) |>
    mutate(coordinate_key = coordinate_key(latitude, longitude))
  conus_cells <- read_grid_centres(conus_grid_file) |>
    select(cellid, coordinate_key)
  
  # set up run
  months_run <- seq(1, 12)
  min_frac <- 1e-6 # minimum fraction to be included in final table
  
  # # Test run: process five source cells for July; remove these lines for all cells/months.
  # source_cells <- source_cells[1:5,]
  # months_run <- 7

  if (anyDuplicated(source_cells$coordinate_key)) {
    stop("Distinct polygon cellids map to the same source-cell centre.")
  }
  if (anyDuplicated(conus_cells$coordinate_key)) {
    stop("CONUS shapefile contains duplicate grid-cell centres.")
  }

  # Run monthly code.
  for (month in months_run) {
    run_month(month, source_cells, conus_cells, min_frac)
  }
}
