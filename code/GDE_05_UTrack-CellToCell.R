## GDE_05_UTrack-CellToCell.R
# Monthly UTrack source-to-sink fractions for global GDE source cells

source(file.path("code", "paths+packages.R"))

polygon_grid_file <- file.path("data", "ERA5_grid", "ERA5_grid_Global-Land-NoPoles.shp")
etgw_file <- file.path(
  path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_04_ETgwByERA",
  "GDE_04_ETgwByERA-AllCells.csv"
)
utrack_directory <- file.path(path_data, "UTrack")
output_directory <- file.path(
  path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_05_UTrack-CellToCell"
)
flow_encoding <- "encoded-log"
source_lon_chunk_size <- 4L
months_run <- seq(1, 12)
min_frac <- 1e-4

normalize_longitude <- function(longitude) {
  (longitude + 180) %% 360 - 180
}

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
  grid <- st_read(grid_file, quiet = TRUE) |> st_transform(4326)
  if (!("cellid" %in% names(grid))) {
    stop(grid_file, " must contain a cellid field.")
  }
  cellid <- grid$cellid
  if (anyNA(cellid) || any(cellid != as.integer(cellid)) || anyDuplicated(cellid) ||
      any(cellid == -9999L)) {
    stop("Grid cellids must be unique, nonmissing integers other than -9999.")
  }
  centres <- t(vapply(st_geometry(grid), function(cell) {
    bbox <- st_bbox(cell)
    c(
      latitude = (bbox[["ymin"]] + bbox[["ymax"]]) / 2,
      longitude = (bbox[["xmin"]] + bbox[["xmax"]]) / 2
    )
  }, numeric(2)))
  cells <- tibble(
    cellid = as.integer(cellid),
    latitude = centres[, "latitude"],
    longitude = centres[, "longitude"],
    coordinate_key = coordinate_key(centres[, "latitude"], centres[, "longitude"])
  )
  if (anyDuplicated(cells$coordinate_key)) {
    stop("Distinct polygon cellids map to the same grid-cell centre.")
  }
  cells
}

read_source_cells <- function(etgw_file, grid_cells, month) {
  if (length(month) != 1L || !is.finite(month) || !month %in% seq_len(12)) {
    stop("month must be an integer from 1 through 12.")
  }
  if (!file.exists(etgw_file)) {
    stop("ETgw input file not found: ", etgw_file)
  }
  etgw <- read_csv(etgw_file, show_col_types = FALSE)
  monthly_column <- sprintf("ETgw_m3_%02d", month)
  if (!all(c("cellid", monthly_column) %in% names(etgw))) {
    stop("ETgw input must contain cellid and ", monthly_column, ".")
  }
  source_ids <- etgw |>
    filter(!is.na(.data[[monthly_column]]), .data[[monthly_column]] > 0) |>
    pull(cellid) |>
    unique()
  if (length(source_ids) == 0) {
    stop("No source cells have positive ETgw in month ", sprintf("%02d", month), ".")
  }
  if (anyNA(source_ids) || any(!source_ids %in% grid_cells$cellid)) {
    stop(
      "Some positive-ETgw source cellids are missing from the polygon grid. ",
      "Rebuild GDE_04_ETgwByERA_AllCells.csv from current tile outputs only; ",
      "exclude AllCells itself from the GDE_04b input file list."
    )
  }
  grid_cells |> filter(cellid %in% source_ids)
}

decode_flow <- function(moisture_flow, encoding) {
  if (encoding == "fraction") {
    return(ifelse(moisture_flow > 0, moisture_flow, 0))
  }
  if (encoding == "encoded-log") {
    return(ifelse(moisture_flow == 0, 0, exp(-0.1 * moisture_flow)))
  }
  stop("flow_encoding must be 'encoded-log' or 'fraction'.")
}

run_month <- function(month, source_cells, sink_cells, min_frac) {
  utrack_file <- file.path(
    utrack_directory, sprintf("utrack_climatology_0.5_%02d.nc", month)
  )
  if (!file.exists(utrack_file)) {
    stop("UTrack input file not found: ", utrack_file)
  }
  if (!flow_encoding %in% c("encoded-log", "fraction")) {
    stop("flow_encoding must be 'encoded-log' or 'fraction'.")
  }
  if (length(source_lon_chunk_size) != 1L || !is.finite(source_lon_chunk_size) ||
      source_lon_chunk_size < 1 || source_lon_chunk_size != as.integer(source_lon_chunk_size)) {
    stop("source_lon_chunk_size must be a positive integer.")
  }
  if (length(min_frac) != 1L || !is.finite(min_frac) || min_frac < 0 || min_frac >= 1) {
    stop("min_frac must be between zero (inclusive) and one (exclusive).")
  }
  if (nrow(source_cells) == 0) {
    stop("No source cells supplied.")
  }

  ut <- tidync::tidync(utrack_file)
  source_chunks <- source_cells |>
    mutate(
      longitude_utrack = longitude %% 360,
      longitude_band = floor(longitude_utrack / (0.5 * source_lon_chunk_size))
    ) |>
    group_by(latitude, longitude_band) |>
    group_split()
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    output_directory, sprintf("GDE_05_UTrack-CellToCell_%02d.csv", month)
  )
  write_csv(
    tibble(
      cellid_source = integer(), cellid_sink = integer(),
      fraction_sourceEvapToSink = double()
    ),
    output_file
  )
  source_ids_found <- integer()
  sinks_validated <- FALSE

  for (chunk_i in seq_along(source_chunks)) {
    chunk_sources <- source_chunks[[chunk_i]]
    source_latitude <- chunk_sources$latitude[[1]]
    source_lon_min <- min(chunk_sources$longitude_utrack)
    source_lon_max <- max(chunk_sources$longitude_utrack)
    flows <- ut |>
      tidync::hyper_filter(
        sourcelat = abs(sourcelat - source_latitude) < 1e-6,
        sourcelon = sourcelon >= source_lon_min - 1e-6 &
          sourcelon <= source_lon_max + 1e-6
      ) |>
      tidync::hyper_tibble(na.rm = FALSE) |>
      mutate(
        source_key = coordinate_key(sourcelat, sourcelon),
        sink_key = coordinate_key(targetlat, targetlon)
      ) |>
      filter(source_key %in% chunk_sources$coordinate_key)

    if (!sinks_validated) {
      missing_sinks <- !sink_cells$coordinate_key %in% unique(flows$sink_key)
      if (any(missing_sinks)) {
        stop(sum(missing_sinks), " polygon sink centres do not match the UTrack grid.")
      }
      sinks_validated <- TRUE
    }

    chunk_result <- flows |>
      transmute(
        cellid_source = chunk_sources$cellid[
          match(source_key, chunk_sources$coordinate_key)
        ],
        cellid_sink = replace_na(
          sink_cells$cellid[match(sink_key, sink_cells$coordinate_key)], -9999L
        ),
        fraction_sourceEvapToSink = decode_flow(moisture_flow, flow_encoding)
      ) |>
      filter(is.finite(fraction_sourceEvapToSink), fraction_sourceEvapToSink > 0) |>
      group_by(cellid_source, cellid_sink) |>
      summarize(
        fraction_sourceEvapToSink = sum(fraction_sourceEvapToSink),
        .groups = "drop"
      ) |>
      group_by(cellid_source) |>
      mutate(fraction_sourceEvapToSink =
               fraction_sourceEvapToSink / sum(fraction_sourceEvapToSink)) |>
      ungroup()

    missing_sources <- setdiff(chunk_sources$cellid, chunk_result$cellid_source)
    if (length(missing_sources) > 0) {
      stop("Source cells have no positive UTrack flow: ", paste(missing_sources, collapse = ", "))
    }
    source_ids_found <- union(source_ids_found, chunk_result$cellid_source)
    chunk_result |>
      filter(fraction_sourceEvapToSink > min_frac) |>
      write_csv(output_file, append = TRUE, col_names = FALSE)
    message(
      "Month ", sprintf("%02d", month), ": chunk ", chunk_i, "/", length(source_chunks)
    )
  }
  if (!setequal(source_ids_found, source_cells$cellid)) {
    stop("Some polygon source cells were not processed.")
  }
  message("Month ", sprintf("%02d", month), ": wrote ", output_file)
  invisible(output_file)
}

if (sys.nframe() == 0) {
  grid_cells <- read_grid_centres(polygon_grid_file)
  
  sink_cells <- grid_cells ## TESTING: RUN WESTERN HEMISPHERE ONLY |> filter(normalize_longitude(longitude) < -60 & normalize_longitude(longitude) > -130)
  for (month in months_run) {
    source_cells <- read_source_cells(etgw_file, grid_cells, month)
    message(
      "Month ", sprintf("%02d", month), ": selected ", nrow(source_cells),
      " sources and ", nrow(sink_cells), " polygon sinks."
    )
    run_month(month, source_cells, sink_cells, min_frac)
  }
}
