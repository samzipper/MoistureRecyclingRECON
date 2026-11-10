## GDE_01a_CheckMonthlyETandP.R
# Check tile-based climatologies from GDE_01_ExportMonthlyETandP.js.
# No yearly stacks are required. Means cannot be independently recomputed without them.
# Numerical checks use native grids; only plot rasters are aggregated.
#
# Author: Sam Zipper + Copilot
# Date: 11/9/2026

# set up workspace
source(file.path("code", "paths+packages.R"))
theme_set(theme_scz() + theme(text = element_text(family = "sans")))

## inputs and outputs
data_directory <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis")
input_directory <- file.path(data_directory, "GEEexport_ET+Pdata")
output_directory <- file.path(input_directory, "GDE_01a_CheckMonthlyETandP")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
sample_size <- 2000
map_max_dimension <- 160

# Process each tile and period separately; subset this list for selected test tiles.
manifest_files <- list.files(input_directory,
                             pattern = "^GDE_monthly_climate_manifest_[0-9]{4}_[0-9]{4}_[ns][0-9]+[ew][0-9]+\\.csv$",
                             full.names = TRUE)
if (length(manifest_files) == 0) stop("No tile-based climate manifests found.")

## helpers for reading stacks and summarizing native-grid pixels
find_stack_files <- function(prefix) {
  list.files(input_directory,
             pattern = paste0("^", prefix, "(-[0-9]+-[0-9]+)?\\.tif$"),
             full.names = TRUE)
}

read_climate_stack <- function(files, expected_names, sf_region, no_data) {
  r <- if (length(files) == 1) rast(files) else merge(sprc(lapply(files, rast)))
  if (nlyr(r) != 12 || !identical(names(r), expected_names)) {
    stop("Unexpected monthly band names or layer count in ", basename(files[1]))
  }
  r <- subst(r, from = no_data, to = NA)
  region <- project(vect(sf_region), crs(r))
  mask(crop(r, region), region)
}

summarize_climatology <- function(r, coverage, year_count, sf_region) {
  area <- cellSize(r[[1]], unit = "m") |>
    mask(project(vect(sf_region), crs(r)))
  valid_area <- global((!is.na(r)) * area, "sum", na.rm = TRUE)[, 1]
  region_area <- global(area, "sum", na.rm = TRUE)[1, 1]
  df <- tibble(
    month = 1:12,
    mean_mm = global(r * area, "sum", na.rm = TRUE)[, 1] / valid_area,
    min_mm = global(r, "min", na.rm = TRUE)[, 1],
    max_mm = global(r, "max", na.rm = TRUE)[, 1],
    valid_fraction = valid_area / region_area,
    negative_cells = global(r < 0, "sum", na.rm = TRUE)[, 1],
    coverage_available = !is.null(coverage),
    expected_years = year_count,
    coverage_min = NA_real_, coverage_max = NA_real_,
    complete_fraction = NA_real_, invalid_count_cells = NA_real_,
    missing_count_cells = NA_real_, mask_mismatch_cells = NA_real_
  )
  if (!is.null(coverage)) {
    compareGeom(r, coverage, stopOnError = TRUE)
    complete <- ifel(is.na(coverage), FALSE, coverage == year_count)
    invalid <- coverage < 0 | coverage > year_count | coverage != round(coverage)
    df$coverage_min <- global(coverage, "min", na.rm = TRUE)[, 1]
    df$coverage_max <- global(coverage, "max", na.rm = TRUE)[, 1]
    df$complete_fraction <- global(complete * area, "sum", na.rm = TRUE)[, 1] / region_area
    df$invalid_count_cells <- global(invalid, "sum", na.rm = TRUE)[, 1]
    df$missing_count_cells <- global(is.na(coverage) * (area > 0), "sum", na.rm = TRUE)[, 1]
    df$mask_mismatch_cells <- global(xor(!is.na(r), complete) * (area > 0), "sum", na.rm = TRUE)[, 1]
  }
  df |>
    mutate(across(c(mean_mm, min_mm, max_mm, coverage_min, coverage_max),
                  ~replace(.x, !is.finite(.x), NA_real_)))
}

## load each tile's manifest, boundary, climatologies, and optional valid-year rasters
for (manifest_file in manifest_files) {
  df_manifest <- read_csv(manifest_file, show_col_types = FALSE)
  tile_names <- unique(df_manifest$tile_name)
  if (length(tile_names) != 1 || is.na(tile_names) || !nzchar(tile_names)) {
    stop("Expected one tile_name in ", basename(manifest_file))
  }
  tile_name <- tile_names[1]
  start_year <- as.integer(format(as.Date(df_manifest$start_date[1]), "%Y"))
  end_year <- as.integer(format(as.Date(df_manifest$end_date_exclusive[1]), "%Y")) - 1L
  year_count <- end_year - start_year + 1L
  period <- paste(start_year, end_year, sep = "_")
  tile_label <- paste(tile_name, period, sep = " | ")
  tile_directory <- file.path(output_directory, paste(tile_name, period, sep = "_"))
  dir.create(tile_directory, recursive = TRUE, showWarnings = FALSE)
  message("Checking ", tile_label)

  sf_region <- st_read(df_manifest$region_geojson[1], quiet = TRUE) |>
    st_make_valid()
  if (any(st_geometry_type(sf_region) == "GEOMETRYCOLLECTION")) {
    sf_region <- st_collection_extract(sf_region, "POLYGON", warn = FALSE)
  }
  if (is.na(st_crs(sf_region))) st_crs(sf_region) <- 4326

  summary_list <- list()
  sample_list <- list()
  map_list <- list()
  grid_list <- list()
  inventory_list <- list()

  for (variable in c("ET", "P")) {
    df_product <- df_manifest |>
      filter(.data$band == .env$variable, product == "monthly_climatology")
    if (nrow(df_product) != 1) {
      stop("Expected one ", variable, " climatology record in ", basename(manifest_file))
    }
    prefix <- df_product$prefix[1]
    coverage_band <- paste0(variable, "_valid_years")
    coverage_prefix <- paste("GDE", coverage_band, period, tile_name, sep = "_")
    files <- find_stack_files(prefix)
    coverage_files <- find_stack_files(coverage_prefix)
    inventory_list[[variable]] <- tibble(
      tile_name = tile_name, period = period, variable = variable,
      product = c("monthly_climatology", "valid_years"),
      prefix = c(prefix, coverage_prefix),
      file_count = c(length(files), length(coverage_files)),
      status = if_else(file_count > 0, "present", "missing")
    )
    if (length(files) == 0) {
      warning("Missing climatology: ", prefix, "; skipping this variable.")
      next
    }
    r <- read_climate_stack(files, sprintf("%s_%02d", variable, 1:12),
                            sf_region, df_manifest$no_data[1])
    coverage <- NULL
    if (length(coverage_files) > 0) {
      coverage <- read_climate_stack(coverage_files, sprintf("%s_%02d", coverage_band, 1:12),
                                     sf_region, df_manifest$no_data[1])
    } else {
      warning("Missing coverage: ", coverage_prefix, "; coverage checks unavailable.")
    }
    summary_list[[variable]] <- summarize_climatology(r, coverage, year_count, sf_region) |>
      mutate(tile_name = tile_name, period = period, variable = variable, .before = 1)
    grid_list[[variable]] <- tibble(
      tile_name = tile_name, period = period, variable = variable,
      rows = nrow(r), columns = ncol(r), bands = nlyr(r),
      x_resolution = res(r)[1], y_resolution = res(r)[2],
      crs = crs(r), expected_crs = df_product$crs[1]
    )

    # Use a regular spatial sample for boxplots; retain missing months explicitly.
    df_sample <- spatSample(r, size = sample_size, method = "regular",
                            na.rm = FALSE, as.df = TRUE)
    names(df_sample) <- sprintf("%02d", 1:12)
    sample_list[[variable]] <- df_sample |>
      pivot_longer(everything(), names_to = "month", values_to = "value_mm") |>
      mutate(variable = variable, month = as.integer(month))

    # Downsample only the map previews, without aligning ET to the precipitation grid.
    map_factor <- max(1, ceiling(max(nrow(r), ncol(r)) / map_max_dimension))
    r_map <- if (map_factor > 1) aggregate(r, fact = map_factor, fun = "mean", na.rm = TRUE) else r
    names(r_map) <- sprintf("%02d", 1:12)
    map_list[[variable]] <- as.data.frame(r_map, xy = TRUE, na.rm = FALSE) |>
      pivot_longer(-c(x, y), names_to = "month", values_to = "value_mm") |>
      mutate(variable = variable, month = as.integer(month))
  }

  ## save per-tile inventories, grid metadata, and monthly diagnostic checks
  df_inventory <- bind_rows(inventory_list)
  write_csv(df_inventory, file.path(tile_directory, "GDE_01a_FileInventory.csv"))
  if (length(summary_list) == 0) {
    warning("No climatologies available for ", tile_label, "; inventory saved only.")
    next
  }
  df_summary <- bind_rows(summary_list)
  df_samples <- bind_rows(sample_list)
  df_maps <- bind_rows(map_list)
  write_csv(df_summary, file.path(tile_directory, "GDE_01a_ClimatologySummary.csv"))
  write_csv(bind_rows(grid_list), file.path(tile_directory, "GDE_01a_RasterGrids.csv"))
  print(df_summary)
  if (any(df_summary$negative_cells > 0, na.rm = TRUE) ||
      any(df_summary$invalid_count_cells > 0, na.rm = TRUE) ||
      any(df_summary$missing_count_cells > 0, na.rm = TRUE) ||
      any(df_summary$mask_mismatch_cells > 0, na.rm = TRUE) ||
      any(df_summary$valid_fraction == 0)) {
    warning("Some climate checks failed for ", tile_label, ". Inspect the summary CSV.")
  }

  ## monthly climatology means and sampled spatial distributions
  p_means <- ggplot(df_summary, aes(x = month, y = mean_mm, color = variable, group = variable)) +
    geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    scale_color_manual(values = c(ET = col.cat.red, P = col.cat.blu)) +
    labs(x = NULL, y = "Area-weighted Monthly Depth [mm]", color = NULL,
         title = "Monthly ET and Precipitation Climatology", subtitle = tile_label)
  ggsave(file.path(tile_directory, "GDE_01a_MonthlyMeans.png"), plot = p_means,
         width = 180, height = 100, units = "mm")

  p_distributions <- ggplot(df_samples,
                            aes(x = factor(month, levels = 1:12, labels = month.abb), y = value_mm)) +
    geom_boxplot(outlier.size = 0.3, na.rm = TRUE) +
    facet_wrap(~variable, ncol = 1) +
    labs(x = NULL, y = "Monthly Depth [mm]", title = "Climatology Spatial Distributions",
         subtitle = paste(tile_label, "| Regular samples on each native grid"))
  ggsave(file.path(tile_directory, "GDE_01a_MonthlyDistributions.png"), plot = p_distributions,
         width = 180, height = 130, units = "mm")

  ## compare retained climatology area with area having every requested year
  df_coverage <- df_summary |>
    select(variable, month, valid_fraction, complete_fraction) |>
    pivot_longer(c(valid_fraction, complete_fraction), names_to = "measure", values_to = "fraction") |>
    filter(!is.na(fraction))
  p_coverage <- ggplot(df_coverage, aes(x = month, y = fraction, color = variable,
                                       linetype = measure, group = interaction(variable, measure))) +
    geom_line() + geom_point() +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_color_manual(values = c(ET = col.cat.red, P = col.cat.blu)) +
    scale_linetype_manual(values = c(valid_fraction = "solid", complete_fraction = "dashed"),
                          labels = c(complete_fraction = "All years available", valid_fraction = "Climatology retained")) +
    labs(x = NULL, y = "Fraction of Tile Area", color = NULL, linetype = NULL,
         title = "Monthly Climatology Coverage", subtitle = tile_label)
  ggsave(file.path(tile_directory, "GDE_01a_MonthlyCoverage.png"), plot = p_coverage,
         width = 180, height = 110, units = "mm")

  ## monthly climatology maps, using the same red-to-blue scale as GEE
  for (variable in names(map_list)) {
    p_map <- ggplot(filter(df_maps, .data$variable == .env$variable), aes(x = x, y = y, fill = value_mm)) +
      geom_raster() +
      facet_wrap(~factor(month, levels = 1:12, labels = month.abb), ncol = 4) +
      scale_fill_viridis_c(oob = scales::squish, na.value = "transparent",
                           name = "Monthly Depth [mm]") +
      coord_equal(ratio = 1 / cos(mean(st_bbox(st_transform(sf_region, 4326))[c("ymin", "ymax")]) * pi / 180)) +
      labs(x = NULL, y = NULL, title = paste(variable, "Monthly Climatology"), subtitle = tile_label) +
      theme(legend.position = "bottom", axis.text = element_blank(), axis.ticks = element_blank())
    ggsave(file.path(tile_directory, paste0("GDE_01a_", variable, "_MonthlyMaps.png")), plot = p_map,
           width = 190, height = 140, units = "mm")
  }
}
