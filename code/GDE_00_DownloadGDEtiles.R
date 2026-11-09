# GDE 5 degree tiles downloader
# Downloads all GDE tiles from Google Cloud Storage following the instructions
# in: https://storage.googleapis.com/gde-globalv6/5deg_v6_snapped/
#
# Author: Sam + Copilot
# Date: 11/9/2026

library(sf)
library(terra)
library(httr)

# Set up paths
gde_tiles_dir <- "C:\\Users\\s947z036\\OneDrive - University of Kansas\\Research\\MoistureRecycling\\RohdeEtAl_GDEs\\D_Data_GDEtiles_5deg"
sqlite_file <- file.path(gde_tiles_dir, "tile_extents.sqlite")
base_url <- "https://storage.googleapis.com/gde-globalv6/5deg_v6_snapped/"

# Create directory if it doesn't exist
if (!dir.exists(gde_tiles_dir)) {
  dir.create(gde_tiles_dir, recursive = TRUE)
}

# Read the tile extents to get tile names
cat("Reading tile extents from SQLite file...\n")
tile_extents <- st_read(sqlite_file)
tile_names <- unique(tile_extents$name)

cat("Found", length(tile_names), "unique tiles to download.\n\n")

# Download each tile
for (i in seq_along(tile_names)) {
  tile_name <- tile_names[i]
  tile_url <- paste0(base_url, tile_name, ".tif")
  output_file <- file.path(gde_tiles_dir, paste0(tile_name, ".tif"))
  
  # Skip if already downloaded
  if (file.exists(output_file)) {
    cat("[", i, "/", length(tile_names), "] Tile", tile_name, "already exists. Skipping...\n")
    next
  }
  
  cat("[", i, "/", length(tile_names), "] Downloading", tile_name, "...")
  
  tryCatch({
    # Download the file
    response <- GET(tile_url, write_disk(output_file, overwrite = FALSE), progress())
    
    # Check if download was successful
    if (status_code(response) == 200) {
      cat(" ✓ Success\n")
    } else {
      cat(" ✗ Failed (HTTP", status_code(response), ")\n")
      if (file.exists(output_file)) {
        file.remove(output_file)
      }
    }
  }, error = function(e) {
    cat(" ✗ Error:", e$message, "\n")
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
  })
}

cat("\nDownload complete! Tiles saved to:\n")
cat(gde_tiles_dir, "\n")
