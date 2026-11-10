## GDE_03_CalculateETgw.R
# This script is intended to integrate extra ET from the GDE_02_CalculateExtraET.R script
# with the GDE presence dataset to estimate the contribution of groundwater to ET at
# each grid cell.
#
# Author: Sam Zipper
# Date: 11/9/2026

# set up workspace
source(file.path("code", "paths+packages.R"))

# paths
path_extraET <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_02_CalculateExtraET")
path_GDEtiles <- file.path(path_data, "RohdeEtAl_GDEs", "D_Data_GDEtiles_5deg")
path_output <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis", "GDE_03_CalculateETgw")

## function here with following inputs
ETfromGW <- function(tile, path_GDEtiles, path_extraET, path_output){

  
  # load datasets
  r_GDE <- rast(file.path(path_GDEtiles, paste0(tile, ".tif")))  # 0 = not analyzed, 1 = probably GDE, 2 = probably not GDE
  r_extraET_mo <- rast(file.path(path_extraET, paste0("GDE_02_MonthlyET-P_0capPrc_", tile, ".tif")))
  
  ## compare summary approaches
  # r_GDE_avgProb <- project(r_GDE["probability"], r_extraET_mo, method = "mean") # this seems to underestimate GDE area
  # r_GDE_mode <- project(r_GDE["classification"], r_extraET_mo, method = "mode")
  # r_GDE_mat <- as.matrix(r_GDE["classification"], wide = T)
  # r_GDE_mode_mat <- as.matrix(r_GDE_mode, wide = T)
  # r_GDE_avgProb_mat <- as.matrix(r_GDE_avgProb, wide = T)
  # 
  # sum(r_GDE_mat == 1)/225 # to account for different cell area
  # sum(r_GDE_mode_mat == 1, na.rm = T)
  # sum(r_GDE_avgProb_mat > 50, na.rm = T)
  # sum(r_GDE_avgProb_mat > 50, na.rm = T)
  
  
  # implement selected summary
  r_GDE_mask <- 
    project(r_GDE["classification"], r_extraET_mo, method = "mode") |>  # summarize by mode
    subst(2, 0) # replace 2 with 0, so GDE = 1 and not analyzed or not GDE both = 0
  
  # multiple to assume all scaled extra ET in GDE cells is from GW
  r_ETfromGW_mo <- r_extraET_mo*r_GDE_mask
  
  
  terra::writeRaster(r_ETfromGW_mo, 
                     file.path(path_output, paste0("GDE_03_ETgw_prc_", tile, ".tif")), 
                     datatype = "INT4S", 
                     scale = 0.001, 
                     offset = 0,
                     overwrite = TRUE)
  
  return(tile)
}

# get list of tiles
tiles_all <- list.files(path_extraET, pattern = "GDE_02_MonthlyET-P_0capPrc_*")
tiles <- str_sub(sapply(str_split(tiles_all, pattern = "_"), `[`, 5), end = -5)

# run in parallel
future::plan(future::multisession, workers = 8)

future.apply::future_lapply(
  tiles,
  ETfromGW,
  path_GDEtiles = path_GDEtiles,
  path_extraET = path_extraET,
  path_output = path_output,
  future.seed = NULL
)
