## GDE_02_CalculateExtraET.R
# This script is intended to reconcile the scales of the ET and P datasets to 
# calculate the 'extra' ET (ET in excess of precipitation) at the scale of the ET data.
#
# Annual extra ET is calculated as multi-year average (ET-P)
# That is distributed to monthly extra ET based on the % of annual ET-P that occurs in each month.
#
# BEFORE RUNNING THIS SCRIPT:
#  Manually copy the files from "G:/My Drive/GEEexports/GlobalGDEs"
#   to "C:/Users/s947z036/OneDrive - University of Kansas/Research/MoistureRecycling/RohdeEtAl_GDEs/SZ_Analysis"
#
# Author: Sam Zipper
# Date: 11/9/2026

# set up workspace
source(file.path("code", "paths+packages.R"))

## inputs and outputs
data_directory <- file.path(path_data, "RohdeEtAl_GDEs", "SZ_Analysis")
input_directory <- file.path(data_directory, "GEEexport_ET+Pdata")

# Process each tile and period separately; subset this list for selected test tiles.
manifest_files <- list.files(input_directory,
                             pattern = "^GDE_monthly_climate_manifest_[0-9]{4}_[0-9]{4}_[ns][0-9]+[ew][0-9]+\\.csv$",
                             full.names = FALSE)
if (length(manifest_files) == 0) stop("No tile-based climate manifests found.")

for (m in manifest_files){
  # m <- manifest_files[1]
  
  # extract name of tile
  tile <- str_sub(str_split(m, pattern = "_")[[1]][7], end = -5)
  
  # load ET and P
  r_ET <- rast(file.path(input_directory, paste0("GDE_ET_monthly_climatology_2008_2017_", tile, ".tif")))
  r_P <- rast(file.path(input_directory, paste0("GDE_P_monthly_climatology_2008_2017_", tile, ".tif")))
  
  # set names as the month
  names(r_ET) <- str_sub(names(r_ET), start = -2, end = -1)
  names(r_P) <- str_sub(names(r_P), start = -2, end = -1)
  
  # confirm names are matching and in same order
  if(sum(names(r_ET) == names(r_P)) != 12) stop("Monthly raster name mismatch")
  
  # resample P to ET grid
  r_P.rs <- resample(r_P, r_ET, method = "near")  # use nearest neighbor instead of interpolation
  crs(r_P.rs) <- crs(r_ET)
  
  # # confirm geometry match
  # res(r_ET)
  # res(r_P.rs)
  # ext(r_ET)
  # ext(r_P.rs)
  
  # sum to annual
  r_ET_ann <- app(r_ET, fun = "sum")
  r_P.rs_ann <- app(r_P.rs, fun = "sum")
  
  # calculate total annual ET - P
  r_ET.P_ann <- r_ET_ann - r_P.rs_ann
  
  # calculate ET - P that can be attributed to each month
  r_ET.P_mo <- r_ET - r_P.rs
  
  # set anything less than 0 to 0
  r_ET.P_mo_0cap <- r_ET.P_mo
  r_ET.P_mo_0cap[r_ET.P_mo < 0] <- 0
  
  # sum to annual ignoring negatives
  r_ET.P_ann_0cap <- sum(r_ET.P_mo_0cap)
  
  # get percent of each month
  r_ET.P_mo_0cap_prc <- r_ET.P_mo_0cap / r_ET.P_ann_0cap
  
  # save two output products:
  #  annual ET-P ("extra ET")
  #  monthly percent of annual ET-P, 0-capped (attribution of percent of extra ET to each month)
  # limit to 3 decimal places to save space
  terra::writeRaster(r_ET.P_ann, 
                     file.path(data_directory, "GDE_02_CalculateExtraET", paste0("GDE_02_AnnualET-P_", tile, ".tif")), 
                     datatype = "INT4S", 
                     scale = 0.001, 
                     offset = 0,
                     overwrite = TRUE)
  terra::writeRaster(r_ET.P_mo_0cap_prc, 
                     file.path(data_directory, "GDE_02_CalculateExtraET", paste0("GDE_02_MonthlyET-P_0capPrc_", tile, ".tif")), 
                     datatype = "INT4S", 
                     scale = 0.001, 
                     offset = 0,
                     overwrite = TRUE)
  
}
