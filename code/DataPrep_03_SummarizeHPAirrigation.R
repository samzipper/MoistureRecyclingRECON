## DataPrep_03_SummarizeHPAirrigation.R
# This script will summarize the 2008-2017 irrigation data for the LANID region
# to identify grid cells that were, on average, irrigated in the HPA region.

source(file.path("code", "paths+packages.R"))

## load data
gee_out <- 
  read_csv(file.path("data", "ERA5-CONUS_LANID_PolygonCounts_Annual.csv")) |> 
  dplyr::select(-`system:index`, -`.geo`)
sf_era_hpa <- st_read(file.path("data", "ERA5_grid", "ERA5_grid_HPA.shp"))

## summarize: mean annual irrigation fraction
irr_avgAnnual <- 
  gee_out |> 
  subset(cellid %in% sf_era_hpa$cellid) |> 
  group_by(cellid) |> 
  summarize(irr_avgPrc = mean(irr_prc))

## join with grid cells
sf_era_hpa_irr <- left_join(sf_era_hpa, irr_avgAnnual, by = "cellid") |> 
  mutate(irr_prc_cut = cut(irr_avgPrc, breaks = c(0, 0.1, 0.2, 0.4, 1)))

# inspect
ggplot() +
  geom_sf(data = sf_era_hpa_irr,
          aes(fill = irr_prc_cut))

# inspect
ggplot() +
  geom_sf(data = sf_era_hpa_irr,
          aes(fill = state))

head(sf_era_hpa_irr)

# save output
sf_era_hpa_irr |> 
  dplyr::select(cellid, irr_avgPrc, state) |> 
  st_write(file.path("data", "ERA5-HPA_LANID_AverageIrrigation.shp"),
            append = FALSE)
