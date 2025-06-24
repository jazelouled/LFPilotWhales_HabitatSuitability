library(raster)
library(dplyr)
library(lubridate)
library(readr)
library(doParallel)
library(foreach)
library(here)
library(rgeos)
library(sp)

# 1. Load presence-absence data
file_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2025_LFPilotWhales_RdeStephanis/LFPilotWhales_HabitatSuitability/000inputOutput/00output/00tracking/04PresAbs/presAbs_with_env.csv"
presAbs <- read_csv(file_path)
presAbs$date <- ymd_hms(presAbs$date)
presAbs$day <- format(presAbs$date, "%Y-%m-%d")

# 2. Load reference raster to set extent/resolution
tif_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2025_LFPilotWhales_RdeStephanis/LFPilotWhales_HabitatSuitability/000inputOutput/00input/daily_tifs"
tif_files <- list.files(tif_dir, recursive = TRUE, full.names = TRUE)
ref_raster <- raster(tif_files[1])

# 3. Load and process bathymetry
bathy_path <- "~/Dropbox/2022_SouthHemisphere_SDM/bathymetry/GEBCO_2014_2D.nc"
bathy <- raster(bathy_path)
bathy_crop <- crop(bathy, ref_raster)
bathy_crop[bathy_crop >= 0] <- NA
bathy_crop_rsp <- resample(bathy_crop, ref_raster)

# 4. Derive slope
slope_raster <- terrain(bathy_crop_rsp, opt = "slope", unit = "degrees", neighbors = 8)

# 5. Derive distance to coast (from coastline into ocean)
# Land = 1, Sea = NA
land_mask <- bathy_crop_rsp
land_mask[!is.na(land_mask)] <- 1
land_mask[is.na(land_mask)] <- NA

# Coastline as boundary between land and sea
coastline <- boundaries(land_mask, type = "inner", classes = TRUE)

# Distance from coast into the sea
dist_raster <- distance(coastline)
dist_raster[!is.na(land_mask)] <- NA   # keep only ocean distances
dist2coast <- dist_raster / 1000       # convert to km

# 6. Parallel backend setup
ncores <- 6
cl <- makeCluster(ncores)
registerDoParallel(cl)

# 7. Get unique days from tracking data
all_days <- unique(na.omit(presAbs$day))

# 8. Loop through each date
results_list <- foreach(d = all_days, .combine = rbind, .packages = c("raster", "dplyr", "lubridate", "sp")) %dopar% {
  message("Processing date: ", d)
  df_day <- presAbs[presAbs$day == d, ]
  
  # Build individual paths for each environmental variable
  year_str <- unique(year(df_day$date))
  thetao_path <- file.path(tif_dir, year_str, paste0("thetao_", d, ".tif"))
  mlotst_path <- file.path(tif_dir, year_str, paste0("mlotst_", d, ".tif"))
  uo_path     <- file.path(tif_dir, year_str, paste0("uo_", d, ".tif"))
  vo_path     <- file.path(tif_dir, year_str, paste0("vo_", d, ".tif"))
  
  # Check if all files exist
  all_exist <- all(file.exists(thetao_path, mlotst_path, uo_path, vo_path))
  
  if (all_exist) {
    # Load rasters
    thetao_raster <- raster(thetao_path)
    mlotst_raster <- raster(mlotst_path)
    uo_raster     <- raster(uo_path)
    vo_raster     <- raster(vo_path)
    
    # Stack all environmental layers
    env_stack <- stack(
      thetao_raster,   # sst
      mlotst_raster,
      uo_raster,
      vo_raster,
      bathy_crop_rsp,
      slope_raster,
      dist2coast
    )
    names(env_stack) <- c("sst", "mlotst", "uo", "vo", "bathy", "slope", "dist2coast")
    
    # Extract environmental values at points
    coords <- df_day %>% select(lon, lat) %>% as.data.frame()
    sp_points <- SpatialPoints(coords, proj4string = CRS(projection(env_stack)))
    extracted_vals <- raster::extract(env_stack, sp_points, buffer = 15000, fun = mean, na.rm = TRUE)
    
    # Combine extracted values with data
    df_day <- bind_cols(df_day, as.data.frame(extracted_vals))
    
  } else {
    warning("Missing environmental files for date: ", d)
    df_day$sst <- NA
    df_day$mlotst <- NA
    df_day$uo <- NA
    df_day$vo <- NA
    df_day$bathy <- NA
    df_day$slope <- NA
    df_day$dist2coast <- NA
  }
  
  df_day
}

# 9. Stop cluster
stopCluster(cl)

# 10. Save final output
output_path <- here("000inputOutput/00output/00tracking/04PresAbs/presAbs_with_env.csv")
write_csv(results_list, output_path)
