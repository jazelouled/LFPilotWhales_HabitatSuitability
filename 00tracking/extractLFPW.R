library(raster)
library(dplyr)
library(lubridate)
library(readr)
library(doParallel)
library(foreach)

# 1. Load presence-absence data
file_path <- "/Users/jazelouled-cheikhbonan/Dropbox/2025_LFPilotWhales_RdeStephanis/LFPilotWhales_HabitatSuitability/000inputOutput/00output/00tracking/04PresAbs/presAbsCells.csv"
presAbs <- read_csv(file_path)

# 2. Convert date to POSIXct and extract date-only string
presAbs$date <- ymd_hms(presAbs$date)
presAbs$day <- format(presAbs$date, "%Y-%m-%d")

# 3. Prepare rasters
tif_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2025_LFPilotWhales_RdeStephanis/LFPilotWhales_HabitatSuitability/000inputOutput/00input/daily_tifs"
tif_copy <- raster::raster(list.files(tif_dir, recursive = TRUE, full.names = TRUE)[1])

# Bathymetry preprocessing
bathy_path <- "~/Dropbox/2022_SouthHemisphere_SDM/bathymetry/GEBCO_2014_2D.nc"
bathy <- raster(bathy_path)
bathy_crop <- crop(bathy, tif_copy)
bathy_crop[bathy_crop >= 0] <- NA
bathy_crop_rsp <- resample(bathy_crop, tif_copy)

# 4. Prepare parallel backend
ncores <- 6
cl <- makeCluster(ncores)
registerDoParallel(cl)

# 5. Unique days to loop over
all_days <- unique(na.omit(presAbs$day))

# 6. Parallel loop using foreach
results_list <- foreach(d = all_days, .combine = rbind, .packages = c("raster", "dplyr", "lubridate")) %dopar% {
  message("Processing date: ", d)
  
  idx <- which(presAbs$day == d)
  df_day <- presAbs[idx, ]
  
  sst_file <- file.path(tif_dir, unique(year(df_day$date)), paste0("thetao_", d, ".tif"))
  
  if (file.exists(sst_file)) {
    sst_raster <- raster(sst_file)
    
    # Stack SST + Bathymetry
    env_stack <- stack(sst_raster, bathy_crop_rsp)
    names(env_stack) <- c("sst", "bathy")
    
    coords <- df_day %>% select(lon, lat) %>% as.data.frame()    
    points <- SpatialPoints(coords, proj4string = CRS(projection(env_stack)))
    values <- raster::extract(env_stack, points, buffer = 15000, fun = mean, na.rm = TRUE)
    
    df_day$sst <- values[, "sst"]
    df_day$bathy <- values[, "bathy"]
    
  } else {
    warning("Missing SST file for date: ", d)
    df_day$sst <- NA
    df_day$bathy <- NA
  }
  
  return(df_day)
}

# 7. Stop cluster
stopCluster(cl)

# 8. Save final output
write_csv(results_list, here::here("000inputOutput/00output/00tracking/04PresAbs/presAbs_with_env.csv"))
