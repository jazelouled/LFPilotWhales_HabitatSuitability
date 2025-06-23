#-------------------------------------------------------------------------------------
# 03_regularize_ssm     Interpolate tracks into regular time steps using foieGras
#-------------------------------------------------------------------------------------
library(here)
library(doParallel)
library(ggplot2)
library(aniMotum)
library(dplyr)
library(availability)
library(TMB)
library(rlang)
library(Rtools)
library(vctrs)


# Main steps are:
# - Regularize tracks

#---------------------------------------------------------------
# Prepare cluster
#---------------------------------------------------------------
cl <- makeCluster(cores)
registerDoParallel(cl)

#---------------------------------------------------------------
# 1. Set data repository
#---------------------------------------------------------------
here()
indir <- here::here("000inputOutput/00output/00tracking/L1_locations")
outdir <- here::here("000inputOutput/00output/00tracking/L2_locations")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)




#---------------------------------------------------------------
# 2. Import data
#---------------------------------------------------------------

# import all location files
loc_files <- list.files(indir, full.names = TRUE, pattern = "L1_locations.csv")
df <- readTrack(loc_files)
df$Date <- as.POSIXct(df$Date, format = "%Y-%m-%d %H:%M:%S")



#---------------------------------------------------------------
# 2. Select trips to run the SSM
#---------------------------------------------------------------

# summarize data per trip 
trips <- summarizeTrips(df)
tripIDs <- trips$PTT 




#---------------------------------------------------------------
# 3. Regularize each track using a SSM
#---------------------------------------------------------------
# foreach(i=tags, .packages=c("dplyr", "ggplot2", "aniMotum", "stringr")) %dopar% {

  
  for(i in 1:length(tripIDs)){
 
  print(paste("Processing deploymentID", tripIDs[i]))

  # subset data
  # filter by id and selected trips
    data <- dplyr::filter(df, PTT == tripIDs[i], tripIDs[i] %in% tripIDs)
  
  
  ###### State-Space Model
  
  # convert to foieGras format
    if (tag_type == "PTT") {
      if (any(grepl("Murcia", data$PTT, ignore.case = TRUE))) {
        indata <- data %>%
          dplyr::select(PTT, date, lc, Longitude, Latitude) %>%
          dplyr::rename(id = PTT)
      } else {
        indata <- data %>%
          dplyr::select(PTT, date, Quality, Longitude, Latitude) %>%
          dplyr::rename(id = PTT)
      }
    }  
  # fit SSM
  # we turn sdafilter off because we previously filtered data
  # we run the model with multiple trips at once
  colnames(indata) <- c("id", "date", "lc", "lon", "lat")
  
  fit <- fit_ssm(indata,
                 # vmax= filt_vmax,
                 # ang = filt_ang,
                 # distlim = filt_distlim,
                 model = "crw",
                 time.step = 4,
                 control = ssm_control(verbose = 0))
  


  # get fitted locations
  # segments that did not converge were not consider
  data <- data.frame(grab(fit, what = "predicted", as_sf = FALSE))
  data <- data %>% arrange(date)
  # check if points on land
  data$onland <- point_on_land(lat = data$lat, lon = data$lon, land = land)
  data$Date <- data$date
  data$argosfilter <- NA
  
  # export track data into individual folder at output path
  outpath <-  outdir
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)
  outfile <- paste0(outpath, "/", tripIDs[i], "_L2_locations.csv")
  write.csv(data, outfile, row.names = FALSE)
  
  # export convergence status
  convergence <- data.frame(trip = fit$id, converged = fit$converged)
  outfile <- paste0(outpath, "/", tripIDs[i], "_L2_convergence.csv")
  write.csv(convergence, outfile, row.names = FALSE)
  
  # plot figures
  
  data$Date <- data$date
  data$argosfilter <- NA
  p <- mapL1(data = data)
  outfile <- paste0(outpath, "/", tripIDs[i], "_L2_locations.png")
  ggsave(outfile, p, width=30, height=15, units = "cm")
  
 }
  
stopCluster(cl) # Stop cluster

#---------------------------------------------------------------
# 4. Summarize processed data
#---------------------------------------------------------------

# import all location files
loc_files <- list.files(outdir, full.names = TRUE, pattern = "L2_locations.csv")
df <- readTrack(loc_files)

# import convergence files
loc_files <- list.files(outdir, full.names = TRUE, pattern = "L2_convergence.csv")
data_proc <- lapply(loc_files, read.csv) %>% rbindlist

# summarize data per trip
tripstats <- summarizeTrips2(df)

# combine track data summary and convergence status
comb <- merge(tripstats, data_proc, by=c("trip"))

# export table
out_file <- paste0(outdir, "summary_ssm.csv")
write.csv(comb, out_file, row.names = FALSE)


print("Regularization ready")



AA <-df[df$deploymentID == "CARCAR_64702a_UNEXE_20060825"]







AA$trip <- timedif.segment(AA$time, thrs = trip_time_gap)



AA$timedif <- timedif(AA$time)








