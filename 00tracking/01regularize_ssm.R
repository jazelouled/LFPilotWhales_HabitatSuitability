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


install.packages("rlang")
install.packages("Rtools")
install.packages("vctrs", dependencies = TRUE)

remove.packages("vctrs")

install.packages("TMB", type = "source")

install.packages("TMB")

# install from my R-universe repository
install.packages("aniMotum", 
                 repos = c("https://cloud.r-project.org",
                           "https://ianjonsen.r-universe.dev"),
                 dependencies = TRUE)


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
indir <- input_data
outdir <- here::here(output_data, "caret/01tracking/L2_loc_NOinternesting")
outdir <- here::here(output_data, "caret/01tracking/ForagingUpwelling_filtering")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#---------------------------------------------------------------
# 2. Import data
#---------------------------------------------------------------

# import all location files
loc_files <- list.files(indir, full.names = TRUE, pattern = "L1_loc.csv")
df <- readTrack(loc_files)
df$date <- parse_date_time(df$time, orders = c("%d/%m/%Y %H:%M", "%Y-%m-%d %H:%M:%S"))


# Deletion of internesting locations (1 month)

# Calculate the number of days from the first observation for each ID
tracking_data <- df %>%
    group_by(tripID) %>%
  mutate(DaysFromStart = as.integer(as.Date(date) - min(as.Date(date))))

# Filter out the first 30 days of data for each individual
filtered_tracking_data <- tracking_data %>%
  filter(DaysFromStart > 30)

# Remove the temporary "DaysFromStart" column
filtered_tracking_data <- select(filtered_tracking_data, -DaysFromStart)

# Print the resulting data frame
print(filtered_tracking_data)

df <- filtered_tracking_data



#---------------------------------------------------------------
# 2. Select trips to run the SSM
#---------------------------------------------------------------

# summarize data per trip
trips <- summarizeTrips(df)

# filter trips
trips <- dplyr::filter(trips,
                duration_h >= sel_min_dur,
                n_loc >= sel_min_loc,
                distance_km >= sel_min_dist, 
                !tripID %in% sel_exclude)
 

trips <- dplyr::filter(trips,
                       duration_h >= sel_min_dur,
                       n_loc >= sel_min_loc,
                       distance_km >= sel_min_dist)



# tripIDs <- trips$tripID
# deploymentID <- trips$deploymentID

tripIDs <- trips$tripID



#---------------------------------------------------------------
# 3. Regularize each track using a SSM
#---------------------------------------------------------------
# foreach(i=tags, .packages=c("dplyr", "ggplot2", "aniMotum", "stringr")) %dopar% {

  
  for(i in 1:length(tripIDs)){
 
  print(paste("Processing deploymentID", tripIDs[i]))

  # subset data
  # filter by id and selected trips
  
  # data <- filter(df, id == i, trip %in% trips$trip)
  data <- dplyr::filter(df, tripID == tripIDs[i], tripIDs[i] %in% tripIDs)
  
  
  ###### State-Space Model
  
  # convert to foieGras format
  # if(tag_type == "GPS") indata <- data %>% dplyr::select(organismID, deploymentID, tripID, date, argosLC, longitude, latitude) %>% dplyr::rename(id = deploymentID)
  if(tag_type == "GPS") indata <- data %>% dplyr::select(organismID, deploymentID, tripID, date, argosLC, longitude, latitude) 
  indata$id <- indata$tripID
  # if(tag_type == "PTT") indata <- data %>% dplyr::select(deploymentID, date, argosLC, longitude, latitude, smaj, smin, eor) %>% dplyr::rename(id = deploymentID)
  
  # fit SSM
  # we turn sdafilter off because we previously filtered data
  # we run the model with multiple trips at once
  
  colnames(indata) <- c("organismID", "deploymentID", "tripID", "date", "lc", "lon", "lat", "id")
  
  fit <- fit_ssm(indata,
                 vmax= filt_vmax,
                 ang = filt_ang,
                 distlim = filt_distlim,
                 model = "crw",
                 time.step = reg_time_step,
                 control = ssm_control(verbose = 0))
  


  # get fitted locations
  # segments that did not converge were not consider
  data <- data.frame(grab(fit, what = "predicted", as_sf = FALSE))
  # data <- data %>% rename(trip = deploymentID) %>% arrange(date)
  data <- data %>% arrange(date)
  # data <- cbind(id = i, data)
  
  # check if points on land
  land <- NULL
  data$onland <- point_on_land(lat = data$lat, lon = data$lon, land = land)
  data$organismID <- indata$organismID[1]
  data$deploymentID <- indata$deploymentID[1]
  data$tripID <- indata$tripID[1]
  
  
  
  # export track data into individual folder at output path
  outpath <-  paste0(outdir, "/", "filteringEMbc/00L2loc")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)
  outfile <- paste0(outpath, "/", tripIDs[i], "_L2_locations.csv")
  write.csv(data, outfile, row.names = FALSE)
  
  # export convergence status
  convergence <- data.frame(trip = fit$id, converged = fit$converged)
  outfile <- paste0(outpath, "/", tripIDs[i], "_L2_convergence.csv")
  write.csv(convergence, outfile, row.names = FALSE)
  
  # plot figures
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








