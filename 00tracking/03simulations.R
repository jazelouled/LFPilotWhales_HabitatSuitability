#-------------------------------------------------------------------------------------
# simulations        Simulate tracks
#-------------------------------------------------------------------------------------
# This script simulates tracks to generate pseudo absences
#
# Main steps are:
# - Simulate tracks


# TODO: Add time-varying land mask (e.g. sea ice)

library(raster)

library(remotes)
# install_github("AustralianAntarcticDivision/availability")
library(availability)
library(stringr)
library(doParallel)

here()
indir <- here::here("000inputOutput/00output/00tracking/L2_locations")
outdir <- here::here("000inputOutput/00output/00tracking/simulations")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)



#------------------------------------------------------------------------------------
# summarizeTrips2        Sumarize tracking data per trip
#------------------------------------------------------------------------------------
summarizeTrips2 <- function(data){
  # data is a data.frame with all tracking data per species
  
  library(geosphere)
  library(dplyr)
  
  df <- data %>%
    dplyr::arrange(dateTime) %>%  # order by date
    dplyr::group_by(tripID, deploymentID) %>%  # select group info
    dplyr::summarize(date_deploy = first(dateTime),
                     lon_deploy = first(lon),
                     lat_deploy = first(lat),
                     date_last = last(dateTime),
                     time_interval_h = median(as.numeric(difftime(tail(dateTime, -1), head(dateTime, -1), units="hours"))),
                     distance_km = sum(distGeo(p1 = cbind(lon, lat)), na.rm=TRUE)/1000,  # segment distance
                     n_loc = n()) %>%  # get first and last observations
    dplyr::mutate(duration_h = round(difftime(date_last, date_deploy, units="hours")))  # calculate duration of the track
  
  return(df)
}
#------------------------------------------------------------------------------------


#-----------------------------------------------------------------
# mapSimTracks       Plot simulations and observed track
#-----------------------------------------------------------------
mapSimTracks  <- function (simData, obsData, title = NULL){
  
  # Load libraries and dependencies
  library(ggplot2)
  
  # Import world map
  data(countriesHigh, package = "rworldxtra", envir = environment())
  wm <- suppressMessages(fortify(countriesHigh))
  
  ### Define extension for plot
  xl <- extendrange(c(simData$lon, obsData$lon), f = 0)
  yl <- extendrange(c(simData$lat, obsData$lat), f = 0)
  
  ### Plot
  p <- ggplot() +
    geom_polygon(data = wm, aes_string(x = "long", y = "lat", group = "group"),
                 fill = grey(0.3)) +
    coord_quickmap(xlim = xl, ylim = yl, expand = TRUE) +
    xlab("Longitude") +
    ylab("Latitude") +
    geom_path(data = simData,
              aes_string(x = "lon", y = "lat", group = "simid"),
              size=0.5, alpha=0.5, color="grey70") +
    geom_path(data = obsData,
              aes_string(x = "lon", y = "lat", group = "id"),
              size=0.5, alpha=1, color="red1") +
    geom_point(data = data.table::first(obsData),
               aes_string(x = "lon", y = "lat"),
               shape = 21, colour = "red4", fill = "white", size = 2, stroke = 2) +
    #size=3, alpha=1, color="black") + #(shape = 21, colour = "red", fill = "white", size = 4, stroke = 2)
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle(title)
  
  return(p)
}
#-----------------------------------------------------------------














#---------------------------------------------------------------
# 4. Summarize processed data
#---------------------------------------------------------------

# import all location files
loc_files <- list.files(indir, full.names = TRUE, pattern = "L2_locations.csv")
df <- readTrack(loc_files)
df$date <- as.POSIXct(df$date, format = "%d/%m/%Y %H:%M")


# import convergence files
loc_files <- list.files(outdir, full.names = TRUE, pattern = "L2_convergence.csv")
data_proc <- lapply(loc_files, read.csv) %>% rbindlist

# summarize data per trip
tripstats <- summarizeTrips2(df)

# combine track data summary and convergence status
comb <- merge(tripstats, data_proc, by=c("trip"))

# export table
out_file <- paste0(outdir, "/", modelSet, "_summary_ssm.csv")
write.csv(tripstats, out_file, row.names = FALSE)


#---------------------------------------------------------------
# Prepare cluster
#---------------------------------------------------------------
cores <- 4  # detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

#---------------------------------------------------------------
# 1. Set data repository
#---------------------------------------------------------------
indir <- paste0(here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/"), modelSet, "/01cutTripsFinal")
outdir <- paste0(here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/"), modelSet, "/02L2_simulations")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


if (modelSet == "neritic") {
  oceanmask <- raster(here::here("00output/caret/00enviro/00terrain/oceanmaskNeriticFiltered.tif"))
} else {
  oceanmask <- raster(here::here("00output/caret/00enviro/00terrain/oceanmaskOceanicFiltered.tif"))
}

require(raster)
r <- oceanmask

myMask <- function() {
  function(tm, pt) {
    xpt <- pt[1]
    ypt <- pt[2]
    v <- raster::extract(oceanmask, matrix(c(xpt, ypt), nrow = 1, ncol = 2))
    if (v > 0 | is.na(v) == T) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}


# xs = d[,c("lon", "lat")], ts = d[,c("date")]
# 
# xs[1]
# 
# v <- raster::extract(oceanmask, xs[14], nrow = 1, ncol = 2)



#---------------------------------------------------------------
# 2. Create oceanmask
#---------------------------------------------------------------

# Import oceanmask
# oceanmask <- raster(here::here(output_data, "00enviro/00terrain/oceanmask.tif"))
oceanmask <- oceanmask+0  # this makes the raster to be in memory and make simulations faster

#---------------------------------------------------------------
# 3. Select data
#---------------------------------------------------------------

# import summary data
db <- tripstats

# select by number of locations and minimum duration
tags <- unique(db$deploymentID)
if (!is.null(sim_exclude)) tags <- tags[which(!tags %in% sim_exclude)]

#---------------------------------------------------------------
# 4. Simulate track
#---------------------------------------------------------------

# foreach(i=1:length(tags), .packages=c("dplyr", "ggplot2", "availability", "data.table", "raster", "stringr")) %dopar% {
  for (i in 1:length(tags)){
  
  print(paste("Processing tag", tags[i]))
  
  # import data
  loc_file <- here::here(indir,paste0(tags[i],"_", "L2_locations.csv"))
 
  data <- readTrack(loc_file)
  
  # # select trips that converged in the SSM
  trips <- unique(data$id)
  # sel <- which(trips %in% db$trip[db$converged==TRUE])
  # trips <- trips[sel]
  # 
  # # if simulations are conducted for the whole track, overwrite trip data
  # if(sim_by_trip == TRUE){
  #   data$trip <- data$id
  #   trips <- data$id[1]
  # }

  trip_list <- list()
  
  
  sim_n <- 50
  
  
  
  for (j in 1:length(trips)) {
    print(trips[j])
    # select data for selected segment
    d <- dplyr::filter(data, id == id[j])

    # Fit a vector-autoregressive movement model to this filtered track.
    # This model assumes that the x- and y-speeds at time t are a linear function
    # of the speeds at time t-1, plus random noise.
    arfit <- surrogateARModel(d[,c("lon", "lat")])
    
    # Now we can use that fitted model to generate new tracks. 
    # Simulated track are fixed to the same start, but no end points
    # Land mask is applied:
    
    ## generate simulations
    data_list <- list()


    for (s in 1:sim_n){
      print(s)
      # Now we can use that fitted model to generate new tracks. 
      # Simulated track are fixed to the same start always. End points fixed for central-place foragers
      # Land mask is applied:
      # I changed "landmask" for gshhsMask(), which is the landmask provided by the package, as I did not find "landmask"
      # anywhere. It works. 
      # 
      # simu <- surrogateAR(arfit, xs = d[,c("lon", "lat")], ts = d[,c("date")], point.check = gshhsMask(),
      #                     fixed = rep(c(TRUE, FALSE, sim_fix_last), c(1, nrow(d) - 2, 1)),
      #                     partial=FALSE)
      
      # The `point.check` argument accepts a function of the form `function(tm, pt)` that returns `TRUE`
      # if the point should be accepted and `FALSE` if not. Note that this function can accept a time coordinate 
      # as well as a location, and so the mask can be made time-varying if required (e.g. dynamically masking out areas covered by sea ice).
      
      
      simu <- availability::surrogateAR(arfit, xs = d[,c("lon", "lat")], ts = d[,c("date")], point.check = myMask(),
                          fixed = rep(c(TRUE, FALSE, sim_fix_last), c(1, nrow(d) - 2, 1)),
                          partial=FALSE)
      
      if(is.null(simu) | is.na(simu$xs[1])) break
      
      ## convert to our generic format for tracks
      sim_code <- str_pad(s, 3, pad = "0")
      df <- data.frame(id = i, trip = trips[j], nsim = s, simid = paste(trips[j], sim_code, sep="_"), date = simu$ts, lon = simu$xs[,1], lat = simu$xs[,2])
      
      ## check overlap with observed track
      #ov <- SToverlap(real_track = d, sim_track = df, exclude_dur = exclude_dur, sp_thrs = sp_thrs, t_thrs = t_thrs)
      #if(ov == T) next
      
      ## append data.frame into list
      data_list[[s]] <- df
     }

    ## combine simulations into a single data.frame
    simdf <- rbindlist(data_list)
    
    ## append to segment list
    trip_list[[j]] <- simdf
  }
  
  ## combine simulations into a single data.frame
  simdf <- rbindlist(trip_list)
  
  # export track data into individual folder at output path
  out_file <- paste0(outdir, "/", tags[i], "_sim_L2_locations.csv")
  write.csv(simdf, out_file, row.names = FALSE)
  
  ## Plot simulations
  p <- mapSimTracks(simData = simdf, obsData = data, title = trips[j])
  out_file <- paste0(outdir, "/", tags[i], "_sim_L2_locations.png")
  ggsave(out_file, p, width=30, height=15, units = "cm")
}
  

#---------------------------------------------------------------
# 5. Summarize processed data
#---------------------------------------------------------------

# identify all location files
loc_files <- list.files(outdir, full.names = TRUE, pattern = "sim_L2_locations.csv")

# sumarize number of simulations per trip
# some trips may have problems for simulation. So, next steps will filter trip
# data that has been simulated successfuly.
sim_stats <- rbindlist(foreach(i=length(loc_files), .packages=c("dplyr", "data.table")) %dopar% {
  df <- readTrack(loc_files[i])
  simdf <- summarizeSim(df)
  return(simdf)
})

# combine with trip data from the SSM (input data)
comb <- merge(db, sim_stats, by="trip", all.x=TRUE)


# export table
out_file <- here::here(outdir, "summary_sim.csv")
write.csv(comb, out_file, row.names = FALSE)


  
#---------------------------------------------------------------
# Stop cluster
#---------------------------------------------------------------
stopCluster(cl)




print("Simulations ready")        
