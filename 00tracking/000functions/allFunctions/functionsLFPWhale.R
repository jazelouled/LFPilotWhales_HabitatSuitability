# FUNCTIONS

readTrack <- function(csvfiles){
  # Description
  # Reads a standardized animal track data in csv.
  # Returns a data frame with parsed time
  # It allows the combination of multiple files
  # csvfiles: string with the location of 1 or more csv files
  
  library(lubridate)
  library(data.table)
  
  ## create empty list
  dt_list <- list()  
  
  ## process and append files
  for (i in 1:length(csvfiles)){
    data <- read.csv(csvfiles[i], header=TRUE)  # read csv
    # data$date <- parse_date_time(data$date, "Ymd HMS") # parse time
    dt_list[[i]] <- data  # append to list
  }
  
  dt <- rbindlist(dt_list, fill=TRUE)  # combine data.frames
  return(dt)
}


summarizeId <- function(data){
  # data is a data.frame with all tracking data per species
  
  df <- data %>%
    dplyr::arrange(Date) %>%  # order by date
    dplyr::group_by(PTT) %>%  # select group info
    dplyr::summarize(date_deploy = first(Date),
                     lon_deploy = first(Longitude),
                     lat_deploy = first(Latitude),
                     date_last = last(Date),
                     time_interval_min = round(median(as.numeric(difftime(tail(Date, -1), head(Date, -1), units="mins")))),
                     n_loc = n()) %>%  # get first and last observations
    dplyr::mutate(duration_d = round(difftime(date_last, date_deploy, units="days")))  # calculate duration of the track
  
  return(df)
}



summarizeTrips <- function(data){
  # data is a data.frame with all tracking data per species
  
  library(geosphere)
  library(dplyr)
  
  df <- data %>%
    dplyr::arrange(Date) %>%  # order by date
    dplyr::group_by(PTT) %>%  # select group info
    dplyr::summarize(date_deploy = first(Date),
                     lon_deploy = first(Longitude),
                     lat_deploy = first(Latitude),
                     date_last = last(Date),
                     time_interval_h = median(as.numeric(difftime(tail(Date, -1), head(Date, -1), units="hours"))),
                     distance_km = sum(distGeo(p1 = cbind(Longitude, Latitude)), na.rm=TRUE)/1000,  # segment distance
                     n_loc = n()) %>%  # get first and last observations
    dplyr::mutate(duration_h = round(difftime(date_last, date_deploy, units="hours")))  # calculate duration of the track
  
  return(df)
}





filter_dup  <- function (data, step.time = 2/60, step.dist = 0.001){
  # Returns a data.frame with filtered data, and the estimates of vmax and vmax loop from SDLfilter
  # vmax: value of the maximum of velocity using in km/h
  #
  # Previous name of function was: L02L1_sdl
  
  require(SDLfilter)
  

  
  ## Keep original number of locations
  loc0 <- nrow(data)
  
  ## Convert standard format to SDLfilter
  
  # Standardize Location clasess
  data$lc <- as.character(data$Quality)
  data$lc[data$Quality == "A"] <- -1 
  data$lc[data$Quality == "B"] <- -2
  data$lc[data$Quality == "Z"] <- -9
  data$lc[data$Quality == "G"] <- 4
  data$lc <- as.numeric(data$lc)
  data$lon <- as.numeric(data$Longitude)
  data$lat <- as.numeric(data$Latitude)
  data$id <- as.numeric(data$PTT)
  
  
  # Rename columns
  names(data)[names(data)=="Date"] <- "DateTime"
  names(data)[names(data)=="lc"] <- "qi"
  
  ### Remove duplicated locations, based on both time and space criteria
  data <- dupfilter(data.frame(data), step.time=step.time, step.dist=step.dist, conditional = FALSE)
  
  ## Back transform data.frame to standar format
  
  # Rename columns
  names(data)[names(data)=="DateTime"] <- "Date"
  names(data)[names(data)=="qi"] <- "lc"
  
  # Standardize Location clasess
  data$lc[data$lc == -1] <- "A" 
  data$lc[data$lc == -2] <- "B"
  data$lc[data$lc == 4] <- "G" 
  
  ## Prepare output
  return(data)
}




landMask <- function(tm, pt){
  xy <- cbind(pt[[1]],pt[[2]])
  ov <- extract(oceanmask, xy)
  if(is.na(ov)) ov <- 1  # if point is outside study area classify as land (1)
  ov == 0  # returns a logical indicating whether the point is at sea (TRUE) or on land (FALSE)
}






mapL1 <- function (data){
  # data: locations L1 data.frame
  #
  # Required columns for L1
  # date: POSIXct
  # argosfilter:
  # onland:
  # lon:
  # lat:
  #
  # Function return a plot
  
  
  # Load libraries and dependencies
  library(ggplot2)
  library(RColorBrewer)
  #source("R/config.R")  # set your Google Drive data folder here
  #source("R/database_tools.R")
  
  # Import world map
  data(countriesHigh, package = "rworldxtra", envir = environment())
  wm <- suppressMessages(fortify(countriesHigh))
  
  ### Parse date format
  data$time <- as.integer(data$Date)
  
  ### Filter location data
  data <- filter(data, argosfilter == "not" & onland == "FALSE")
  
  ### Get metadata
  sdate <- data$Date[1]
  edate <- data$Date[nrow(data)]
  days <- round(as.numeric(difftime(edate, sdate, units="days")))
  
  ### Define extension for plot
  xl <- extendrange(data$lon, f = 0.5)
  yl <- extendrange(data$lat, f = 0.5)
  b <- c(data$time[1], data$time[nrow(data)])
  blab <- c(as.Date(data$date[1]), as.Date(data$date[nrow(data)]))
  
  ### Plot
  p <- ggplot() +
    geom_polygon(data = wm, aes_string(x = "long", y = "lat", group = "group"),
                 fill = grey(0.3)) +
    coord_quickmap(xlim = xl, ylim = yl, expand = TRUE) +
    xlab("Longitude") +
    ylab("Latitude") +
    geom_point(data = data,
               aes_string(x = "lon", y = "lat", group = NULL, colour = "time"),
               size = 2) +
    #scale_colour_gradientn(colours = terrain.colors(10)) + 
    # scale_colour_gradientn(colours=rev(brewer.pal(10,"RdYlGn")),
    #                        breaks=b,labels=format(blab)) +
    #scale_colour_brewer() +
    geom_path(data = data,
              aes_string(x = "lon", y = "lat", group = "id")) +
    labs(title = paste(data$sp_code[1], "Id:", data$id[1]),
         subtitle = paste("Start:", sdate, "End:", edate, paste0("(",days), "days)")) 
  
  return(p)
  
}






point_on_land <- function(lon, lat, land = NULL){
  
  require(sp)
  require(maptools)
  require(maps)
  
  # If no landmask provided, get continent map
  if(is.null(land)){
    land <- map("world", fill = TRUE, plot = FALSE)
    land <- map2SpatialPolygons(land, IDs=land$names,proj4string=CRS("+proj=longlat +ellps=WGS84"))
  }
  
  # Convert to spatial point
  xy <- cbind(Longitude,Latitude)
  pts <- SpatialPoints(xy, proj4string=land@proj4string)
  
  # Overlay points with continents
  ov <- over(pts, as(land,"SpatialPolygons"))  # returns NA when no overlap, and poly ID when there is an overlap
  onland <- !is.na(ov)  # returns TRUE when there is overalp, FALSE when not
  
  # Return output
  return(onland)
}


