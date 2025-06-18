#-----------------------------------------------------------------------------------------
# 01_preproc_GAZ.R        Pre-process tracking data
#-----------------------------------------------------------------------------------------
# This script pre-processes tracking data. The main goal is to standardize data
# formats
#
# About Antarctic fur seal data:
# Dataset has been prepared by Lluis Cardona (UB). All tags have been previously standardized
# into a custom common format. All location data correspond to Argos
# install.packages("openxlsx")

library(openxlsx)


library(doParallel)
library(here)


#---------------------------------------------------------------
# Prepare cluster
#---------------------------------------------------------------
cores <- 6
cl <- makeCluster(cores)
registerDoParallel(cl)

#---------------------------------------------------------------
# 1. Set data repository
#---------------------------------------------------------------

# input file with batch data
file <- here::here("000inputOutput/00input/datos_calderones_todos_MOD.csv") # Tag_Gm_delpoint.csv only has the start of the trip until there is a big gap (>48)
df <- read.csv(file)

# output directory
outdir <- here::here("000inputOutput/00output/00tracking/L0_locations")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


#---------------------------------------------------------------
# 2. Process metadata
#---------------------------------------------------------------


df$date <- as.POSIXct(df$date, format = "%Y-%m-%d %H:%M:%S")
df <- dplyr::select(df,  DeployID = deployid, PTT = ptt, Date = date, Longitude = longitude, Latitude = latitude, Quality = quality)



#---------------------------------------------------------------
# 3. Process batch data
#---------------------------------------------------------------

# process each tag 
foreach(i=unique(df$PTT), .packages=c("dplyr", "ggplot2", "stringr")) %dopar% {
  
  # subset data and standardize
  id <- i
  dataL0 <- dplyr::filter(df, PTT == i)
  
  # define "trips"
  # at this step we use animal ID. See filter_locs.R for further steps
  # dataL0$trip <- dataL0$DeployID
  dataL0 <- dplyr::select(dataL0, DeployID, everything())
  
  # store into individual folder at output path
  out_file <- paste0(outdir, "/", id, "_L0_locations.csv")
  write.csv(dataL0, out_file, row.names = FALSE)
  
  # plot track
  # Load libraries
  library(ggplot2)
  library(rworldxtra)
  
  # Optional species code (remove if not needed)
  dataL0$sp_code <- "GMEL"
  
  # Load and fortify world map
  data(countriesHigh, package = "rworldxtra")
  wm <- suppressMessages(fortify(countriesHigh))
  
  # Extract metadata
  sdate <- min(dataL0$Date, na.rm = TRUE)
  edate <- max(dataL0$Date, na.rm = TRUE)
  days <- round(as.numeric(difftime(edate, sdate, units = "days")))
  
  # Define map extent
  xl <- extendrange(dataL0$Longitude, f = 0.2)
  yl <- extendrange(dataL0$Latitude, f = 0.2)
  
  # Create plot
  p <- ggplot() +
    geom_polygon(data = wm, aes_string(x = "long", y = "lat", group = "group"),
                 fill = grey(0.3)) +
    coord_quickmap(xlim = xl, ylim = yl, expand = TRUE) +
    xlab("Longitude") +
    ylab("Latitude") +
    geom_path(data = dataL0,
              aes(x = Longitude, y = Latitude),
              size = 0.3,
              color = "grey40") +
    geom_point(data = dataL0,
               aes(x = Longitude, y = Latitude, colour = Quality),
               size = 0.9,
               alpha=0.4) +
    labs(
      title = paste(df$sp_code[1], "Id:", df$DeployID[1]),
      subtitle = paste("Start:", sdate, "End:", edate, "(", days, "days)")
    ) +
    theme_bw()+
    scale_colour_viridis_d(option = "C", name = "Location Class", direction = -1)
  
  # Show plot
  print(p)
  
  out_file <- paste0(outdir, "/", id, "_L0_locations.png")
  ggsave(out_file, p, width=30, height=15, units = "cm")
}

stopCluster(cl)  # Stop cluster


#---------------------------------------------------------------
# 5. Summarize processed data
#---------------------------------------------------------------

# summarize data per animal id
idstats <- summarizeId(df)

# export table
out_file <- paste0(outdir, "/", sp_code, "_summary_id.csv")
write.csv(idstats, out_file, row.names = FALSE)

print("Pre-processing ready")





