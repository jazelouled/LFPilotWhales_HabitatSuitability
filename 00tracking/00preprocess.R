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
cores <- 2
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
df <- dplyr::select(df, IF = DeployID, Tag_ID = Ptt, Tag_type = Instr, date = Date, lon = Longitude, lat = Latitude, lc = Quality)
df <- dplyr::select(df, DeployID = id, Ptt = Tag_ID, Instr = Tag_type, Date = date, Longitude = lon, Latitude = lat, Quality = lc)


# summarize data per id
df <- dplyr::select(df, deployid = DeployID, ptt = PTT, date = Date, lon = Longitude, lat = Latitude, lc = Quality)

df <- dplyr::select(df,  DeployID = deployid, PTT = ptt, Date = date, Longitude = longitude, Latitude = latitude, Quality = lc)


#---------------------------------------------------------------
# 3. Process batch data
#---------------------------------------------------------------

# process each tag 
foreach(i=db$id, .packages=c("dplyr", "ggplot2", "stringr")) %dopar% {
  
  # subset data and standardize
  id <- i
  dataL0 <- dplyr::filter(df, id == i)
  
  # define "trips"
  # at this step we use animal ID. See filter_locs.R for further steps
  dataL0$trip <- dataL0$id
  dataL0 <- dplyr::select(dataL0, id, trip, everything())
  
  # store into individual folder at output path
  out_file <- paste0(outdir, "/", id, "_L0_locations.csv")
  write.csv(dataL0, out_file, row.names = FALSE)
  
  # plot track
  p <- map_argos(dataL0)
  out_file <- paste0(outdir, "/", id, "_L0_locations.png")
  ggsave(out_file, p, width=30, height=15, units = "cm")
}

stopCluster(cl)  # Stop cluster


#---------------------------------------------------------------
# 5. Summarize processed data
#---------------------------------------------------------------

# import all location files
loc_files <- list.files(outdir, full.names = TRUE, pattern = "L0_locations.csv")
df <- readTrack(loc_files)
df$date <- as.POSIXct(df$date, format = "%Y-%m-%d %H:%M")

# summarize data per animal id
idstats <- summarizeId(df)

# export table
out_file <- paste0(outdir, "/", sp_code, "_summary_id.csv")
write.csv(idstats, out_file, row.names = FALSE)


print("Pre-processing ready")





