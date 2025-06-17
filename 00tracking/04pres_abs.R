

# 1. Combine all observed tracks. We avoid placing an absence for individual (A) in an cell (i) and day (j),
# where there is the presence of indivdiaul (B).
# 2. Combine all pseudo-absences.
# 3. Get cell (i) and day (j) for observed tracks. Removed duplicates: considers temporal and spatial autocorrelation of tracking data.
# Note individual variability is not considered in this study.
# 4. Remove pseudo-absences that overlap in time and space with observations.


aa <- read.csv(here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/03PresAbs_1monthThreshold/L2_PresAbsNeritic_1monthThreshold.csv"))
pres_aa <- aa %>% filter(occ == 1)
abs_aa <- aa %>% filter(occ == 0)

plot(abs_aa$lon, abs_aa$lat)
points(pres_aa$lon, pres_aa$lat, col="red")


 
library(dplyr)
library(raster)
library(doParallel)



#---------------------------------------------------------------
# 1. Set data repository
#---------------------------------------------------------------

ssm_data <- paste0(output_data, "/tracking/AGAZ", "/L2_locations")
sim_data <- paste0(output_data, "/tracking/AGAZ", "/L2_simulations")
sim_data <- paste0(output_data, "/tracking/AGAZ", "/L2_simulations_FScpf")


ssm_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/01cutTripsFinal/")
sim_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/02L2_simulations/")
outdir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/03PresAbs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)



ssm_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/01cutTripsFinal/")
sim_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/02L2_simulations/")
outdir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/03PresAbs_1monthThreshold")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)



ssm_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/01cutTripsFinal/")
sim_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/02L2_simulations/")
outdir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/03PresAbs_1monthThreshold")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

 
ssm_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/01cutTripsFinal/")
sim_data <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/02L2_simulations/")
outdir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/03PresAbs_2monthThreshold")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

 
sim_n <- 50
#---------------------------------------------------------------
# 2. Import presence and absences
#---------------------------------------------------------------

# 77856 to 59326 cell absences when doing it by organismID
# 77856 to 59326 cell absences when doing it by tripID

# Presence data (observed state-space models)
loc_files <- list.files(ssm_data, full.names = TRUE, pattern="cutTrack.csv")

pres <- readTrack(loc_files)

pres <- pres %>%
  dplyr::mutate(dateTime = ifelse(grepl("/", timestamp), dmy_hm(timestamp), ymd_hms(timestamp)))

pres$dateTime <- as.POSIXct(pres$dateTime, origin = "1970-01-01", tz = "UTC")
pres$date <- as.Date(pres$dateTime)


# Absence data (simulations)
loc_files <- list.files(sim_data, full.names = TRUE, pattern=".csv")
abs <- readTrack(loc_files)

abs <- abs %>%
  dplyr::mutate(dateTime = ifelse(grepl("/", timestamp), dmy_hm(timestamp), ymd_hms(timestamp)))

abs$dateTime <- as.POSIXct(abs$dateTime, origin = "1970-01-01", tz = "UTC")
abs$date <- as.Date(abs$dateTime)




# Filter data by number of simulations
abs <- dplyr::filter(abs, nsim <= sim_n)


#---------------------------------------------------------------
# 2. Generate oceanmask
#---------------------------------------------------------------
# Create a ocean mask to grid all observations
# It is based on the following parameters:
# res: resolution
# ext: extent estimated from the data
# define bounding box
xmin <- floor(min(min(pres$lon), min(abs$lon)))
xmax <- ceiling(max(max(pres$lon), max(abs$lon)))
ymin <- floor(min(min(pres$lat), min(abs$lat)))
ymax <- ceiling(max(max(pres$lat), max(abs$lat)))
ext <- extent(xmin, xmax, ymin, ymax)

#create grid
res <- 0.25
grid <- raster(ext, res = res, crs = crs("+proj=longlat +datum=WGS84"))


#---------------------------------------------------------------
# 3. Generate presences
#---------------------------------------------------------------

# Extract cell ID from raster for each observation
pres$cell <- cellFromXY(grid, cbind(pres$lon, pres$lat))

# Transform observation to presence by cell and date
cpres <- pres %>%
  dplyr::group_by(organismID, tripID, cell, date) %>%
  dplyr::summarize(occ = 1,
            n = dplyr::n())

colnames(cpres)[colnames(cpres) == 'tripID'] <- 'trip'


# Define the trips to be excluded
trips_to_exclude <- c("CARCAR_60523b_CEAMAR_003", "CARCAR_60523b_CEAMAR_004", "CARCAR_60523b_CEAMAR_006")




# Filter out the specified trips
cpres <- cpres %>% 
  filter(!trip %in% trips_to_exclude)

# mixind <- pres %>% filter(pres$deploymentID == "CARCAR_60523b_CEAMAR_20060513")
# mixind <- pres %>% filter(pres$tripID == "CARCAR_60523b_CEAMAR_001")
# plot(mixind$lon, mixind$lat)



# Get raster coordinates
xy <- xyFromCell(grid, cpres$cell)
cpres$lon <- xy[,"x"]
cpres$lat <- xy[,"y"]


#---------------------------------------------------------------
# 4. Generate absences
#---------------------------------------------------------------

# Extract cell ID from raster for each observation
abs$cell <- cellFromXY(grid, cbind(abs$lon, abs$lat))

# Transform observation to absence by cell and date

abs <- abs %>%
  mutate(organismID = sub("_[^_]+$", "", trip))



cabs <- abs %>%
  dplyr::group_by(organismID, trip, cell, date) %>%
  dplyr::summarize(occ = 0,
            n = dplyr::n())

# Get raster coordinates
xy <- xyFromCell(grid, cabs$cell)
cabs$lon <- xy[,"x"]
cabs$lat <- xy[,"y"]

points(cpres$lon, cpres$lat, col="green")
plot(cpres$lon, cpres$lat)
points(cpres$lon, cpres$lat, col="red")

points(cabs$lon, cabs$lat, col="green")
plot(cabs$lon, cabs$lat)
points(cpres$lon, cpres$lat, col="red")


#--------------------------------------------------------------------------
# FILTER OUT ABSENCES BY SPATIAL AND TEMPORAL CRITERIA
# We overlap absence with the presence bins. If there was a presence, we remove such absence
# separate non-overlapping absence and presence locations
#--------------------------------------------------------------------------
detectCores()
cores <- 7
# list unique ids
id_list <- unique(cabs$trip)

# Register number of cores to use in parallel
cl <- parallel::makeCluster(cores) # 10 cores work at around 63% CPU (no major problem with RAM)
registerDoParallel(cl)


#create empty list
cabs_list <- list()

temporal_thrs <- 31 # sequential processing for each tag

for(j in 1:length(id_list)){
  
  print(paste("Processing tag", j, "from", length(id_list)))
  
  # selected absences for a given animal id
  jcabs <- dplyr::filter(cabs, organismID == id_list[j])
  
  # for each absence, check if there is a presence in adjacent cells within the temporal period defined
  # if there is a match, remove absence. if not, keep it.
  # Note: This part of the code is computer intensive and time consuming. Using parallel works fine.
  keep <- foreach(i=1:nrow(jcabs), .packages=c("dplyr", "raster")) %dopar% {
    spt_overlap(abs_cell = jcabs$cell[i], abs_date = jcabs$date[i],
                pres_df = cpres, temporal_thrs, grid)
  }
  
  # filter out absences that match presences
  jcabs$keep <- unlist(keep)
  jcabs <- dplyr::filter(jcabs, keep == TRUE) %>%
    dplyr::select(-keep)
  
  # append
  cabs_list[[j]] <- jcabs
}

#combine absences
cabs_all <- rbindlist(cabs_list)

# Stop cluster
parallel::stopCluster(cl)




#--------------------------------------------------------------------------
# EXPORT DATA
#--------------------------------------------------------------------------

# combine presence and absence into a single data.frame and save
comb <- rbind(data.frame(cpres), data.frame(cabs_all))
comb$sp_code <- sp_code
comb <- dplyr::select(comb, organismID, trip, cell, lon, lat, date, occ)

# export to species folder
outfile <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/03PresAbs/L2_PresAbsNeritic.csv")
write.csv(comb, outfile, row.names = FALSE)

outfile <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/03PresAbs/L2_PresAbsOceanic.csv")
write.csv(comb, outfile, row.names = FALSE)

outfile <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/03PresAbs_1monthThreshold/L2_PresAbsOceanic_1monthThreshold.csv")
write.csv(comb, outfile, row.names = FALSE)

outfile <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/03PresAbs_2monthThreshold/L2_PresAbsOceanic.csv")
write.csv(comb, outfile, row.names = FALSE)

outfile <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/03PresAbs_1monthThreshold/L2_PresAbsNeritic_1monthThreshold.csv")
write.csv(comb, outfile, row.names = FALSE)



# ff <- read.csv("/Users/jazelouled-cheikhbonan/Dropbox/2023_LoggerheadWestAfrica_cc/Loggerhead_Rproject/00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/03PresAbs/L2_PresAbsNeritic.csv")
# fff <- subset(ff, ff$date == "2011-11-30")
# plot(ff$occ)
# 
# library(ggplot2)
# 
# # Assuming 'df' is your data frame
# ggplot(fff, aes(x = lon, y = lat, color = factor(occ))) +
#   geom_point() +
#   scale_color_manual(values = c("blue", "red"), labels = c("Absence", "Presence")) +
#   labs(title = "Presence and Absence Data",
#        x = "Longitude",
#        y = "Latitude",
#        color = "Occurrence") +
#   theme_minimal()






