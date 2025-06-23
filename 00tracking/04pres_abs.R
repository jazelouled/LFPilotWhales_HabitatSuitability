

# 1. Combine all observed tracks. We avoid placing an absence for individual (A) in an cell (i) and day (j),
# where there is the presence of indivdiaul (B).
# 2. Combine all pseudo-absences.
# 3. Get cell (i) and day (j) for observed tracks. Removed duplicates: considers temporal and spatial autocorrelation of tracking data.
# Note individual variability is not considered in this study.
# 4. Remove pseudo-absences that overlap in time and space with observations.
 
library(dplyr)
library(raster)
library(doParallel)



#---------------------------------------------------------------
# 1. Set data repository
#---------------------------------------------------------------
ssm_data <- here::here("000inputOutput/00output/00tracking/L2_locations")
sim_data <- here::here("000inputOutput/00output/00tracking/simulations")
outdir <- here::here("000inputOutput/00output/00tracking/03PresAbs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#---------------------------------------------------------------
# 2. Import presence and absences
#---------------------------------------------------------------

# 77856 to 59326 cell absences when doing it by organismID
# 77856 to 59326 cell absences when doing it by tripID

# Presence data (observed state-space models)
loc_files <- list.files(ssm_data, full.names = TRUE, pattern="L2_locations.csv")
pres <- readTrack(loc_files)
pres$dateTime <- as.POSIXct(pres$date, format = "%Y-%m-%d %H:%M:%S")


# Absence data (simulations)
loc_files <- list.files(sim_data, full.names = TRUE, pattern="L2_locations.csv")
abs <- readTrack(loc_files)
abs$dateTime <- as.POSIXct(abs$date, format = "%Y-%m-%d %H:%M:%S")


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
res <- 0.083
grid <- raster(ext, res = res, crs = crs("+proj=longlat +datum=WGS84"))


#---------------------------------------------------------------
# 3. Generate presences
#---------------------------------------------------------------

# Extract cell ID from raster for each observation
pres$cell <- cellFromXY(grid, cbind(pres$lon, pres$lat))

# Transform observation to presence by cell and date
cpres <- pres %>%
  dplyr::group_by(id, cell, date) %>%
  dplyr::summarize(occ = 1,
            n = dplyr::n())

colnames(cpres)[colnames(cpres) == 'id'] <- 'trip'

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

temporal_thrs <- 2 # sequential processing for each tag
 


cpres$date <- as.POSIXct(cpres$date, format = "%Y-%m-%d %H:%M:%S")
cabs$date <- as.POSIXct(cabs$date, format = "%Y-%m-%d %H:%M:%S")



for(j in 1:length(id_list)){
  
  print(paste("Processing tag", j, "from", length(id_list)))
  
  # selected absences for a given animal id
  jcabs <- dplyr::filter(cabs, trip == id_list[j])
  
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
cabs_all <- cabs_all %>% dplyr::select(-organismID)


# Stop cluster
parallel::stopCluster(cl)




#--------------------------------------------------------------------------
# EXPORT DATA
#--------------------------------------------------------------------------

# combine presence and absence into a single data.frame and save
comb <- rbind(data.frame(cpres), data.frame(cabs_all))
comb <- dplyr::select(comb, organismID, trip, cell, lon, lat, date, occ)

# export to species folder
outfile <- here::here("000inputOutput/00output/00tracking/03PresAbs/presAbsCells.csv")
write.csv(comb, outfile, row.names = FALSE)


