library(lubridate)
library(data.table)
library(raster)
library(terra)
library(here)
library(sf)
library(doParallel)
`%notlike%` <- Negate(`%like%`)
library(doParallel)

# Directories
# indir <- here::here("00output/caret/01tracking/PresAbs/")
# outdir <- here::here("00output/caret/01tracking/04extract")
# if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)  # create output directory if does not exist
# 
# indir <- here::here("00output/caret/01tracking/PresAbs_NOinternesting/")
# outdir <- here::here("00output/caret/01tracking/extract_NOinternesting")
# if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)  # create output directory if does not exist

indir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/03PresAbs")
outdir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/neritic/04extract/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)  # create output directory if does not exist


indir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/03PresAbs")
outdir <- here::here("00output/caret/01tracking/01filteringForagingUpwelling/00L2_loc_cut/cutTracks/oceanic/04extract/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)  # create output directory if does not exist



# Get stack paths
# stacks <- dir(here::here("00output/00enviro/03EnviroLayersFit_MeanSSP/"), recursive = T, full.names = T)

stacks <- dir("/Volumes/phalaropus/2023_LoggerheadWestAfrica_cc/ForcingData/04fitLayers/", recursive = T, full.names = T)


# Set buffer  
env_buffer <- 25000
# Read presence absence data
file_presabs <- read.csv(paste0(indir, "/L2_PresAbsNeritic.csv"))
file_presabs <- read.csv(paste0(indir, "/L2_PresAbsOceanic.csv"))
file_presabs$date <- as.Date(file_presabs$date) 
file_presabs$datestring <-  gsub("-", "", file_presabs$date)

# datestring <- (file_presabs$datestring)[1]
# env_stack <-  stacks[stacks %like% ".grd" & stacks %like% datestring]
# env_stack_ <- stack(env_stack)
# mask_env <- env_stack_$NPP/env_stack_$NPP
# plot(mask_env)

# polygonDeletePoints <- rasterToPolygons(mask_env, dissolve = T)
# polygonDeletePoints_ <- st_as_sf(polygonDeletePoints)
# 
# coordinatesAndAttributes <- sf::st_as_sf(file_presabs, coords = c("lon","lat"))
# xysf <- st_as_sf(as.data.frame(file_presabs), coords = c("lon","lat"), crs = 4326)
# xy_intersect <- st_intersection(polygonDeletePoints_, xysf)
# plot(xy_intersect$geometry, add=T)
# 
# 
# xy_intersect_ <- tidyr::extract(xy_intersect, geometry, into = c('Lon', 'Lat'), '\\((.*),(.*)\\)', conv = T)
# 
# 



extract_WestAfrica <- list() 

extract_SSP <- list()


# Prepare cluster
cl <- parallel::makeCluster(4)
registerDoParallel(cl)



dynamic_data <- dir("/Volumes/phalaropus/2023_LoggerheadWestAfrica_cc/ForcingData/03extractedNetCDFs/", recursive = T, full.names = T)
dynamic_data_ <- dynamic_data[dynamic_data %notlike% "regridded" & dynamic_data %notlike% "cropped"]
datestring <- sub(".*([0-9]{8})\\.nc$", "\\1", dynamic_data_)
datestring_ <- unique(datestring)
dates <- as.Date(datestring_, format = "%Y%m%d")

# outdir <- "/Volumes/phalaropus/2023_LoggerheadWestAfrica_cc/ForcingData/04stacksCMIP6"
cmip6_files_histor <- dir("/Volumes/phalaropus/2023_LoggerheadWestAfrica_cc/ForcingData/03extractedNetCDFs/historical/", recursive = T)
cmip6_files_histor_ <- cmip6_files_histor[cmip6_files_histor %notlike% "regridded" & cmip6_files_histor %notlike% "cropped"]
datestring_histor <- unique(sub(".*([0-9]{8})\\.nc$", "\\1", cmip6_files_histor_))

cmip6_files_ssp <- dir("/Volumes/phalaropus/2023_LoggerheadWestAfrica_cc/ForcingData/03extractedNetCDFs/ssp126/", recursive = T)
cmip6_files_ssp_ <- cmip6_files_ssp[cmip6_files_ssp %notlike% "regridded" & cmip6_files_ssp %notlike% "cropped"]
datestring_ssp <- unique(sub(".*([0-9]{8})\\.nc$", "\\1", cmip6_files_ssp_))


min(file_presabs$date)
max(file_presabs$date)




ssp <- c("historical", "ssp126", "ssp585")

cl <- parallel::makeCluster(4)
registerDoParallel(cl)



for (j in 1:length(ssp)){
          print(j)
          ssp_ <- ssp[j]
          
          if(ssp_ == "historical"){
            datestring <- datestring_histor
          }   else if (ssp_ %like% "ssp") {
            datestring <- datestring_ssp
          }
          
          dates <- as.Date(datestring, format = "%Y%m%d")
          
aa <- foreach(h=1:length(file_presabs$cell), .packages=c("raster", "data.table")) %dopar% {
# for (h in 65800:length(file_presabs$cell)){
    print(h)
  # if(h == 13453) next
    datestring <- (file_presabs$datestring)[h]
    env_stack <-  stacks[stacks %like% ".grd" & stacks %like% paste0("/", substr(datestring, 1, 6))]
    env_stack_ <- stack(env_stack)
    extracted_env <- raster::extract(env_stack_, cbind(file_presabs$lon[h], file_presabs$lat[h]), buffer = env_buffer, fun=mean, na.rm=TRUE)  # extract data
    extracted_env_ <- as.data.frame(extracted_env)
    extracted_env_xy <- cbind(file_presabs$lon[h], file_presabs$lat[h], file_presabs$trip[h], file_presabs$cell[h], file_presabs$date[h], file_presabs$occ[h], extracted_env_)
    colnames(extracted_env_xy)[1] <- "Longitude"
    colnames(extracted_env_xy)[2] <- "Latitude"
    colnames(extracted_env_xy)[3] <- "TripID"
    colnames(extracted_env_xy)[4] <- "Cell"
    colnames(extracted_env_xy)[5] <- "Date"
    colnames(extracted_env_xy)[6] <- "Occurrence"
    extract_WestAfrica[[h]] <- extracted_env_xy

}

extract_SSP[[j]] <- aa


}
  
extract_WestAfrica_ <- do.call(rbind, aa)


# num_cols <- sapply(extract_WestAfrica, function(df) ncol(df))
# if (length(unique(num_cols)) > 1) {
#   cat("There are dataframes with different numbers of columns:\n")
#   print(names(num_cols[which(table(num_cols) > 1)]))
# } else {
#   cat("All dataframes have the same number of columns.\n")
# }
# 
# 
# # Check if there are dataframes with different numbers of columns
# num_cols <- sapply(extract_WestAfrica, function(df) ncol(df))
# if (length(unique(num_cols)) > 1) {
#   cat("There are dataframes with different numbers of columns:\n")
#   different_cols <- names(num_cols[which(table(num_cols) > 1)])
#   print(different_cols)
#   
#   # Find dataframes with more than 13 columns
#   more_than_13 <- names(num_cols[num_cols > 13])
#   cat("\nDataframes with more than 13 columns:\n")
#   print(more_than_13)
# } else {
#   cat("All dataframes have the same number of columns.\n")
# }
# 
# 
# 
# positions_more_than_13 <- which(sapply(extract_WestAfrica, function(df) ncol(df) > 13))
# 
# if (length(positions_more_than_13) > 0) {
#   cat("Dataframes with more than 13 columns are at positions:\n")
#   print(positions_more_than_13)
# } else {
#   cat("No dataframes with more than 13 columns found.\n")
# }
# 



product_folder <- outdir
if (!dir.exists(product_folder)) dir.create(product_folder, recursive = TRUE)  # create output directory if does not exist
outfile <- "extract_WestAfrica_neriticAWIhist126.csv"
write.csv(extract_WestAfrica_, paste0(product_folder, "/", outfile))
  

product_folder <- outdir
if (!dir.exists(product_folder)) dir.create(product_folder, recursive = TRUE)  # create output directory if does not exist
outfile <- "extract_WestAfrica_oceanicAWIhist126.csv"
write.csv(extract_WestAfrica_, paste0(product_folder, "/", outfile))











