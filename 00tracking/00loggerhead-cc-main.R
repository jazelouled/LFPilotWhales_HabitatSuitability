#----------------------------------------------------------------------------------
# loggerhead-cc-main.R            Main script for processing satellite tracking data
#----------------------------------------------------------------------------------
`%notlike%` <- Negate(`%like%`)

cores <- 5
input_data <- here::here("00input/L1_loc_CV")
output_data <- here::here("00output")


cl <- makeCluster(cores)
registerDoParallel(cl)


aaa
#---------------------------------------------------------------
# 1. Set parameters
#---------------------------------------------------------------
sp_code <- "LFPW"  # species code
tag_type <- "PTT"

# Trip definition
trip_type <- "time"  # haul: trim track by haul-out locations; time: trim track by time gaps
trip_time_gap <- 7 * 24  # (used if trip_type == time) Tracks with data gaps in excess of [seg_time_gap] hours were broken up for separate modeling
# trip_time_gap <- 200 * 24  # (used if trip_type == time) Tracks with data gaps in excess of [seg_time_gap] hours were broken up for separate modeling

# Track selection
sel_min_loc <- 10  # minimum number of locations
sel_min_dur <- 12 # minimum duration of track, in hours
sel_exclude <- NULL # custom selection of tags based on exploration of data
sel_min_dist <- 15 # minimum distance of tracks, in km

# Track filtering
filt_step_time <- 2/60  # time difference to consider duplicated positions, in hours
filt_step_dist <- 0/1000  # spatial distance to consider duplicated poisitions, in km
filt_land <- FALSE  # remove locations on land
filt_vmax <- 3  # value of the maximum of velocity using in sdafilter, in m/s
filt_ang <- c(15, 25) # value of the angle using in sdafilter, no spikes are removed if ang=-1
filt_distlim <- c(2500, 5000) # value of the limite distance using in sdafilter, no spikes are removed if ang=-1

# Track regularization
 reg_time_step <- 6  # time step to interpolate positions, in hours

# create ocean mask
mcp_expand <- 5  # expand the minimum convex polygon, in degrees.

# Simulations for CPF FS
# sim_n <- 50  # number of simulations
# sim_fix_last <- TRUE  # fix last track location
# sim_exclude <- NULL # remove individuals from simulations
# sim_by_trip <- TRUE  # generate simulation by trip rather than full track

# Simulations for non-CPF FS, MW, FW
sim_n <- 50  # number of simulations
sim_fix_last <- FALSE  # fix last track location
sim_exclude <- NULL # remove individuals from simulations
sim_by_trip <- FALSE  # generate simulation by trip rather than full track


# Extract environment
env_buffer <- 15000  # radius of buffer to average environmental data around each location, in meters. (15000)
# env_buffer <- 100000  # radius of buffer to average environmental data around each location, in meters. (100000)
# all_vars <- c("BAT", "SLP", "SDIST", "SST", "SSTg", "SAL", "SALg", "SSH", "EKE", "CHL", "SIC", "SIT", "MLD", "EDGE", "D2COL")
# env_max_date <- as.Date("2019-09-16")
all_vars <- c("BAT", "SLP", "SDIST", "SST", "CHL")


#---------------------------------------------------------------
# 2. Set data paths and import libraries
#---------------------------------------------------------------

# Load dependencies
source("setup.R")
source("scr/fun_track_reading.R")  # read multiple tracking data formats
source("scr/fun_track_plot.R")  # plot tracking data
source("scr/fun_track_proc.R")  # miscellanea of processing functions
source("scr/fun_enviro.R")  # function to process environmental data
source("scr/fun_survey.R")


#---------------------------------------------------------------
# 3. Processing workflow
#---------------------------------------------------------------

# Set number of cores for parallel processing
cores <- 10

#detectCores()-2

# Step 1. Pre-process data and standardize data
# Transforms different data sources into a common format
# Output data can be found: main_dir/output/species/sp_code/L0_locations
source("analysis/tracking/scr/preproc_GAZ.R")



# Step 2. Filter location data
# Filtering is based on selected parameters from above
source("analysis/tracking/scr/filter_locs.R")

# Step 3. Regularize location data (PROBLEMS WITH SERVER - USE PC)
# Uses correlated random walk state-space model from Jonsen et al. 2019 doi:10.1002/ecy.2566
source("analysis/tracking/scr/regularize_ssm.R")

# Step 4. Generate oceanmask
# Combines MCP from all L2 locations with bathymetry and ice extent
source("analysis/tracking/scr/oceanmask.R")

# Step 5. Generate pseudo-absences using simulations for habitat model
# We select individuals that remain within the study area.
source("analysis/tracking/scr/simulations.R")

# Step 6. Generate Presence-Absence data eliminating overlap between observations
cores <- 3
# res <- 0.1  # size of spatial bin, in decimal degrees
res <- 0.25 # size of spatial bin, in decimal degrees
temporal_thrs <- 2  # length of temporal bin, in days
sim_n <- 50  # select number of simulations to use
source("analysis/tracking/scr/pres_abs.R")
 
# Step 7. Extract environmental data
# stack_repo <- paste0(output_data, "/stack_monthly")
stack_repo <- "D:/ISIMIP_netCDF/output/montly_stacks_antarctica_ISIMIP_"

source("analysis/tracking/scr/extract_stack_presabs.R")
