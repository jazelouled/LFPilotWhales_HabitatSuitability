# FUNCTIONS


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