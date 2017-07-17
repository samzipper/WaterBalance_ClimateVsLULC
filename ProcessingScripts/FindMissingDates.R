## FindMissingDates.R
#
# Goal: Take a vector of dates (in POSIXct format) and return a vector of dates missing between the start & end
#
# Required inputs:
#  datesIn = vector of dates in POSIXct format

FindMissingDates <- function(datesIn){
  # get start & end of input vector
  start <- datesIn[1]
  end   <- datesIn[length(datesIn)]
  
  # make a vector that has every date between the two
  datesAll <- seq(start, end, by="day")
  
  # figure out difference between the two
  MissingDates <- datesAll[!(datesAll %in% datesIn)]
  
  # return MissingDates
  return(MissingDates)
}
