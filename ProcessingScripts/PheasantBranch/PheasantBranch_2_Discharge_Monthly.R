## PheasantBranch_2_Discharge_Monthly.R
# This reads in Pheasant Branch USGS daily data (downloaded with script PheasantBranch_1_Discharge_DownloadDaily.R) and
# summarizes it to monthly.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(lubridate)
require(dplyr)
require(rgdal)
require(raster)
source(paste0(git.dir, "ProcessingScripts/FindMissingDates.R"))

setwd(paste0(git.dir, "Data/PheasantBranch/"))

# read in daily data
df <- read.csv("PheasantBranch_Discharge_Daily.csv", stringsAsFactors=F)

# reclassify date
df$Date <- ymd(df$Date)

# data starts on 1930-09-10, trim to start at 1930-10-01
df <- subset(df, Date >= ymd("1930-10-01"))

# find missing dates
missing <- FindMissingDates(df$Date)

# make a data frame with all dates
df.all <- data.frame(Date = seq(min(df$Date), max(df$Date), by="day"))
df <- merge(df.all, df, by="Date", all.x=T)

# get year and month columns
df$year <- year(df$Date)
df$month <- month(df$Date)

# summarize to monthly
df.mo <- dplyr::summarize(group_by(df, year, month),
                   discharge.cfs.mean = mean(discharge.cfs, na.rm=T),
                   missing = sum(is.na(discharge.cfs)))

# for each month, figure out how many days are in that month
df.mo$ndaypm <- days_in_month(ymd(paste0(df.mo$year, "-", df.mo$month, "-01")))

# get rid of months with >10% of days missing
df.mo <- subset(df.mo, missing/ndaypm<0.10)  # this is NO MONTHS! 

# convert discharge in cfs to total discharge in cubic meters
df.mo$discharge.m3 <- df.mo$discharge.cfs.mean*86400*df.mo$ndaypm*(0.3048^3)

## convert to a depth
# read in shapefile of upstream watershed
area.upstream <- (17.08)*(5280*5280)*(0.3048*0.3048)  # [m] contributing area, from USGS: https://waterdata.usgs.gov/nwis/inventory/?site_no=05427948&agency_cd=USGS

# convert discharge to mm
df.mo$discharge.mm <- 1000*df.mo$discharge.m3/area.upstream

# save parts of interest
write.csv(df.mo[,c("year", "month", "discharge.mm")], "PheasantBranch_Discharge_Monthly.csv", row.names=F)
