## PheasantBranch_1_Discharge_DownloadDaily.R
# This script downloads daily USGS streamflow data for the Pheasant Branch gauging station.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(dataRetrieval)
require(lubridate)

# get USGS data
dv <- readNWISdv(siteNumber="05427948",                          # site code (can be a vector of multiple sites)
                 parameterCd="00060",                            # parameter code: "00060" is cubic ft/sec
                 startDate="1900-01-01",endDate="2016-12-31",    # start & end dates (YYYY-MM-DD format)
                 statCd = "00003")                               # statistic code: "00003" is daily mean (default)

# do some conversions and renaming
colnames(dv) <- c("agency_cd", "site_no", "Date", "discharge.cfs", "QA.code")

# trim to only interesting columns
dv <- dv[,c("Date", "discharge.cfs", "QA.code")]

# save as output
write.csv(dv, paste0(git.dir, "Data/PheasantBranch/PheasantBranch_Discharge_Daily.csv"),
          row.names=F)
