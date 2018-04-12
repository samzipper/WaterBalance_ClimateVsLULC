## PheasantBranch_2_DailyToMonthly.R
# This reads in Pheasant Branch USGS daily data (downloaded with script PheasantBranch_1_Discharge_DownloadDaily.R) and
# summarizes it to monthly.
#
# Also summarizes monthly GHCN data and corrects SW radiation based on Arlington weather data.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(lubridate)
require(dplyr)
require(rgdal)
require(raster)
require(EcoHydRology)
require(ggplot2)
require(hydroGOF)
source(paste0(git.dir, "ProcessingScripts/FindMissingDates.R"))
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))
source(paste0(git.dir, "ProcessingScripts/PenmanMonteith_Daily.R"))

setwd(paste0(git.dir, "Data/PheasantBranch/"))

## First: process streamflow data
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
area.upstream <- (17.08)*(5280*5280)*(0.3048*0.3048)  # [m] contributing area, from USGS: https://waterdata.usgs.gov/nwis/inventory/?site_no=05427948&agency_cd=USGS

# convert discharge to mm
df.mo$discharge.mm <- 1000*df.mo$discharge.m3/area.upstream

# save parts of interest
write.csv(df.mo[,c("year", "month", "discharge.mm")], "PheasantBranch_Discharge_Monthly.csv", row.names=F)

## repeat for baseflow, which was separated using the WHAT online tool
df.Q.d <- read.csv("PheasantBranch_BaseflowSeparation_Daily.csv")   # this is from the WHAT online baseflow separation filter for Pheasant branch

# summarize runoff to mean monthly mm
df.Q.d$Date <- mdy(df.Q.d$Date)
df.Q.d$year <- year(df.Q.d$Date)
df.Q.d$month <- month(df.Q.d$Date)
df.Q.mo <- dplyr::summarize(group_by(df.Q.d, year, month),
                            discharge.cfs.mean = mean(discharge.cfs, na.rm=T),
                            runoff.cfs.mean = mean(runoff.cfs, na.rm=T),
                            baseflow.cfs.mean = mean(baseflow.cfs, na.rm=T),
                            missing = sum(is.na(discharge.cfs)))

# for each month, figure out how many days are in that month
df.Q.mo$ndaypm <- days_in_month(ymd(paste0(df.Q.mo$year, "-", df.Q.mo$month, "-01")))

# convert discharge in cfs to total discharge in cubic meters
df.Q.mo$discharge.m3 <- df.Q.mo$discharge.cfs.mean*86400*df.Q.mo$ndaypm*(0.3048^3)
df.Q.mo$runoff.m3 <- df.Q.mo$runoff.cfs.mean*86400*df.Q.mo$ndaypm*(0.3048^3)
df.Q.mo$baseflow.m3 <- df.Q.mo$baseflow.cfs.mean*86400*df.Q.mo$ndaypm*(0.3048^3)

# convert discharge to mm
df.Q.mo$discharge.mm <- 1000*df.Q.mo$discharge.m3/area.upstream
df.Q.mo$runoff.mm <- 1000*df.Q.mo$runoff.m3/area.upstream
df.Q.mo$baseflow.mm <- 1000*df.Q.mo$baseflow.m3/area.upstream

# save parts of interest
write.csv(df.Q.mo[,c("year", "month", "discharge.mm", "runoff.mm", "baseflow.mm")], "PheasantBranch_BaseflowSeparation_Monthly.csv", row.names=F)

## Second: process meteorological data
# load GHCN data
df.GHCN <- read.csv("USW00014837_GapFilled.csv", stringsAsFactors=F)
df.GHCN$date <- mdy(df.GHCN$DATE)
df.GHCN$year <- year(df.GHCN$date)
df.GHCN$DOY <- yday(df.GHCN$date)

# load Arlington data, which will be used to correct GHCN data
df.ARL <- read.csv(paste0(git.dir, "Data/PheasantBranch/ArlingtonWI_Metdata.csv"), stringsAsFactors=F)
df.ARL$date <- ymd(df.ARL$date)

# set arlington NaNs
df.ARL$DToPcpn[df.ARL$DToPcpn==-6999 | df.ARL$DToPcpn==6999] <- NaN
df.ARL$DAvSol[df.ARL$DAvSol==-6999 | df.ARL$DAvSol==6999 | df.ARL$DAvSol==0 | df.ARL$DAvSol>400] <- NaN
df.ARL$DAvTAir[df.ARL$DAvTAir==-6999 | df.ARL$DAvTAir==6999] <- NaN
df.ARL$DMxTAir[df.ARL$DMxTAir==-6999 | df.ARL$DMxTAir==6999] <- NaN
df.ARL$DMnTAir[df.ARL$DMnTAir==-6999 | df.ARL$DMnTAir==6999] <- NaN
df.ARL$DAvRHum[df.ARL$DAvRHum==-6999 | df.ARL$DAvRHum==6999] <- NaN
df.ARL$DAvWind[df.ARL$DAvWind==-6999 | df.ARL$DAvWind==6999] <- NaN

# get rid of repeated dates
df.ARL <- unique(df.ARL)

# merge Arlington with GHCN
df.GHCN <- merge(df.GHCN, df.ARL, by=c("date"), all.x=T)

# station information from GHCN website (https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND:USW00014837/detail)
GHCN.lat <- 43.1405
GHCN.lon <- -89.3452
GHCN.elev <- 264 # [m]

# estimate solar radiation at edge of atmosphere based on latitude and DOY
df.GHCN$solar.toa <- PotentialSolar(lat=GHCN.lat*pi/180, Jday=df.GHCN$DOY)   # [kJ/m2/day]
df.GHCN$solar.toa.W_m2 <- df.GHCN$solar.toa*1000/86400

## calculate top-of-canopy radiation
# figure out max transmissivity (A parameter)
#p.trans.A <- 
#  ggplot(df.GHCN, aes(x=DAvSol/solar.toa.W_m2)) +
#  geom_histogram(binwidth=0.01) +
#  scale_x_continuous() +
#  geom_vline(xintercept=0.75, color="red")
trans.A <- 0.75  # the default, 0.75, seems pretty close based on the plot

# tune the C parameter
df.tune <- subset(df.GHCN, is.finite(solar.toa) & is.finite(DAvSol))
#qplot(df.tune$DAvSol, binwidth=2)
df.fit <- data.frame(trans.C = seq(0.1,5,0.025),
                     NSE = NaN,
                     RMSE = NaN)
for (trans.C in df.fit$trans.C){
  i <- which(df.fit$trans.C == trans.C) 
  df.tune$solar.canopy.W_m2 <- df.tune$solar.toa*(1000/86400)*
    transmissivity(Tx=df.tune$TMAX, Tn=df.tune$TMIN, A=trans.A, C=trans.C, opt="1day", JD=df.tune$DOY)
  df.fit$NSE[i] <- NashSutcliffe(df.tune$solar.canopy.W_m2, df.tune$DAvSol)
  df.fit$RMSE[i] <- RMSE(df.tune$solar.canopy.W_m2, df.tune$DAvSol)
}

# choose best trans.C
#qplot(RMSE, trans.C, data=df.fit)
trans.C <- df.fit$trans.C[which.min(df.fit$RMSE)]

df.tune$solar.canopy.W_m2 <- df.tune$solar.toa*(1000/86400)*
  transmissivity(Tx=df.tune$TMAX, Tn=df.tune$TMIN, A=trans.A, C=trans.C, opt="1day", JD=df.tune$DOY)

# linear bias adjustment
#qplot(DAvSol, solar.canopy.W_m2, data=df.tune) + geom_abline(intercept=0, slope=1, color="red") + stat_smooth(method="lm")
fit.tune <- lm(solar.canopy.W_m2 ~ DAvSol, data=df.tune)
df.tune$solar.canopy.W_m2.adj <- (df.tune$solar.canopy.W_m2-coef(fit.tune)[1])/coef(fit.tune)[2]
#qplot(DAvSol, solar.canopy.W_m2.adj, data=df.tune) + geom_abline(intercept=0, slope=1, color="red") + stat_smooth(method="lm")

# calculate solar radiation
df.GHCN$solar.canopy <- df.GHCN$solar.toa*transmissivity(Tx=df.GHCN$TMAX, Tn=df.GHCN$TMIN, A=trans.A, C=trans.C, opt="1day", JD=df.tune$DOY)
df.GHCN$rads <- df.GHCN$solar.canopy*1000/86400  # convert to incoming SW radiation [kJ/m2/day] to [W/m2]=[J/m2/s]
df.GHCN$rads <- (df.GHCN$rads-coef(fit.tune)[1])/coef(fit.tune)[2]
df.GHCN$rads[df.GHCN$rads<=0] <- 0
#qplot(DAvSol, rads, data=df.GHCN) + geom_abline(intercept=0, slope=1, color="red") + stat_smooth(method="lm")

# estimate daily actual vapor pressure based on the saturation vapor pressure at Tmin
# (assumes that air is ~saturated at dawn for dew formation), from FAO56 Eq. 48
df.GHCN$ea <- SatVaporPressure(df.GHCN$TMIN)

# estimate wind speed based on DOY
df.GHCN.longterm <- read.csv(paste0(git.dir, "Data/PheasantBranch/MadisonAirport_1939-2015.csv"))
df.GHCN.longterm$DATE <- ymd(df.GHCN.longterm$DATE)
df.GHCN.longterm$AWND[df.GHCN.longterm$AWND==-9999] <- NaN
df.GHCN.longterm$AWND <- df.GHCN.longterm$AWND*0.44704  # convert from mph to m/s (accidentally downloaded GHCN longterm in standard units)

# summarize wind speed to DOY
df.GHCN.longterm$DOY <- yday(df.GHCN.longterm$DATE)
df.GHCN.longterm.DOY <- dplyr::summarize(group_by(df.GHCN.longterm, DOY),
                                         wind = mean(AWND, na.rm=T))

# join with overall data frame
df.GHCN <- left_join(df.GHCN, df.GHCN.longterm.DOY[,c("DOY","wind")], by="DOY")

# make your data frame to plug in to Penman Monteith
df.PM <- data.frame(DATE = df.GHCN$date,
                    year = df.GHCN$year,
                    month = month(df.GHCN$date),
                    DOY = df.GHCN$DOY,
                    srad = df.GHCN$rads,
                    TMAX = df.GHCN$TMAX,
                    TMIN = df.GHCN$TMIN,
                    es_mean = (SatVaporPressure(df.GHCN$TMAX)+df.GHCN$ea)/2,
                    ea_mean = df.GHCN$ea,
                    prec = df.GHCN$PRCP,
                    wind = df.GHCN$wind,
                    LAT = GHCN.lat,
                    ELEV = GHCN.elev)

# calculate RET
ETo <- PenmanMonteith(df.PM)
df.PM$RET <- ETo$ETo

# save output
write.csv(df.PM, "USW00014837_GHCN_Daily.csv", row.names=F)

# summarize to monthly
df.PM.mo <- dplyr::summarize(group_by(df.PM, year, month),
                             prec.days = sum(prec>0),
                             prec.max = max(prec),
                             prec.gt.12 = sum(prec>12.7),
                             prec.gt.25 = sum(prec>25.4),
                             prec.gt.50 = sum(prec>50.8),
                             prec.gt.76 = sum(prec>76.2),
                             tmin = mean(TMIN),
                             tmax = mean(TMAX),
                             rads = mean(srad),
                             relh = mean(100*ea_mean/es_mean),
                             prec = sum(prec),
                             wspd = mean(wind),
                             VPD = mean(es_mean-ea_mean),
                             es = mean(es_mean),
                             ea = mean(ea_mean),
                             RET = sum(RET))

# save output
write.csv(df.PM.mo, "USW00014837_GHCN_Monthly.csv", row.names=F)

## Third: process Arlington GHCN meteorological data
# load GHCN data
df.GHCN.ARL <- read.csv("USC00470308_GapFilled.csv", stringsAsFactors=F)
df.GHCN.ARL$date <- ymd(df.GHCN.ARL$DATE)
df.GHCN.ARL$year <- year(df.GHCN.ARL$date)
df.GHCN.ARL$DOY <- yday(df.GHCN.ARL$date)

# make sure no TMAX<TMIN
df.GHCN.ARL$TMAX[df.GHCN.ARL$TMAX<df.GHCN.ARL$TMIN] <- df.GHCN.ARL$TMIN[df.GHCN.ARL$TMAX<df.GHCN.ARL$TMIN]

# merge Arlington with GHCN
df.GHCN.ARL <- merge(df.GHCN.ARL, df.ARL, by=c("date"), all.x=T)

# station information from GHCN website (https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND:USW00014837/detail)
GHCN.lat <- 43.3
GHCN.lon <- -89.3426
GHCN.elev <- 329 # [m]

# estimate solar radiation at edge of atmosphere based on latitude and DOYdf.GHCN.ARL$solar.toa <- PotentialSolar(lat=GHCN.lat*pi/180, Jday=df.GHCN.ARL$DOY)   # [kJ/m2/day]
df.GHCN.ARL$solar.toa <- PotentialSolar(lat=GHCN.lat*pi/180, Jday=df.GHCN.ARL$DOY)   # [kJ/m2/day]
df.GHCN.ARL$solar.toa.W_m2 <- df.GHCN.ARL$solar.toa*1000/86400

## calculate top-of-canopy radiation
# figure out max transmissivity (A parameter)
#p.trans.A <- 
#  ggplot(df.GHCN.ARL, aes(x=DAvSol/solar.toa.W_m2)) +
#  geom_histogram(binwidth=0.01) +
#  scale_x_continuous() +
#  geom_vline(xintercept=0.75, color="red")
trans.A <- 0.75  # the default, 0.75, seems pretty close based on the plot

# tune the C parameter
df.tune <- subset(df.GHCN.ARL, is.finite(solar.toa) & is.finite(DAvSol))
#qplot(df.tune$DAvSol, binwidth=2)
df.fit <- data.frame(trans.C = seq(0.1,5,0.025),
                     NSE = NaN,
                     RMSE = NaN)
for (trans.C in df.fit$trans.C){
  i <- which(df.fit$trans.C == trans.C) 
  df.tune$solar.canopy.W_m2 <- df.tune$solar.toa*(1000/86400)*
    transmissivity(Tx=df.tune$TMAX, Tn=df.tune$TMIN, A=trans.A, C=trans.C, opt="1day", JD=df.tune$DOY)
  df.fit$NSE[i] <- NashSutcliffe(df.tune$solar.canopy.W_m2, df.tune$DAvSol)
  df.fit$RMSE[i] <- RMSE(df.tune$solar.canopy.W_m2, df.tune$DAvSol)
}

# choose best trans.C
#qplot(RMSE, trans.C, data=df.fit)
trans.C <- df.fit$trans.C[which.min(df.fit$RMSE)]

df.tune$solar.canopy.W_m2 <- df.tune$solar.toa*(1000/86400)*
  transmissivity(Tx=df.tune$TMAX, Tn=df.tune$TMIN, A=trans.A, C=trans.C, opt="1day", JD=df.tune$DOY)

# linear bias adjustment
#qplot(DAvSol, solar.canopy.W_m2, data=df.tune) + geom_abline(intercept=0, slope=1, color="red") + stat_smooth(method="lm")
fit.tune <- lm(solar.canopy.W_m2 ~ DAvSol, data=df.tune)
df.tune$solar.canopy.W_m2.adj <- (df.tune$solar.canopy.W_m2-coef(fit.tune)[1])/coef(fit.tune)[2]
#qplot(DAvSol, solar.canopy.W_m2.adj, data=df.tune) + geom_abline(intercept=0, slope=1, color="red") + stat_smooth(method="lm")

# calculate solar radiation
df.GHCN.ARL$solar.canopy <- df.GHCN.ARL$solar.toa*transmissivity(Tx=df.GHCN.ARL$TMAX, Tn=df.GHCN.ARL$TMIN, A=trans.A, C=trans.C, opt="1day", JD=df.tune$DOY)
df.GHCN.ARL$rads <- df.GHCN.ARL$solar.canopy*1000/86400  # convert to incoming SW radiation [kJ/m2/day] to [W/m2]=[J/m2/s]
df.GHCN.ARL$rads <- (df.GHCN.ARL$rads-coef(fit.tune)[1])/coef(fit.tune)[2]
df.GHCN.ARL$rads[df.GHCN.ARL$rads<=0] <- 0
#qplot(DAvSol, rads, data=df.GHCN.ARL) + geom_abline(intercept=0, slope=1, color="red") + stat_smooth(method="lm")

# estimate daily actual vapor pressure based on the saturation vapor pressure at Tmin
# (assumes that air is ~saturated at dawn for dew formation), from FAO56 Eq. 48
df.GHCN.ARL$ea <- SatVaporPressure(df.GHCN.ARL$TMIN)

# join with overall data frame
df.GHCN.ARL <- left_join(df.GHCN.ARL, df.GHCN.longterm.DOY[,c("DOY","wind")], by="DOY")

# make your data frame to plug in to Penman Monteith
df.PM.ARL <- data.frame(DATE = df.GHCN.ARL$date,
                        year = df.GHCN.ARL$year,
                        month = month(df.GHCN.ARL$date),
                        DOY = df.GHCN.ARL$DOY,
                        srad = df.GHCN.ARL$rads,
                        TMAX = df.GHCN.ARL$TMAX,
                        TMIN = df.GHCN.ARL$TMIN,
                        es_mean = (SatVaporPressure(df.GHCN.ARL$TMAX)+df.GHCN.ARL$ea)/2,
                        ea_mean = df.GHCN.ARL$ea,
                        prec = df.GHCN.ARL$PRCP,
                        wind = df.GHCN.ARL$wind,
                        LAT = GHCN.lat,
                        ELEV = GHCN.elev)

# calculate RET
ETo <- PenmanMonteith(df.PM.ARL)
df.PM.ARL$RET <- ETo$ETo

# save output
write.csv(df.PM.ARL, "USC00470308_GHCN_Daily.csv", row.names=F)

# summarize to monthly
df.PM.ARL.mo <- dplyr::summarize(group_by(df.PM.ARL, year, month),
                                 prec.days = sum(prec>0),
                                 prec.max = max(prec),
                                 prec.gt.12 = sum(prec>12.7),
                                 prec.gt.25 = sum(prec>25.4),
                                 prec.gt.50 = sum(prec>50.8),
                                 prec.gt.76 = sum(prec>76.2),
                                 tmin = mean(TMIN),
                                 tmax = mean(TMAX),
                                 rads = mean(srad),
                                 relh = mean(100*ea_mean/es_mean),
                                 prec = sum(prec),
                                 wspd = mean(wind),
                                 VPD = mean(es_mean-ea_mean),
                                 es = mean(es_mean),
                                 ea = mean(ea_mean),
                                 RET = sum(RET))

# save output
write.csv(df.PM.ARL.mo, "USC00470308_GHCN_Monthly.csv", row.names=F)
