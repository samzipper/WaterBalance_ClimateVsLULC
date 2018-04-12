## PheasantBranch_7_SensitivityInputDataPLS.R
# This script is intended to quantify the sensitivity of the results for the Pheasant Branch watershed
# to the input data.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

# load packages
source(paste0(git.dir, "ProcessingScripts/PCA_Functions.R"))
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))
source(paste0(git.dir, "ProcessingScripts/CalculateClimatePLS.R"))
require(ggplot2)
require(dplyr)
require(tidyr)
require(gridExtra)
require(lubridate)
require(pls)
require(reshape2)

# function used for scaling data
scale.baseline <- function(x,year,yr.start,yr.end){
  # this function is intended to take in a vector x of data and associated 
  # year timestamps, and and scale x based on the mean and standard deviation only 
  # within the baseline period, defined by yr.start and yr.end
  
  # find measurements within baseline period
  i.baseline <- (year >= yr.start & year <= yr.end)
  
  # find mean and standard deviation of baseline period
  mean.baseline <- mean(x[i.baseline], na.rm=T)
  sd.baseline <- sd(x[i.baseline], na.rm=T)
  
  # scale whole record based on baseline period
  x.scale <- (x-mean.baseline)/sd.baseline
  
  return(x.scale)
}

# path to folder with output data
setwd(paste0(git.dir, "Data/PheasantBranch/"))

# directory to save plots
plot.dir <- paste0(git.dir, "Plots/PheasantBranch/")

# which flux to analyze?
flux.name <- "discharge.mm"

# set significance threshold for pruning
p.thres <- 0.10
cum.var <- 0.8

## define "baseline" period (inclusive)
yr.baseline.start <- 1974
yr.baseline.end <- 1995

# number of permutations to run
n.perm <- 250

# read in discharge and met data frames
df.met.GHCN <- read.csv(paste0(git.dir, "Data/PheasantBranch/USW00014837_GHCN_Monthly.csv"))  # Madison airport met data
df.met.ARL <- read.csv(paste0(git.dir, "Data/PheasantBranch/USC00470308_GHCN_Monthly.csv"))  # Arlington met data

# process to make sure it is the same columns and years
df.met.GHCN <- df.met.GHCN[, intersect(names(df.met.GHCN), names(df.met.ARL))]
df.met.ARL <- df.met.ARL[, intersect(names(df.met.GHCN), names(df.met.ARL))]

df.met.GHCN <- subset(df.met.GHCN, year %in% intersect(unique(df.met.GHCN$year), unique(df.met.ARL$year)))
df.met.ARL <- subset(df.met.ARL, year %in% intersect(unique(df.met.GHCN$year), unique(df.met.ARL$year)))

# read in discharge data
df.Q <- read.csv("PheasantBranch_BaseflowSeparation_Monthly.csv")   # this is from the WHAT online baseflow separation filter for Pheasant branch

## scroll through different met datasets
for (met.data in c("GHCN", "ARL")){
  # merge data frames, with df.met defining the temporal extent
  if (met.data=="GHCN") df <- merge(df.met.GHCN, df.Q, all.x=T)
  if (met.data=="ARL") df <- merge(df.met.ARL, df.Q, all.x=T)
  
  # make column for output variable
  df$discharge.est <- NaN
  
  ## calculate a bunch of meteorological variables
  # calculate deficit variable
  df$defc <- df$prec - df$RET
  
  ## calculate a bunch of variables that just keep track of the preceding months
  df$prec.1mo <- NaN  # cumulative precipitation for current + preceding month
  df$prec.2mo <- NaN  # cumulative precipitation for current + preceding 2 months
  df$prec.3mo <- NaN  # cumulative precipitation for current + preceding 3 months
  df$prec.6mo <- NaN  # cumulative precipitation for current + preceding 6 months
  df$prec.12mo <- NaN  # cumulative precipitation for current + preceding 6 months
  df$tmin.1mo <- NaN  # average tmin for current + preceding month
  df$tmin.2mo <- NaN  # average
  df$tmin.3mo <- NaN  # average
  df$tmin.6mo <- NaN  # average
  df$tmin.12mo <- NaN  # average
  df$tmax.1mo <- NaN  # average
  df$tmax.2mo <- NaN  # average
  df$tmax.3mo <- NaN  # average
  df$tmax.6mo <- NaN  # average
  df$tmax.12mo <- NaN  # average
  df$relh.1mo <- NaN  # average
  df$relh.2mo <- NaN  # average
  df$relh.3mo <- NaN  # average
  df$relh.6mo <- NaN  # average
  df$relh.12mo <- NaN  # average
  df$rads.1mo <- NaN  # average
  df$rads.2mo <- NaN  # average
  df$rads.3mo <- NaN  # average
  df$rads.6mo <- NaN  # average
  df$rads.12mo <- NaN  # average
  df$wspd.1mo <- NaN  # average
  df$wspd.2mo <- NaN  # average
  df$wspd.3mo <- NaN  # average
  df$wspd.6mo <- NaN  # average
  df$wspd.12mo <- NaN  # average
  df$ea.1mo <- NaN    # average
  df$ea.2mo <- NaN    # average
  df$ea.3mo <- NaN    # average
  df$ea.6mo <- NaN    # average
  df$ea.12mo <- NaN    # average
  df$es.1mo <- NaN    # average
  df$es.2mo <- NaN    # average
  df$es.3mo <- NaN    # average
  df$es.6mo <- NaN    # average
  df$es.12mo <- NaN    # average
  df$VPD.1mo <- NaN   # average
  df$VPD.2mo <- NaN   # average
  df$VPD.3mo <- NaN   # average
  df$VPD.6mo <- NaN   # average
  df$VPD.12mo <- NaN   # average
  df$RET.1mo <- NaN   # cumulative
  df$RET.2mo <- NaN   # cumulative
  df$RET.3mo <- NaN   # cumulative
  df$RET.6mo <- NaN   # cumulative
  df$RET.12mo <- NaN   # cumulative
  df$defc.1mo <- NaN  # cumulative
  df$defc.2mo <- NaN  # cumulative
  df$defc.3mo <- NaN  # cumulative
  df$defc.6mo <- NaN  # cumulative
  df$defc.12mo <- NaN  # cumulative
  for (j in 1:dim(df)[1]){
    if (j >= 2){
      df$prec.1mo[j] <- sum(df$prec[(j-1):j])
      df$tmin.1mo[j] <- mean(df$tmin[(j-1):j])
      df$tmax.1mo[j] <- mean(df$tmax[(j-1):j])
      df$relh.1mo[j] <- mean(df$relh[(j-1):j])
      df$rads.1mo[j] <- mean(df$rads[(j-1):j])
      df$wspd.1mo[j] <- mean(df$wspd[(j-1):j])
      df$ea.1mo[j]   <- mean(df$ea[(j-1):j])
      df$es.1mo[j]   <- mean(df$es[(j-1):j])
      df$VPD.1mo[j]  <- mean(df$VPD[(j-1):j])
      df$RET.1mo[j]  <- sum(df$RET[(j-1):j])
      df$defc.1mo[j] <- sum(df$defc[(j-1):j])
    }
    if (j >= 3){
      df$prec.2mo[j] <- sum(df$prec[(j-2):j])
      df$tmin.2mo[j] <- mean(df$tmin[(j-2):j])
      df$tmax.2mo[j] <- mean(df$tmax[(j-2):j])
      df$relh.2mo[j] <- mean(df$relh[(j-2):j])
      df$rads.2mo[j] <- mean(df$rads[(j-2):j])
      df$wspd.2mo[j] <- mean(df$wspd[(j-2):j])
      df$ea.2mo[j]   <- mean(df$ea[(j-2):j])
      df$es.2mo[j]   <- mean(df$es[(j-2):j])
      df$VPD.2mo[j]  <- mean(df$VPD[(j-2):j])
      df$RET.2mo[j]  <- sum(df$RET[(j-2):j])
      df$defc.2mo[j] <- sum(df$defc[(j-2):j])
    }
    if (j >= 4){
      df$prec.3mo[j] <- sum(df$prec[(j-3):j])
      df$tmin.3mo[j] <- mean(df$tmin[(j-3):j])
      df$tmax.3mo[j] <- mean(df$tmax[(j-3):j])
      df$relh.3mo[j] <- mean(df$relh[(j-3):j])
      df$rads.3mo[j] <- mean(df$rads[(j-3):j])
      df$wspd.3mo[j] <- mean(df$wspd[(j-3):j])
      df$ea.3mo[j]   <- mean(df$ea[(j-3):j])
      df$es.3mo[j]   <- mean(df$es[(j-3):j])
      df$VPD.3mo[j]  <- mean(df$VPD[(j-3):j])
      df$RET.3mo[j]  <- sum(df$RET[(j-3):j])
      df$defc.3mo[j] <- sum(df$defc[(j-3):j])
    }
    if (j >= 7){
      df$prec.6mo[j] <- sum(df$prec[(j-6):j])
      df$tmin.6mo[j] <- mean(df$tmin[(j-6):j])
      df$tmax.6mo[j] <- mean(df$tmax[(j-6):j])
      df$relh.6mo[j] <- mean(df$relh[(j-6):j])
      df$rads.6mo[j] <- mean(df$rads[(j-6):j])
      df$wspd.6mo[j] <- mean(df$wspd[(j-6):j])
      df$ea.6mo[j]   <- mean(df$ea[(j-6):j])
      df$es.6mo[j]   <- mean(df$es[(j-6):j])
      df$VPD.6mo[j]  <- mean(df$VPD[(j-6):j])
      df$RET.6mo[j]  <- sum(df$RET[(j-6):j])
      df$defc.6mo[j] <- sum(df$defc[(j-6):j])
    }
    if (j >= 13){
      df$prec.12mo[j] <- sum(df$prec[(j-12):j])
      df$tmin.12mo[j] <- mean(df$tmin[(j-12):j])
      df$tmax.12mo[j] <- mean(df$tmax[(j-12):j])
      df$relh.12mo[j] <- mean(df$relh[(j-12):j])
      df$rads.12mo[j] <- mean(df$rads[(j-12):j])
      df$wspd.12mo[j] <- mean(df$wspd[(j-12):j])
      df$ea.12mo[j]   <- mean(df$ea[(j-12):j])
      df$es.12mo[j]   <- mean(df$es[(j-12):j])
      df$VPD.12mo[j]  <- mean(df$VPD[(j-12):j])
      df$RET.12mo[j]  <- sum(df$RET[(j-12):j])
      df$defc.12mo[j] <- sum(df$defc[(j-12):j])
    }
  }
  
  # some other candidate variables
  df$prec.sq <- df$prec^2
  df$defc.sq <- df$defc^2
  df$defc.1mo.sq <- df$defc.1mo^2
  df$defc.2mo.sq <- df$defc.2mo^2
  df$defc.3mo.sq <- df$defc.3mo^2
  df$defc.6mo.sq <- df$defc.6mo^2
  df$defc.12mo.sq <- df$defc.12mo^2
  
  # get list of all possible predictor variables
  var.options <- colnames(df)[!(colnames(df) %in% c("year","month","discharge.mm", "runoff.mm", "baseflow.mm", "flux.est", "prec.gt.101"))]
  
  # list of months
  mo.list <- seq(1,12)
  # number of years to save for validation
  n.val.yr <- floor(length(seq(yr.baseline.start, yr.baseline.end))*0.25)
  for (mo in mo.list){
    
    # run PLS script
    df.out.mo <- CalculateClimatePLS(df=df, mo=mo, flux.name=flux.name, var.options=var.options,
                                     yr.baseline.start=yr.baseline.start, yr.baseline.end=yr.baseline.end,
                                     n.perm=n.perm, n.val.yr=n.val.yr, p.thres=p.thres, 
                                     neg.allowed=F, write.vars.keep=F, write.perm=F)
    df.out.mo$met.data <- met.data

    # combined with other data
    if (exists("df.out")){
      df.out <- rbind(df.out, df.out.mo)
    } else {
      df.out <- df.out.mo
    }
    
    # status update
    print(paste0(flux.name, " met data ", met.data, " month ", mo, " complete"))
    
  }
}

# save output csv
write.csv(df.out, paste0(flux.name, "/GHCN_SensitivityInputData_PLS_OutputAll.csv"), row.names=F, quote=F)

## calculate fit metrics by yr.baseline.end
df.baseline <- dplyr::summarize(group_by(subset(df.out, group != "prediction"), met.data, group),
                                RMSE = RMSE(PLS, flux),
                                NSE = NashSutcliffe(PLS, flux))
df.baseline.melt <- melt(df.baseline, id=c("group", "met.data"))

## calculate baseline period discharge
df.year.perm.baseline <-
  dplyr::summarize(group_by(df.out, year, perm, met.data),
                   obs.sum = sum(flux),
                   PLS.sum = sum(PLS))
df.year.baseline <-
  dplyr::summarize(group_by(df.year.perm.baseline, year, met.data),
                   obs.mean = mean(obs.sum),
                   PLS.mean = mean(PLS.sum))

df.Q.baseline <-
  dplyr::summarize(group_by(df.year.baseline, met.data),
                   obs.baseline.mean = mean(obs.mean[year<=yr.baseline.end]))

## separate LULC vs climate
df.LULCvClimate <- merge(df.year.perm.baseline, df.Q.baseline, by=c("met.data"), all.x=T)
df.LULCvClimate$change.overall <- df.LULCvClimate$obs.sum - df.LULCvClimate$obs.baseline.mean
df.LULCvClimate$change.climate <- df.LULCvClimate$PLS.sum - df.LULCvClimate$obs.baseline.mean
df.LULCvClimate$change.LULC <- df.LULCvClimate$change.overall - df.LULCvClimate$change.climate

df.LULCvClimate.baseline <-
  dplyr::summarize(group_by(df.LULCvClimate, year, met.data),
                   obs = mean(obs.sum),
                   PLS.mean = mean(PLS.sum),
                   PLS.sd = sd(PLS.sum),
                   change.overall.mean = mean(change.overall),
                   change.overall.sd = sd(change.overall),
                   change.climate.mean = mean(change.climate),
                   change.climate.sd = sd(change.climate),
                   change.LULC.mean = mean(change.LULC),
                   change.LULC.sd = sd(change.LULC))

df.LULCvClimate.prediction <- subset(df.LULCvClimate.baseline, year>yr.baseline.end)

df.LULCvClimate.baseline.melt <- melt(df.LULCvClimate.baseline[,c("year", "met.data", "change.overall.mean", "change.climate.mean", "change.LULC.mean")],
                                      id=c("year", "met.data"))

# calculate some min/max columns to use for ribbons
df.LULCvClimate.baseline$change.overall.min <- df.LULCvClimate.baseline$change.overall.mean - df.LULCvClimate.baseline$change.overall.sd
df.LULCvClimate.baseline$change.overall.max <- df.LULCvClimate.baseline$change.overall.mean + df.LULCvClimate.baseline$change.overall.sd
df.LULCvClimate.baseline$change.climate.min <- df.LULCvClimate.baseline$change.climate.mean - df.LULCvClimate.baseline$change.climate.sd
df.LULCvClimate.baseline$change.climate.max <- df.LULCvClimate.baseline$change.climate.mean + df.LULCvClimate.baseline$change.climate.sd
df.LULCvClimate.baseline$change.LULC.min <- df.LULCvClimate.baseline$change.LULC.mean - df.LULCvClimate.baseline$change.LULC.sd
df.LULCvClimate.baseline$change.LULC.max <- df.LULCvClimate.baseline$change.LULC.mean + df.LULCvClimate.baseline$change.LULC.sd

# save output csv
write.csv(df.LULCvClimate.baseline, paste0(flux.name, "/GHCN_SensitivityInputData_PLS_ClimateVsLULC.csv"), row.names=F, quote=F)

p.change.climate <- 
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_ribbon(data=subset(df.LULCvClimate.baseline, met.data=="GHCN"), 
              aes(x=year, ymin=change.climate.min, ymax=change.climate.max), fill="red", alpha=0.5) +
  geom_ribbon(data=subset(df.LULCvClimate.baseline, met.data=="ARL"), 
              aes(x=year, ymin=change.climate.min, ymax=change.climate.max), fill="blue", alpha=0.5) +
  geom_line(data=subset(df.LULCvClimate.baseline, met.data=="ARL"),
            aes(x=year, y=change.overall.mean), color="black") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))
