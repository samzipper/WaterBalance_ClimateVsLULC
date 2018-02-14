## PheasantBranch_3c_MonthlyRegressionsPLS.R
# This script is intended to read in monthly discharge and meteorological data for the Pheasant Branch USGS gauge
# and fit a partial least squares (PLS) regression. The PLS is calculated using only data from a historical 
# baseline period (1974-1995) and applied to the entire timeseries (1974-2016). The PLS is fit 250 times using
# different calibration years while retaining independent validation data.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

# load packages
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))
source(paste0(git.dir, "ProcessingScripts/CalculateClimatePLS.R"))   # this is the script that does most of the work...
require(ggplot2)
require(tidyr)
require(WVPlots)  # devtools::install_github('WinVector/WVPlots',build_vignettes=TRUE)
require(gridExtra)
require(lubridate)

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

# set significance threshold for pruning
p.thres <- 0.10

## decide "baseline" period (inclusive)
# logical end options:
#   1992: Gebert et al. (2012) identify breakpoint in Pheasant Branch streamflow starting in 1993
#   1995: splits data in half (21 years before, 21 years after) - this is the one i'll use
#   1999: 25 years
#   2002: 28 years (same as model)
#   2004: 30 years
#   
yr.baseline.start <- 1974
yr.baseline.end.all <- c(1995)

## which flux to analyze?
flux.name.all <- c("discharge.mm")
#flux.name.all <- c("discharge.mm", "runoff.mm", "baseflow.mm")

# read in discharge and met data frames
df.met <- read.csv(paste0(git.dir, "Data/PheasantBranch/USW00014837_GHCN_Monthly.csv"))  # use same met data as McFarland
df.Q <- read.csv("PheasantBranch_BaseflowSeparation_Monthly.csv")   # this is from the WHAT online baseflow separation filter for Pheasant branch

# merge data frames, with df.met defining the temporal extent
df <- merge(df.met, df.Q, all.x=T)

# make column for output variable
df$flux.est <- NaN

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

## data frame to store overall output
df.summary <- data.frame(month = numeric(0),
                         R2.cal = numeric(0),
                         R2.val = numeric(0),
                         RMSE.cal = numeric(0),
                         RMSE.val = numeric(0),
                         NRMSE.cal = numeric(0),
                         NRMSE.val = numeric(0))

# list of months
mo.list <- seq(1,12)
for (yr.baseline.end in yr.baseline.end.all){
  n.val.yr <- floor(length(seq(yr.baseline.start, yr.baseline.end))*0.25)
  for (flux.name in flux.name.all){
    for (mo in mo.list){
      
      # figure out number of predictors to retain, based on PC.keep
      PCs.keep <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_PC.keep_", sprintf("%02d", mo), ".csv"))
      min.PCs.keep <- length(PCs.keep$PC)
      
      # run regressions
      df.out.mo <- CalculateClimatePLS(df=df, mo=mo, flux.name=flux.name, var.options=var.options,
                                       yr.baseline.start=yr.baseline.start, yr.baseline.end=yr.baseline.end, n.val.yr=n.val.yr,
                                       p.thres=p.thres, neg.allowed=F, min.PCs.keep=min.PCs.keep, write.vars.keep=T, write.perm=F)
      
      # combine with other data
      if (exists("df.out")){
        df.out <- rbind(df.out, df.out.mo)
      } else {
        df.out <- df.out.mo
      }
      
      # status update
      print(paste0("flux ", flux.name, " year ", yr.baseline.end, " month ", mo, " complete"))
      
    }
    
    # move CSV output from CalculateClimatePLS script
    files <- list.files(paste0(git.dir, "Data/PheasantBranch/"), pattern="GHCN_MonthlyRegressions_PLS_.")
    file.rename(paste0(git.dir, "Data/PheasantBranch/", files), paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/", files))
    
    # save output csv
    write.csv(df.out, paste0(flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_PLS_OutputAll.csv"), row.names=F, quote=F)
    
    # fit metric
    RMSE(subset(df.out, group=="val")$PLS, subset(df.out, group=="val")$flux)
    NRMSE(subset(df.out, group=="val")$PLS, subset(df.out, group=="val")$flux)
    NashSutcliffe(subset(df.out, group=="val")$PLS, subset(df.out, group=="val")$flux)
    R2(subset(df.out, group=="val")$PLS, subset(df.out, group=="val")$flux)
    
    #rm(df.out)
  }
}