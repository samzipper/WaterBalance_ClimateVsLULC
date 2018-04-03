## Yahara_1b_MonthlyRegressions-MLR.R
# This script is intended to load monthly model output data (generated using the script
# Scenarios_SummarizeMonthlyWaterBalance.R) and calculate a separate statistical relationship 
# between meteorological data and a water balance output variable for each month using 
# multiple linear regression. A permutation-based approach is used to provide independent
# estimates of model fit and uncertainty associated with estimates.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

# load packages
source("C:/Users/Sam/Dropbox/Work/R/lmp.R")
source("C:/Users/Sam/Dropbox/Work/R/FitMetrics.R")
source(paste0(git.dir, "ProcessingScripts/CalculateClimateMLR.R"))   # this is the script that does most of the work...
require(ggplot2)
require(tidyr)
require(WVPlots)  # devtools::install_github('WinVector/WVPlots',build_vignettes=TRUE)
require(gridExtra)
require(lubridate)

setwd(paste0(git.dir, "Data/Yahara-AgroIBIS/"))

# directory to save plots
plot.dir <- paste0(git.dir, "Plots/Yahara-AgroIBIS/")

# use HistLULC_XXclim simulations?
hist.LULC <- F   # T or F

# set significance threshold for pruning
p.thres <- 0.10

# variance for selecting options
cum.var <- 0.8
min.var <- 0.01

# baseline period
yr.baseline.start <- 1986
yr.baseline.end <- 2013
n.baseline.yrs <- length(seq(yr.baseline.start, yr.baseline.end))

# number of permutations
n.perm <- 250

# which fluxes to analyze?
flux.name.all <- c("aet", "srunoff", "drainage")

# which months to analyze?
mo.all <- seq(1,12)

# which scenarios to analyze?
scen.climate.all <- c("AI", "AR", "CC", "NW")
scen.LULC.all <- c("AI", "AR", "CC", "NW")

# which output to write?
write.vars.keep <- T
write.PC.keep <- T
write.perm <- F
PC.plot <- F

# initialize logical to determine whether to write historical data as CSV
hist.csv <- T   # T means write, F means don't

## first: preprocess HistLULC_XXclim datasets
# use the new HistLULC_XXclim dataset that is based on historical land use for the entire scenario
df.AI <- read.csv("input/HistLULC_AIclim_waterbalance_monthly.csv")
df.AR <- read.csv("input/HistLULC_ARClim_waterbalance_monthly.csv")
df.CC <- read.csv("input/HistLULC_CCclim_waterbalance_monthly.csv")
df.NW <- read.csv("input/HistLULC_NWclim_waterbalance_monthly.csv")

# convert daily mean values to cumulative monthly totals for depth variables
df.AI$ndaypm <- days_in_month(ymd(paste0(df.AI$year, "-", df.AI$month, "-01")))
df.AR$ndaypm <- days_in_month(ymd(paste0(df.AR$year, "-", df.AR$month, "-01")))
df.CC$ndaypm <- days_in_month(ymd(paste0(df.CC$year, "-", df.CC$month, "-01")))
df.NW$ndaypm <- days_in_month(ymd(paste0(df.NW$year, "-", df.NW$month, "-01")))

df.AI$aet <- df.AI$aet*df.AI$ndaypm
df.AR$aet <- df.AR$aet*df.AR$ndaypm
df.CC$aet <- df.CC$aet*df.CC$ndaypm
df.NW$aet <- df.NW$aet*df.NW$ndaypm

df.AI$srunoff <- df.AI$srunoff*df.AI$ndaypm
df.AR$srunoff <- df.AR$srunoff*df.AR$ndaypm
df.CC$srunoff <- df.CC$srunoff*df.CC$ndaypm
df.NW$srunoff <- df.NW$srunoff*df.NW$ndaypm

df.AI$drainage <- df.AI$drainage*df.AI$ndaypm
df.AR$drainage <- df.AR$drainage*df.AR$ndaypm
df.CC$drainage <- df.CC$drainage*df.CC$ndaypm
df.NW$drainage <- df.NW$drainage*df.NW$ndaypm

# calculate deficit variable
df.AI$defc <- df.AI$prec - df.AI$RET
df.AR$defc <- df.AR$prec - df.AR$RET
df.CC$defc <- df.CC$prec - df.CC$RET
df.NW$defc <- df.NW$prec - df.NW$RET

## calculate a bunch of variables that just keep track of the preceding months
df.AI$prec.1mo <- NaN  # cumulative precipitation for current + preceding month
df.AI$prec.2mo <- NaN  # cumulative precipitation for current + preceding 2 months
df.AI$prec.3mo <- NaN  # cumulative precipitation for current + preceding 3 months
df.AI$prec.6mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.AI$prec.12mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.AI$tmin.1mo <- NaN  # average tmin for current + preceding month
df.AI$tmin.2mo <- NaN  # average
df.AI$tmin.3mo <- NaN  # average
df.AI$tmin.6mo <- NaN  # average
df.AI$tmin.12mo <- NaN  # average
df.AI$tmax.1mo <- NaN  # average
df.AI$tmax.2mo <- NaN  # average
df.AI$tmax.3mo <- NaN  # average
df.AI$tmax.6mo <- NaN  # average
df.AI$tmax.12mo <- NaN  # average
df.AI$relh.1mo <- NaN  # average
df.AI$relh.2mo <- NaN  # average
df.AI$relh.3mo <- NaN  # average
df.AI$relh.6mo <- NaN  # average
df.AI$relh.12mo <- NaN  # average
df.AI$rads.1mo <- NaN  # average
df.AI$rads.2mo <- NaN  # average
df.AI$rads.3mo <- NaN  # average
df.AI$rads.6mo <- NaN  # average
df.AI$rads.12mo <- NaN  # average
df.AI$wspd.1mo <- NaN  # average
df.AI$wspd.2mo <- NaN  # average
df.AI$wspd.3mo <- NaN  # average
df.AI$wspd.6mo <- NaN  # average
df.AI$wspd.12mo <- NaN  # average
df.AI$ea.1mo <- NaN    # average
df.AI$ea.2mo <- NaN    # average
df.AI$ea.3mo <- NaN    # average
df.AI$ea.6mo <- NaN    # average
df.AI$ea.12mo <- NaN    # average
df.AI$es.1mo <- NaN    # average
df.AI$es.2mo <- NaN    # average
df.AI$es.3mo <- NaN    # average
df.AI$es.6mo <- NaN    # average
df.AI$es.12mo <- NaN    # average
df.AI$VPD.1mo <- NaN   # average
df.AI$VPD.2mo <- NaN   # average
df.AI$VPD.3mo <- NaN   # average
df.AI$VPD.6mo <- NaN   # average
df.AI$VPD.12mo <- NaN   # average
df.AI$RET.1mo <- NaN   # cumulative
df.AI$RET.2mo <- NaN   # cumulative
df.AI$RET.3mo <- NaN   # cumulative
df.AI$RET.6mo <- NaN   # cumulative
df.AI$RET.12mo <- NaN   # cumulative
df.AI$defc.1mo <- NaN  # cumulative
df.AI$defc.2mo <- NaN  # cumulative
df.AI$defc.3mo <- NaN  # cumulative
df.AI$defc.6mo <- NaN  # cumulative
df.AI$defc.12mo <- NaN  # cumulative
df.AR$prec.1mo <- NaN  # cumulative precipitation for current + preceding month
df.AR$prec.2mo <- NaN  # cumulative precipitation for current + preceding 2 months
df.AR$prec.3mo <- NaN  # cumulative precipitation for current + preceding 3 months
df.AR$prec.6mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.AR$prec.12mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.AR$tmin.1mo <- NaN  # average tmin for current + preceding month
df.AR$tmin.2mo <- NaN  # average
df.AR$tmin.3mo <- NaN  # average
df.AR$tmin.6mo <- NaN  # average
df.AR$tmin.12mo <- NaN  # average
df.AR$tmax.1mo <- NaN  # average
df.AR$tmax.2mo <- NaN  # average
df.AR$tmax.3mo <- NaN  # average
df.AR$tmax.6mo <- NaN  # average
df.AR$tmax.12mo <- NaN  # average
df.AR$relh.1mo <- NaN  # average
df.AR$relh.2mo <- NaN  # average
df.AR$relh.3mo <- NaN  # average
df.AR$relh.6mo <- NaN  # average
df.AR$relh.12mo <- NaN  # average
df.AR$rads.1mo <- NaN  # average
df.AR$rads.2mo <- NaN  # average
df.AR$rads.3mo <- NaN  # average
df.AR$rads.6mo <- NaN  # average
df.AR$rads.12mo <- NaN  # average
df.AR$wspd.1mo <- NaN  # average
df.AR$wspd.2mo <- NaN  # average
df.AR$wspd.3mo <- NaN  # average
df.AR$wspd.6mo <- NaN  # average
df.AR$wspd.12mo <- NaN  # average
df.AR$ea.1mo <- NaN    # average
df.AR$ea.2mo <- NaN    # average
df.AR$ea.3mo <- NaN    # average
df.AR$ea.6mo <- NaN    # average
df.AR$ea.12mo <- NaN    # average
df.AR$es.1mo <- NaN    # average
df.AR$es.2mo <- NaN    # average
df.AR$es.3mo <- NaN    # average
df.AR$es.6mo <- NaN    # average
df.AR$es.12mo <- NaN    # average
df.AR$VPD.1mo <- NaN   # average
df.AR$VPD.2mo <- NaN   # average
df.AR$VPD.3mo <- NaN   # average
df.AR$VPD.6mo <- NaN   # average
df.AR$VPD.12mo <- NaN   # average
df.AR$RET.1mo <- NaN   # cumulative
df.AR$RET.2mo <- NaN   # cumulative
df.AR$RET.3mo <- NaN   # cumulative
df.AR$RET.6mo <- NaN   # cumulative
df.AR$RET.12mo <- NaN   # cumulative
df.AR$defc.1mo <- NaN  # cumulative
df.AR$defc.2mo <- NaN  # cumulative
df.AR$defc.3mo <- NaN  # cumulative
df.AR$defc.6mo <- NaN  # cumulative
df.AR$defc.12mo <- NaN  # cumulative
df.CC$prec.1mo <- NaN  # cumulative precipitation for current + preceding month
df.CC$prec.2mo <- NaN  # cumulative precipitation for current + preceding 2 months
df.CC$prec.3mo <- NaN  # cumulative precipitation for current + preceding 3 months
df.CC$prec.6mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.CC$prec.12mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.CC$tmin.1mo <- NaN  # average tmin for current + preceding month
df.CC$tmin.2mo <- NaN  # average
df.CC$tmin.3mo <- NaN  # average
df.CC$tmin.6mo <- NaN  # average
df.CC$tmin.12mo <- NaN  # average
df.CC$tmax.1mo <- NaN  # average
df.CC$tmax.2mo <- NaN  # average
df.CC$tmax.3mo <- NaN  # average
df.CC$tmax.6mo <- NaN  # average
df.CC$tmax.12mo <- NaN  # average
df.CC$relh.1mo <- NaN  # average
df.CC$relh.2mo <- NaN  # average
df.CC$relh.3mo <- NaN  # average
df.CC$relh.6mo <- NaN  # average
df.CC$relh.12mo <- NaN  # average
df.CC$rads.1mo <- NaN  # average
df.CC$rads.2mo <- NaN  # average
df.CC$rads.3mo <- NaN  # average
df.CC$rads.6mo <- NaN  # average
df.CC$rads.12mo <- NaN  # average
df.CC$wspd.1mo <- NaN  # average
df.CC$wspd.2mo <- NaN  # average
df.CC$wspd.3mo <- NaN  # average
df.CC$wspd.6mo <- NaN  # average
df.CC$wspd.12mo <- NaN  # average
df.CC$ea.1mo <- NaN    # average
df.CC$ea.2mo <- NaN    # average
df.CC$ea.3mo <- NaN    # average
df.CC$ea.6mo <- NaN    # average
df.CC$ea.12mo <- NaN    # average
df.CC$es.1mo <- NaN    # average
df.CC$es.2mo <- NaN    # average
df.CC$es.3mo <- NaN    # average
df.CC$es.6mo <- NaN    # average
df.CC$es.12mo <- NaN    # average
df.CC$VPD.1mo <- NaN   # average
df.CC$VPD.2mo <- NaN   # average
df.CC$VPD.3mo <- NaN   # average
df.CC$VPD.6mo <- NaN   # average
df.CC$VPD.12mo <- NaN   # average
df.CC$RET.1mo <- NaN   # cumulative
df.CC$RET.2mo <- NaN   # cumulative
df.CC$RET.3mo <- NaN   # cumulative
df.CC$RET.6mo <- NaN   # cumulative
df.CC$RET.12mo <- NaN   # cumulative
df.CC$defc.1mo <- NaN  # cumulative
df.CC$defc.2mo <- NaN  # cumulative
df.CC$defc.3mo <- NaN  # cumulative
df.CC$defc.6mo <- NaN  # cumulative
df.CC$defc.12mo <- NaN  # cumulative
df.NW$prec.1mo <- NaN  # cumulative precipitation for current + preceding month
df.NW$prec.2mo <- NaN  # cumulative precipitation for current + preceding 2 months
df.NW$prec.3mo <- NaN  # cumulative precipitation for current + preceding 3 months
df.NW$prec.6mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.NW$prec.12mo <- NaN  # cumulative precipitation for current + preceding 6 months
df.NW$tmin.1mo <- NaN  # average tmin for current + preceding month
df.NW$tmin.2mo <- NaN  # average
df.NW$tmin.3mo <- NaN  # average
df.NW$tmin.6mo <- NaN  # average
df.NW$tmin.12mo <- NaN  # average
df.NW$tmax.1mo <- NaN  # average
df.NW$tmax.2mo <- NaN  # average
df.NW$tmax.3mo <- NaN  # average
df.NW$tmax.6mo <- NaN  # average
df.NW$tmax.12mo <- NaN  # average
df.NW$relh.1mo <- NaN  # average
df.NW$relh.2mo <- NaN  # average
df.NW$relh.3mo <- NaN  # average
df.NW$relh.6mo <- NaN  # average
df.NW$relh.12mo <- NaN  # average
df.NW$rads.1mo <- NaN  # average
df.NW$rads.2mo <- NaN  # average
df.NW$rads.3mo <- NaN  # average
df.NW$rads.6mo <- NaN  # average
df.NW$rads.12mo <- NaN  # average
df.NW$wspd.1mo <- NaN  # average
df.NW$wspd.2mo <- NaN  # average
df.NW$wspd.3mo <- NaN  # average
df.NW$wspd.6mo <- NaN  # average
df.NW$wspd.12mo <- NaN  # average
df.NW$ea.1mo <- NaN    # average
df.NW$ea.2mo <- NaN    # average
df.NW$ea.3mo <- NaN    # average
df.NW$ea.6mo <- NaN    # average
df.NW$ea.12mo <- NaN    # average
df.NW$es.1mo <- NaN    # average
df.NW$es.2mo <- NaN    # average
df.NW$es.3mo <- NaN    # average
df.NW$es.6mo <- NaN    # average
df.NW$es.12mo <- NaN    # average
df.NW$VPD.1mo <- NaN   # average
df.NW$VPD.2mo <- NaN   # average
df.NW$VPD.3mo <- NaN   # average
df.NW$VPD.6mo <- NaN   # average
df.NW$VPD.12mo <- NaN   # average
df.NW$RET.1mo <- NaN   # cumulative
df.NW$RET.2mo <- NaN   # cumulative
df.NW$RET.3mo <- NaN   # cumulative
df.NW$RET.6mo <- NaN   # cumulative
df.NW$RET.12mo <- NaN   # cumulative
df.NW$defc.1mo <- NaN  # cumulative
df.NW$defc.2mo <- NaN  # cumulative
df.NW$defc.3mo <- NaN  # cumulative
df.NW$defc.6mo <- NaN  # cumulative
df.NW$defc.12mo <- NaN  # cumulative
for (j in 1:dim(df.AI)[1]){
  if (j >= 2){
    df.AI$prec.1mo[j] <- sum(df.AI$prec[(j-1):j])
    df.AI$tmin.1mo[j] <- mean(df.AI$tmin[(j-1):j])
    df.AI$tmax.1mo[j] <- mean(df.AI$tmax[(j-1):j])
    df.AI$relh.1mo[j] <- mean(df.AI$relh[(j-1):j])
    df.AI$rads.1mo[j] <- mean(df.AI$rads[(j-1):j])
    df.AI$wspd.1mo[j] <- mean(df.AI$wspd[(j-1):j])
    df.AI$ea.1mo[j]   <- mean(df.AI$ea[(j-1):j])
    df.AI$es.1mo[j]   <- mean(df.AI$es[(j-1):j])
    df.AI$VPD.1mo[j]  <- mean(df.AI$VPD[(j-1):j])
    df.AI$RET.1mo[j]  <- sum(df.AI$RET[(j-1):j])
    df.AI$defc.1mo[j] <- sum(df.AI$defc[(j-1):j])
    df.AR$prec.1mo[j] <- sum(df.AR$prec[(j-1):j])
    df.AR$tmin.1mo[j] <- mean(df.AR$tmin[(j-1):j])
    df.AR$tmax.1mo[j] <- mean(df.AR$tmax[(j-1):j])
    df.AR$relh.1mo[j] <- mean(df.AR$relh[(j-1):j])
    df.AR$rads.1mo[j] <- mean(df.AR$rads[(j-1):j])
    df.AR$wspd.1mo[j] <- mean(df.AR$wspd[(j-1):j])
    df.AR$ea.1mo[j]   <- mean(df.AR$ea[(j-1):j])
    df.AR$es.1mo[j]   <- mean(df.AR$es[(j-1):j])
    df.AR$VPD.1mo[j]  <- mean(df.AR$VPD[(j-1):j])
    df.AR$RET.1mo[j]  <- sum(df.AR$RET[(j-1):j])
    df.AR$defc.1mo[j] <- sum(df.AR$defc[(j-1):j])
    df.CC$prec.1mo[j] <- sum(df.CC$prec[(j-1):j])
    df.CC$tmin.1mo[j] <- mean(df.CC$tmin[(j-1):j])
    df.CC$tmax.1mo[j] <- mean(df.CC$tmax[(j-1):j])
    df.CC$relh.1mo[j] <- mean(df.CC$relh[(j-1):j])
    df.CC$rads.1mo[j] <- mean(df.CC$rads[(j-1):j])
    df.CC$wspd.1mo[j] <- mean(df.CC$wspd[(j-1):j])
    df.CC$ea.1mo[j]   <- mean(df.CC$ea[(j-1):j])
    df.CC$es.1mo[j]   <- mean(df.CC$es[(j-1):j])
    df.CC$VPD.1mo[j]  <- mean(df.CC$VPD[(j-1):j])
    df.CC$RET.1mo[j]  <- sum(df.CC$RET[(j-1):j])
    df.CC$defc.1mo[j] <- sum(df.CC$defc[(j-1):j])
    df.NW$prec.1mo[j] <- sum(df.NW$prec[(j-1):j])
    df.NW$tmin.1mo[j] <- mean(df.NW$tmin[(j-1):j])
    df.NW$tmax.1mo[j] <- mean(df.NW$tmax[(j-1):j])
    df.NW$relh.1mo[j] <- mean(df.NW$relh[(j-1):j])
    df.NW$rads.1mo[j] <- mean(df.NW$rads[(j-1):j])
    df.NW$wspd.1mo[j] <- mean(df.NW$wspd[(j-1):j])
    df.NW$ea.1mo[j]   <- mean(df.NW$ea[(j-1):j])
    df.NW$es.1mo[j]   <- mean(df.NW$es[(j-1):j])
    df.NW$VPD.1mo[j]  <- mean(df.NW$VPD[(j-1):j])
    df.NW$RET.1mo[j]  <- sum(df.NW$RET[(j-1):j])
    df.NW$defc.1mo[j] <- sum(df.NW$defc[(j-1):j])
  }
  if (j >= 3){
    df.AI$prec.2mo[j] <- sum(df.AI$prec[(j-2):j])
    df.AI$tmin.2mo[j] <- mean(df.AI$tmin[(j-2):j])
    df.AI$tmax.2mo[j] <- mean(df.AI$tmax[(j-2):j])
    df.AI$relh.2mo[j] <- mean(df.AI$relh[(j-2):j])
    df.AI$rads.2mo[j] <- mean(df.AI$rads[(j-2):j])
    df.AI$wspd.2mo[j] <- mean(df.AI$wspd[(j-2):j])
    df.AI$ea.2mo[j]   <- mean(df.AI$ea[(j-2):j])
    df.AI$es.2mo[j]   <- mean(df.AI$es[(j-2):j])
    df.AI$VPD.2mo[j]  <- mean(df.AI$VPD[(j-2):j])
    df.AI$RET.2mo[j]  <- sum(df.AI$RET[(j-2):j])
    df.AI$defc.2mo[j] <- sum(df.AI$defc[(j-2):j])
    df.AR$prec.2mo[j] <- sum(df.AR$prec[(j-2):j])
    df.AR$tmin.2mo[j] <- mean(df.AR$tmin[(j-2):j])
    df.AR$tmax.2mo[j] <- mean(df.AR$tmax[(j-2):j])
    df.AR$relh.2mo[j] <- mean(df.AR$relh[(j-2):j])
    df.AR$rads.2mo[j] <- mean(df.AR$rads[(j-2):j])
    df.AR$wspd.2mo[j] <- mean(df.AR$wspd[(j-2):j])
    df.AR$ea.2mo[j]   <- mean(df.AR$ea[(j-2):j])
    df.AR$es.2mo[j]   <- mean(df.AR$es[(j-2):j])
    df.AR$VPD.2mo[j]  <- mean(df.AR$VPD[(j-2):j])
    df.AR$RET.2mo[j]  <- sum(df.AR$RET[(j-2):j])
    df.AR$defc.2mo[j] <- sum(df.AR$defc[(j-2):j])
    df.CC$prec.2mo[j] <- sum(df.CC$prec[(j-2):j])
    df.CC$tmin.2mo[j] <- mean(df.CC$tmin[(j-2):j])
    df.CC$tmax.2mo[j] <- mean(df.CC$tmax[(j-2):j])
    df.CC$relh.2mo[j] <- mean(df.CC$relh[(j-2):j])
    df.CC$rads.2mo[j] <- mean(df.CC$rads[(j-2):j])
    df.CC$wspd.2mo[j] <- mean(df.CC$wspd[(j-2):j])
    df.CC$ea.2mo[j]   <- mean(df.CC$ea[(j-2):j])
    df.CC$es.2mo[j]   <- mean(df.CC$es[(j-2):j])
    df.CC$VPD.2mo[j]  <- mean(df.CC$VPD[(j-2):j])
    df.CC$RET.2mo[j]  <- sum(df.CC$RET[(j-2):j])
    df.CC$defc.2mo[j] <- sum(df.CC$defc[(j-2):j])
    df.NW$prec.2mo[j] <- sum(df.NW$prec[(j-2):j])
    df.NW$tmin.2mo[j] <- mean(df.NW$tmin[(j-2):j])
    df.NW$tmax.2mo[j] <- mean(df.NW$tmax[(j-2):j])
    df.NW$relh.2mo[j] <- mean(df.NW$relh[(j-2):j])
    df.NW$rads.2mo[j] <- mean(df.NW$rads[(j-2):j])
    df.NW$wspd.2mo[j] <- mean(df.NW$wspd[(j-2):j])
    df.NW$ea.2mo[j]   <- mean(df.NW$ea[(j-2):j])
    df.NW$es.2mo[j]   <- mean(df.NW$es[(j-2):j])
    df.NW$VPD.2mo[j]  <- mean(df.NW$VPD[(j-2):j])
    df.NW$RET.2mo[j]  <- sum(df.NW$RET[(j-2):j])
    df.NW$defc.2mo[j] <- sum(df.NW$defc[(j-2):j])
  }
  if (j >= 4){
    df.AI$prec.3mo[j] <- sum(df.AI$prec[(j-3):j])
    df.AI$tmin.3mo[j] <- mean(df.AI$tmin[(j-3):j])
    df.AI$tmax.3mo[j] <- mean(df.AI$tmax[(j-3):j])
    df.AI$relh.3mo[j] <- mean(df.AI$relh[(j-3):j])
    df.AI$rads.3mo[j] <- mean(df.AI$rads[(j-3):j])
    df.AI$wspd.3mo[j] <- mean(df.AI$wspd[(j-3):j])
    df.AI$ea.3mo[j]   <- mean(df.AI$ea[(j-3):j])
    df.AI$es.3mo[j]   <- mean(df.AI$es[(j-3):j])
    df.AI$VPD.3mo[j]  <- mean(df.AI$VPD[(j-3):j])
    df.AI$RET.3mo[j]  <- sum(df.AI$RET[(j-3):j])
    df.AI$defc.3mo[j] <- sum(df.AI$defc[(j-3):j])
    df.AR$prec.3mo[j] <- sum(df.AR$prec[(j-3):j])
    df.AR$tmin.3mo[j] <- mean(df.AR$tmin[(j-3):j])
    df.AR$tmax.3mo[j] <- mean(df.AR$tmax[(j-3):j])
    df.AR$relh.3mo[j] <- mean(df.AR$relh[(j-3):j])
    df.AR$rads.3mo[j] <- mean(df.AR$rads[(j-3):j])
    df.AR$wspd.3mo[j] <- mean(df.AR$wspd[(j-3):j])
    df.AR$ea.3mo[j]   <- mean(df.AR$ea[(j-3):j])
    df.AR$es.3mo[j]   <- mean(df.AR$es[(j-3):j])
    df.AR$VPD.3mo[j]  <- mean(df.AR$VPD[(j-3):j])
    df.AR$RET.3mo[j]  <- sum(df.AR$RET[(j-3):j])
    df.AR$defc.3mo[j] <- sum(df.AR$defc[(j-3):j])
    df.CC$prec.3mo[j] <- sum(df.CC$prec[(j-3):j])
    df.CC$tmin.3mo[j] <- mean(df.CC$tmin[(j-3):j])
    df.CC$tmax.3mo[j] <- mean(df.CC$tmax[(j-3):j])
    df.CC$relh.3mo[j] <- mean(df.CC$relh[(j-3):j])
    df.CC$rads.3mo[j] <- mean(df.CC$rads[(j-3):j])
    df.CC$wspd.3mo[j] <- mean(df.CC$wspd[(j-3):j])
    df.CC$ea.3mo[j]   <- mean(df.CC$ea[(j-3):j])
    df.CC$es.3mo[j]   <- mean(df.CC$es[(j-3):j])
    df.CC$VPD.3mo[j]  <- mean(df.CC$VPD[(j-3):j])
    df.CC$RET.3mo[j]  <- sum(df.CC$RET[(j-3):j])
    df.CC$defc.3mo[j] <- sum(df.CC$defc[(j-3):j])
    df.NW$prec.3mo[j] <- sum(df.NW$prec[(j-3):j])
    df.NW$tmin.3mo[j] <- mean(df.NW$tmin[(j-3):j])
    df.NW$tmax.3mo[j] <- mean(df.NW$tmax[(j-3):j])
    df.NW$relh.3mo[j] <- mean(df.NW$relh[(j-3):j])
    df.NW$rads.3mo[j] <- mean(df.NW$rads[(j-3):j])
    df.NW$wspd.3mo[j] <- mean(df.NW$wspd[(j-3):j])
    df.NW$ea.3mo[j]   <- mean(df.NW$ea[(j-3):j])
    df.NW$es.3mo[j]   <- mean(df.NW$es[(j-3):j])
    df.NW$VPD.3mo[j]  <- mean(df.NW$VPD[(j-3):j])
    df.NW$RET.3mo[j]  <- sum(df.NW$RET[(j-3):j])
    df.NW$defc.3mo[j] <- sum(df.NW$defc[(j-3):j])
  }
  if (j >= 7){
    df.AI$prec.6mo[j] <- sum(df.AI$prec[(j-6):j])
    df.AI$tmin.6mo[j] <- mean(df.AI$tmin[(j-6):j])
    df.AI$tmax.6mo[j] <- mean(df.AI$tmax[(j-6):j])
    df.AI$relh.6mo[j] <- mean(df.AI$relh[(j-6):j])
    df.AI$rads.6mo[j] <- mean(df.AI$rads[(j-6):j])
    df.AI$wspd.6mo[j] <- mean(df.AI$wspd[(j-6):j])
    df.AI$ea.6mo[j]   <- mean(df.AI$ea[(j-6):j])
    df.AI$es.6mo[j]   <- mean(df.AI$es[(j-6):j])
    df.AI$VPD.6mo[j]  <- mean(df.AI$VPD[(j-6):j])
    df.AI$RET.6mo[j]  <- sum(df.AI$RET[(j-6):j])
    df.AI$defc.6mo[j] <- sum(df.AI$defc[(j-6):j])
    df.AR$prec.6mo[j] <- sum(df.AR$prec[(j-6):j])
    df.AR$tmin.6mo[j] <- mean(df.AR$tmin[(j-6):j])
    df.AR$tmax.6mo[j] <- mean(df.AR$tmax[(j-6):j])
    df.AR$relh.6mo[j] <- mean(df.AR$relh[(j-6):j])
    df.AR$rads.6mo[j] <- mean(df.AR$rads[(j-6):j])
    df.AR$wspd.6mo[j] <- mean(df.AR$wspd[(j-6):j])
    df.AR$ea.6mo[j]   <- mean(df.AR$ea[(j-6):j])
    df.AR$es.6mo[j]   <- mean(df.AR$es[(j-6):j])
    df.AR$VPD.6mo[j]  <- mean(df.AR$VPD[(j-6):j])
    df.AR$RET.6mo[j]  <- sum(df.AR$RET[(j-6):j])
    df.AR$defc.6mo[j] <- sum(df.AR$defc[(j-6):j])
    df.CC$prec.6mo[j] <- sum(df.CC$prec[(j-6):j])
    df.CC$tmin.6mo[j] <- mean(df.CC$tmin[(j-6):j])
    df.CC$tmax.6mo[j] <- mean(df.CC$tmax[(j-6):j])
    df.CC$relh.6mo[j] <- mean(df.CC$relh[(j-6):j])
    df.CC$rads.6mo[j] <- mean(df.CC$rads[(j-6):j])
    df.CC$wspd.6mo[j] <- mean(df.CC$wspd[(j-6):j])
    df.CC$ea.6mo[j]   <- mean(df.CC$ea[(j-6):j])
    df.CC$es.6mo[j]   <- mean(df.CC$es[(j-6):j])
    df.CC$VPD.6mo[j]  <- mean(df.CC$VPD[(j-6):j])
    df.CC$RET.6mo[j]  <- sum(df.CC$RET[(j-6):j])
    df.CC$defc.6mo[j] <- sum(df.CC$defc[(j-6):j])
    df.NW$prec.6mo[j] <- sum(df.NW$prec[(j-6):j])
    df.NW$tmin.6mo[j] <- mean(df.NW$tmin[(j-6):j])
    df.NW$tmax.6mo[j] <- mean(df.NW$tmax[(j-6):j])
    df.NW$relh.6mo[j] <- mean(df.NW$relh[(j-6):j])
    df.NW$rads.6mo[j] <- mean(df.NW$rads[(j-6):j])
    df.NW$wspd.6mo[j] <- mean(df.NW$wspd[(j-6):j])
    df.NW$ea.6mo[j]   <- mean(df.NW$ea[(j-6):j])
    df.NW$es.6mo[j]   <- mean(df.NW$es[(j-6):j])
    df.NW$VPD.6mo[j]  <- mean(df.NW$VPD[(j-6):j])
    df.NW$RET.6mo[j]  <- sum(df.NW$RET[(j-6):j])
    df.NW$defc.6mo[j] <- sum(df.NW$defc[(j-6):j])
  }
  if (j >= 13){
    df.AI$prec.12mo[j] <- sum(df.AI$prec[(j-12):j])
    df.AI$tmin.12mo[j] <- mean(df.AI$tmin[(j-12):j])
    df.AI$tmax.12mo[j] <- mean(df.AI$tmax[(j-12):j])
    df.AI$relh.12mo[j] <- mean(df.AI$relh[(j-12):j])
    df.AI$rads.12mo[j] <- mean(df.AI$rads[(j-12):j])
    df.AI$wspd.12mo[j] <- mean(df.AI$wspd[(j-12):j])
    df.AI$ea.12mo[j]   <- mean(df.AI$ea[(j-12):j])
    df.AI$es.12mo[j]   <- mean(df.AI$es[(j-12):j])
    df.AI$VPD.12mo[j]  <- mean(df.AI$VPD[(j-12):j])
    df.AI$RET.12mo[j]  <- sum(df.AI$RET[(j-12):j])
    df.AI$defc.12mo[j] <- sum(df.AI$defc[(j-12):j])
    df.AR$prec.12mo[j] <- sum(df.AR$prec[(j-12):j])
    df.AR$tmin.12mo[j] <- mean(df.AR$tmin[(j-12):j])
    df.AR$tmax.12mo[j] <- mean(df.AR$tmax[(j-12):j])
    df.AR$relh.12mo[j] <- mean(df.AR$relh[(j-12):j])
    df.AR$rads.12mo[j] <- mean(df.AR$rads[(j-12):j])
    df.AR$wspd.12mo[j] <- mean(df.AR$wspd[(j-12):j])
    df.AR$ea.12mo[j]   <- mean(df.AR$ea[(j-12):j])
    df.AR$es.12mo[j]   <- mean(df.AR$es[(j-12):j])
    df.AR$VPD.12mo[j]  <- mean(df.AR$VPD[(j-12):j])
    df.AR$RET.12mo[j]  <- sum(df.AR$RET[(j-12):j])
    df.AR$defc.12mo[j] <- sum(df.AR$defc[(j-12):j])
    df.CC$prec.12mo[j] <- sum(df.CC$prec[(j-12):j])
    df.CC$tmin.12mo[j] <- mean(df.CC$tmin[(j-12):j])
    df.CC$tmax.12mo[j] <- mean(df.CC$tmax[(j-12):j])
    df.CC$relh.12mo[j] <- mean(df.CC$relh[(j-12):j])
    df.CC$rads.12mo[j] <- mean(df.CC$rads[(j-12):j])
    df.CC$wspd.12mo[j] <- mean(df.CC$wspd[(j-12):j])
    df.CC$ea.12mo[j]   <- mean(df.CC$ea[(j-12):j])
    df.CC$es.12mo[j]   <- mean(df.CC$es[(j-12):j])
    df.CC$VPD.12mo[j]  <- mean(df.CC$VPD[(j-12):j])
    df.CC$RET.12mo[j]  <- sum(df.CC$RET[(j-12):j])
    df.CC$defc.12mo[j] <- sum(df.CC$defc[(j-12):j])
    df.NW$prec.12mo[j] <- sum(df.NW$prec[(j-12):j])
    df.NW$tmin.12mo[j] <- mean(df.NW$tmin[(j-12):j])
    df.NW$tmax.12mo[j] <- mean(df.NW$tmax[(j-12):j])
    df.NW$relh.12mo[j] <- mean(df.NW$relh[(j-12):j])
    df.NW$rads.12mo[j] <- mean(df.NW$rads[(j-12):j])
    df.NW$wspd.12mo[j] <- mean(df.NW$wspd[(j-12):j])
    df.NW$ea.12mo[j]   <- mean(df.NW$ea[(j-12):j])
    df.NW$es.12mo[j]   <- mean(df.NW$es[(j-12):j])
    df.NW$VPD.12mo[j]  <- mean(df.NW$VPD[(j-12):j])
    df.NW$RET.12mo[j]  <- sum(df.NW$RET[(j-12):j])
    df.NW$defc.12mo[j] <- sum(df.NW$defc[(j-12):j])
  }
}

# some other candidate variables
df.AI$prec.sq <- df.AI$prec^2
df.AI$defc.sq <- df.AI$defc^2
df.AI$defc.1mo.sq <- df.AI$defc.1mo^2
df.AI$defc.2mo.sq <- df.AI$defc.2mo^2
df.AI$defc.3mo.sq <- df.AI$defc.3mo^2
df.AI$defc.6mo.sq <- df.AI$defc.6mo^2
df.AI$defc.12mo.sq <- df.AI$defc.12mo^2
df.AR$prec.sq <- df.AR$prec^2
df.AR$defc.sq <- df.AR$defc^2
df.AR$defc.1mo.sq <- df.AR$defc.1mo^2
df.AR$defc.2mo.sq <- df.AR$defc.2mo^2
df.AR$defc.3mo.sq <- df.AR$defc.3mo^2
df.AR$defc.6mo.sq <- df.AR$defc.6mo^2
df.AR$defc.12mo.sq <- df.AR$defc.12mo^2
df.CC$prec.sq <- df.CC$prec^2
df.CC$defc.sq <- df.CC$defc^2
df.CC$defc.1mo.sq <- df.CC$defc.1mo^2
df.CC$defc.2mo.sq <- df.CC$defc.2mo^2
df.CC$defc.3mo.sq <- df.CC$defc.3mo^2
df.CC$defc.6mo.sq <- df.CC$defc.6mo^2
df.CC$defc.12mo.sq <- df.CC$defc.12mo^2
df.NW$prec.sq <- df.NW$prec^2
df.NW$defc.sq <- df.NW$defc^2
df.NW$defc.1mo.sq <- df.NW$defc.1mo^2
df.NW$defc.2mo.sq <- df.NW$defc.2mo^2
df.NW$defc.3mo.sq <- df.NW$defc.3mo^2
df.NW$defc.6mo.sq <- df.NW$defc.6mo^2
df.NW$defc.12mo.sq <- df.NW$defc.12mo^2

# get rid of baseline period to avoid repeating;
# it will be loaded in with the actual scenarios
# and is the same for all
df.AI <- subset(df.AI, year>yr.baseline.end)
df.AR <- subset(df.AR, year>yr.baseline.end)
df.CC <- subset(df.CC, year>yr.baseline.end)
df.NW <- subset(df.NW, year>yr.baseline.end)

# bind into single data frame
df.hist <- rbind(df.AI, df.AR, df.CC, df.NW)

# all of the "future" years (2014-2070) should be set to the past so they are used to
# build MLR relationships; set equal to yr.baseline.start-1
df.hist$year <- yr.baseline.start-1

# get list of all possible predictor variables
var.options <- colnames(df.hist)[!(colnames(df.hist) %in% c("year","month","aet","srunoff","drainage", "ndaypm", "climate", "prec.gt.101"))]

## data frame to store overall output
df.summary <- data.frame(var = character(0),
                         month = numeric(0),
                         R2.cal = numeric(0),
                         R2.val = numeric(0),
                         RMSE.cal = numeric(0),
                         RMSE.val = numeric(0),
                         NRMSE.cal = numeric(0),
                         NRMSE.val = numeric(0))

## 1. Scroll through each LULC and climate scenario
for (LULC in scen.LULC.all){
  for (climate in scen.climate.all){
    
    ## 2. Load in data.
    df <- read.csv(paste0("input/", LULC, climate, "_waterbalance_monthly.csv"), stringsAsFactors=F)
    
    # convert daily mean values to cumulative monthly totals for depth variables
    df$ndaypm <- days_in_month(ymd(paste0(df$year, "-", df$month, "-01")))
    df$aet <- df$aet*df$ndaypm
    df$srunoff <- df$srunoff*df$ndaypm
    df$drainage <- df$drainage*df$ndaypm
    
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
    
    # if you are using HistLULC_XXclim, bind data frames together
    n.baseline.yrs <- length(seq(yr.baseline.start, yr.baseline.end))
    
    if (hist.LULC){
      df <- rbind(df.hist, df)
      yr.baseline.start <- yr.baseline.start-1
      n.baseline.yrs <- n.baseline.yrs + length(seq(yr.baseline.end+1,2070))*4
    }
    
    # number of validation years
    n.val.yr <- round(n.baseline.yrs*0.25)
    
    # scroll through fluxes
    for (flux.name in flux.name.all){
      for (mo in mo.all){
        if (flux.name %in% c("aet", "drainage")){
          neg.allowed=T
        } else {
          neg.allowed=F
        }
        
        # figure out number of predictors to retain, based on PC.keep
        PCs.keep <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/", flux.name, "_GHCN_MonthlyRegressions_PC.keep_", sprintf("%02d", mo), ".csv"))
        n.vars.keep <- length(PCs.keep$PC)
        
        # run regressions
        df.out.mo <- CalculateClimateMLR(df=df, mo=mo, flux.name=flux.name, var.options=var.options,
                                         yr.baseline.start=yr.baseline.start, yr.baseline.end=yr.baseline.end, n.val.yr=n.val.yr,
                                         p.thres=p.thres, neg.allowed=neg.allowed, n.vars.keep=n.vars.keep, write.vars.keep=T, write.perm=F)
        
        # add flux name
        df.out.mo$flux.name <- flux.name
        
        # combine with other data
        if (exists("df.out")){
          df.out <- rbind(df.out, df.out.mo)
        } else {
          df.out <- df.out.mo
        }
        
        # status update
        print(paste0(LULC, climate, " ", flux.name, " month ", mo, " complete"))
        
      }  # end of mo loop
      
    }  # end of flux.name loop
    
    # divide into cal/val and prediction periods
    df.out.hist <- subset(df.out, year<=yr.baseline.end)
    df.out.pred <- subset(df.out, year>yr.baseline.end)
    
    # save output csv
    if (hist.LULC){
      write.csv(df.out.pred, paste0("HistLULC_AllClim_", LULC, climate, "_MonthlyRegressions_MLR_OutputAll.csv"), row.names=F, quote=F)
      if (hist.csv){
        write.csv(df.out.hist, "HistLULC_AllClim_CalVal_MonthlyRegressions_MLR_OutputAll.csv", row.names=F, quote=F)
        hist.csv <- F   # disable
        write.vars.keep <- F
        write.PC.keep <- F
      }
    } else {
      write.csv(df.out.pred, paste0(LULC, climate, "_MonthlyRegressions_MLR_OutputAll.csv"), row.names=F, quote=F)
      if (hist.csv){
        write.csv(df.out.hist, "CalVal_MonthlyRegressions_MLR_OutputAll.csv", row.names=F, quote=F)
        hist.csv <- F   # disable
        write.vars.keep <- F
        write.PC.keep <- F
      }
    }
    
    rm(df.out)
    
  }  # end of climate loop
}  # end of LULC loop