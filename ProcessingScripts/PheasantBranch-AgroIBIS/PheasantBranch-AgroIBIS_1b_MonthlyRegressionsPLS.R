## PheasantBranch-AgroIBIS_1_MonthlyRegressions.R
# This script is intended to load monthly model output data (generated using the script
# Scenarios_SummarizeMonthlyWaterBalance.R) and calculate a separate statistical relationship 
# between meteorological data and a water balance output variable for each month using 
# principal components regression. A permutation-based approach is used to provide independent
# estimates of model fit and uncertainty associated with estimates.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

# load packages
source(paste0(git.dir, "ProcessingScripts/PCA_Functions.R"))
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))
source(paste0(git.dir, "ProcessingScripts/CalculateClimatePLS.R"))   # this is the script that does most of the work...
require(ggplot2)
require(tidyr)
require(WVPlots)  # devtools::install_github('WinVector/WVPlots',build_vignettes=TRUE)
require(gridExtra)
require(lubridate)

setwd(paste0(git.dir, "Data/PheasantBranch-AgroIBIS/"))

# directory to save plots
plot.dir <- paste0(git.dir, "Plots/PheasantBranch-AgroIBIS/")

# set significance threshold for pruning
p.thres <- 0.10

# variance for selecting options
cum.var <- 0.8

## decide "baseline" period (inclusive)
# 1986-2013 = 28 years, 1986-1999 = 14 years
# 1974-2016 = 43 years, 1974-1995 = 22 years
yr.baseline.start <- 1974
yr.baseline.end <- 1995
n.baseline.yrs <- length(seq(yr.baseline.start, yr.baseline.end))

# number of permutations
n.perm <- 250

# which fluxes to analyze?
flux.name.all <- c("aet", "srunoff", "drainage")

# which months to analyze?
mo.all <- seq(1,12)

# which scenario to analyze? only need one, because looking at historical data
LULC <- "AI"
climate <- "AI"

# which output to write?
write.vars.keep <- T
write.PC.keep <- T
write.perm <- F
PC.plot <- F

## data frame to store overall output
df.summary <- data.frame(var = character(0),
                         month = numeric(0),
                         R2.cal = numeric(0),
                         R2.val = numeric(0),
                         RMSE.cal = numeric(0),
                         RMSE.val = numeric(0),
                         NRMSE.cal = numeric(0),
                         NRMSE.val = numeric(0))

## 2. Load in data.
df <- read.csv(paste0("PheasantBranch_", LULC, climate, "_waterbalance_monthly.csv"), stringsAsFactors=F)

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

# get list of all possible predictor variables
var.options <- colnames(df)[!(colnames(df) %in% c("year", "year.met", "month","aet","srunoff","drainage", "ndaypm", "climate", "prec.gt.101"))]

# number of validation years
n.val.yr <- floor(n.baseline.yrs*0.25)

# scroll through fluxes
for (flux.name in flux.name.all){
  for (mo in mo.all){
    if (flux.name %in% c("aet", "drainage")){
      neg.allowed=T
    } else {
      neg.allowed=F
    }
    
    # run PLS script
    df.out.mo <- CalculateClimatePLS(df=df, mo=mo, flux.name=flux.name, var.options=var.options,
                                     yr.baseline.start=yr.baseline.start, yr.baseline.end=yr.baseline.end, n.val.yr=n.val.yr,
                                     p.thres=p.thres, neg.allowed=neg.allowed, min.PCs.keep=1, write.vars.keep=F, write.perm=F)
    
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
  
  # rename output files from CalculateClimatePLS script
  files <- list.files(paste0(git.dir, "Data/PheasantBranch-AgroIBIS/"), pattern="GHCN_MonthlyRegressions_.")
  file.rename(paste0(git.dir, "Data/PheasantBranch-AgroIBIS/", files), paste0(git.dir, "Data/PheasantBranch-AgroIBIS/PLS_Output/", flux.name, "_", files))
  
  
}  # end of flux.name loop

# save output csv
write.csv(df.out, paste0(LULC, climate, "_MonthlyRegressions_PLS_OutputAll.csv"), row.names=F, quote=F)
