## Yahara_2a_ClimateVsLULC-PCR.R
# This script is intended to take the output of Yahara_1_MonthlyRegressions.R and
# separate out the climate and LULC contributions to change since the baseline period
# at an annual scale.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(ggplot2)
require(dplyr)
require(lubridate)
require(reshape2)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))
# moving average function
avg <- function(x, avg.yr){
  as.numeric(stats::filter(x, filter=rep(1/avg.yr, avg.yr), sides=1))
}

# path to folder with output data
setwd(paste0(git.dir, "Data/Yahara-AgroIBIS/"))

# baseline period
yr.baseline.start <- 1986
yr.baseline.end <- 2013

# which fluxes to analyze?
flux.name.all <- c("aet", "srunoff", "drainage")

# which scenarios to analyze?
scen.climate.all <- c("AI", "AR", "CC", "NW")
scen.LULC.all <- c("AI", "AR", "CC", "NW")

# years for moving average
avg.yr <- 11

# use HistLULC_XXclim simulations?
hist.LULC <- T   # T or F

# read in cal/val to fill in the baseline values
if (hist.LULC){
  df.baseline <- read.csv("HistLULC_AllClim_CalVal_MonthlyRegressions_OutputAll.csv", stringsAsFactors=F)
} else {
  df.baseline <- read.csv("CalVal_MonthlyRegressions_OutputAll.csv", stringsAsFactors=F)
}
df.baseline <- subset(df.baseline, year >= yr.baseline.start & year<=yr.baseline.end)     # subset to get rid of HistLULC_XXclim results (which all have year 1985)
df.baseline.mo <- dplyr::summarize(group_by(df.baseline, month, flux.name),
                                   flux.baseline.mean = mean(flux),
                                   PCR.baseline.mean= mean(PCR))

for (LULC in scen.LULC.all){
  for (climate in scen.climate.all){
    
    # read in data
    if (hist.LULC){
      df <- read.csv(paste0("HistLULC_AllClim_", LULC, climate, "_MonthlyRegressions_OutputAll.csv"), stringsAsFactors=F)
    } else {
      df <- read.csv(paste0(LULC, climate, "_MonthlyRegressions_OutputAll.csv"), stringsAsFactors=F)
    }
    
    # join with baseline
    df <- rbind(df.baseline, df)
    df <- left_join(df, df.baseline.mo, by=c("month", "flux.name"))
    
    # date column
    df$date <- ymd(paste0(df$year, "-", df$month, "-15"))
    
    # sum to annual
    df.perm.ann <- dplyr::summarize(group_by(df, flux.name, year, perm),
                                    flux.mm = sum(flux),
                                    PCR.mm = sum(PCR),
                                    flux.baseline.mm = sum(flux.baseline.mean),
                                    PCR.baseline.mm = sum(PCR.baseline.mean))
    
    # calculate overall change, climate, and LULC contribution
    df.perm.ann$change.overall <- df.perm.ann$flux.mm - df.perm.ann$flux.baseline.mm
    df.perm.ann$change.climate <- df.perm.ann$PCR.mm - df.perm.ann$flux.baseline.mm
    df.perm.ann$change.LULC <- df.perm.ann$change.overall - df.perm.ann$change.climate
    
    # calculate mean and standard deviation
    df.ann <- dplyr::summarize(group_by(df.perm.ann, flux.name, year),
                               LULC = LULC,
                               climate = climate,
                               obs.baseline = mean(flux.baseline.mm),
                               obs = mean(flux.mm),
                               PCR.mean = mean(PCR.mm),
                               PCR.sd = sd(PCR.mm),
                               change.overall.mean = mean(change.overall),
                               change.overall.sd = sd(change.overall),
                               change.climate.mean = mean(change.climate),
                               change.climate.sd = sd(change.climate),
                               change.LULC.mean = mean(change.LULC),
                               change.LULC.sd = sd(change.LULC))
    
    # calculate some min/max columns to use for ribbons
    df.ann$change.overall.min <- df.ann$change.overall.mean - df.ann$change.overall.sd
    df.ann$change.overall.max <- df.ann$change.overall.mean + df.ann$change.overall.sd
    df.ann$change.climate.min <- df.ann$change.climate.mean - df.ann$change.climate.sd
    df.ann$change.climate.max <- df.ann$change.climate.mean + df.ann$change.climate.sd
    df.ann$change.LULC.min <- df.ann$change.LULC.mean - df.ann$change.LULC.sd
    df.ann$change.LULC.max <- df.ann$change.LULC.mean + df.ann$change.LULC.sd
    
    # moving average, separately for each flux
    for (flux in unique(df.ann$flux.name)){
      df.flux <- subset(df.ann, flux.name==flux)
      df.flux.avg <- df.flux
      df.flux.avg <- df.flux.avg[order(df.flux.avg$year), ]
      df.flux.avg[, !(colnames(df.flux.avg) %in% c("flux.name", "year", "LULC", "climate", "obs.baseline"))] <- 
        apply(df.flux.avg[, !(colnames(df.flux.avg) %in% c("flux.name", "year", "LULC", "climate", "obs.baseline"))], 2, avg, avg.yr=avg.yr)
      colnames(df.flux.avg)[!(colnames(df.flux.avg) %in% c("flux.name", "year", "LULC", "climate", "obs.baseline"))] <- paste0("avg.", colnames(df.flux.avg)[!(colnames(df.flux.avg) %in% c("flux.name", "year", "LULC", "climate", "obs.baseline"))])
      
      if (exists("df.ann.avg")){
        df.ann.avg <- rbind(df.ann.avg, df.flux.avg)
      } else {
        df.ann.avg <- df.flux.avg
      }
    }
    
    # add averaged to overall data frame
    df.ann <- left_join(df.ann, df.ann.avg, by = c("flux.name", "year", "LULC", "climate", "obs.baseline"))
    
    rm(df.ann.avg)
    
    # add to overall data frame
    if (exists("df.out")){
      df.out <- rbind(df.out, df.ann)
    } else {
      df.out <- df.ann
    }
    
    # status update
    print(paste0(LULC, climate, " complete"))
  }
}

# calculate percent change
df.out$change.overall.mean.prc <- 100*df.out$change.overall.mean/df.out$obs.baseline
df.out$change.overall.sd.prc <- 100*df.out$change.overall.sd/df.out$obs.baseline
df.out$change.overall.min.prc <- 100*df.out$change.overall.min/df.out$obs.baseline
df.out$change.overall.max.prc <- 100*df.out$change.overall.max/df.out$obs.baseline

df.out$avg.change.overall.mean.prc <- 100*df.out$avg.change.overall.mean/df.out$obs.baseline
df.out$avg.change.overall.sd.prc <- 100*df.out$avg.change.overall.sd/df.out$obs.baseline
df.out$avg.change.overall.min.prc <- 100*df.out$avg.change.overall.min/df.out$obs.baseline
df.out$avg.change.overall.max.prc <- 100*df.out$avg.change.overall.max/df.out$obs.baseline

df.out$change.LULC.mean.prc <- 100*df.out$change.LULC.mean/df.out$obs.baseline
df.out$change.LULC.sd.prc <- 100*df.out$change.LULC.sd/df.out$obs.baseline
df.out$change.LULC.min.prc <- 100*df.out$change.LULC.min/df.out$obs.baseline
df.out$change.LULC.max.prc <- 100*df.out$change.LULC.max/df.out$obs.baseline

df.out$avg.change.LULC.mean.prc <- 100*df.out$avg.change.LULC.mean/df.out$obs.baseline
df.out$avg.change.LULC.sd.prc <- 100*df.out$avg.change.LULC.sd/df.out$obs.baseline
df.out$avg.change.LULC.min.prc <- 100*df.out$avg.change.LULC.min/df.out$obs.baseline
df.out$avg.change.LULC.max.prc <- 100*df.out$avg.change.LULC.max/df.out$obs.baseline

df.out$change.climate.mean.prc <- 100*df.out$change.climate.mean/df.out$obs.baseline
df.out$change.climate.sd.prc <- 100*df.out$change.climate.sd/df.out$obs.baseline
df.out$change.climate.min.prc <- 100*df.out$change.climate.min/df.out$obs.baseline
df.out$change.climate.max.prc <- 100*df.out$change.climate.max/df.out$obs.baseline

df.out$avg.change.climate.mean.prc <- 100*df.out$avg.change.climate.mean/df.out$obs.baseline
df.out$avg.change.climate.sd.prc <- 100*df.out$avg.change.climate.sd/df.out$obs.baseline
df.out$avg.change.climate.min.prc <- 100*df.out$avg.change.climate.min/df.out$obs.baseline
df.out$avg.change.climate.max.prc <- 100*df.out$avg.change.climate.max/df.out$obs.baseline

# write output
if (hist.LULC){
  write.csv(df.out, "HistLULC_AllClim_ClimateVsLULC.csv", row.names=F, quote=F)
} else {
  write.csv(df.out, "ClimateVsLULC.csv", row.names=F, quote=F)
}
