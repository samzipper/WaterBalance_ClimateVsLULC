## PheasantBranch_3_SensitivityBaseline.R
# This script is intended to quantify the sensitivity of the results for the Pheasant Branch watershed
# to the selection of the baseline period.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

# load packages
source(paste0(git.dir, "ProcessingScripts/PCA_Functions.R"))
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))
source(paste0(git.dir, "ProcessingScripts/CalculateClimatePCR.R"))
require(ggplot2)
require(dplyr)
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

# directory to save plots
plot.dir <- paste0(git.dir, "Plots/PheasantBranch/")

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
yr.baseline.end.all <- seq(1988,1996)  # 9 years: 1992 +/- 4

# number of permutations to run
n.perm <- 250

# read in discharge and met data frames
df.met <- read.csv(paste0(git.dir, "Data/PheasantBranch/USW00014837_GHCN_MetData_Monthly.csv"))  # use same met data as McFarland
df.Q <- read.csv("PheasantBranch_Discharge_Monthly.csv")

# merge data frames, with df.met defining the temporal extent
df <- merge(df.met, df.Q, all.x=T)

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
var.options <- colnames(df)[!(colnames(df) %in% c("year","month","discharge.mm", "discharge.est", "prec.gt.101"))]

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
  for (mo in mo.list){
    
    # run PCR script
    df.out.mo <- CalculateClimatePCR(df=df, mo=mo, flux.name="discharge.mm", var.predictors=var.options,
                                     yr.baseline.start=yr.baseline.start, yr.baseline.end=yr.baseline.end,
                                     n.perm=n.perm, n.val.yr=5, p.thres=p.thres, cum.var=0.95, min.var=0.005, 
                                     neg.allowed=F, write.vars.keep=F, write.PC.keep=F, PC.plot=F, write.perm=F)
    df.out.mo$yr.baseline.start <- yr.baseline.start
    df.out.mo$yr.baseline.end <- yr.baseline.end
    
    # combined with other data
    if (exists("df.out")){
      df.out <- rbind(df.out, df.out.mo)
    } else {
      df.out <- df.out.mo
    }
    
    # status update
    print(paste0("baseline end ", yr.baseline.end, " month ", mo, " complete"))
    
  }
}

# save output csv
write.csv(df.out, "GHCN_SensitivityBaseline_OutputAll.csv", row.names=F, quote=F)

## calculate fit metrics by yr.baseline.end
df.baseline <- dplyr::summarize(group_by(subset(df.out, group != "prediction"), yr.baseline.end, group),
                                RMSE = RMSE(PCR, flux),
                                NSE = NashSutcliffe(PCR, flux))
df.baseline.melt <- melt(df.baseline, id=c("group", "yr.baseline.end"))

## calculate baseline period discharge
df.year.perm.baseline <-
  dplyr::summarize(group_by(df.out, year, perm, yr.baseline.end),
                   obs.sum = sum(flux),
                   PCR.sum = sum(PCR))
df.year.baseline <-
  dplyr::summarize(group_by(df.year.perm.baseline, year, yr.baseline.end),
                   obs.mean = mean(obs.sum),
                   PCR.mean = mean(PCR.sum))

df.Q.baseline <-
  dplyr::summarize(group_by(df.year.baseline, yr.baseline.end),
                   obs.baseline.mean = mean(obs.mean[year<=yr.baseline.end]))

## separate LULC vs climate
df.LULCvClimate <- merge(df.year.perm.baseline, df.Q.baseline, by=c("yr.baseline.end"), all.x=T)
df.LULCvClimate$change.overall <- df.LULCvClimate$obs.sum - df.LULCvClimate$obs.baseline.mean
df.LULCvClimate$change.climate <- df.LULCvClimate$PCR.sum - df.LULCvClimate$obs.baseline.mean
df.LULCvClimate$change.LULC <- df.LULCvClimate$change.overall - df.LULCvClimate$change.climate

df.LULCvClimate.baseline <-
  dplyr::summarize(group_by(df.LULCvClimate, year, yr.baseline.end),
                   obs = mean(obs.sum),
                   PCR.mean = mean(PCR.sum),
                   PCR.sd = sd(PCR.sum),
                   change.overall.mean = mean(change.overall),
                   change.overall.sd = sd(change.overall),
                   change.climate.mean = mean(change.climate),
                   change.climate.sd = sd(change.climate),
                   change.LULC.mean = mean(change.LULC),
                   change.LULC.sd = sd(change.LULC))

df.LULCvClimate.prediction <- subset(df.LULCvClimate.baseline, year>max(df.LULCvClimate$yr.baseline.end))

# calculate some min/max columns to use for ribbons
df.LULCvClimate.baseline$change.overall.min <- df.LULCvClimate.baseline$change.overall.mean - df.LULCvClimate.baseline$change.overall.sd
df.LULCvClimate.baseline$change.overall.max <- df.LULCvClimate.baseline$change.overall.mean + df.LULCvClimate.baseline$change.overall.sd
df.LULCvClimate.baseline$change.climate.min <- df.LULCvClimate.baseline$change.climate.mean - df.LULCvClimate.baseline$change.climate.sd
df.LULCvClimate.baseline$change.climate.max <- df.LULCvClimate.baseline$change.climate.mean + df.LULCvClimate.baseline$change.climate.sd
df.LULCvClimate.baseline$change.LULC.min <- df.LULCvClimate.baseline$change.LULC.mean - df.LULCvClimate.baseline$change.LULC.sd
df.LULCvClimate.baseline$change.LULC.max <- df.LULCvClimate.baseline$change.LULC.mean + df.LULCvClimate.baseline$change.LULC.sd

## make plots
# barplot of fit metrics by group and yr.baseline.end
p.fit.group.baseline <-
  ggplot(df.baseline.melt, aes(x=yr.baseline.end, y=value, fill=group)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity", position="dodge") +
  facet_grid(variable~., scales="free_y") +
  scale_x_continuous(name="End of Baseline Period [year]", 
                     breaks=seq(min(df.baseline$yr.baseline.end), max(df.baseline$yr.baseline.end))) +
  scale_y_continuous(name="PCR Fit") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(plot.dir, "GHCN_PheasantBranch_SensitivityBaseline_p.fit.group.baseline.png"),
       p.fit.group.baseline, width=8, height=6, units="in")

# line plot of observed vs PCR facet by yr.baseline.end
p.facet.obs.PCR <-
  ggplot(melt(df.year.baseline, id=c("year", "yr.baseline.end")), aes(x=year, y=value, color=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() +
  facet_wrap(~yr.baseline.end) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Annual Discharge [mm]") +
  scale_color_discrete(name="Data", labels=c("obs.mean"="Observed", "PCR.mean"="PCR")) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")
ggsave(paste0(plot.dir, "GHCN_PheasantBranch_SensitivityBaseline_p.facet.obs.PCR.png"),
       p.facet.obs.PCR, width=8, height=6, units="in")

# line plot of baseline period discharge as a function of yr.baseline.end
p.baseline.discharge <-
  ggplot(df.Q.baseline, aes(x=yr.baseline.end, y=obs.baseline.mean)) +
  geom_line() + geom_point() +
  scale_x_continuous(name="End of Baseline Period [year]", 
                     breaks=seq(min(df.Q.baseline$yr.baseline.end), max(df.Q.baseline$yr.baseline.end))) +
  scale_y_continuous(name="Baseline Period Mean Discharge [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(plot.dir, "GHCN_PheasantBranch_SensitivityBaseline_p.baseline.discharge.png"),
       p.baseline.discharge, width=8, height=6, units="in")

# density plot of LULC and climate effects for prediction period
p.dens.LULC.baseline <-
  ggplot(melt(subset(df.LULCvClimate.prediction, select=c("year", "yr.baseline.end", "change.overall.mean", "change.climate.mean", "change.LULC.mean")), 
              id=c("year", "yr.baseline.end")), 
         aes(x=value, color=factor(yr.baseline.end))) +
  geom_vline(xintercept=0, color="gray65") +
  geom_density() +
  facet_grid(.~variable, scales="free", 
             labeller=as_labeller(c("change.overall.mean"="Overall", "change.climate.mean"="Climate", "change.LULC.mean"="LULC"))) +
  scale_x_continuous(name="Discharge Change from Baseline Period [mm]") +
  scale_y_continuous(name="Density") +
  scale_color_discrete(name="End of Baseline Period [year]") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")
ggsave(paste0(plot.dir, "GHCN_PheasantBranch_SensitivityBaseline_p.dens.LULC.baseline.png"),
       p.dens.LULC.baseline, width=8, height=6, units="in")

# ribbon plot of LULC and climate effects facet by yr.baseline.end
p.ribbon <- 
  ggplot(df.LULCvClimate.baseline, aes(x=year)) +
  geom_line(aes(y=change.overall.mean), color="black") +
  geom_ribbon(aes(ymin=change.climate.min, ymax=change.climate.max), fill="red", alpha=0.25) +
  geom_ribbon(aes(ymin=change.LULC.min, ymax=change.LULC.max), fill="green", alpha=0.25) +
  geom_hline(yintercept=0, color="gray65") +
  facet_wrap(~yr.baseline.end, scales="free_y") +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(plot.dir, "GHCN_PheasantBranch_SensitivityBaseline_p.ribbon.png"),
       p.ribbon, width=8, height=6, units="in")

ks.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1989)$change.overall.mean, 
        subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.overall.mean)$p.value

ks.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1989)$change.LULC.mean, 
        subset(df.LULCvClimate.prediction, yr.baseline.end==1990)$change.LULC.mean)$p.value

ks.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1989)$change.climate.mean, 
        subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.climate.mean)$p.value
