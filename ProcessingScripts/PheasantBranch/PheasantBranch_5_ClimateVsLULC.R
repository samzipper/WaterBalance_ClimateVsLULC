## PheasantBranch_5_ClimateVsLULC.R
# This script is intended to take the output of PheasantBranch_3_MonthlyRegressions.R and
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

# path to folder with output data
setwd(paste0(git.dir, "Data/PheasantBranch/"))

## which baseline year to use?
yr.baseline.end.all <- c(1992, 1995)

## which flux to analyze?
flux.name.all <- c("discharge.mm", "runoff.mm", "baseflow.mm")

for (flux.name in flux.name.all){
  for (yr.baseline.end in yr.baseline.end.all){
    
    # read in data
    df <- read.csv(paste0(flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_OutputAll.csv"))
    df.Q <- read.csv("PheasantBranch_BaseflowSeparation_Monthly.csv")
    
    # directory to save plots
    plot.dir <- paste0(git.dir, "Plots/PheasantBranch/")
    
    # rename column
    colnames(df.Q)[colnames(df.Q)==flux.name] <- "flux.mm"
    
    # date column
    df$date <- ymd(paste0(df$year, "-", df$month, "-15"))
    
    # figure out first/last year of baseline period
    yr.baseline.start <- min(subset(df, group != "prediction")$year)
    
    # subset discharge data to common years and sum to annual
    df.Q <- subset(df.Q, year %in% df$year)
    df.Q.ann <- dplyr::summarize(group_by(df.Q, year),
                                 flux.mm = sum(flux.mm))
    
    # mean annual discharge for baseline period
    Q.ann.baseline.mean <- mean(subset(df.Q.ann, year>= yr.baseline.start & year<= yr.baseline.end)$flux.mm)
    Q.ann.baseline.sd <- sd(subset(df.Q.ann, year>= yr.baseline.start & year<= yr.baseline.end)$flux.mm)
    
    # for each permutation, sum to annual total
    df.perm.ann <- dplyr::summarize(group_by(df, year, perm),
                                    obs.sum = sum(flux),
                                    PCR.sum = sum(PCR))
    
    # merge baseline with df.perm.ann
    df.perm.ann$obs.baseline.static <- Q.ann.baseline.mean
    
    # calculate overall change, climate, and LULC contribution
    df.perm.ann$change.overall.static <- df.perm.ann$obs.sum - df.perm.ann$obs.baseline.static
    df.perm.ann$change.climate.static <- df.perm.ann$PCR.sum - df.perm.ann$obs.baseline.static
    df.perm.ann$change.LULC.static <- df.perm.ann$change.overall.static - df.perm.ann$change.climate.static
    
    # calculate mean and standard deviation
    df.ann <- dplyr::summarize(group_by(df.perm.ann, year),
                               obs = mean(obs.sum),
                               PCR.mean = mean(PCR.sum),
                               PCR.sd = sd(PCR.sum),
                               change.overall.static.mean = mean(change.overall.static),
                               change.overall.static.sd = sd(change.overall.static),
                               change.climate.static.mean = mean(change.climate.static),
                               change.climate.static.sd = sd(change.climate.static),
                               change.LULC.static.mean = mean(change.LULC.static),
                               change.LULC.static.sd = sd(change.LULC.static))
    
    # calculate some min/max columns to use for ribbons
    df.ann$change.overall.static.min <- df.ann$change.overall.static.mean - df.ann$change.overall.static.sd
    df.ann$change.overall.static.max <- df.ann$change.overall.static.mean + df.ann$change.overall.static.sd
    df.ann$change.climate.static.min <- df.ann$change.climate.static.mean - df.ann$change.climate.static.sd
    df.ann$change.climate.static.max <- df.ann$change.climate.static.mean + df.ann$change.climate.static.sd
    df.ann$change.LULC.static.min <- df.ann$change.LULC.static.mean - df.ann$change.LULC.static.sd
    df.ann$change.LULC.static.max <- df.ann$change.LULC.static.mean + df.ann$change.LULC.static.sd
    
    # melt for histogram
    df.ann.melt <- melt(df.ann[,c("year", "change.overall.static.mean", "change.climate.static.mean", "change.LULC.static.mean")], id="year")
    df.ann.melt <- subset(df.ann.melt, year>yr.baseline.end)
    
    # histogram
    p.climate.LULC.hist <-
      ggplot(subset(df.ann.melt, variable != "change.overall.static.mean")) +
      geom_vline(xintercept=0, color="gray65") +
      geom_density(aes(x=value, fill=variable), alpha=0.5, color=NA) +
      geom_density(data=subset(df.ann.melt, variable=="change.overall.static.mean"), aes(x=value), color="black", fill=NA) +
      ggtitle(flux.name) +
      scale_x_continuous(name="Change in Annual Depth [mm]") +
      scale_y_continuous(name="Density") +
      scale_fill_manual(name="Driver: ", 
                        values=c("change.climate.static.mean"="red", "change.LULC.static.mean"="forestgreen"), 
                        labels=c("change.climate.static.mean"="Climate", "change.LULC.static.mean"="LULC")) +
      theme_bw() +
      theme(panel.grid=element_blank(),
            legend.position="bottom")
    
    # make some ribbon plots
    p.ribbon.static <- 
      ggplot(df.ann, aes(x=year)) +
      geom_line(aes(y=change.overall.static.mean), color="black") +
      geom_ribbon(aes(ymin=change.climate.static.min, ymax=change.climate.static.max), fill="red", alpha=0.25) +
      geom_ribbon(aes(ymin=change.LULC.static.min, ymax=change.LULC.static.max), fill="green", alpha=0.25) +
      geom_hline(yintercept=0, color="gray65") +
      ggtitle(flux.name) +
      scale_x_continuous(name="Year", expand=c(0,0)) +
      scale_y_continuous(name="Change from Baseline Period [mm]") +
      theme_bw() +
      theme(panel.grid=element_blank())
    ggsave(paste0(plot.dir, "GHCN_PheasantBranch_ClimateVsLULC_AnnualStatic_", flux.name, "_", yr.baseline.end, ".png"),
           arrangeGrob(p.ribbon.static, p.climate.LULC.hist, ncol=2), width=12, height=6, units="in")
    
  }
}