## PheasantBranch_4_CalculateFitMetrics.R
# This script is intended to take the output of PheasantBranch_3_MonthlyRegressions.R and
# calculate fit metrics and make plots of fit.

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

# directory to save plots
plot.dir <- paste0(git.dir, "Plots/PheasantBranch/")

## which baseline year to use?
yr.baseline.end.all <- c(1992, 1995)

## which flux to analyze?
flux.name.all <- c("discharge.mm", "runoff.mm", "baseflow.mm")

for (flux.name in flux.name.all){
  for (yr.baseline.end in yr.baseline.end.all){
    
    # read in data
    df <- read.csv(paste0(flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_OutputAll.csv"))
    
    # date column
    df$date <- ymd(paste0(df$year, "-", df$month, "-15"))
    
    # subset to validation data only and make a plot
    df.val <- subset(df, group=="val")
    df.val.melt <- melt(df.val, id=c("year", "month", "perm", "group", "date"))
    df.val.mean <- dplyr::summarize(group_by(df.val, year, month),
                                    date = mean(date),
                                    flux.mean = mean(flux),
                                    PCR.mean = mean(PCR),
                                    PCR.sd = sd(PCR))
    df.val.mean.melt <- melt(df.val.mean[,c("year", "month", "date", "flux.mean", "PCR.mean")], id=c("year", "month", "date"))
    df.val.mean.yr <- dplyr::summarize(group_by(df.val.mean, year),
                                       flux = sum(flux.mean),
                                       PCR = sum(PCR.mean))
    
    # statistics
    val.R2 <- R2(df.val$PCR, df.val$flux)
    val.RMSE <- RMSE(df.val$PCR, df.val$flux)
    val.NRMSE <- NRMSE(df.val$PCR, df.val$flux)
    val.NSE <- NashSutcliffe(df.val$PCR, df.val$flux)
    val.NSE.mean <- NashSutcliffe(df.val.mean$PCR.mean, df.val.mean$flux.mean)
    
    val.yr.R2 <- R2(df.val.mean.yr$PCR, df.val.mean.yr$flux)
    val.yr.RMSE <- RMSE(df.val.mean.yr$PCR, df.val.mean.yr$flux)
    val.yr.NRMSE <- NRMSE(df.val.mean.yr$PCR, df.val.mean.yr$flux)
    val.yr.NSE <- NashSutcliffe(df.val.mean.yr$PCR, df.val.mean.yr$flux)
    
    df.val.fit <- data.frame(time=c("mo", "yr"),
                             R2=round(c(val.R2, val.yr.R2), 3),
                             RMSE=round(c(val.RMSE, val.yr.RMSE), 3),
                             NRMSE=round(c(val.NRMSE, val.yr.NRMSE), 3),
                             NSE=round(c(val.NSE, val.yr.NSE), 3),
                             NSE.mean=round(c(val.NSE.mean, NaN), 3))
    
    ## monthly plots
    # simple x/y plot
    p.val <- 
      ggplot(df.val, aes(x=PCR, y=flux)) +
      geom_point(aes(color=factor(month))) +
      geom_abline(intercept=0, slope=1, color="gray65") +
      stat_smooth(method="lm") +
      theme_bw() +
      theme(panel.grid=element_blank(),
            legend.position="bottom")
    
    # timeseries
    p.time <- 
      ggplot(df.val.mean.melt, aes(x=date, y=value, color=variable)) +
      geom_line() +
      scale_x_date(name="Date", expand=c(0,0)) +
      scale_y_continuous(name="Discharge [mm]") +
      scale_color_discrete(name="Source", labels=c("flux.mean"="Obs.", "PCR.mean"="PCR")) +
      theme_bw() +
      theme(panel.grid=element_blank(),
            legend.position="bottom")
    
    # monthly boxplots
    p.box <- 
      ggplot(df.val.melt, aes(x=factor(month), y=value, fill=variable)) +
      geom_boxplot() +
      scale_x_discrete(name="Month") +
      scale_y_continuous(name="Discharge [mm]") +
      scale_fill_discrete(name="Source", labels=c("flux"="Obs.", "PCR"="PCR")) +
      theme_bw() +
      theme(panel.grid=element_blank(),
            legend.position="bottom")
    
    ggsave(paste0(plot.dir, "GHCN_PheasantBranch_Validation_", flux.name, "_", yr.baseline.end, ".png"),
           arrangeGrob(arrangeGrob(p.val, p.time, p.box, ncol=3), tableGrob(df.val.fit), ncol=1, heights=c(1,0.3)), 
           width=12, height=5, units="in")
  }
}