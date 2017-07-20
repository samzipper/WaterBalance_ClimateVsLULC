## Figure_PheasantBranch-AgroIBIS_Runoff.R
# requires output from script PheasantBranch-AgroIBIS_1_MonthlyRegressions.R

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(ggplot2)
require(dplyr)
require(lubridate)
require(reshape2)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path to save figure output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# name of flux and baseline year to plot
flux <- "srunoff"
yr.baseline.end <- 1995

# read in data
df <- read.csv(paste0(git.dir, "Data/PheasantBranch-AgroIBIS/AIAI_MonthlyRegressions_OutputAll.csv"), stringsAsFactors=F)

# subset to flux of interest
df <- df[df$flux.name==flux, ]

# date column
df$date <- ymd(paste0(df$year, "-", df$month, "-15"))

# subset to validation data only and make a plot
df.val <- subset(df, group=="val")
df.val.melt <- melt(df.val, id=c("year", "month", "perm", "group", "date", "flux.name"))
df.val.mean <- dplyr::summarize(group_by(df.val, year, month),
                                date = mean(date),
                                obs.mean = mean(flux),
                                PCR.mean = mean(PCR),
                                PCR.sd = sd(PCR))
df.val.mean.melt <- melt(df.val.mean[,c("year", "month", "date", "obs.mean", "PCR.mean")], id=c("year", "month", "date"))
df.val.mean.yr <- dplyr::summarize(group_by(df.val.mean, year),
                                   obs = sum(obs.mean),
                                   PCR = sum(PCR.mean))

## now process climate vs LULC data
df.ann <- read.csv(paste0(git.dir, "Data/PheasantBranch-AgroIBIS/ClimateVsLULC.csv"))

# subset to flux of interest
df.ann <- df.ann[df.ann$flux.name==flux, ]

# figure out first/last year of baseline period
yr.baseline.start <- min(subset(df, group != "prediction")$year)
yr.baseline.end <- max(subset(df, group != "prediction")$year)

# melt for histogram
df.ann.melt <- melt(df.ann[,c("year", "change.overall.mean", "change.climate.mean", "change.LULC.mean")], id="year")
df.ann.melt <- subset(df.ann.melt, year>yr.baseline.end)

#### MAKE PLOTS
## validation plots
# statistics
val.R2 <- R2(df.val$PCR, df.val$flux)
val.RMSE <- RMSE(df.val$PCR, df.val$flux)
val.NRMSE <- NRMSE(df.val$PCR, df.val$flux)
val.NSE <- NashSutcliffe(df.val$PCR, df.val$flux)
val.NSE.mean <- NashSutcliffe(df.val.mean$PCR.mean, df.val.mean$obs.mean)

val.yr.R2 <- R2(df.val.mean.yr$PCR, df.val.mean.yr$obs)
val.yr.RMSE <- RMSE(df.val.mean.yr$PCR, df.val.mean.yr$obs)
val.yr.NRMSE <- NRMSE(df.val.mean.yr$PCR, df.val.mean.yr$obs)
val.yr.NSE <- NashSutcliffe(df.val.mean.yr$PCR, df.val.mean.yr$obs)

# timeseries
p.val.time <- 
  ggplot(df.val.mean.melt, aes(x=date, y=value, color=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() +
  scale_x_date(name="Date", expand=c(0,0)) +
  scale_y_continuous(name="Runoff [mm]", limits=c(0,66.2), breaks=seq(0,60,15)) +
  scale_color_manual(name="Source", labels=c("obs.mean"="Obs.", "PCR.mean"="PCR"),
                     values=c("obs.mean"="black", "PCR.mean"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

# monthly boxplots
p.val.box <- 
  ggplot(df.val.melt, aes(x=factor(month), y=value, fill=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_boxplot(outlier.shape=1, outlier.fill=NULL) +
  scale_x_discrete(name="Month") +
  scale_y_continuous(name="Discharge [mm]", limits=c(0,66.2), breaks=seq(0,60,15)) +
  scale_fill_manual(name="Source", labels=c("obs"="Obs.", "PCR"="PCR"), 
                    values=c("obs"="white", "PCR"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

## results plots
# statistics
mean(subset(df.ann, year>=1996 & year<=2013)$change.overall.mean)     # overall mean
sd(subset(df.ann, year>=1996 & year<=2013)$change.overall.mean)
t.test(subset(df.ann, year>=1996 & year<=2013)$change.overall.mean)
mean(subset(df.ann, year>=1996 & year<=2013)$change.climate.mean)     # climate mean
sd(subset(df.ann, year>=1996 & year<=2013)$change.climate.mean)
t.test(subset(df.ann, year>=1996 & year<=2013)$change.climate.mean)
mean(subset(df.ann, year>=1996 & year<=2013)$change.LULC.mean)        # LULC mean
sd(subset(df.ann, year>=1996 & year<=2013)$change.LULC.mean)
t.test(subset(df.ann, year>=1996 & year<=2013)$change.LULC.mean)

mean(subset(df.ann, year<1996)$change.overall.mean)     # overall mean
sd(subset(df.ann, year<1996)$change.overall.mean)
t.test(subset(df.ann, year<1996)$change.overall.mean)
mean(subset(df.ann, year<1996)$change.climate.mean)     # climate mean
sd(subset(df.ann, year<1996)$change.climate.mean)
t.test(subset(df.ann, year<1996)$change.climate.mean)
mean(subset(df.ann, year<1996)$change.LULC.mean)        # LULC mean
sd(subset(df.ann, year<1996)$change.LULC.mean)
t.test(subset(df.ann, year<1996)$change.LULC.mean)

mean(subset(df.ann, year>=1996 & year<=2013 & year<=2013)$change.climate.mean)/mean(subset(df.ann, year>=1996 & year<=2013)$change.overall.mean)
mean(subset(df.ann, year>=1996 & year<=2013)$change.LULC.mean)/mean(subset(df.ann, year>=1996 & year<=2013)$change.overall.mean)

sum(subset(df.ann, year>=1996 & year<=2013)$change.LULC.mean>0)        # LULC positive effect number of years
sum(subset(df.ann, year>=1996 & year<=2013)$change.climate.mean>0)


summary(lm(change.overall.mean ~ year, data=subset(df.ann, year>=1996 & year<=2013 & year<=2013)))   # overall trend
summary(lm(change.climate.mean ~ year, data=subset(df.ann, year>=1996 & year<=2013 & year<=2013)))   # climate trend
summary(lm(change.LULC.mean ~ year, data=subset(df.ann, year>=1996 & year<=2013 & year<=2013)))      # LULC trend


summary(lm(change.overall.mean ~ year, data=subset(df.ann, year<1996)))   # overall trend
summary(lm(change.climate.mean ~ year, data=subset(df.ann, year<1996)))   # climate trend
summary(lm(change.LULC.mean ~ year, data=subset(df.ann, year<1996)))      # LULC trend

coef(lm(change.LULC.mean ~ year, data=subset(df.ann, year>=1996 & year<=2013)))[2]/coef(lm(change.overall.mean ~ year, data=subset(df.ann, year>=1996 & year<=2013)))[2]

# plots
p.ribbon <- 
  ggplot(df.ann, aes(x=year)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=change.overall.mean), color="black") +
  geom_ribbon(aes(ymin=change.climate.min, ymax=change.climate.max), fill="#D01D1D", alpha=0.5) +
  geom_ribbon(aes(ymin=change.LULC.min, ymax=change.LULC.max), fill="#18A718", alpha=0.5) +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.climate.LULC.hist <-
  ggplot(subset(df.ann.melt, variable != "change.overall.mean")) +
  geom_vline(xintercept=0, color="gray65") +
  geom_density(aes(x=value, fill=variable), alpha=0.5, color=NA) +
  geom_density(data=subset(df.ann.melt, variable=="change.overall.mean"), aes(x=value), color="black", fill=NA) +
  scale_x_continuous(name="Change in Annual Runoff Depth [mm]", expand=c(0,0), breaks=seq(-30,90,30)) +
  scale_y_continuous(name="Density", expand=c(0,0)) +
  scale_fill_manual(name="Driver: ", 
                    values=c("change.climate.mean"="#D01D1D", "change.LULC.mean"="#18A718"), 
                    labels=c("change.climate.mean"="Climate", "change.LULC.mean"="LULC"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

mean(subset(df.ann, year<=yr.baseline.end)$change.climate.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.overall.mean)

## save plots
pdf(file=paste0(path.fig, "Figure_PheasantBranch-AgroIBIS_Runoff_NoText.pdf"), width=(181/25.4), height=(120/25.4))
grid.arrange(
  arrangeGrob(p.val.time+theme(text=element_blank(), plot.margin=unit(c(0,6,4,0), "mm")),
              p.val.box+theme(text=element_blank(), plot.margin=unit(c(4,6,0,0), "mm")),
              ncol=1),
  arrangeGrob(p.ribbon+theme(text=element_blank(), plot.margin=unit(c(0,0,4,6), "mm")),
              p.climate.LULC.hist+theme(text=element_blank(), plot.margin=unit(c(4,0,0,6), "mm")),
              ncol=1),
  ncol=2)
dev.off()
