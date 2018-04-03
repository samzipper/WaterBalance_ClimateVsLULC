## Figure_PheasantBranch_CalVal_Runoff.R
# requires output from script PheasantBranch_3_MonthlyRegressions.R

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(ggplot2)
require(dplyr)
require(lubridate)
require(reshape2)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path to save figure output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# name of flux and baseline year to plot
flux.name <- "runoff.mm"
yr.baseline.start <- 1975
yr.baseline.end <- 1995

# which method to use? options are PCR, PLS, MLR
method <- "PLS"
if (method=="PCR"){
  df <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_OutputAll.csv"))
  colnames(df)[colnames(df)=="PCR"] <- "stat"
} else if (method=="MLR") {
  df <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_MLR_OutputAll.csv"))
  colnames(df)[colnames(df)=="MLR"] <- "stat"
} else if (method=="PLS") {
  df <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_PLS_OutputAll.csv"))
  colnames(df)[colnames(df)=="PLS"] <- "stat"
}

# date column
df$date <- ymd(paste0(df$year, "-", df$month, "-15"))

# subset to validation data only and make a plot
df.val <- subset(df, group=="val")
df.val.melt <- melt(df.val, id=c("year", "month", "perm", "group", "date"))
df.val.mean <- dplyr::summarize(group_by(df.val, year, month),
                                date = mean(date),
                                obs.mean = mean(flux),
                                stat.mean = mean(stat),
                                stat.sd = sd(stat))
df.val.mean.melt <- melt(df.val.mean[,c("year", "month", "date", "obs.mean", "stat.mean")], id=c("year", "month", "date"))
df.val.mean.yr <- dplyr::summarize(group_by(df.val.mean, year),
                                   obs = sum(obs.mean),
                                   stat = sum(stat.mean))

## now process climate vs LULC data
# number of years for moving average
avg.years <- 11

# read in data
df.Q <- read.csv(paste0(git.dir, "Data/PheasantBranch/PheasantBranch_BaseflowSeparation_Monthly.csv"))

# make column for flux of interest
colnames(df.Q)[colnames(df.Q)==flux.name] <- "flux.mm"

# figure out first/last year of baseline period
yr.baseline.start <- min(subset(df, group != "prediction")$year)
yr.baseline.end <- max(subset(df, group != "prediction")$year)

# subset discharge data to common years and sum to annual
df.Q <- subset(df.Q, year %in% df$year)
df.Q.ann <- dplyr::summarize(group_by(df.Q, year),
                             flux.mm = sum(flux.mm))

# mean annual discharge for baseline period
Q.ann.baseline.mean <- mean(subset(df.Q.ann, year >= yr.baseline.start & year <= yr.baseline.end)$flux.mm)
Q.ann.baseline.sd <- sd(subset(df.Q.ann, year >= yr.baseline.start & year <= yr.baseline.end)$flux.mm)

# for each permutation, sum to annual total
df.perm.ann <- dplyr::summarize(group_by(df, year, perm),
                                obs.sum = sum(flux),
                                stat.sum = sum(stat))
df.perm.ann$obs.baseline.static <- Q.ann.baseline.mean

df.perm.ann$change.overall.static <- df.perm.ann$obs.sum - df.perm.ann$obs.baseline.static
df.perm.ann$change.climate.static <- df.perm.ann$stat.sum - df.perm.ann$obs.baseline.static
df.perm.ann$change.LULC.static <- df.perm.ann$change.overall.static - df.perm.ann$change.climate.static

# calculate mean and standard deviation
df.ann <- dplyr::summarize(group_by(df.perm.ann, year),
                           obs = mean(obs.sum),
                           stat.mean = mean(stat.sum),
                           stat.sd = sd(stat.sum),
                           change.overall.static.mean = mean(change.overall.static),
                           change.overall.static.sd = sd(change.overall.static),
                           change.climate.static.mean = mean(change.climate.static),
                           change.climate.static.sd = sd(change.climate.static),
                           change.LULC.static.mean = mean(change.LULC.static),
                           change.LULC.static.sd = sd(change.LULC.static))

# calculate some min/max columns to use for ribbons
df.ann$stat.min <- df.ann$stat.mean - df.ann$stat.sd
df.ann$stat.max <- df.ann$stat.mean + df.ann$stat.sd

df.ann$change.overall.static.min <- df.ann$change.overall.static.mean - df.ann$change.overall.static.sd
df.ann$change.overall.static.max <- df.ann$change.overall.static.mean + df.ann$change.overall.static.sd
df.ann$change.climate.static.min <- df.ann$change.climate.static.mean - df.ann$change.climate.static.sd
df.ann$change.climate.static.max <- df.ann$change.climate.static.mean + df.ann$change.climate.static.sd
df.ann$change.LULC.static.min <- df.ann$change.LULC.static.mean - df.ann$change.LULC.static.sd
df.ann$change.LULC.static.max <- df.ann$change.LULC.static.mean + df.ann$change.LULC.static.sd

# melt for histogram
df.ann.melt <- melt(df.ann[,c("year", "change.overall.static.mean", "change.climate.static.mean", "change.LULC.static.mean")], id="year")
df.ann.melt <- subset(df.ann.melt, year>yr.baseline.end)

df.perm.ann.melt <- melt(df.perm.ann[,c("year", "change.overall.static", "change.climate.static", "change.LULC.static")], id="year")
df.perm.ann.melt <- subset(df.perm.ann.melt, year>yr.baseline.end)

# statistics
val.R2 <- R2(df.val$stat, df.val$flux)
val.RMSE <- RMSE(df.val$stat, df.val$flux)
val.NRMSE <- NRMSE(df.val$stat, df.val$flux)
val.NSE <- round(NashSutcliffe(df.val$stat, df.val$flux), 3)
val.NSE.mean <- NashSutcliffe(df.val.mean$stat.mean, df.val.mean$obs.mean)

val.yr.R2 <- R2(df.val.mean.yr$stat, df.val.mean.yr$obs)
val.yr.RMSE <- RMSE(df.val.mean.yr$stat, df.val.mean.yr$obs)
val.yr.NRMSE <- NRMSE(df.val.mean.yr$stat, df.val.mean.yr$obs)
val.yr.NSE <- round(NashSutcliffe(df.val.mean.yr$stat, df.val.mean.yr$obs), 3)

# timeseries
p.val.time <- 
  ggplot(df.val.mean.melt, aes(x=date, y=value, color=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() +
  scale_x_date(name="Date", expand=c(0,0)) +
  scale_y_continuous(name="Runoff [mm]", limits=c(min(df.val.melt$value), max(df.val.melt$value))) +
  scale_color_manual(name="Source", labels=c("obs.mean"="Obs.", "stat.mean"="stat"),
                     values=c("obs.mean"="black", "stat.mean"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

# monthly boxplots
p.val.box <- 
  ggplot(df.val.melt, aes(x=factor(month), y=value, fill=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_boxplot(outlier.shape=1, outlier.fill=NULL) +
  scale_x_discrete(name="Month") +
  scale_y_continuous(name="Runoff [mm]", limits=c(min(df.val.melt$value), max(df.val.melt$value))) +
  scale_fill_manual(name="Source", labels=c("obs"="Obs.", "stat"="stat"), 
                    values=c("obs"="white", "stat"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

## results plots
# statistics
mean(subset(df.ann, year>yr.baseline.end)$change.overall.static.mean)     # overall mean
sd(subset(df.ann, year>yr.baseline.end)$change.overall.static.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.overall.static.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.climate.static.mean)     # climate mean
sd(subset(df.ann, year>yr.baseline.end)$change.climate.static.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.climate.static.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.static.mean)        # LULC mean
sd(subset(df.ann, year>yr.baseline.end)$change.LULC.static.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.LULC.static.mean)

mean(subset(df.ann, year<=yr.baseline.end)$change.overall.static.mean)     # overall mean
sd(subset(df.ann, year<=yr.baseline.end)$change.overall.static.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.overall.static.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.climate.static.mean)     # climate mean
sd(subset(df.ann, year<=yr.baseline.end)$change.climate.static.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.climate.static.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.static.mean)        # LULC mean
sd(subset(df.ann, year<=yr.baseline.end)$change.LULC.static.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.LULC.static.mean)

mean(subset(df.ann, year>yr.baseline.end)$change.climate.static.mean)/mean(subset(df.ann, year>yr.baseline.end)$change.overall.static.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.static.mean)/mean(subset(df.ann, year>yr.baseline.end)$change.overall.static.mean)

sum(subset(df.ann, year>yr.baseline.end)$change.LULC.static.mean>0)        # LULC positive effect number of years
sum(subset(df.ann, year>yr.baseline.end)$change.climate.static.mean>0)


summary(lm(change.overall.static.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))   # overall trend
summary(lm(change.climate.static.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))   # climate trend
summary(lm(change.LULC.static.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))      # LULC trend


summary(lm(change.overall.static.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))   # overall trend
summary(lm(change.climate.static.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))   # climate trend
summary(lm(change.LULC.static.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))      # LULC trend

coef(lm(change.LULC.static.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]/coef(lm(change.overall.static.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]

# plots
p.ribbon.static <- 
  ggplot(df.ann, aes(x=year)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=change.overall.static.mean), color="black") +
  geom_ribbon(aes(ymin=change.climate.static.min, ymax=change.climate.static.max), fill="#D01D1D", alpha=0.5) +
  geom_ribbon(aes(ymin=change.LULC.static.min, ymax=change.LULC.static.max), fill="#18A718", alpha=0.5) +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]", breaks=seq(-50,150,50)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.climate.LULC.hist <-
  ggplot(subset(df.ann.melt, variable != "change.overall.static.mean")) +
  geom_vline(xintercept=0, color="gray65") +
  geom_density(aes(x=value, fill=variable), alpha=0.5, color=NA) +
  geom_density(data=subset(df.ann.melt, variable=="change.overall.static.mean"), aes(x=value), color="black", fill=NA) +
  scale_x_continuous(name="Change in Annual Runoff Depth [mm]", expand=c(0,0)) +
  scale_y_continuous(name="Density", expand=c(0,0)) +
  scale_fill_manual(name="Driver: ", 
                    values=c("change.climate.static.mean"="#D01D1D", "change.LULC.static.mean"="#18A718"), 
                    labels=c("change.climate.static.mean"="Climate", "change.LULC.static.mean"="LULC"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

mean(subset(df.ann, year<=yr.baseline.end)$change.climate.static.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.static.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.overall.static.mean)

## urbanization through time as % of watershed
# data from file PBS_LULCchange.xlsx, received from Eric Booth via email in June/July 2017
# only considers low, medium, and high intensity urban land use - does not include open space
df.urb <- data.frame(year=c(1992, 2001, 2006, 2011, 2012),
                     urb.prc = c(12.16, 20.85, 28.80, 29.70, 39.06),
                     source = c("WISCLAND", "NLCD", "NLCD", "NLCD", "WISCLAND"))
p.urb <- 
  ggplot(df.urb, aes(x=year, y=urb.prc)) +
  geom_line(color="black") +
  geom_point(aes(shape=source), fill="black") +
  scale_y_continuous(name="% of Watershed Urbanized", limits=c(0, max(df.urb$urb.prc))) +
  scale_x_continuous(name="Year", limits=c(min(df.ann$year), max(df.ann$year)), expand=c(0,0)) +
  scale_shape_manual(values=c("WISCLAND"=24, "NLCD"=22), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

summary(lm(urb.prc ~ year, data=df.urb))  # urban LULC % trend
coef(lm(change.LULC.static.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]/coef(lm(urb.prc ~ year, data=df.urb))[2]

## save plots
pdf(file=paste0(path.fig, "Figure_PheasantBranch_CalVal_Runoff_NoText.pdf"), width=(181/25.4), height=(120/25.4))
grid.arrange(
  arrangeGrob(p.val.time+theme(text=element_blank(), plot.margin=unit(c(0,6,4,0), "mm")),
              p.val.box+theme(text=element_blank(), plot.margin=unit(c(4,6,0,0), "mm")),
              ncol=1),
  arrangeGrob(p.urb+theme(text=element_blank(), plot.margin=unit(c(0,0,8,6), "mm")),
              p.ribbon.static+theme(text=element_blank(), plot.margin=unit(c(0,0,4,6), "mm")),
              p.climate.LULC.hist+theme(text=element_blank(), plot.margin=unit(c(4,0,0,6), "mm")),
              ncol=1),
  ncol=2)
dev.off()