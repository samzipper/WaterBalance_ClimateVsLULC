## Figure_PheasantBranch_ComparePCR-MLR-PLS.R
# requires output from script PheasantBranch_3a_MonthlyRegressions-PCR.R,
# PheasantBranch_3b_MonthlyRegressions-MLR.R, and PheasantBranch_3c_MonthlyRegressions-PLS.R

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
flux.name <- "discharge.mm"
yr.baseline.end <- 1995

# read in data
df.PCR <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_OutputAll.csv"))
df.MLR <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_MLR_OutputAll.csv"))
df.PLS <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_PLS_OutputAll.csv"))

# combine
df <- left_join(df.PCR, df.MLR, by=c("year", "month", "perm", "group", "flux"))
df <- left_join(df, df.PLS, by=c("year", "month", "perm", "group", "flux"))

# date column
df$date <- ymd(paste0(df$year, "-", df$month, "-15"))

# subset to validation data only and make a plot
df.val <- subset(df, group=="val")
df.val.melt <- melt(df.val, id=c("year", "month", "perm", "group", "date"))
df.val.mean <- dplyr::summarize(group_by(df.val, year, month),
                                date = mean(date),
                                obs.mean = mean(flux),
                                PCR.mean = mean(PCR),
                                PCR.sd = sd(PCR),
                                MLR.mean = mean(MLR),
                                MLR.sd = sd(MLR),
                                PLS.mean = mean(PLS),
                                PLS.sd = sd(PLS))
df.val.mean.melt <- melt(df.val.mean[,c("year", "month", "date", "obs.mean", "PCR.mean", "MLR.mean", "PLS.mean")], id=c("year", "month", "date"))
df.val.mean.yr <- dplyr::summarize(group_by(df.val.mean, year),
                                   obs = sum(obs.mean),
                                   PCR = sum(PCR.mean),
                                   MLR = sum(MLR.mean),
                                   PLS = sum(PLS.mean))

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
                                PCR.sum = sum(PCR),
                                MLR.sum = sum(MLR),
                                PLS.sum = sum(PLS))
df.perm.ann$obs.baseline <- Q.ann.baseline.mean

df.perm.ann$change.overall <- df.perm.ann$obs.sum - df.perm.ann$obs.baseline
df.perm.ann$change.climate.PCR <- df.perm.ann$PCR.sum - df.perm.ann$obs.baseline
df.perm.ann$change.climate.MLR <- df.perm.ann$MLR.sum - df.perm.ann$obs.baseline
df.perm.ann$change.climate.PLS <- df.perm.ann$PLS.sum - df.perm.ann$obs.baseline
df.perm.ann$change.LULC.PCR <- df.perm.ann$change.overall - df.perm.ann$change.climate.PCR
df.perm.ann$change.LULC.MLR <- df.perm.ann$change.overall - df.perm.ann$change.climate.MLR
df.perm.ann$change.LULC.PLS <- df.perm.ann$change.overall - df.perm.ann$change.climate.PLS

# calculate mean and standard deviation
df.ann <- dplyr::summarize(group_by(df.perm.ann, year),
                           obs = mean(obs.sum),
                           PCR.mean = mean(PCR.sum),
                           PCR.sd = sd(PCR.sum),
                           MLR.mean = mean(MLR.sum),
                           MLR.sd = sd(MLR.sum),
                           PLS.mean = mean(PLS.sum),
                           PLS.sd = sd(PLS.sum),
                           change.overall.mean = mean(change.overall),
                           change.overall.sd = sd(change.overall),
                           change.climate.PCR.mean = mean(change.climate.PCR),
                           change.climate.PCR.sd = sd(change.climate.PCR),
                           change.LULC.PCR.mean = mean(change.LULC.PCR),
                           change.LULC.PCR.sd = sd(change.LULC.PCR),
                           change.climate.MLR.mean = mean(change.climate.MLR),
                           change.climate.MLR.sd = sd(change.climate.MLR),
                           change.LULC.MLR.mean = mean(change.LULC.MLR),
                           change.LULC.MLR.sd = sd(change.LULC.MLR),
                           change.climate.PLS.mean = mean(change.climate.PLS),
                           change.climate.PLS.sd = sd(change.climate.PLS),
                           change.LULC.PLS.mean = mean(change.LULC.PLS),
                           change.LULC.PLS.sd = sd(change.LULC.PLS))

# calculate some min/max columns to use for ribbons
df.ann$PCR.min <- df.ann$PCR.mean - df.ann$PCR.sd
df.ann$PCR.max <- df.ann$PCR.mean + df.ann$PCR.sd

df.ann$change.overall.min <- df.ann$change.overall.mean - df.ann$change.overall.sd
df.ann$change.overall.max <- df.ann$change.overall.mean + df.ann$change.overall.sd
df.ann$change.climate.PCR.min <- df.ann$change.climate.PCR.mean - df.ann$change.climate.PCR.sd
df.ann$change.climate.PCR.max <- df.ann$change.climate.PCR.mean + df.ann$change.climate.PCR.sd
df.ann$change.LULC.PCR.min <- df.ann$change.LULC.PCR.mean - df.ann$change.LULC.PCR.sd
df.ann$change.LULC.PCR.max <- df.ann$change.LULC.PCR.mean + df.ann$change.LULC.PCR.sd
df.ann$change.climate.MLR.min <- df.ann$change.climate.MLR.mean - df.ann$change.climate.MLR.sd
df.ann$change.climate.MLR.max <- df.ann$change.climate.MLR.mean + df.ann$change.climate.MLR.sd
df.ann$change.LULC.MLR.min <- df.ann$change.LULC.MLR.mean - df.ann$change.LULC.MLR.sd
df.ann$change.LULC.MLR.max <- df.ann$change.LULC.MLR.mean + df.ann$change.LULC.MLR.sd
df.ann$change.climate.PLS.min <- df.ann$change.climate.PLS.mean - df.ann$change.climate.PLS.sd
df.ann$change.climate.PLS.max <- df.ann$change.climate.PLS.mean + df.ann$change.climate.PLS.sd
df.ann$change.LULC.PLS.min <- df.ann$change.LULC.PLS.mean - df.ann$change.LULC.PLS.sd
df.ann$change.LULC.PLS.max <- df.ann$change.LULC.PLS.mean + df.ann$change.LULC.PLS.sd

# melt for histogram
df.ann.melt <- melt(df.ann[,c("year", "change.overall.mean", 
                              "change.climate.PCR.mean", "change.LULC.PCR.mean", 
                              "change.climate.MLR.mean", "change.LULC.MLR.mean",
                              "change.climate.PLS.mean", "change.LULC.PLS.mean")], id="year")
df.ann.melt <- subset(df.ann.melt, year>yr.baseline.end)

df.perm.ann.melt <- melt(df.perm.ann[,c("year", "change.overall", "change.climate.PCR", "change.LULC.PCR", 
                                        "change.climate.MLR", "change.LULC.MLR", 
                                        "change.climate.PLS", "change.LULC.PLS")], id="year")
df.perm.ann.melt <- subset(df.perm.ann.melt, year>yr.baseline.end)

#### MAKE PLOTS
## validation plots
# statistics
val.PCR.R2 <- R2(df.val$PCR, df.val$flux)
val.PCR.RMSE <- RMSE(df.val$PCR, df.val$flux)
val.PCR.NRMSE <- NRMSE(df.val$PCR, df.val$flux)
val.PCR.NSE <- round(NashSutcliffe(df.val$PCR, df.val$flux), 3)
val.PCR.NSE.mean <- NashSutcliffe(df.val.mean$PCR.mean, df.val.mean$obs.mean)

val.MLR.R2 <- R2(df.val$MLR, df.val$flux)
val.MLR.RMSE <- RMSE(df.val$MLR, df.val$flux)
val.MLR.NRMSE <- NRMSE(df.val$MLR, df.val$flux)
val.MLR.NSE <- round(NashSutcliffe(df.val$MLR, df.val$flux), 3)
val.MLR.NSE.mean <- NashSutcliffe(df.val.mean$MLR.mean, df.val.mean$obs.mean)

val.PLS.R2 <- R2(df.val$PLS, df.val$flux)
val.PLS.RMSE <- RMSE(df.val$PLS, df.val$flux)
val.PLS.NRMSE <- NRMSE(df.val$PLS, df.val$flux)
val.PLS.NSE <- round(NashSutcliffe(df.val$PLS, df.val$flux), 3)
val.PLS.NSE.mean <- NashSutcliffe(df.val.mean$PLS.mean, df.val.mean$obs.mean)

val.PCR.yr.R2 <- R2(df.val.mean.yr$PCR, df.val.mean.yr$obs)
val.PCR.yr.RMSE <- RMSE(df.val.mean.yr$PCR, df.val.mean.yr$obs)
val.PCR.yr.NRMSE <- NRMSE(df.val.mean.yr$PCR, df.val.mean.yr$obs)
val.PCR.yr.NSE <- round(NashSutcliffe(df.val.mean.yr$PCR, df.val.mean.yr$obs), 3)

val.MLR.yr.R2 <- R2(df.val.mean.yr$MLR, df.val.mean.yr$obs)
val.MLR.yr.RMSE <- RMSE(df.val.mean.yr$MLR, df.val.mean.yr$obs)
val.MLR.yr.NRMSE <- NRMSE(df.val.mean.yr$MLR, df.val.mean.yr$obs)
val.MLR.yr.NSE <- round(NashSutcliffe(df.val.mean.yr$MLR, df.val.mean.yr$obs), 3)

val.PLS.yr.R2 <- R2(df.val.mean.yr$PLS, df.val.mean.yr$obs)
val.PLS.yr.RMSE <- RMSE(df.val.mean.yr$PLS, df.val.mean.yr$obs)
val.PLS.yr.NRMSE <- NRMSE(df.val.mean.yr$PLS, df.val.mean.yr$obs)
val.PLS.yr.NSE <- round(NashSutcliffe(df.val.mean.yr$PLS, df.val.mean.yr$obs), 3)

# set factor order
df.val.mean.melt$variable <- factor(df.val.mean.melt$variable, levels=c("obs.mean", "PLS.mean", "PCR.mean", "MLR.mean"))
df.val.melt$variable <- factor(df.val.melt$variable, levels=c("flux", "PLS", "PCR", "MLR"))

# timeseries
p.val.time.PLS <- 
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.val.mean.melt, variable=="obs.mean"), aes(x=date, y=value), color="black") +
  geom_line(data=subset(df.val.mean.melt, variable=="PLS.mean"), aes(x=date, y=value), color="#127D7D") +
  scale_x_date(name="Date", expand=c(0,0)) +
  scale_y_continuous(name="Discharge [mm]", limits=c(0,60), breaks=seq(0,60,20)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.val.time.PCR <- 
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.val.mean.melt, variable=="obs.mean"), aes(x=date, y=value), color="black") +
  geom_line(data=subset(df.val.mean.melt, variable=="PCR.mean"), aes(x=date, y=value), color="#D01E1E") +
  scale_x_date(name="Date", expand=c(0,0)) +
  scale_y_continuous(name="Discharge [mm]", limits=c(0,60), breaks=seq(0,60,20)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.val.time.MLR <- 
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.val.mean.melt, variable=="obs.mean"), aes(x=date, y=value), color="black") +
  geom_line(data=subset(df.val.mean.melt, variable=="MLR.mean"), aes(x=date, y=value), color="#D0951E") +
  scale_x_date(name="Date", expand=c(0,0)) +
  scale_y_continuous(name="Discharge [mm]", limits=c(0,60), breaks=seq(0,60,20)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

# monthly boxplots
p.val.box <- 
  ggplot(df.val.melt, aes(x=factor(month), y=value, fill=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_boxplot(outlier.shape=1, outlier.fill=NULL) +
  scale_x_discrete(name="Month") +
  scale_y_continuous(name="Discharge [mm]", limits=c(0,77.5), breaks=seq(0,75,15)) +
  scale_fill_manual(name="Source", labels=c("flux"="Obs.", "PCR"="PCR", "MLR"="MLR", "PLS"="PLS"), 
                      values=c("flux"="white", "PLS"="#127D7D", "MLR"="#D0951E", "PCR"="#D01E1E"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

## results plots
# plots
p.ribbon.static.PLS <- 
  ggplot(df.ann, aes(x=year)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=change.overall.mean), color="black") +
  geom_ribbon(aes(ymin=change.climate.PLS.min, ymax=change.climate.PLS.max), fill="#127D7D", alpha=0.5) +
  geom_line(aes(y=change.climate.PLS.mean), color="#127D7D", alpha=0.5) +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]", breaks=seq(-50,150,50)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.ribbon.static.PCR <- 
  ggplot(df.ann, aes(x=year)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=change.overall.mean), color="black") +
  geom_ribbon(aes(ymin=change.climate.PCR.min, ymax=change.climate.PCR.max), fill="#D01E1E", alpha=0.5) +
  geom_line(aes(y=change.climate.PCR.mean), color="#D01E1E", alpha=0.5) +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]", breaks=seq(-50,150,50)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.ribbon.static.MLR <- 
  ggplot(df.ann, aes(x=year)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=change.overall.mean), color="black") +
  geom_ribbon(aes(ymin=change.climate.MLR.min, ymax=change.climate.MLR.max), fill="#D0951E", alpha=0.5) +
  geom_line(aes(y=change.climate.MLR.mean), color="#D0951E", alpha=0.5) +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]", breaks=seq(-50,150,50)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

## save plots
ggsave(paste0(path.fig, "Figure_PheasantBranch_ComparePCR-MLR-PLS_NoText.pdf"), 
       grid.arrange(
         arrangeGrob(
           arrangeGrob(p.val.time.PLS+theme(text=element_blank(), plot.margin=unit(c(0,6,0.5,0), "mm")),
                       p.val.time.PCR+theme(text=element_blank(), plot.margin=unit(c(0,6,0.5,0), "mm")),
                       p.val.time.MLR+theme(text=element_blank(), plot.margin=unit(c(0,6,0.5,0), "mm")),
                       ncol=1),
           arrangeGrob(p.ribbon.static.PLS+theme(text=element_blank(), plot.margin=unit(c(0,0,0.5,6), "mm")),
                       p.ribbon.static.PCR+theme(text=element_blank(), plot.margin=unit(c(0,0,0.5,6), "mm")),
                       p.ribbon.static.MLR+theme(text=element_blank(), plot.margin=unit(c(0,0,0.5,6), "mm")),
                       ncol=1),
           ncol=2),
         p.val.box+theme(text=element_blank(), plot.margin=unit(c(8,0,0,0), "mm")), ncol=1, heights=c(1, 0.6)), 
       width=(181/25.4), height=(180/25.4), units="in")


## statistics
# prediction period
mean(subset(df.ann, year>yr.baseline.end)$change.overall.mean)     # overall mean
sd(subset(df.ann, year>yr.baseline.end)$change.overall.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.overall.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.climate.PCR.mean)     # PCR climate mean
sd(subset(df.ann, year>yr.baseline.end)$change.climate.PCR.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.climate.PCR.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.PCR.mean)        # PCR LULC mean
sd(subset(df.ann, year>yr.baseline.end)$change.LULC.PCR.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.LULC.PCR.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.climate.MLR.mean)     # MLR climate mean
sd(subset(df.ann, year>yr.baseline.end)$change.climate.MLR.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.climate.MLR.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.MLR.mean)        # MLR LULC mean
sd(subset(df.ann, year>yr.baseline.end)$change.LULC.MLR.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.LULC.MLR.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.climate.PLS.mean)     # PLS climate mean
sd(subset(df.ann, year>yr.baseline.end)$change.climate.PLS.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.climate.PLS.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.PLS.mean)        # PLS LULC mean
sd(subset(df.ann, year>yr.baseline.end)$change.LULC.PLS.mean)
t.test(subset(df.ann, year>yr.baseline.end)$change.LULC.PLS.mean)

# calibration period
mean(subset(df.ann, year<=yr.baseline.end)$change.overall.mean)     # overall mean
sd(subset(df.ann, year<=yr.baseline.end)$change.overall.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.overall.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.climate.PCR.mean)     # PCR climate mean
sd(subset(df.ann, year<=yr.baseline.end)$change.climate.PCR.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.climate.PCR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.PCR.mean)        # PCR LULC mean
sd(subset(df.ann, year<=yr.baseline.end)$change.LULC.PCR.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.LULC.PCR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.climate.MLR.mean)     # MLR climate mean
sd(subset(df.ann, year<=yr.baseline.end)$change.climate.MLR.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.climate.MLR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.MLR.mean)        # MLR LULC mean
sd(subset(df.ann, year<=yr.baseline.end)$change.LULC.MLR.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.LULC.MLR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.climate.PLS.mean)     # PLS climate mean
sd(subset(df.ann, year<=yr.baseline.end)$change.climate.PLS.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.climate.PLS.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.PLS.mean)        # PLS LULC mean
sd(subset(df.ann, year<=yr.baseline.end)$change.LULC.PLS.mean)
t.test(subset(df.ann, year<=yr.baseline.end)$change.LULC.PLS.mean)

mean(subset(df.ann, year>yr.baseline.end)$change.climate.mean)/mean(subset(df.ann, year>yr.baseline.end)$change.overall.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.PCR.mean)/mean(subset(df.ann, year>yr.baseline.end)$change.overall.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.MLR.mean)/mean(subset(df.ann, year>yr.baseline.end)$change.overall.mean)
mean(subset(df.ann, year>yr.baseline.end)$change.LULC.PLS.mean)/mean(subset(df.ann, year>yr.baseline.end)$change.overall.mean)

sum(subset(df.ann, year>yr.baseline.end)$change.LULC.PCR.mean>0)        # PCR LULC positive effect number of years
sum(subset(df.ann, year>yr.baseline.end)$change.climate.PCR.mean>0)
sum(subset(df.ann, year>yr.baseline.end)$change.LULC.MLR.mean>0)        # MLR LULC positive effect number of years
sum(subset(df.ann, year>yr.baseline.end)$change.climate.MLR.mean>0)
sum(subset(df.ann, year>yr.baseline.end)$change.LULC.PLS.mean>0)        # PLS LULC positive effect number of years
sum(subset(df.ann, year>yr.baseline.end)$change.climate.PLS.mean>0)


summary(lm(change.overall.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))   # overall trend
summary(lm(change.climate.PCR.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))   # PCR climate trend
summary(lm(change.LULC.PCR.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))      # PCR LULC trend
summary(lm(change.climate.MLR.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))   # MLR climate trend
summary(lm(change.LULC.MLR.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))      # MLR LULC trend
summary(lm(change.climate.PLS.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))   # PLS climate trend
summary(lm(change.LULC.PLS.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))      # PLS LULC trend


summary(lm(change.overall.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))   # overall trend
summary(lm(change.climate.PCR.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))   # PCR climate trend
summary(lm(change.LULC.PCR.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))      # PCR LULC trend
summary(lm(change.climate.MLR.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))   # MLR climate trend
summary(lm(change.LULC.MLR.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))      # MLR LULC trend
summary(lm(change.climate.PLS.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))   # PLS climate trend
summary(lm(change.LULC.PLS.mean ~ year, data=subset(df.ann, year<=yr.baseline.end)))      # PLS LULC trend

coef(lm(change.LULC.PCR.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]/coef(lm(change.overall.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]
coef(lm(change.LULC.MLR.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]/coef(lm(change.overall.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]
coef(lm(change.LULC.PLS.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]/coef(lm(change.overall.mean ~ year, data=subset(df.ann, year>yr.baseline.end)))[2]


p.climate.LULC.hist <-
  ggplot(subset(df.ann.melt, variable %in% c("change.climate.MLR.mean", "change.climate.PCR.mean", "change.climate.PLS.mean"))) +
  geom_vline(xintercept=0, color="gray65") +
  geom_density(aes(x=value, fill=variable), alpha=0.5, color=NA) +
  geom_density(data=subset(df.ann.melt, variable=="change.overall.mean"), aes(x=value), color="black", fill=NA) +
  scale_x_continuous(name="Change in Annual Runoff Depth [mm]", expand=c(0,0)) +
  scale_y_continuous(name="Density", expand=c(0,0)) +
  scale_fill_manual(name="Driver: ", 
                    values=c("change.climate.PLS.mean"="#127D7D", "change.climate.MLR.mean"="#D0951E", "change.climate.PCR.mean"="#D01E1E"), 
                    labels=c("change.climate.PCR.mean"="PCR", "change.climate.MLR.mean"="MLR", "change.climate.PLS.mean"="PLS"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

mean(subset(df.ann, year<=yr.baseline.end)$change.climate.PCR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.PCR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.climate.MLR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.MLR.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.climate.PLS.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.LULC.PLS.mean)
mean(subset(df.ann, year<=yr.baseline.end)$change.overall.mean)
