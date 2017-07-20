## PheasantBranch-AgroIBIS_3_CompareToObservations.R
# This script is intended to take the output of PheasantBranch-AgroIBIS_1_MonthlyRegressions.R and
# compares three different runoff values: observed (separated with baseflow separation), AgroIBIS
# simulated, and PCR predicted.

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
setwd(paste0(git.dir, "Data/PheasantBranch-AgroIBIS/"))

# period for comparison
yr.baseline.start <- 1986
yr.baseline.end <- 2013

# which scenarios to analyze?
climate <- "AI"
LULC <- "AI"

## get AgroIBIS and PCR results
# read in output from PheasantBranch-AgroIBIS_1_MonthlyRegressions.R
df <- read.csv(paste0(LULC, climate, "_MonthlyRegressions_OutputAll.csv"), stringsAsFactors=F)

# subset to only runoff and period of interest
df <- subset(df, flux.name=="srunoff" & year >= yr.baseline.start & year <= yr.baseline.end)

# summarize to monthly (mean of all permutations)
df.mo <- dplyr::summarize(group_by(df, year, month),
                          AgroIBIS = mean(flux),
                          PCR = mean(PCR))

## get observations
df.Q <- read.csv(paste0(git.dir, "Data/PheasantBranch/PheasantBranch_BaseflowSeparation_Monthly.csv"))

# subset to period of interest
df.Q <- subset(df.Q, year >= yr.baseline.start & year <= yr.baseline.end)

# get only columns of interest
df.Q <- data.frame(year = df.Q$year, 
                   month = df.Q$month, 
                   Obs = df.Q$runoff.mm)

## merge data frames
df.all <- left_join(df.mo, df.Q, by=c("year", "month"))

# PCR should only be plotted through 1995 (the cal/val period)
df.all$PCR[df.all$year>1995] <- NaN

# date column
df.all$date <- ymd(paste0(df.all$year, "-", df.all$month, "-15"))

# sum to annual
df.ann <- dplyr::summarize(group_by(df.all, year),
                           AgroIBIS = sum(AgroIBIS),
                           PCR = sum(PCR),
                           Obs = sum(Obs))

# melt for plotting
df.melt <- melt(df.all, id=c("year", "month", "date"))
df.ann.melt <- melt(df.ann, id=c("year"))

## make plots
# timeseries
p.time.mo <-
  ggplot(subset(df.melt, variable != "PCR"), aes(x=date, y=value, color=variable)) +
  geom_line() +
  scale_x_date(expand=c(0,0)) +
  theme_bw()

p.time.yr <-
  ggplot(subset(df.ann.melt, variable != "PCR"), aes(x=year, y=value, color=variable)) +
  geom_line() +
  scale_x_continuous(expand=c(0,0)) +
  theme_bw()

# scatters
p.Obs.AgroIBIS.mo <-
  ggplot(df.all, aes(x=AgroIBIS, y=Obs)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()

p.Obs.AgroIBIS.yr <-
  ggplot(df.all.ann, aes(x=AgroIBIS, y=Obs)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()

## statistics
# monthly
NashSutcliffe(df.all$AgroIBIS, df.all$Obs)
RMSE(df.all$AgroIBIS, df.all$Obs)
NRMSE(df.all$AgroIBIS, df.all$Obs)
R2(df.all$AgroIBIS, df.all$Obs)

# yearly
NashSutcliffe(df.all.ann$AgroIBIS, df.all.ann$Obs)
RMSE(df.all.ann$AgroIBIS, df.all.ann$Obs)
NRMSE(df.all.ann$AgroIBIS, df.all.ann$Obs)
R2(df.all.ann$AgroIBIS, df.all.ann$Obs)
