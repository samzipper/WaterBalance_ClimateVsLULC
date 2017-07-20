## Figure_PheasantBranch-AgroIBIS_CompareToObs.R
# This script is intended to take the output of PheasantBranch-AgroIBIS_1_MonthlyRegressions.R and
# compares observed runoff (at gauging station, separated with baseflow separation) to AgroIBIS 
# simulated runoff (average of all grid cells within Pheasant Branch subwatershed).

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

# set factor levels
df.melt$variable <- factor(df.melt$variable, levels=c("Obs", "PCR", "AgroIBIS"))
df.ann.melt$variable <- factor(df.ann.melt$variable, levels=c("Obs", "PCR", "AgroIBIS"))

## make plots - don't include PCR
# timeseries
p.time.mo <-
  ggplot(subset(df.melt, variable != "PCR"), aes(x=date, y=value, color=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() +
  scale_x_date(expand=c(0,0)) +
  scale_color_manual(name="Source", labels=c("Obs"="Obs.", "AgroIBIS"="AgroIBIS"),
                     values=c("Obs"="black", "AgroIBIS"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.time.yr <-
  ggplot(subset(df.ann.melt, variable != "PCR"), aes(x=year, y=value, color=variable)) +
  geom_line() +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(name="Source", labels=c("Obs"="Obs.", "AgroIBIS"="AgroIBIS"),
                     values=c("Obs"="black", "AgroIBIS"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

pdf(file=paste0(path.fig, "Figure_PheasantBranch-AgroIBIS_CompareToObs_NoText.pdf"), width=(81/25.4), height=(120/25.4))
grid.arrange(p.time.mo+theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm")), 
             p.time.yr+theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm")),
             ncol=1)
dev.off()


## statistics
# monthly
NashSutcliffe(df.all$AgroIBIS, df.all$Obs)
RMSE(df.all$AgroIBIS, df.all$Obs)
NRMSE(df.all$AgroIBIS, df.all$Obs)
R2(df.all$AgroIBIS, df.all$Obs)

# yearly
NashSutcliffe(df.ann$AgroIBIS, df.ann$Obs)
RMSE(df.ann$AgroIBIS, df.ann$Obs)
NRMSE(df.ann$AgroIBIS, df.ann$Obs)
R2(df.ann$AgroIBIS, df.ann$Obs)
