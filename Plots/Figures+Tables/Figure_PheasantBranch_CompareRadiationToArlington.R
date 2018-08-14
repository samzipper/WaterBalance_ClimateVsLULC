## Figure_PheasantBranch_CompareRadiationToArlington.R
#' This script is intended to compare estimated incoming SW radiation
#' estimated from GHCN meteorological data to measurements from the
#' Arlington WI agricultural weather station (1985-2016).
#' 
#' Arlington data source: http://agwx.soils.wisc.edu/uwex_agwx/awon

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(ggplot2)
require(dplyr)
require(lubridate)
require(gridExtra)
require(reshape2)
require(hydroGOF)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path to save figure output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

## load data
# arlington data
df.arl.d <- read.csv(paste0(git.dir, "Data/PheasantBranch/ArlingtonWI_Metdata.csv"), stringsAsFactors=F)
df.arl.d <- unique(df.arl.d)

# GHCN data
df.GHCN.d <- read.csv(paste0(git.dir, "Data/PheasantBranch/USW00014837_GHCN_Daily.csv"), stringsAsFactors=F)

## convert date column
df.arl.d$date <- ymd(df.arl.d$date)
df.GHCN.d$DATE <- ymd(df.GHCN.d$DATE)

## data frame for comparison
df.arl <- data.frame(date = df.arl.d$date,
                     srad.ARL = df.arl.d$DAvSol)

df.GHCN <- data.frame(date = df.GHCN.d$DATE,
                      srad.GHCN = df.GHCN.d$srad)

## merge into single data frame
df <- merge(df.arl, df.GHCN, by=c("date"))

# get rid of nodata and a weird artefact at srad=0 for arlington
df$srad.ARL[df$srad.ARL %in% c(-6999, 6999, 0)] <- NaN
df$srad.ARL[df$srad.ARL > 400] <- NaN

# summarize to monthly
df$year <- year(df$date)
df$month <- month(df$date)

df.mo <- dplyr::summarize(group_by(df, year, month),
                          srad.ARL.mean = mean(srad.ARL, na.rm=F),
                          srad.GHCN.mean = mean(srad.GHCN, na.rm=F))
df.mo <- df.mo[complete.cases(df.mo), ]

## plot
p.d <- 
  ggplot(df, aes(x=srad.ARL, y=srad.GHCN)) +
  geom_point(shape=21, alpha=0.1) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  stat_smooth(method="lm") +
  scale_x_continuous(expand=c(0,0), limits=c(0,max(c(max(df$srad.ARL, na.rm=T), max(df$srad.GHCN, na.rm=T))))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,max(c(max(df$srad.ARL, na.rm=T), max(df$srad.GHCN, na.rm=T))))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.mo <-
  ggplot(df.mo, aes(x=srad.ARL.mean, y=srad.GHCN.mean)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  stat_smooth(method="lm") +
  scale_x_continuous(expand=c(0,0), limits=c(0,max(c(max(df.mo$srad.ARL.mean, na.rm=T), max(df.mo$srad.GHCN.mean, na.rm=T)))+30)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,max(c(max(df.mo$srad.ARL.mean, na.rm=T), max(df.mo$srad.GHCN.mean, na.rm=T)))+30)) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

pdf(file=paste0(path.fig, "Figure_PheasantBranch_CompareRadiationToArlington_NoText.pdf"), width=(120/25.4), height=(55/25.4))
grid.arrange(p.d + theme(text=element_blank(), plot.margin=unit(c(0.5,6,0,0.5), "mm")),
             p.mo + theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,0,6), "mm")),
  ncol=2)
dev.off()


## fit metrics
round(NSE(df$srad.GHCN, df$srad.ARL, na.rm=T), 3)
RMSE(subset(df, is.finite(srad.GHCN) & is.finite(srad.ARL))$srad.GHCN, subset(df, is.finite(srad.GHCN) & is.finite(srad.ARL))$srad.ARL)
pbias(df$srad.GHCN, df$srad.ARL, na.rm=T)

round(NSE(df.mo$srad.GHCN.mean, df.mo$srad.ARL.mean, na.rm=T), 3)
RMSE(df.mo$srad.GHCN.mean, df.mo$srad.ARL.mean)
pbias(df.mo$srad.GHCN.mean, df.mo$srad.ARL.mean, na.rm=T)
