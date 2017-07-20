## Figure_Scenarios_LULC+Climate.R
# This makes a plot showing climate and LULC for different scenarios.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(dplyr)
require(ggplot2)
require(reshape2)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path to save figure output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# path to model data
path.data <-paste0(git.dir, "Data/Yahara-AgroIBIS/input/")

# length of moving average
yrs.avg <- 11

### First: LULC data
## open LULC data
df.LULC.hist <- read.csv(paste0(path.data, "LULCsums_1786-2014.csv"), stringsAsFactors=F)
df.LULC.AI <- read.csv(paste0(path.data, "LULCsums_2010_2070_AI.csv"), stringsAsFactors=F)
df.LULC.AR <- read.csv(paste0(path.data, "LULCsums_2010_2070_AR.csv"), stringsAsFactors=F)
df.LULC.CC <- read.csv(paste0(path.data, "LULCsums_2010_2070_CC.csv"), stringsAsFactors=F)
df.LULC.NW <- read.csv(paste0(path.data, "LULCsums_2010_2070_NW.csv"), stringsAsFactors=F)

# make a total urban category
df.LULC.AI$urb.total <- df.LULC.AI$Developed.High.Intensity + df.LULC.AI$Developed.Medium.Intensity + df.LULC.AI$Developed.Low.Intensity
df.LULC.AR$urb.total <- df.LULC.AR$Developed.High.Intensity + df.LULC.AR$Developed.Medium.Intensity + df.LULC.AR$Developed.Low.Intensity
df.LULC.CC$urb.total <- df.LULC.CC$Developed.High.Intensity + df.LULC.CC$Developed.Medium.Intensity + df.LULC.CC$Developed.Low.Intensity
df.LULC.NW$urb.total <- df.LULC.NW$Developed.High.Intensity + df.LULC.NW$Developed.Medium.Intensity + df.LULC.NW$Developed.Low.Intensity
df.LULC.hist$urb.total <- df.LULC.hist$Developed.High.Intensity + df.LULC.hist$Developed.Medium.Intensity + df.LULC.hist$Developed.Low.Intensity
  
df.LULC.AI$wetland.total <- df.LULC.AI$Herbaceous.Wetlands + df.LULC.AI$Woody.Wetlands
df.LULC.AR$wetland.total <- df.LULC.AR$Herbaceous.Wetlands + df.LULC.AR$Woody.Wetlands
df.LULC.CC$wetland.total <- df.LULC.CC$Herbaceous.Wetlands + df.LULC.CC$Woody.Wetlands
df.LULC.NW$wetland.total <- df.LULC.NW$Herbaceous.Wetlands + df.LULC.NW$Woody.Wetlands
df.LULC.hist$wetland.total <- df.LULC.hist$Herbaceous.Wetlands + df.LULC.hist$Woody.Wetlands

# subset historical data to 1986-2013
df.LULC.hist <- melt(subset(df.LULC.hist, Year>=1986 & Year <= 2013), id=c("Year"),
                     value.name="area", variable.name="LULC")

# add scenario column
df.LULC.AI$scenario <- "AI"
df.LULC.AR$scenario <- "AR"
df.LULC.CC$scenario <- "CC"
df.LULC.NW$scenario <- "NW"

# bind scenarios and melt
df.LULC <- melt(rbind(df.LULC.AI, df.LULC.AR, df.LULC.CC, df.LULC.NW), 
                id=c("Year", "scenario"), value.name="area", variable.name="LULC")

# subset scenario data to 2014-2070
df.LULC <- subset(df.LULC, Year >= 2014)

# convert from number of pixels to area (square meters)
df.LULC.hist$area <- df.LULC.hist$area*(219.456^2)   # [m2]
df.LULC$area <- df.LULC$area*(219.456^2)   # [m2]

# figure out the total area using a random year
area.total <- sum(subset(df.LULC, Year==2021 & scenario=="AI" & LULC != "urb.total" & LULC != "wetland.total")$area)

# calculate area as % of watershed
df.LULC$area.prc <- 100*df.LULC$area/area.total
df.LULC.hist$area.prc <- 100*df.LULC.hist$area/area.total

# for historical data, figure out some summary statistics for each LULC
df.LULC.hist.summary <- dplyr::summarize(group_by(df.LULC.hist, LULC),
                                         area.prc.min = min(area.prc),
                                         area.prc.max = max(area.prc),
                                         area.prc.mean = mean(area.prc),
                                         area.prc.std = sd(area.prc))

### Second: Climate data
## load data
df.clim.AI <- read.csv(paste0(path.data, "AI_met_yearly.csv"), stringsAsFactors=F)
df.clim.AR <- read.csv(paste0(path.data, "AR_met_yearly.csv"), stringsAsFactors=F)
df.clim.CC <- read.csv(paste0(path.data, "NW_met_yearly.csv"), stringsAsFactors=F)
df.clim.NW <- read.csv(paste0(path.data, "CC_met_yearly.csv"), stringsAsFactors=F)

## read in daily and calculate days with > 1" precip
df.clim.AI.d <- read.csv(paste0(path.data, "AI_met_daily.csv"), stringsAsFactors=F)
df.clim.AR.d <- read.csv(paste0(path.data, "AR_met_daily.csv"), stringsAsFactors=F)
df.clim.CC.d <- read.csv(paste0(path.data, "NW_met_daily.csv"), stringsAsFactors=F)
df.clim.NW.d <- read.csv(paste0(path.data, "CC_met_daily.csv"), stringsAsFactors=F)

df.clim.AI.d.y <- dplyr::summarize(group_by(df.clim.AI.d, year),
                                   prec.gt.25 = sum(prec>25.4))
df.clim.AR.d.y <- dplyr::summarize(group_by(df.clim.AR.d, year),
                                   prec.gt.25 = sum(prec>25.4))
df.clim.CC.d.y <- dplyr::summarize(group_by(df.clim.CC.d, year),
                                   prec.gt.25 = sum(prec>25.4))
df.clim.NW.d.y <- dplyr::summarize(group_by(df.clim.NW.d, year),
                                   prec.gt.25 = sum(prec>25.4))

df.clim.AI <- merge(df.clim.AI, df.clim.AI.d.y, by="year")
df.clim.AR <- merge(df.clim.AR, df.clim.AR.d.y, by="year")
df.clim.CC <- merge(df.clim.CC, df.clim.CC.d.y, by="year")
df.clim.NW <- merge(df.clim.NW, df.clim.NW.d.y, by="year")

# add deficit column
df.clim.AI$defc <- df.clim.AI$RET - df.clim.AI$prec
df.clim.AR$defc <- df.clim.AR$RET - df.clim.AR$prec
df.clim.CC$defc <- df.clim.CC$RET - df.clim.CC$prec
df.clim.NW$defc <- df.clim.NW$RET - df.clim.NW$prec

# get historical data
df.clim.hist <- subset(df.clim.AI, year <= 2013)

# add scenario column
df.clim.AI$scenario <- "AI"
df.clim.AR$scenario <- "AR"
df.clim.CC$scenario <- "CC"
df.clim.NW$scenario <- "NW"

# calculate moving averages
df.clim.AI$prec.avg <- as.numeric(stats::filter(df.clim.AI$prec, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AI$prec.gt.25.avg <- as.numeric(stats::filter(df.clim.AI$prec.gt.25, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AI$tmax.avg <- as.numeric(stats::filter(df.clim.AI$tmax, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AI$rads.avg <- as.numeric(stats::filter(df.clim.AI$rads, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AI$RET.avg <- as.numeric(stats::filter(df.clim.AI$RET, rep(1/yrs.avg, yrs.avg), sides=1))

df.clim.AR$prec.avg <- as.numeric(stats::filter(df.clim.AR$prec, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AR$prec.gt.25.avg <- as.numeric(stats::filter(df.clim.AR$prec.gt.25, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AR$tmax.avg <- as.numeric(stats::filter(df.clim.AR$tmax, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AR$rads.avg <- as.numeric(stats::filter(df.clim.AR$rads, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.AR$RET.avg <- as.numeric(stats::filter(df.clim.AR$RET, rep(1/yrs.avg, yrs.avg), sides=1))

df.clim.CC$prec.avg <- as.numeric(stats::filter(df.clim.CC$prec, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.CC$prec.gt.25.avg <- as.numeric(stats::filter(df.clim.CC$prec.gt.25, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.CC$tmax.avg <- as.numeric(stats::filter(df.clim.CC$tmax, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.CC$rads.avg <- as.numeric(stats::filter(df.clim.CC$rads, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.CC$RET.avg <- as.numeric(stats::filter(df.clim.CC$RET, rep(1/yrs.avg, yrs.avg), sides=1))

df.clim.NW$prec.avg <- as.numeric(stats::filter(df.clim.NW$prec, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.NW$prec.gt.25.avg <- as.numeric(stats::filter(df.clim.NW$prec.gt.25, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.NW$tmax.avg <- as.numeric(stats::filter(df.clim.NW$tmax, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.NW$rads.avg <- as.numeric(stats::filter(df.clim.NW$rads, rep(1/yrs.avg, yrs.avg), sides=1))
df.clim.NW$RET.avg <- as.numeric(stats::filter(df.clim.NW$RET, rep(1/yrs.avg, yrs.avg), sides=1))

# bind scenarios and melt
df.clim <- rbind(subset(df.clim.AI, year>=2014), subset(df.clim.AR, year>=2014), 
                 subset(df.clim.CC, year>=2014), subset(df.clim.NW, year>=2014))

df.clim <- melt(df.clim, 
                id=c("year", "scenario"))

# melt historical data and get summary statistics
df.clim.hist <- melt(df.clim.hist, id=c("year"))
df.clim.hist.summary <- dplyr::summarize(group_by(df.clim.hist, variable),
                                         var.min = min(value),
                                         var.max = max(value),
                                         var.mean = mean(value),
                                         var.std = sd(value))

### Last: Make plots

# climate plots - not smoothed
p.clim.prec <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="prec"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="prec"]) +
  geom_line(data=subset(df.clim, variable=="prec"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="prec") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.clim.tmax <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="tmax"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="tmax"]) +
  geom_line(data=subset(df.clim, variable=="tmax"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="tmax") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.clim.rads <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="rads"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="rads"]) +
  geom_line(data=subset(df.clim, variable=="rads"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="rads") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.clim.RET <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="RET"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="RET"]) +
  geom_line(data=subset(df.clim, variable=="RET"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="RET") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

# climate plots - smoothed
p.clim.prec <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="prec"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="prec"]) +
  geom_line(data=subset(df.clim, variable=="prec.avg"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="prec") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.clim.prec.gt.25 <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="prec.gt.25"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="prec.gt.25"]) +
  geom_line(data=subset(df.clim, variable=="prec.gt.25.avg"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="prec.gt.25", limits=c(0,13)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.clim.tmax <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="tmax"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="tmax"]) +
  geom_line(data=subset(df.clim, variable=="tmax.avg"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="tmax") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.clim.rads <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="rads"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="rads"]) +
  geom_line(data=subset(df.clim, variable=="rads.avg"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="rads") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.clim.RET <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.clim.hist.summary$var.min[df.clim.hist.summary$variable=="RET"],
           ymax=df.clim.hist.summary$var.max[df.clim.hist.summary$variable=="RET"]) +
  geom_line(data=subset(df.clim, variable=="RET.avg"), aes(x=year, y=value, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="RET", breaks=seq(650,950,100)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

# LULC plots
p.LULC.Corn <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.LULC.hist.summary$area.prc.min[df.LULC.hist.summary$LULC=="Corn"],
           ymax=df.LULC.hist.summary$area.prc.max[df.LULC.hist.summary$LULC=="Corn"]) +
  geom_line(data=subset(df.LULC, LULC=="Corn"), aes(x=Year, y=area.prc, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="% Corn", limits=c(0,42), breaks=seq(0,40,10)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.LULC.Soybeans <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.LULC.hist.summary$area.prc.min[df.LULC.hist.summary$LULC=="Soybeans"],
           ymax=df.LULC.hist.summary$area.prc.max[df.LULC.hist.summary$LULC=="Soybeans"]) +
  geom_line(data=subset(df.LULC, LULC=="Soybeans"), aes(x=Year, y=area.prc, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="% Soybeans", limits=c(0,42), breaks=seq(0,40,10)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.LULC.urb.total <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.LULC.hist.summary$area.prc.min[df.LULC.hist.summary$LULC=="urb.total"],
           ymax=df.LULC.hist.summary$area.prc.max[df.LULC.hist.summary$LULC=="urb.total"]) +
  geom_line(data=subset(df.LULC, LULC=="urb.total"), aes(x=Year, y=area.prc, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="% urb.total", limits=c(0,42), breaks=seq(0,40,10)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.LULC.Alfalfa <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.LULC.hist.summary$area.prc.min[df.LULC.hist.summary$LULC=="Alfalfa"],
           ymax=df.LULC.hist.summary$area.prc.max[df.LULC.hist.summary$LULC=="Alfalfa"]) +
  geom_line(data=subset(df.LULC, LULC=="Alfalfa"), aes(x=Year, y=area.prc, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="% Alfalfa", limits=c(0,42), breaks=seq(0,40,10)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.LULC.Deciduous.Forest <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.LULC.hist.summary$area.prc.min[df.LULC.hist.summary$LULC=="Deciduous.Forest"],
           ymax=df.LULC.hist.summary$area.prc.max[df.LULC.hist.summary$LULC=="Deciduous.Forest"]) +
  geom_line(data=subset(df.LULC, LULC=="Deciduous.Forest"), aes(x=Year, y=area.prc, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="% Deciduous.Forest", limits=c(0,42), breaks=seq(0,40,10)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.LULC.wetland.total <- 
  ggplot() +
  annotate("rect", alpha=0.25, xmin=-Inf, xmax=Inf, 
           ymin=df.LULC.hist.summary$area.prc.min[df.LULC.hist.summary$LULC=="wetland.total"],
           ymax=df.LULC.hist.summary$area.prc.max[df.LULC.hist.summary$LULC=="wetland.total"]) +
  geom_line(data=subset(df.LULC, LULC=="wetland.total"), aes(x=Year, y=area.prc, color=scenario)) +
  scale_color_manual(name="Scenario", values=c("AR"="#8cc63f", "AI"="#ff931e", "CC"="#0071bc", "NW"="#ff1d25"), guide=F) +
  scale_x_continuous(name="Year", limits=c(2014,2070), breaks=seq(2020,2070, 10), expand=c(0,0)) +
  scale_y_continuous(name="% wetland.total", limits=c(0,42), breaks=seq(0,40,10)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

pdf(file=paste0(path.fig, "Figure_Yahara_Scenarios_LULC+Climate_NoText.pdf"), width=(85/25.4), height=(220/25.4))
grid.arrange(p.clim.prec+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")),
             p.clim.prec.gt.25+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")), 
             p.clim.tmax+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")),
             p.clim.RET+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")), 
             p.LULC.Corn+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")),
             p.LULC.urb.total+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")), 
             p.LULC.Deciduous.Forest+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")), 
             p.LULC.wetland.total+theme(text=element_blank(), plot.margin=unit(c(0.5,1,0,0), "mm")), 
             ncol=1)
dev.off()

