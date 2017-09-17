## Figure_PheasantBranch_GroundwaterLevels.R
#' This script is intended to plot groundwater levels, baseflow, and total discharge
#' through time for the PBS.
#' 
#' Groundwater data source: https://nwis.waterdata.usgs.gov/nwis/gwlevels/?site_no=430638089353101

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(dataRetrieval)
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
# groundwater data
dv.gw <- readNWISgwl(siteNumber="430638089353101",                          # site code (can be a vector of multiple sites)
                    startDate="1900-01-01",endDate="2016-12-31")    # start & end dates (YYYY-MM-DD format)

df.Q <- read.csv(paste0(git.dir, "Data/PheasantBranch/PheasantBranch_BaseflowSeparation_Monthly.csv"))


## organize data
# groundwater
df.gw <- data.frame(date = dv.gw$lev_dt,
                    value = -dv.gw$lev_va*0.3048,
                    variable="WTD.m")

# discharge
df.Q$date <- ymd(paste0(df.Q$year, "-", df.Q$month, "-", round(days_in_month(df.Q$month)/2)))

# combine
df.Q.melt <- melt(subset(df.Q, select=c("date", "baseflow.mm")), id=c("date"))
df <- rbind(df.gw, df.Q.melt)
df$variable <- factor(df$variable, levels=c("WTD.m", "baseflow.mm"))

# format for scatterplot
df.gw$year <- year(df.gw$date)
df.gw$month <- month(df.gw$date)

df.scatter <- left_join(df.gw[,c("year", "month", "value")], df.Q[,c("year", "month", "baseflow.mm")])

## plot
p.time <- 
  ggplot(df, aes(x=date, y=value)) +
  geom_line() +
  facet_wrap(~variable, ncol=1, scales="free_y") +
  scale_x_date(expand=c(0,0)) +
  stat_smooth(method="lm") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.scatter <-
  ggplot(df.scatter, aes(x=-value, y=baseflow.mm)) + 
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point(shape=21) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))
  
pdf(file=paste0(path.fig, "Figure_PheasantBranch_GroundwaterLevels_NoText.pdf"), width=(182/25.4), height=(70/25.4))
grid.arrange(p.time + theme(text=element_blank(), plot.margin=unit(c(0.5,8,0.5,0.5), "mm")),
             p.scatter + theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm")),
             ncol=2, widths=c(1,0.75))
dev.off()

# statistics
coef(lm(value ~ date, data=subset(df, variable=="WTD.m")))[2]*365.25*10
summary(lm(value ~ date, data=subset(df, variable=="WTD.m")))

coef(lm(value ~ date, data=subset(df, variable=="baseflow.mm")))[2]*365.25*10
summary(lm(value ~ date, data=subset(df, variable=="baseflow.mm")))
