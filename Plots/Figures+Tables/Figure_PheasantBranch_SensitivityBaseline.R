## Figure_PheasantBranch_SensitivityBaseline.R
# requires output from script PheasantBranch_6_SensitivityBaseline.R

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

# years to plot
yrs.plot <- seq(1989,1995,1)

# read in data
df.out <- read.csv(paste0(git.dir, "Data/PheasantBranch/GHCN_SensitivityBaseline_OutputAll.csv"), stringsAsFactors=F)

# subset to plot years
df.out <- subset(df.out, yr.baseline.end %in% yrs.plot)

## calculate fit metrics by yr.baseline.end
df.baseline <- dplyr::summarize(group_by(subset(df.out, group != "prediction"), yr.baseline.end, group),
                                RMSE = RMSE(PCR, flux))
df.baseline.melt <- melt(df.baseline, id=c("group", "yr.baseline.end"))

## calculate baseline period discharge
df.year.perm.baseline <-
  dplyr::summarize(group_by(df.out, year, perm, yr.baseline.end),
                   obs.sum = sum(flux),
                   PCR.sum = sum(PCR))
df.year.baseline <-
  dplyr::summarize(group_by(df.year.perm.baseline, year, yr.baseline.end),
                   obs.mean = mean(obs.sum),
                   PCR.mean = mean(PCR.sum))

df.Q.baseline <-
  dplyr::summarize(group_by(df.year.baseline, yr.baseline.end),
                   obs.baseline.mean = mean(obs.mean[year<=yr.baseline.end]))

## separate LULC vs climate
df.LULCvClimate <- merge(df.year.perm.baseline, df.Q.baseline, by=c("yr.baseline.end"), all.x=T)
df.LULCvClimate$change.overall <- df.LULCvClimate$obs.sum - df.LULCvClimate$obs.baseline.mean
df.LULCvClimate$change.climate <- df.LULCvClimate$PCR.sum - df.LULCvClimate$obs.baseline.mean
df.LULCvClimate$change.LULC <- df.LULCvClimate$change.overall - df.LULCvClimate$change.climate

df.LULCvClimate.baseline <-
  dplyr::summarize(group_by(df.LULCvClimate, year, yr.baseline.end),
                   obs = mean(obs.sum),
                   PCR.mean = mean(PCR.sum),
                   PCR.sd = sd(PCR.sum),
                   change.overall.mean = mean(change.overall),
                   change.overall.sd = sd(change.overall),
                   change.climate.mean = mean(change.climate),
                   change.climate.sd = sd(change.climate),
                   change.LULC.mean = mean(change.LULC),
                   change.LULC.sd = sd(change.LULC))

df.LULCvClimate.prediction <- subset(df.LULCvClimate.baseline, year>max(df.LULCvClimate$yr.baseline.end))

# density plot of LULC and climate effects for prediction period
p.dens.LULC.baseline <-
  ggplot(melt(subset(df.LULCvClimate.prediction, select=c("year", "yr.baseline.end", "change.overall.mean", "change.climate.mean", "change.LULC.mean")), 
              id=c("year", "yr.baseline.end")), 
         aes(x=value, color=factor(yr.baseline.end))) +
  geom_vline(xintercept=0, color="gray65") +
  geom_density() +
  facet_grid(.~variable, scales="free", 
             labeller=as_labeller(c("change.overall.mean"="Overall", "change.climate.mean"="Climate", "change.LULC.mean"="LULC"))) +
  scale_x_continuous(name="Discharge Change from Baseline Period [mm]") +
  scale_y_continuous(name="Density") +
  scale_color_discrete(name="End of Baseline Period [year]") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

# barplot of fit metrics by group and yr.baseline.end
p.fit.group.baseline <-
  ggplot(df.baseline.melt, aes(x=yr.baseline.end, y=value, fill=group)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity", position="dodge") +
  scale_x_continuous(name="End of Baseline Period [year]", 
                     breaks=seq(min(df.baseline$yr.baseline.end), max(df.baseline$yr.baseline.end))) +
  scale_y_continuous(name="RMSE [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(plot.dir, "GHCN_PheasantBranch_SensitivityBaseline_p.fit.group.baseline.png"),
       p.fit.group.baseline, width=8, height=6, units="in")





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

## save plots
pdf(file=paste0(path.fig, "Figure_PheasantBranch_NoText.pdf"), width=(181/25.4), height=(120/25.4))
grid.arrange(p.val.time+theme(text=element_blank(), plot.margin=unit(c(0.5,6,4,0), "mm")),
             p.val.box+theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,4,6), "mm")),
             p.ribbon.static+theme(text=element_blank(), plot.margin=unit(c(4,6,0,0), "mm")),
             p.climate.LULC.hist+theme(text=element_blank(), plot.margin=unit(c(4,0.5,0,6), "mm")),
             ncol=2)
dev.off()

