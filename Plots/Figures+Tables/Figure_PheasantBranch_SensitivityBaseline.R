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
yrs.plot <- seq(1992,1998,1)  # my year +/- 3 years

# name of flux to plot
flux.name <- "discharge.mm"

# read in data
df.out <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/GHCN_SensitivityBaseline_OutputAll.csv"), stringsAsFactors=F)

# subset to plot years
df.out <- subset(df.out, yr.baseline.end %in% yrs.plot)

## calculate fit metrics by yr.baseline.end
df.baseline <- dplyr::summarize(group_by(subset(df.out, group != "prediction"), yr.baseline.end, group),
                                NSE = NashSutcliffe(PCR, flux))
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
                   change.LULC.mean = mean(change.LULC),
                   change.LULC.sd = sd(change.LULC),
                   change.climate.mean = mean(change.climate),
                   change.climate.sd = sd(change.climate))

df.LULCvClimate.prediction <- subset(df.LULCvClimate.baseline, year>max(df.LULCvClimate$yr.baseline.end))

df.LULCvClimate.baseline.melt <- melt(df.LULCvClimate.baseline[,c("year", "yr.baseline.end", "change.overall.mean", "change.LULC.mean", "change.climate.mean")],
                                      id=c("year", "yr.baseline.end"))

# barplot of fit metrics by group and yr.baseline.end
p.fit.group.baseline <-
  ggplot(df.baseline.melt, aes(x=yr.baseline.end, y=value, fill=group)) +
  geom_bar(stat="identity", position="dodge", color="white") +
  geom_hline(yintercept=0, color="gray65") +
  scale_x_continuous(name="End of Baseline Period [year]", 
                     breaks=seq(min(df.baseline$yr.baseline.end), max(df.baseline$yr.baseline.end))) +
  scale_y_continuous(name="NSE [mm]") +
  scale_fill_manual(values=c("cal"="#ff1d25", "val"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.change.baseline <- 
  ggplot(subset(df.LULCvClimate.baseline.melt, variable != "change.overall.mean"), aes(x=year, y=value, color=factor(yr.baseline.end))) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() +
  geom_line(data=subset(df.LULCvClimate.baseline.melt, variable != "change.overall.mean" & yr.baseline.end==1995), aes(x=year, y=value), color="black") +
  facet_grid(variable~.,
             labeller=as_labeller(c("change.climate.mean"="Climate", "change.LULC.mean"="LULC"))) +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]") +
  scale_color_manual(guide=F, values=c("1992"="#8181ffff", "1993"="#5656abff", "1994"="#2d2d59ff", 
                                         "1995"="white", 
                                         "1996"="#582a2aff", "1997"="#a50000ff", "1998"="#ff7a7aff")) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom",
        panel.border=element_rect(color="black"))

## save plots
pdf(file=paste0(path.fig, "Figure_PheasantBranch_SensitivityBaseline_NoText.pdf"), width=(82/25.4), height=(120/25.4))
grid.arrange(p.fit.group.baseline+theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,4,0), "mm")),
             p.change.baseline+theme(text=element_blank(), plot.margin=unit(c(4,0.5,0,0), "mm")),
             ncol=1, heights=c(0.25,1))
dev.off()


# comparison of mean LULC and climate effects depending on baseline period end
df.sig.LULC <- data.frame(year=seq(1992,1998),
                     yr.1992=NaN,
                     yr.1993=NaN,
                     yr.1994=NaN,
                     yr.1995=NaN,
                     yr.1996=NaN,
                     yr.1997=NaN,
                     yr.1998=NaN)
df.sig.LULC$yr.1992[df.sig$year==1993] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.LULC.mean)$p.value
df.sig.LULC$yr.1992[df.sig$year==1994] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.LULC.mean)$p.value
df.sig.LULC$yr.1992[df.sig$year==1995] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.LULC.mean)$p.value
df.sig.LULC$yr.1992[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.LULC.mean)$p.value
df.sig.LULC$yr.1992[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.LULC.mean)$p.value
df.sig.LULC$yr.1992[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.LULC.mean)$p.value

df.sig.LULC$yr.1993[df.sig$year==1994] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.LULC.mean)$p.value
df.sig.LULC$yr.1993[df.sig$year==1995] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.LULC.mean)$p.value
df.sig.LULC$yr.1993[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.LULC.mean)$p.value
df.sig.LULC$yr.1993[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.LULC.mean)$p.value
df.sig.LULC$yr.1993[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.LULC.mean)$p.value

df.sig.LULC$yr.1994[df.sig$year==1995] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.LULC.mean)$p.value
df.sig.LULC$yr.1994[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.LULC.mean)$p.value
df.sig.LULC$yr.1994[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.LULC.mean)$p.value
df.sig.LULC$yr.1994[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.LULC.mean)$p.value

df.sig.LULC$yr.1995[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.LULC.mean)$p.value
df.sig.LULC$yr.1995[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.LULC.mean)$p.value
df.sig.LULC$yr.1995[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.LULC.mean)$p.value

df.sig.LULC$yr.1996[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.LULC.mean)$p.value
df.sig.LULC$yr.1996[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.LULC.mean)$p.value

df.sig.LULC$yr.1997[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.LULC.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.LULC.mean)$p.value

df.sig.climate <- data.frame(year=seq(1992,1998),
                             yr.1992=NaN,
                             yr.1993=NaN,
                             yr.1994=NaN,
                             yr.1995=NaN,
                             yr.1996=NaN,
                             yr.1997=NaN,
                             yr.1998=NaN)
df.sig.climate$yr.1992[df.sig$year==1993] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.climate.mean)$p.value
df.sig.climate$yr.1992[df.sig$year==1994] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.climate.mean)$p.value
df.sig.climate$yr.1992[df.sig$year==1995] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.climate.mean)$p.value
df.sig.climate$yr.1992[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.climate.mean)$p.value
df.sig.climate$yr.1992[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.climate.mean)$p.value
df.sig.climate$yr.1992[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1992)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.climate.mean)$p.value

df.sig.climate$yr.1993[df.sig$year==1994] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.climate.mean)$p.value
df.sig.climate$yr.1993[df.sig$year==1995] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.climate.mean)$p.value
df.sig.climate$yr.1993[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.climate.mean)$p.value
df.sig.climate$yr.1993[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.climate.mean)$p.value
df.sig.climate$yr.1993[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1993)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.climate.mean)$p.value

df.sig.climate$yr.1994[df.sig$year==1995] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.climate.mean)$p.value
df.sig.climate$yr.1994[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.climate.mean)$p.value
df.sig.climate$yr.1994[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.climate.mean)$p.value
df.sig.climate$yr.1994[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1994)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.climate.mean)$p.value

df.sig.climate$yr.1995[df.sig$year==1996] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.climate.mean)$p.value
df.sig.climate$yr.1995[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.climate.mean)$p.value
df.sig.climate$yr.1995[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1995)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.climate.mean)$p.value

df.sig.climate$yr.1996[df.sig$year==1997] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.climate.mean)$p.value
df.sig.climate$yr.1996[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1996)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.climate.mean)$p.value

df.sig.climate$yr.1997[df.sig$year==1998] <- t.test(subset(df.LULCvClimate.prediction, yr.baseline.end==1997)$change.climate.mean, subset(df.LULCvClimate.prediction, yr.baseline.end==1998)$change.climate.mean)$p.value

df.sig.climate.logical <- df.sig.climate<0.05
df.sig.LULC.logical <- df.sig.LULC<0.05
df.sig.climate.logical[,1] <- df.sig.climate$year
df.sig.LULC.logical[,1] <- df.sig.LULC$year
