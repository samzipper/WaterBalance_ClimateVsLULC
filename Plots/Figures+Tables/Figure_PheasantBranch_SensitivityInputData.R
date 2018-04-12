## Figure_PheasantBranch_SensitivityInputData.R
# requires output from script PheasantBranch_7_SensitivityInputData.R

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

# name of flux to plot
flux.name <- "discharge.mm"

yr.baseline.end <- 1995

# which method to use? options are PCR, PLS, MLR
method <- "PLS"
if (method=="PLS") {
  df <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/GHCN_SensitivityInputData_PLS_ClimateVsLULC.csv"), stringsAsFactors=F)
}

p.change.climate <- 
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_ribbon(data=subset(df, met.data=="GHCN"), 
              aes(x=year, ymin=change.climate.min, ymax=change.climate.max), fill="#8181ffff", alpha=0.75) +
  geom_ribbon(data=subset(df, met.data=="ARL"), 
              aes(x=year, ymin=change.climate.min, ymax=change.climate.max), fill="#ff7a7aff", alpha=0.75) +
  geom_line(data=subset(df, met.data=="ARL"),
            aes(x=year, y=change.overall.mean), color="black") +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.change.LULC <- 
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_ribbon(data=subset(df, met.data=="GHCN"), 
              aes(x=year, ymin=change.LULC.min, ymax=change.LULC.max), fill="#8181ffff", alpha=0.75) +
  geom_ribbon(data=subset(df, met.data=="ARL"), 
              aes(x=year, ymin=change.LULC.min, ymax=change.LULC.max), fill="#ff7a7aff", alpha=0.75) +
  geom_line(data=subset(df, met.data=="ARL"),
            aes(x=year, y=change.overall.mean), color="black") +
  scale_x_continuous(name="Year", expand=c(0,0)) +
  scale_y_continuous(name="Change from Baseline Period [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

pdf(file=paste0(path.fig, "Figure_PheasantBranch_SensitivityInputData_NoText.pdf"), width=(82/25.4), height=(60/25.4))
grid.arrange(p.change.LULC+theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,0.5,0), "mm")),
             p.change.climate+theme(text=element_blank(), plot.margin=unit(c(0.5,0.5,0.5,0), "mm")),
             ncol=1)
dev.off()

# calculate fit statistics
df.all <- read.csv(paste0(flux.name, "/GHCN_SensitivityInputData_PLS_OutputAll.csv"))
df.fit <- 
  df.all %>% 
  subset(group != "prediction") %>% 
  group_by(group, met.data) %>% 
  summarize(NSE = NashSutcliffe(PLS, flux),
            RMSE = RMSE(PLS, flux),
            NRMSE = NRMSE(PLS, flux))

df.fit.ann <- 
  df.all %>% 
  subset(group != "prediction") %>% 
  group_by(group, met.data, year, month) %>% 
  summarize(flux.mean = mean(flux),
            PLS.mean = mean(PLS)) %>% 
  group_by(group, met.data, year) %>% 
  summarize(flux.sum = sum(flux.mean),
            PLS.sum = sum(PLS.mean)) %>% 
  group_by(group, met.data) %>% 
  summarize(NSE = NashSutcliffe(PLS.sum, flux.sum),
            RMSE = RMSE(PLS.sum, flux.sum),
            NRMSE = NRMSE(PLS.sum, flux.sum))

df.wide <- 
  subset(df, select=c("year", "met.data", "change.climate.mean")) %>% 
  dcast(year ~ met.data, value.var="change.climate.mean")

ggplot(df.wide, aes(x=GHCN, y=ARL)) + 
  geom_abline(intercept=0, slope=1) + 
  geom_point() +
  stat_smooth(method="lm")

summary(lm(ARL ~ GHCN, data=df.wide))
cor(df.wide$ARL, df.wide$GHCN, method="pearson")

sum(df.wide$GHCN[df.wide$year>yr.baseline.end] > 0 & df.wide$ARL[df.wide$year>yr.baseline.end] < 0) + 
  sum(df.wide$GHCN[df.wide$year>yr.baseline.end] < 0 & df.wide$ARL[df.wide$year>yr.baseline.end] > 0)

t.test(df.wide$GHCN[df.wide$year>yr.baseline.end], df.wide$ARL[df.wide$year>yr.baseline.end])
ks.test(df.wide$GHCN[df.wide$year>yr.baseline.end], df.wide$ARL[df.wide$year>yr.baseline.end])


RMSE(df.wide$ARL[df.wide$year>yr.baseline.end], df.wide$GHCN[df.wide$year>yr.baseline.end])
NRMSE(df.wide$ARL[df.wide$year>yr.baseline.end], df.wide$GHCN[df.wide$year>yr.baseline.end])
