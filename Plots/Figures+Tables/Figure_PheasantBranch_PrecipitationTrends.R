## Figure_PheasantBranch_PrecipitationTrends.R
# requires output from script PheasantBranch_2_DailyToMonthly.R

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

## load output
df <- read.csv(file.path(git.dir, "Data", "PheasantBranch", "USW00014837_GHCN_Monthly.csv"))

# summarize to annual
df.ann <- 
  df %>% 
  group_by(year) %>% 
  summarize(prec.mm = sum(prec),
            prec.max = max(prec.max),
            prec.gt.50 = sum(prec.gt.50),
            prec.gt.76 = sum(prec.gt.76))

# plot
p.prec <-
  ggplot(df.ann, aes(x=year, y=prec.mm)) +
  geom_point() + 
  geom_line() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Year", limits=c(min(df.ann$year), max(df.ann$year)), expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))


summary(lm(prec.mm ~ year, data=df.ann))

## save plots
pdf(file=paste0(path.fig, "Figure_PheasantBranch_PrecipitationTrend_NoText.pdf"), width=(85/25.4), height=(30/25.4))
p.prec+theme(text=element_blank(), plot.margin=unit(c(0,0,0,0), "mm"))
dev.off()