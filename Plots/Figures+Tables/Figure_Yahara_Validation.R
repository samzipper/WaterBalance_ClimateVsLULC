## Figure_Yahara_Validation.R
# This makes a plot showing only validation data.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(dplyr)
require(ggplot2)
require(reshape2)
require(lubridate)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path to save figure output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# use HistLULC_XXclim simulations?
hist.LULC <- F   # T or F

# which method to use? options are PCR, PLS
method <- "PLS"
if (method=="PCR"){
  if (hist.LULC){
    df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/HistLULC_AllClim_CalVal_MonthlyRegressions_OutputAll.csv"), stringsAsFactors=F)
  } else {
    df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/CalVal_MonthlyRegressions_OutputAll.csv"), stringsAsFactors=F)
  }
  colnames(df)[colnames(df)=="PCR"] <- "stat"
  
} else if (method=="PLS") {
  if (hist.LULC){
    df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/HistLULC_AllClim_CalVal_MonthlyRegressions_PLS_OutputAll.csv"), stringsAsFactors=F)
  } else {
    df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/CalVal_MonthlyRegressions_PLS_OutputAll.csv"), stringsAsFactors=F)
  }
  colnames(df)[colnames(df)=="PLS"] <- "stat"
}

# subset to validation only
df.val <- subset(df, group=="val")

# melt for boxplot
df.val.melt <- melt(df.val, id=c("year","month","flux.name","perm", "group"))

# summarize for timeseries
df.val.mo <- dplyr::summarize(group_by(df.val, year, month, flux.name),
                              AI.mean = mean(flux), 
                              stat.mean = mean(stat))
df.val.mo$date <- ymd(paste0(df.val.mo$year, "-", df.val.mo$month, "-15"))

df.val.mo.melt <- melt(df.val.mo, id=c("year", "month", "date", "flux.name"))

# summarize yearly for fit statistics
df.val.yr <- dplyr::summarize(group_by(df.val, year, flux.name),
                              AI.mean = mean(flux), 
                              stat.mean = mean(stat))

### fit statistics
NSE.ET.mo <- NashSutcliffe(subset(df.val, flux.name=="aet")$stat, subset(df.val, flux.name=="aet")$flux)
RMSE.ET.mo <- RMSE(subset(df.val, flux.name=="aet")$stat, subset(df.val, flux.name=="aet")$flux)
NRMSE.ET.mo <- NRMSE(subset(df.val, flux.name=="aet")$stat, subset(df.val, flux.name=="aet")$flux)

NSE.ET.yr <- NashSutcliffe(subset(df.val.yr, flux.name=="aet")$stat.mean, subset(df.val.yr, flux.name=="aet")$AI.mean)
RMSE.ET.yr <- RMSE(subset(df.val.yr, flux.name=="aet")$stat.mean, subset(df.val.yr, flux.name=="aet")$AI.mean)
NRMSE.ET.yr <- NRMSE(subset(df.val.yr, flux.name=="aet")$stat.mean, subset(df.val.yr, flux.name=="aet")$AI.mean)

NSE.drainage.mo <- NashSutcliffe(subset(df.val, flux.name=="drainage")$stat, subset(df.val, flux.name=="drainage")$flux)
RMSE.drainage.mo <- RMSE(subset(df.val, flux.name=="drainage")$stat, subset(df.val, flux.name=="drainage")$flux)
NRMSE.drainage.mo <- NRMSE(subset(df.val, flux.name=="drainage")$stat, subset(df.val, flux.name=="drainage")$flux)

NSE.drainage.yr <- NashSutcliffe(subset(df.val.yr, flux.name=="drainage")$stat.mean, subset(df.val.yr, flux.name=="drainage")$AI.mean)
RMSE.drainage.yr <- RMSE(subset(df.val.yr, flux.name=="drainage")$stat.mean, subset(df.val.yr, flux.name=="drainage")$AI.mean)
NRMSE.drainage.yr <- NRMSE(subset(df.val.yr, flux.name=="drainage")$stat.mean, subset(df.val.yr, flux.name=="drainage")$AI.mean)

NSE.srunoff.mo <- NashSutcliffe(subset(df.val, flux.name=="srunoff")$stat, subset(df.val, flux.name=="srunoff")$flux)
RMSE.srunoff.mo <- RMSE(subset(df.val, flux.name=="srunoff")$stat, subset(df.val, flux.name=="srunoff")$flux)
NRMSE.srunoff.mo <- NRMSE(subset(df.val, flux.name=="srunoff")$stat, subset(df.val, flux.name=="srunoff")$flux)

NSE.srunoff.yr <- NashSutcliffe(subset(df.val.yr, flux.name=="srunoff")$stat.mean, subset(df.val.yr, flux.name=="srunoff")$AI.mean)
RMSE.srunoff.yr <- RMSE(subset(df.val.yr, flux.name=="srunoff")$stat.mean, subset(df.val.yr, flux.name=="srunoff")$AI.mean)
NRMSE.srunoff.yr <- NRMSE(subset(df.val.yr, flux.name=="srunoff")$stat.mean, subset(df.val.yr, flux.name=="srunoff")$AI.mean)

# names for facet variables
var_names <- c("aet"="ET", "drainage"="Drainage", "srunoff"="Runoff")

# make plot
p.val.scatter <-
  ggplot(df.val, aes(y=flux, x=stat)) +
  geom_point(shape=21, alpha=0.1) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  stat_smooth(method="lm") +
  facet_wrap(~flux.name, scales="free", labeller=as_labeller(var_names), ncol=1) +
  scale_x_continuous(name=paste0(method, " Estimated Monthly Flux [mm]")) +
  scale_y_continuous(name="AgroIBIS Modeled Monthly Flux [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

p.val.box <- 
  ggplot(df.val.melt, aes(x=factor(month), y=value, fill=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_boxplot(outlier.shape=1, outlier.fill=NULL, width=1) +
  facet_wrap(~flux.name, scales="free", labeller=as_labeller(var_names), ncol=1) +
  scale_x_discrete(name="Month") +
  scale_y_continuous(name="Monthly Flux [mm]") +
  scale_fill_manual(name="Source", labels=c("AI"="AI", "stat"=method), 
                    values=c("AI"="white", "stat"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

tiff(file=paste0(path.fig, "Figure_Yahara_Validation_NoText.tif"), width=(181/25.4), height=(200/25.4), units="in", res=300)
grid.arrange(p.val.scatter+theme(text=element_blank(), plot.margin=unit(c(0,6,0,0), "mm"), panel.spacing=unit(9, "mm")),
             p.val.box+theme(text=element_blank(), plot.margin=unit(c(0,0,0,6), "mm"), panel.spacing=unit(9, "mm")),
             ncol=2)
dev.off()