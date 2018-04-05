## Table_Yahara_FitMetrics.R
# This summarizes the calibration results for each month and overall.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(dplyr)
require(magrittr)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path to save table output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# load validation data
# which method to use? options are PCR, PLS, MLR
method <- "PLS"
if (method=="PCR"){
  df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/CalVal_MonthlyRegressions_OutputAll.csv"))
  colnames(df)[colnames(df)=="PCR"] <- "stat"
} else if (method=="MLR") {
  df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/CalVal_MonthlyRegressions_MLR_OutputAll.csv"))
  colnames(df)[colnames(df)=="MLR"] <- "stat"
} else if (method=="PLS") {
  df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/CalVal_MonthlyRegressions_PLS_OutputAll.csv"))
  colnames(df)[colnames(df)=="PLS"] <- "stat"
}

# subset to validation only
df.val <- subset(df, group=="val")

# summarize to monthly
df.val.mo <- dplyr::summarize(group_by(df.val, flux.name, month),
                              R2 = R2(stat, flux),
                              RMSE = RMSE(stat, flux),
                              NRMSE = NRMSE(stat, flux),
                              NSE = NashSutcliffe(stat, flux))
df.val.mo$month <- as.character(df.val.mo$month)

df.val.wide <- cbind(subset(df.val.mo, flux.name=="aet")[,c("NSE", "RMSE", "NRMSE")], 
                     subset(df.val.mo, flux.name=="drainage")[,c("NSE", "RMSE", "NRMSE")], 
                     subset(df.val.mo, flux.name=="srunoff")[,c("NSE", "RMSE", "NRMSE")]) %>% 
  set_colnames(c("aet.NSE", "aet.RMSE", "aet.NRMSE", "drainage.NSE", "drainage.RMSE", "drainage.NRMSE", 
                 "srunoff.NSE", "srunoff.RMSE", "srunoff.NRMSE"))

# overall output
df.val.flux.name <- dplyr::summarize(group_by(df.val, flux.name),
                                     month = "overall",
                                     R2 = R2(stat, flux),
                                     RMSE = RMSE(stat, flux),
                                     NRMSE = NRMSE(stat, flux),
                                     NSE = NashSutcliffe(stat, flux))

df.val.overall.wide <- cbind(subset(df.val.flux.name, flux.name=="aet")[,c("NSE", "RMSE", "NRMSE")], 
                             subset(df.val.flux.name, flux.name=="drainage")[,c("NSE", "RMSE", "NRMSE")], 
                             subset(df.val.flux.name, flux.name=="srunoff")[,c("NSE", "RMSE", "NRMSE")]) %>% 
  set_colnames(c("aet.NSE", "aet.RMSE", "aet.NRMSE", "drainage.NSE", "drainage.RMSE", "drainage.NRMSE", 
                 "srunoff.NSE", "srunoff.RMSE", "srunoff.NRMSE"))

# combine and reorganize
df.val.out <- rbind(df.val.wide, df.val.overall.wide)
df.val.out <- df.val.out[,c("aet.NSE", "drainage.NSE", "srunoff.NSE",
                            "aet.RMSE", "drainage.RMSE", "srunoff.RMSE",
                            "aet.NRMSE", "drainage.NRMSE", "srunoff.NRMSE")]

# write as CSV
write.csv(df.val.out, paste0(path.fig, "Table_Yahara_FitMetrics.csv"), row.names=F)
