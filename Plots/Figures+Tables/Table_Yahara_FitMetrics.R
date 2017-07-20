## Table_Yahara_FitMetrics.R
# This summarizes the calibration results for each month and overall.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(dplyr)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path to save table output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# load validation data (produced by script 1_WaterBalance_MonthlyRegressions_PastOnly.R)
df <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/HistLULC_AllClim_CalVal_MonthlyRegressions_OutputAll.csv"), stringsAsFactors=F)

# subset to validation only
df.val <- subset(df, group=="val")

# summarize to monthly
df.val.mo <- dplyr::summarize(group_by(df.val, flux.name, month),
                              R2 = R2(PCR, flux),
                              RMSE = RMSE(PCR, flux),
                              NRMSE = NRMSE(PCR, flux),
                              NSE = NashSutcliffe(PCR, flux))
df.val.mo$month <- as.character(df.val.mo$month)

# overall output
df.val.flux.name <- dplyr::summarize(group_by(df.val, flux.name),
                                     month = "overall",
                                     R2 = R2(PCR, flux),
                                     RMSE = RMSE(PCR, flux),
                                     NRMSE = NRMSE(PCR, flux),
                                     NSE = NashSutcliffe(PCR, flux))

# combine data frames
df.val.out <- df.val.mo
df.val.out[(dim(df.val.out)[1]+1):(dim(df.val.out)[1]+3), ] <- df.val.flux.name

# write as CSV
write.csv(df.val.out, paste0(path.fig, "Table_Yahara_FitMetrics.csv"), row.names=F)
