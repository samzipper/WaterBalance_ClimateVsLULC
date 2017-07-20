## Figure_PCR_VarsKeep.R
#' This script is intended to make a figure showing the variables that are
#' retained for each PCR analysis.

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(ggplot2)
require(dplyr)
require(reshape2)

# path to save figure output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# list of months
mo.list <- seq(1,12)

## vector of all variable options
var.options <- c("prec", "prec.sq", "prec.1mo", "prec.2mo", "prec.3mo", "prec.6mo", "prec.12mo",
                 "prec.days", "prec.gt.12", "prec.gt.25", "prec.gt.50", "prec.gt.76", "prec.max",
                 "RET", "RET.1mo", "RET.2mo", "RET.3mo", "RET.6mo", "RET.12mo",
                 "defc", "defc.1mo", "defc.2mo", "defc.3mo", "defc.6mo", "defc.12mo",
                 "defc.sq", "defc.1mo.sq", "defc.2mo.sq", "defc.3mo.sq", "defc.6mo.sq", "defc.12mo.sq",
                 "rads", "rads.1mo", "rads.2mo", "rads.3mo", "rads.6mo", "rads.12mo",
                 "tmin", "tmin.1mo", "tmin.2mo", "tmin.3mo", "tmin.6mo", "tmin.12mo",
                 "tmax", "tmax.1mo", "tmax.2mo", "tmax.3mo", "tmax.6mo", "tmax.12mo",
                 "ea", "ea.1mo", "ea.2mo", "ea.3mo", "ea.6mo", "ea.12mo",
                 "es", "es.1mo", "es.2mo", "es.3mo", "es.6mo", "es.12mo",
                 "VPD", "VPD.1mo", "VPD.2mo", "VPD.3mo", "VPD.6mo", "VPD.12mo",
                 "wspd", "wspd.1mo", "wspd.2mo", "wspd.3mo", "wspd.6mo", "wspd.12mo",
                 "relh", "relh.1mo", "relh.2mo", "relh.3mo", "relh.6mo", "relh.12mo")

## first: Pheasant Branch
yr.baseline.end <- 1995
flux.name.all <- c("discharge.mm", "baseflow.mm", "runoff.mm")

for (flux.name in flux.name.all){
  for (mo in mo.list){
    vars <- read.csv(paste0(git.dir, "Data/PheasantBranch/", flux.name, "/", yr.baseline.end, "/GHCN_MonthlyRegressions_vars.keep_", sprintf("%02d", mo), ".csv"),
                     stringsAsFactors=F)
    vars$month <- mo
    vars$flux.name <- flux.name
    
    if (exists("vars.PB")){
      vars.PB <- rbind(vars.PB, vars)
    } else {
      vars.PB <- vars
    }
  }
}

# make plot
vars.PB$var <- factor(vars.PB$var, levels=rev(var.options))
vars.PB$flux.name <- factor(vars.PB$flux.name, levels=c("discharge.mm", "baseflow.mm", "runoff.mm"))

# var.options which were never selected
var.missing <- var.options[!(var.options %in% unique(vars.PB$var))]
vars.PB <- rbind(vars.PB, data.frame(var=var.missing, month=NaN, flux.name="discharge.mm"))

p.PB <- 
  ggplot(vars.PB, aes(x=month, y=var)) +
  geom_point(color="#127D7D") +
  facet_grid(.~flux.name) +
  scale_x_continuous(name="Month", breaks=seq(1,12,1), labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme_bw() +
  theme(panel.border=element_rect(color="black"),
        panel.grid.minor=element_blank())

pdf(file=paste0(path.fig, "Figure_PCR_VarsKeep_PheasantBranch_NoText.pdf"), width=(100/25.4), height=(220/25.4))
p.PB + theme(text=element_blank())
dev.off()

## second: Yahara
## first: Pheasant Branch
flux.name.all <- c("aet", "drainage", "srunoff")

for (flux.name in flux.name.all){
  for (mo in mo.list){
    vars <- read.csv(paste0(git.dir, "Data/Yahara-AgroIBIS/", flux.name, "_GHCN_MonthlyRegressions_vars.keep_", sprintf("%02d", mo), ".csv"),
                     stringsAsFactors=F)
    vars$month <- mo
    vars$flux.name <- flux.name
    
    if (exists("vars.YW")){
      vars.YW <- rbind(vars.YW, vars)
    } else {
      vars.YW <- vars
    }
  }
}

# make plot
vars.YW$var <- factor(vars.YW$var, levels=rev(var.options))
vars.YW$flux.name <- factor(vars.YW$flux.name, levels=c("aet", "drainage", "srunoff"))

# var.options which were never selected
var.missing <- var.options[!(var.options %in% unique(vars.YW$var))]
vars.YW <- rbind(vars.YW, data.frame(var=var.missing, month=NaN, flux.name="aet"))

p.YW <- 
  ggplot(vars.YW, aes(x=month, y=var)) +
  geom_point(color="#127D7D") +
  facet_grid(.~flux.name) +
  scale_x_continuous(name="Month", breaks=seq(1,12,1), labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  theme_bw() +
  theme(panel.border=element_rect(color="black"),
        panel.grid.minor=element_blank())

pdf(file=paste0(path.fig, "Figure_PCR_VarsKeep_Yahara_NoText.pdf"), width=(100/25.4), height=(220/25.4))
p.YW + theme(text=element_blank())
dev.off()
