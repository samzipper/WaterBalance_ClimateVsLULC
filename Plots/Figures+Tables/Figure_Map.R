## Figure_Map.R
# This script is intended to plot the proportion of LULC in different 
# classes for the Yahara Watershed as a whole and the Pheasant Branch
# subwatershed.
#
# Data received by email from Eric Booth on 7/25/2017

rm(list=ls())

# path to local git repository
git.dir <- "C:/Users/Sam/WorkGits/WaterBalance_ClimateVsLULC/"

require(ggplot2)
require(reshape2)

# path to save figure output
path.fig <- paste0(git.dir, "Plots/Figures+Tables/")

# read in data
df <- read.csv(paste0(git.dir, "Data/LULC_2014_Yah_PhBr.csv"), stringsAsFactors=F)
df$LULC <- factor(df$LULC, levels=c("Agriculture", "Urban/Barren", "Natural", "Water"))

# melt
df.melt <- melt(df, id=c("LULC"))

p.prc <- 
  ggplot(df.melt, aes(x=LULC, y=100*value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0, color="gray65") +
  scale_fill_manual(values=c("Yahara"="#ff1d25", "PhBr"="#127D7D"), guide=F) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black"))

pdf(file=paste0(path.fig, "Figure_Map_Prc_NoText.pdf"), width=(82/25.4), height=(30/25.4))
p.prc + theme(text=element_blank(), plot.margin=unit(c(0,0,0,0), "mm"))
dev.off()
