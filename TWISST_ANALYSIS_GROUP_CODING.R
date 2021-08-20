rm(list=ls())
library(tidyverse)
library(ggpubr)
library(viridis)
library(reshape2)
library(wesanderson)
library(scales)
library(ggbeeswarm)
library(WVPlots)
setwd("~/Desktop/PhD_Labbook/Comparative_analysis/virilis_tree/")
source("plot_twisst.R")
setwd("~/Desktop/PhD_Labbook/Comparative_analysis/virilis_tree/TWISST/")


################# MONTANA LACICOLA BOREALIS ########################

####################################################################
############# BY CHROMOSOME ##########

########### X_chromosome

weights_file<-"Chrom_X.weights.csv.gz"
window_data_file <- "Chrom_X.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')
data.frame(smooth_twisst_melted)

smooth_twisst_melted %>%
  group_by(variable) %>%
  summarise(mean(value))

twisst_lm<-lm(value ~ variable, data = smooth_twisst_melted)
anova(twisst_lm)

twisst_aov<-aov(value ~ variable, data = smooth_twisst_melted)
TukeyHSD(twisst_aov)

twisst_X_bar <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_viridis(discrete = T, option='D')+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_X_weightings_violin.png', width=3, height=3, units='in', res=320)

twisst_X_bar
dev.off()

png('~/Desktop/PhD_Labbook/Comparative_analysis/virilis_tree/twisst_X_weightings_windows.png', width=6, height=2, units='in', res=320)
twisst_sliding_window_X<-ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_minimal()+theme(legend.position = 'none')+
  scale_colour_viridis(discrete = T, option='D')+ylab('Topology weight')+
  xlab('Genomic position (bp)')+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=12.5,face="bold"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))
twisst_sliding_window_X
dev.off()






########## Chromosome 2




weights_file<-"Chrom_2.weights.csv.gz"
window_data_file <- "Chrom_2.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_2_bar <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_viridis(discrete = T, option='D')+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_2_weightings_violin.png', width=3, height=3, units='in', res=320)
twisst_2_bar
dev.off()



png('twisst_2_weightings_windows.png', width=6, height=3, units='in', res=320)
twisst_sliding_window_2<-ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_pubr()+
  scale_colour_viridis(discrete = T, option='D')+ylab('Topology weight')+
  xlab('Genomic position (Mb)')+theme(legend.position = 'none')+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))
dev.off()
twisst_sliding_window_2







########## Chromosome 3



weights_file<-"Chrom_3.weights.csv.gz"
window_data_file <- "Chrom_3.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')

colMeans(smooth_twisst_df)

twisst_3_bar <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_viridis(discrete = T, option='D')+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_3_weightings_violin.png', width=3, height=3, units='in', res=320)
twisst_3_bar
dev.off()

png('twisst_3_weightings_windows.png', width=6, height=3, units='in', res=320)
twisst_sliding_window_3<-ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_pubr()+
  scale_colour_viridis(discrete = T, option='D')+ylab('Topology weight')+
  xlab('Genomic position (Mb)')+theme(legend.position = 'none')+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))
dev.off()
twisst_sliding_window_3







########## Chromosome 4




weights_file<-"Chrom_4.weights.csv.gz"
window_data_file <- "Chrom_4.phyml.data.tsv"

twisst_data_4 <- import.twisst(weights_files=weights_file,
                               window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data_4, lwd=3, cex=0.7)
twisst_data_smooth_4 <- smooth.twisst(twisst_data_4, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth_4$pos[[1]], twisst_data_smooth_4$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')

png('twisst_4_weightings_windows.png', width=6, height=3, units='in', res=320)
twisst_sliding_window_4<-ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_pubr()+
  scale_colour_viridis(discrete = T, option='D')+ylab('Topology weight')+
  xlab('Genomic position (Mb)')+theme(legend.position = 'none')+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))
dev.off()
twisst_sliding_window_4
twisst_4_bar <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_viridis(discrete = T, option='D')+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_4_weightings_violin.png', width=3, height=3, units='in', res=320)
twisst_4_bar
dev.off()








########## Chromosome 5
weights_file<-"Chrom_5.weights.csv.gz"
window_data_file <- "Chrom_5.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)

twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')

png('twisst_5_weightings_windows.png', width=6, height=3, units='in', res=320)
twisst_sliding_window_5<-ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_pubr()+
  scale_colour_viridis(discrete = T, option='D')+ylab('Topology weight')+
  xlab('Genomic position (Mb)')+theme(legend.position = 'none')+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))
dev.off()
twisst_sliding_window_5
twisst_5_bar <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_viridis(discrete = T, option='D')+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_5_weightings_violin.png', width=3, height=3, units='in', res=320)
twisst_5_bar
dev.off()

png('twisst_weightings_windows.png', width=12, height=8, units='in', res=320)
ggarrange(twisst_sliding_window_2, twisst_sliding_window_3, 
          twisst_sliding_window_4, twisst_sliding_window_5,
          twisst_sliding_window_X, nrow=5)
dev.off()



##############################################################
##############################################################
##############################################################
##############################################################
######################### NEXT CLADE #########################
################## EZOANA KANEKOI LITTORALIS ########################



######## CHROMOSOME 2
weights_file<-"Chrom_2_EKL.weights.csv.gz"
window_data_file <- "Chrom_2_EKL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_2_EKL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest1"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
#png('twisst_2_weightings_violin.png', width=3, height=3, units='in', res=320)
png('twisst_EKL_2.png', width=3, height=3, units='in', res=320)
twisst_2_EKL+theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()


####### CHROMOSOME 3

weights_file<-"Chrom_3_EKL.weights.csv.gz"
window_data_file <- "Chrom_3_EKL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_3_EKL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest1"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_EKL_3.png', width=3, height=3, units='in', res=320)
twisst_3_EKL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()


########### CHROMOSOME 4

weights_file<-"Chrom_4_EKL.weights.csv.gz"
window_data_file <- "Chrom_4_EKL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')

ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_pubr()+
  scale_colour_manual(values = wes_palette("GrandBudapest1"))+ylab('Topology weight')+
  xlab('Genomic position (Mb)')+theme(legend.position = 'none')+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))


twisst_4_EKL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest1"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_EKL_4.png', width=3, height=3, units='in', res=320)
twisst_4_EKL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()


########### CHROMOSOME 5

weights_file<-"Chrom_5_EKL.weights.csv.gz"
window_data_file <- "Chrom_5_EKL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')

as_tibble(phy_x)


twisst_5_EKL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest1"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_EKL_5.png', width=3, height=3, units='in', res=320)
twisst_5_EKL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()
########### X CHROMOSOME


weights_file<-"TWISST/Chrom_X_EKL.weights.csv.gz"
window_data_file <- "TWISST/Chrom_X_EKL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_X_EKL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest1"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_EKL_X.png', width=3, height=3, units='in', res=320)
twisst_X_EKL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))


png('twisst_ekl_X_sliding_window.png', width=6, height=2, units='in', res=320)
ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_minimal()+theme(legend.position = 'none')+
  scale_colour_manual(values = wes_palette("GrandBudapest1"))+ylab('Topology weight')+
  xlab('Genomic position (bp)')+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=12.5,face="bold"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))
dev.off()

###########################################################################
###########################################################################
################# NOVAMEXICANA AMERICANA LUMMEI ############################



######## CHROMOSOME 2
weights_file<-"Chrom_2_NAL.weights.csv.gz"
window_data_file <- "Chrom_2_NAL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_2_bar_NAL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest2"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_nal_2.png', width=3, height=3, units='in', res=320)
twisst_2_bar_NAL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()


####### CHROMOSOME 3

weights_file<-"Chrom_3_NAL.weights.csv.gz"
window_data_file <- "Chrom_3_NAL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_3_bar_NAL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest2"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_NAL_3.png', width=3, height=3, units='in', res=320)
twisst_3_bar_NAL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()


########### CHROMOSOME 4

weights_file<-"Chrom_4_NAL.weights.csv.gz"
window_data_file <- "Chrom_4_NAL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_4_bar_NAL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest2"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_NAL_4.png', width=3, height=3, units='in', res=320)
twisst_4_bar_NAL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()


########### CHROMOSOME 5

weights_file<-"Chrom_5_NAL.weights.csv.gz"
window_data_file <- "Chrom_5_NAL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_5_bar_NAL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest2"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_NAL_5.png', width=3, height=3, units='in', res=320)
twisst_5_bar_NAL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()
########### X CHROMOSOME


weights_file<-"TWISST/Chrom_X_NAL.weights.csv.gz"
window_data_file <- "TWISST/Chrom_X_NAL.phyml.data.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

weights <- weights / apply(weights, 1, sum)

#melt(windows_twisst, id.vars='rn')

plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, spacing = 10000)

smooth_twisst_df<-cbind(twisst_data_smooth$pos[[1]], twisst_data_smooth$weights[[1]])
names(smooth_twisst_df)[1]<-"position"
smooth_twisst_melted<-melt(smooth_twisst_df, id.vars='position')



twisst_X_bar_NAL <- ggplot(smooth_twisst_melted, aes(x=variable, y=value, fill=variable))+
  geom_violin()+scale_fill_manual(values = wes_palette("GrandBudapest2"))+theme_pubr()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  ylab('Topology weight')+xlab('Topology')
png('twisst_NAL_X.png', width=3, height=3, units='in', res=320)
twisst_X_bar_NAL+  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
dev.off()

png('twisst_NAL_X_sliding_window.png', width=6, height=2, units='in', res=320)
ggplot(smooth_twisst_melted, aes(x=position, y=value, colour=variable))+
  geom_line(size=1)+theme_minimal()+
  scale_colour_manual(values = wes_palette("GrandBudapest2"))+ylab('Topology weight')+
  xlab('Genomic position (bp)')+theme(legend.position = 'none')+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25))
dev.off()
