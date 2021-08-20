library(tidyverse)
library(ggpubr)
library(viridis)
library(reshape2)
library(wesanderson)
library(scales)
library(ggExtra)
library(cowplot)
library(ggridges)
library(qqman)
library(ggrepel)
library(GGally)
devtools::install_github("eclarke/ggbeeswarm")
install.packages(c("googleway", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
  library("sf")
library("rnaturalearth")
library("rnaturalearthdata")



source("plot_twisst.R")

vp=c('#440154FF','#29AF7FFF', '#95D840FF', '#FDE725FF')
dxy_colours=c('gray10', 'orange1', 'orangered2')
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

setwd("~/Desktop/PhD_Labbook/Comparative_analysis/virilis_tree")
#Get all the files read in
Montana_Lacicola_genic <- read.csv("species_pairs_popgen/mont_laci.genic.out.windows.h5.window_metrics.csv")
Novamexicana_Americana_genic <- read.csv("species_pairs_popgen/ameri_nov.genic.out.windows.h5.window_metrics.csv")
Kanekoi_Ezoana_genic <- read.csv('species_pairs_popgen/ezoana_kanekoi.genic.out.windows.h5.window_metrics.csv')
Borealis_Flavomontana_genic <- read.csv('species_pairs_popgen/bor_flavo.genic.out.windows.h5.window_metrics.csv')
chrom<-read.table("monCan3F9_chromosomes.txt",header=T)

###intergenic bits
Montana_Lacicola_intergenic <- read.csv("species_pairs_popgen/mont_laci.intergenic.out.windows.h5.window_metrics.csv")
Novamexicana_Americana_intergenic <- read.csv("species_pairs_popgen/ameri_nov.intergenic.out.windows.h5.window_metrics.csv")
Kanekoi_Ezoana_intergenic <- read.csv('species_pairs_popgen/ezoana_kanekoi.intergenic.out.windows.h5.window_metrics.csv')
Borealis_Flavomontana_intergenic <- read.csv('species_pairs_popgen/bor_flavo.intergenic.out.windows.h5.window_metrics.csv')



#Make sure chromosome information is attached to each for intergenic
Montana_Lacicola_intergenic$Chromosome <- chrom$Chromosome[ match(Montana_Lacicola_intergenic$sequence_id, chrom$Contig)]
Novamexicana_Americana_intergenic$Chromosome <- chrom$Chromosome[ match(Novamexicana_Americana_intergenic$sequence_id, chrom$Contig)]
Kanekoi_Ezoana_intergenic$Chromosome <- chrom$Chromosome[ match(Kanekoi_Ezoana_intergenic$sequence_id, chrom$Contig)]
Borealis_Flavomontana_intergenic$Chromosome <- chrom$Chromosome[ match(Borealis_Flavomontana_intergenic$sequence_id, chrom$Contig)]


###genic
Montana_Lacicola_genic$Chromosome <- chrom$Chromosome[ match(Montana_Lacicola_genic$sequence_id, chrom$Contig)]
Novamexicana_Americana_genic$Chromosome <- chrom$Chromosome[ match(Novamexicana_Americana_genic$sequence_id, chrom$Contig)]
Kanekoi_Ezoana_genic$Chromosome <- chrom$Chromosome[ match(Kanekoi_Ezoana_genic$sequence_id, chrom$Contig)]
Borealis_Flavomontana_genic$Chromosome <- chrom$Chromosome[ match(Borealis_Flavomontana_genic$sequence_id, chrom$Contig)]



#omit things that don't have chromosome data
Montana_Lacicola_genic<-na.omit(Montana_Lacicola_genic)
Novamexicana_Americana_genic<-na.omit(Novamexicana_Americana_genic)
Kanekoi_Ezoana_genic<-na.omit(Kanekoi_Ezoana_genic)
Borealis_Flavomontana_genic<-na.omit(Borealis_Flavomontana_genic)

#do the same for intergenic

Montana_Lacicola_intergenic<-na.omit(Montana_Lacicola_intergenic)
Novamexicana_Americana_intergenic<-na.omit(Novamexicana_Americana_intergenic)
Kanekoi_Ezoana_intergenic<-na.omit(Kanekoi_Ezoana_intergenic)
Borealis_Flavomontana_intergenic<-na.omit(Borealis_Flavomontana_intergenic)




###look at some the basic differences across chromosomes using boxplot
BOX_ML<-ggplot(Montana_Lacicola_genic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (M-L)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_NA<-ggplot(Novamexicana_Americana_genic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (N-A)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_KE<-ggplot(Kanekoi_Ezoana_genic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (K-E)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_BF<-ggplot(Borealis_Flavomontana_genic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (B-F)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')
png('dxy_chromosomes_ml_bf.png', width=5, height=2.5, units='in', res=320)
ggarrange(BOX_ML, BOX_BF, ncol = 2)
dev.off()
png('dxy_chromosomes_ke_na.png', width=5, height=2.5, units='in', res=320)
ggarrange(BOX_KE, BOX_NA, ncol = 2)
dev.off()
#Montana-Lacicola
ml_dxy_lm<-lm(dxy ~ Chromosome, data = Montana_Lacicola_genic)
ml_dxy_aov<-aov(dxy ~ Chromosome, data = Montana_Lacicola_genic)

anova(ml_dxy_lm)
TukeyHSD(ml_dxy_aov)

#Borealis_Flavomontana
bf_dxy_lm<-lm(dxy ~ Chromosome, data = Borealis_Flavomontana_genic)
anova(bf_dxy_lm)
bf_dxy_aov<-aov(dxy ~ Chromosome, data = Borealis_Flavomontana_genic)
TukeyHSD(bf_dxy_aov)


#Kanekoi_Ezoana
ke_dxy_lm<-lm(dxy ~ Chromosome, data = Kanekoi_Ezoana_genic)
anova(ke_dxy_lm)
ke_dxy_aov<-aov(dxy ~ Chromosome, data = Kanekoi_Ezoana_genic)
TukeyHSD(ke_dxy_aov)

#Novamexicana_americana
na_dxy_lm<-lm(dxy ~ Chromosome, data = Novamexicana_Americana_genic)
anova(na_dxy_lm)
na_dxy_aov<-aov(dxy ~ Chromosome, data = Novamexicana_Americana_genic)
TukeyHSD(na_dxy_aov)




#######Now is that the case for intergenic regions?
#Montana-Lacicola
ml_inter_dxy_lm<-lm(dxy ~ Chromosome, data = Montana_Lacicola_intergenic)
ml_inter_dxy_aov<-aov(dxy ~ Chromosome, data = Montana_Lacicola_intergenic)

anova(ml_inter_dxy_lm)
TukeyHSD(ml_inter_dxy_aov)

#Borealis_Flavomontana
bf_dxy_lm_inter<-lm(dxy ~ Chromosome, data = Borealis_Flavomontana_intergenic)
anova(bf_dxy_lm_inter)
bf_dxy_aov_inter<-aov(dxy ~ Chromosome, data = Borealis_Flavomontana_intergenic)
TukeyHSD(bf_dxy_aov_inter)


#Kanekoi_Ezoana
ke_dxy_lm<-lm(dxy ~ Chromosome, data = Kanekoi_Ezoana_intergenic)
anova(ke_dxy_lm)
ke_dxy_aov<-aov(dxy ~ Chromosome, data = Kanekoi_Ezoana_intergenic)
TukeyHSD(ke_dxy_aov)

#Novamexicana_americana
na_dxy_lm<-lm(dxy ~ Chromosome, data = Novamexicana_Americana_intergenic)
anova(na_dxy_lm)
na_dxy_aov<-aov(dxy ~ Chromosome, data = Novamexicana_Americana_intergenic)
TukeyHSD(na_dxy_aov)









#now plot by genomic position to detect outliers:
#Montana_Lacicola
dxyML_piM<-ggplot(Montana_Lacicola_genic, aes(x=pi_Montana, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+labs(y=expression('d'[XY]))+xlab('π (Montana)')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+stat_cor(method = "pearson", label.x = 0.008,label.y=0.005)

dxyML_piL<-ggplot(Montana_Lacicola_genic, aes(x=pi_Lacicola, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(y=expression('d'[XY]))+xlab('π (Lacicola)')+
  stat_cor(method = "pearson", label.x = 0.008, label.y=0.005)

#Borealis_Flavomontana
dxyBF_piB<-ggplot(Borealis_Flavomontana_genic, aes(x=pi_Borealis, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+
  theme(legend.position='none')+labs(y=expression('d'[XY]))+xlab('π (Borealis)')+
  stat_cor(method = "pearson", label.x = 0.01, label.y=0.01)+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

dxyBF_piF<-ggplot(Borealis_Flavomontana_genic, aes(x=pi_Flavomontana, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(y=expression('d'[XY]))+xlab('π (Flavomontana)')+
  stat_cor(method = "pearson", label.x = 0.008, label.y=0.01)+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

#Novamexicana_Americana
dxyNA_piN<-ggplot(Novamexicana_Americana_genic, aes(x=pi_Novamexicana, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(y=expression('d'[XY]))+xlab('π (Novamexicana)')+
  stat_cor(method = "pearson", label.x = 0.0035, label.y=0.005)+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

dxyNA_piA<-ggplot(Novamexicana_Americana_genic, aes(x=pi_Americana, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(y=expression('d'[XY]))+xlab('π (Americana)')+
  stat_cor(method = "pearson", label.x = 0.00275, label.y=0.005)+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

#Kanekoi_Ezoana
dxyKE_piK <-ggplot(Kanekoi_Ezoana_genic, aes(x=pi_Kanekoi, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(y=expression('d'[XY]))+xlab('π (Kanekoi)')+
  stat_cor(method = "pearson", label.x = 0.005, label.y=0.01)

dxyKE_piE<-ggplot(Kanekoi_Ezoana_genic, aes(x=pi_Ezoana, y=dxy))+
  geom_point(size=3, alpha=0.7)+theme_pubr()+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(y=expression('d'[XY]))+xlab('π (Ezoana)')+
  stat_cor(method = "pearson", label.x = 0.005, label.y=0.01)+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')


ggarrange(dxyML_piM, dxyML_piL, dxyBF_piB, dxyBF_piF ,nrow = 2, ncol=2)
ggarrange(dxyKE_piE, dxyKE_piK, dxyNA_piA, dxyNA_piN ,nrow = 2, ncol=2)


########################## Pi Correlation ##############


PI_ML<-ggplot(Montana_Lacicola_genic, aes(x=pi_Montana, y=pi_Lacicola))+
  geom_jitter(size=4, alpha=0.7)+ylab('π (Lacicola)')+xlab('π (Montana)')+
  scale_colour_viridis(discrete=T, option='D')+theme_pubr()+
  stat_cor(method = "pearson", label.x = 0.008, label.y=0.017)+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

PI_ML


xplot <- ggplot(Montana_Lacicola_genic, aes(x=pi_Montana, colour=Chromosome))+geom_density(size=1.5)+
  scale_colour_viridis(option='D', discrete=T)+theme_pubr()+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+xlab('π (Montana)')

yplot <- ggplot(Montana_Lacicola_genic, aes(x=pi_Lacicola, colour=Chromosome))+geom_density(size=1.5)+
  scale_colour_viridis(option='D', discrete=T)+theme_pubr()+rotate()+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+xlab('π (Lacicola)')


plot_grid(xplot, NULL, PI_ML, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))


ggplot(Borealis_Flavomontana_genic, aes(x=pi_Flavomontana, y=pi_Borealis, colour=Chromosome))+
  geom_jitter(size=4, alpha=0.7)+
  scale_colour_viridis(discrete=T, option='D')+theme_pubr()

ggplot(Novamexicana_Americana_genic, aes(x=pi_Americana, y=pi_Novamexicana, colour=Chromosome))+
  geom_jitter(size=4, alpha=0.7)+
  scale_colour_viridis(discrete=T, option='D')+theme_bw()

picorr_KE<-ggplot(Kanekoi_Ezoana_genic, aes(x=pi_Kanekoi, y=pi_Ezoana, colour=Chromosome))+
  geom_jitter(size=2.5)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  scale_colour_viridis(discrete=T, option='D')+theme_bw()+xlab('π (Kanekoi)')+
  ylab('π (Ezoana)')+theme(axis.title = element_text(face='bold'))+theme(legend.position = 'none')

picorr_ML<-ggplot(Montana_Lacicola_genic, aes(x=pi_Montana, y=pi_Lacicola, colour=Chromosome))+
  geom_jitter(size=2.5)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  scale_colour_viridis(discrete=T, option='D')+theme_bw()+xlab('π (Montana)')+
  ylab('π (Lacicola)')+theme(axis.title = element_text(face='bold'))+theme(legend.position = 'none')

picorr_NA<-ggplot(Novamexicana_Americana_genic, aes(x=pi_Novamexicana, y=pi_Americana, colour=Chromosome))+
  geom_jitter(size=2.5)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  scale_colour_viridis(discrete=T, option='D')+theme_bw()+xlab('π (Novamexicana)')+
  ylab('π (Americana)')+theme(axis.title = element_text(face='bold'))+theme(legend.position = 'none')

picorr_BF<-ggplot(Borealis_Flavomontana_genic, aes(x=pi_Borealis, y=pi_Flavomontana, colour=Chromosome))+
  geom_jitter(size=2.5)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  scale_colour_viridis(discrete=T, option='D')+theme_bw()+xlab('π (Borealis)')+
  ylab('π (Flavomontana)')+theme(axis.title = element_text(face='bold'))+theme(legend.title = element_text(face='bold'))

png('pi_corr_altogether.png', width=6, height=5, units='in', res=320)
ggarrange(picorr_ML, picorr_BF, picorr_KE, picorr_NA, ncol=2, nrow=2)
dev.off()

##################DXY CODING GENOME-WIDE #################

## MONTANA LACICOLA

meanData_ML = Montana_Lacicola_genic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_ML

quants <- quantile(Montana_Lacicola_genic$dxy, c(0.95, 0.99))
Montana_Lacicola_genic$quant  <- with(Montana_Lacicola_genic, factor(ifelse(dxy < quants[1], 0, 
                                      ifelse(dxy < quants[2], 1, 2))))



png('genomewide_dxy_ml.png', width=6, height=2, units='in', res=320)
ggplot(Montana_Lacicola_genic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant))+
  geom_hline(data=meanData_ML, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (M-L)))+
  theme_minimal_hgrid()+ylim(0.0, 0.045)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()


#### BOREALIS flavomontana

meanData_bf = Borealis_Flavomontana_genic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_bf

quants_bf <- quantile(Borealis_Flavomontana_genic$dxy, c(0.95, 0.99))
Borealis_Flavomontana_genic$quant_bf  <- with(Borealis_Flavomontana_genic, factor(ifelse(dxy < quants_bf[1], 0, 
                                                                            ifelse(dxy < quants_bf[2], 1, 2))))

png('genomewide_dxy_bf.png', width=6, height=2, units='in', res=320)
ggplot(Borealis_Flavomontana_genic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_bf))+
  geom_hline(data=meanData_bf, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (B-F)))+
  theme_minimal_hgrid()+ylim(0.0, 0.045)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()

###### kanekoi and ezoana

meanData_ke = Kanekoi_Ezoana_genic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_ke

quants_ke <- quantile(Kanekoi_Ezoana_genic$dxy, c(0.95, 0.99))
Kanekoi_Ezoana_genic$quant_ke  <- with(Kanekoi_Ezoana_genic, factor(ifelse(dxy < quants_ke[1], 0, 
                                                                                         ifelse(dxy < quants_ke[2], 1, 2))))

png('genomewide_dxy_ke.png', width=6, height=2, units='in', res=320)
ggplot(Kanekoi_Ezoana_genic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_ke))+
  geom_hline(data=meanData_ke, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (K-E)))+
  theme_minimal_hgrid()+ylim(0.02, 0.07)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()


###### NOVAMEXICANA AMERICANA

meanData_na = Novamexicana_Americana_genic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_na

quants_na <- quantile(Novamexicana_Americana_genic$dxy, c(0.95, 0.99))
Novamexicana_Americana_genic$quant_na <- with(Novamexicana_Americana_genic, factor(ifelse(dxy < quants_na[1], 0, 
                                                                           ifelse(dxy < quants_na[2], 1, 2))))

png('genomewide_dxy_na.png', width=6, height=2, units='in', res=320)
ggplot(Novamexicana_Americana_genic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_na))+
  geom_hline(data=meanData_na, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (N-A)))+
  theme_minimal_hgrid()+ylim(0.0, 0.04)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()



###################### OUTLIER DXY BY CHROMOSOME -- GENIC #################

ML_dxy_outlier <- filter(Montana_Lacicola_genic, quant==1 | quant==2 )
BF_dxy_outlier <- filter(Borealis_Flavomontana_genic, quant_bf==1 | quant_bf==2 )
KE_dxy_outlier <- filter(Kanekoi_Ezoana_genic, quant_ke==1 | quant_ke==2 )
NA_dxy_outlier <- filter(Novamexicana_Americana_genic, quant_na==1 | quant_na==2 )

levels(ML_dxy_outlier$Chromosome) <- c(levels(ML_dxy_outlier$Chromosome), "A") 
ML_dxy_outlier$Chromosome[ML_dxy_outlier$Chromosome==2]  <- "A"

levels(BF_dxy_outlier$Chromosome) <- c(levels(BF_dxy_outlier$Chromosome), "A") 
BF_dxy_outlier$Chromosome[BF_dxy_outlier$Chromosome==5]  <- "A"

levels(KE_dxy_outlier$Chromosome) <- c(levels(KE_dxy_outlier$Chromosome), "A") 
KE_dxy_outlier$Chromosome[KE_dxy_outlier$Chromosome==2]  <- "A"

levels(NA_dxy_outlier$Chromosome) <- c(levels(NA_dxy_outlier$Chromosome), "A") 
NA_dxy_outlier$Chromosome[NA_dxy_outlier$Chromosome==2]  <- "A"

ggplot(ML_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+
  labs(y=expression('d'[XY] (M-L)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(BF_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+
  labs(y=expression('d'[XY] (B-F)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(KE_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+
  labs(y=expression('d'[XY] (K-E)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(NA_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+
  labs(y=expression('d'[XY] (N-A)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')



kruskal.test(dxy ~ Chromosome, data=ML_dxy_outlier)
kruskal.test(dxy ~ Chromosome, data=BF_dxy_outlier)
kruskal.test(dxy ~ Chromosome, data=KE_dxy_outlier)
kruskal.test(dxy ~ Chromosome, data=NA_dxy_outlier)


##################DXY INTERGENIC GENOME-WIDE #################

## MONTANA LACICOLA


BOX_INTER_ML<-ggplot(Montana_Lacicola_intergenic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (M-L)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_INTER_BF<-ggplot(Borealis_Flavomontana_intergenic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (B-F)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_INTER_KE<-dxy_inter_box_ke<-ggplot(Kanekoi_Ezoana_intergenic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (K-E)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_INTER_NA<-ggplot(Novamexicana_Americana_intergenic, aes(x=Chromosome, y=dxy, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (N-A)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

png('dxy_INTERBOX_ml_bf.png', width=5, height=2.5, units='in', res=320)
ggarrange(BOX_INTER_ML, BOX_INTER_BF, ncol = 2)
dev.off()

png('dxy_INTERBOX_ke_na.png', width=5, height=2.5, units='in', res=320)
ggarrange(BOX_INTER_KE, BOX_INTER_NA, ncol = 2)
dev.off()




############### GENOMEWIDE ###############
meanData_ML_intergenic = Montana_Lacicola_intergenic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_ML

quants_ml_inter <- quantile(Montana_Lacicola_intergenic$dxy, c(0.95, 0.99))
Montana_Lacicola_intergenic$quant_ml_inter  <- with(Montana_Lacicola_intergenic, factor(ifelse(dxy < quants_ml_inter[1], 0, 
                                                                            ifelse(dxy < quants_ml_inter[2], 1, 2))))

png('genomewide_dxy_ml_inter.png', width=6, height=2, units='in', res=320)
dxy_gw_ml_inter<-ggplot(Montana_Lacicola_intergenic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_ml_inter))+
  geom_hline(data=meanData_ML_intergenic, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (M-L)))+
  theme_minimal_hgrid()+ylim(0.0, 0.045)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dxy_gw_ml_inter
dev.off()

#### BOREALIS flavomontana

meanData_bf_inter = Borealis_Flavomontana_intergenic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_bf_inter

quants_bf_inter <- quantile(Borealis_Flavomontana_intergenic$dxy, c(0.95, 0.99))
Borealis_Flavomontana_intergenic$quant_bf_inter  <- with(Borealis_Flavomontana_intergenic, factor(ifelse(dxy < quants_bf_inter[1], 0, 
                                                                                         ifelse(dxy < quants_bf_inter[2], 1, 2))))

png('genomewide_dxy_bf_inter.png', width=6, height=2, units='in', res=320)
dxy_gw_inter_bf<-ggplot(Borealis_Flavomontana_intergenic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_bf_inter))+
  geom_hline(data=meanData_bf_inter, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (B-F)))+
  theme_minimal_hgrid()+ylim(0.0, 0.05)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dxy_gw_inter_bf
dev.off()

###### kanekoi and ezoana

meanData_ke_inter = Kanekoi_Ezoana_intergenic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_ke_inter

quants_ke_inter <- quantile(Kanekoi_Ezoana_intergenic$dxy, c(0.95, 0.99))
Kanekoi_Ezoana_intergenic$quant_ke_inter  <- with(Kanekoi_Ezoana_intergenic, factor(ifelse(dxy < quants_ke_inter[1], 0, 
                                                                           ifelse(dxy < quants_ke_inter[2], 1, 2))))

png('genomewide_dxy_ke_inter.png', width=6, height=2, units='in', res=320)
dxy_gw_kl_inter<-ggplot(Kanekoi_Ezoana_intergenic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_ke_inter))+
  geom_hline(data=meanData_ke_inter, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (K-E)))+
  theme_minimal_hgrid()+ylim(0.01, 0.07)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dxy_gw_kl_inter
dev.off()


###### NOVAMEXICANA AMERICANA

meanData_na_inter = Novamexicana_Americana_intergenic %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy)) 

meanData_na_inter

quants_na_inter <- quantile(Novamexicana_Americana_intergenic$dxy, c(0.95, 0.99))
Novamexicana_Americana_intergenic$quant_na_inter <- with(Novamexicana_Americana_intergenic, factor(ifelse(dxy < quants_na_inter[1], 0, 
                                                                                          ifelse(dxy < quants_na_inter[2], 1, 2))))

png('genomewide_dxy_na_inter.png', width=6, height=2, units='in', res=320)
dxy_na_inter<-ggplot(Novamexicana_Americana_intergenic, aes(x=sequence_id, y=dxy))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_na_inter))+
  geom_hline(data=meanData_na_inter, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (N-A)))+
  theme_minimal_hgrid()+ylim(0.0, 0.04)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dxy_na_inter
dev.off()


png('genomewide_dxy_inter.png', width=7, height=9, units='in', res=320)
ggarrange(dxy_gw_ml_inter,dxy_gw_inter_bf,dxy_gw_kl_inter, dxy_na_inter, nrow=4)
dev.off()
############# OUTLIER DXY REGIONS FOR INTERGENIC WINDOWS ########

ML_inter_dxy_outlier <- filter(Montana_Lacicola_intergenic, quant_ml_inter==1 | quant_ml_inter==2 )
BF_inter_dxy_outlier <- filter(Borealis_Flavomontana_intergenic, quant_bf_inter==1 | quant_bf_inter==2 )
KE_inter_dxy_outlier <- filter(Kanekoi_Ezoana_intergenic, quant_ke_inter==1 | quant_ke_inter==2 )
NA_inter_dxy_outlier <- filter(Novamexicana_Americana_intergenic, quant_na_inter==1 | quant_na_inter==2 )



levels(ML_inter_dxy_outlier$Chromosome) <- c(levels(ML_inter_dxy_outlier$Chromosome), "A") 
ML_inter_dxy_outlier$Chromosome[ML_inter_dxy_outlier$Chromosome==5]  <- "A"

levels(BF_inter_dxy_outlier$Chromosome) <- c(levels(BF_inter_dxy_outlier$Chromosome), "A") 
BF_inter_dxy_outlier$Chromosome[BF_inter_dxy_outlier$Chromosome==2]  <- "A"

levels(KE_inter_dxy_outlier$Chromosome) <- c(levels(KE_inter_dxy_outlier$Chromosome), "A") 
KE_inter_dxy_outlier$Chromosome[KE_inter_dxy_outlier$Chromosome==5]  <- "A"

levels(NA_inter_dxy_outlier$Chromosome) <- c(levels(NA_inter_dxy_outlier$Chromosome), "A") 
NA_inter_dxy_outlier$Chromosome[NA_inter_dxy_outlier$Chromosome==5]  <- "A"




ggplot(ML_inter_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (M-L)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(BF_inter_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (B-F)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(KE_inter_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (K-E)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(NA_inter_dxy_outlier, aes(x=Chromosome, y=dxy, fill=Chromosome))+theme_pubr()+
  geom_boxplot(width=.3, position='dodge', coef = 6, alpha=0.9)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (N-A)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')


kruskal.test(dxy ~ Chromosome, data=ML_inter_dxy_outlier)
kruskal.test(dxy ~ Chromosome, data=BF_inter_dxy_outlier)
kruskal.test(dxy ~ Chromosome, data=KE_inter_dxy_outlier)
kruskal.test(dxy ~ Chromosome, data=NA_inter_dxy_outlier)


##################### Dmin analyses #############
#Borealis_Flavo_Montana
BFM_fdm<-read.table("BF_M_100_200_fdm.txt",header=T)
BFM_fdm$Chromosome <- chrom$Chromosome[ match(BFM_fdm$chr, chrom$Contig)]
BFM_fdm <- na.omit(BFM_fdm)

BFM_fdm_plot<-ggplot(BFM_fdm, aes(x=f_dM, y=Chromosome, fill=Chromosome))+
  geom_density_ridges2(alpha=0.7, scale=1.7)+
  theme_bw()+geom_vline(xintercept = 0.0, linetype='dashed')+
  coord_cartesian(clip = "off")+
  scale_y_discrete(expand = expansion(add = c(0.2, 2)))+
  scale_fill_viridis(discrete=T, option='D')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(x=expression('f'[dm]))

###Montana_Lacicola_Flavomontana
MLF_fdm<-read.table("ML_F_100_200.txt",header=T)
MLF_fdm$Chromosome <- chrom$Chromosome[ match(MLF_fdm$chr, chrom$Contig)]
MLF_fdm <- na.omit(MLF_fdm)

MLF_fdm_plot<-ggplot(MLF_fdm, aes(x=f_dM, y=Chromosome, fill=Chromosome))+
  geom_density_ridges2(alpha=0.7, scale=1.7)+
  theme_bw()+geom_vline(xintercept = 0.0, linetype='dashed')+
  coord_cartesian(clip = "off")+
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5)))+
  scale_fill_viridis(discrete=T, option='D')+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(x=expression('f'[dm]))

MLF_fdm_plot

#Kanekoi_Ezoana_Littoralis
KEL_fdm<-read.table("KE_L_100_200.txt",header=T)
KEL_fdm$Chromosome <- chrom$Chromosome[ match(KEL_fdm$chr, chrom$Contig)]
KEL_fdm <- na.omit(KEL_fdm)

kel_fdm_plot<-ggplot(KEL_fdm, aes(x=f_dM, y=Chromosome, fill=Chromosome))+
  geom_density_ridges2(alpha=0.7, scale=1.7)+
  theme_bw()+geom_vline(xintercept = 0.0, linetype='dashed')+
  coord_cartesian(clip = "off")+
  scale_y_discrete(expand = expansion(add = c(0.2, 2)))+
  scale_fill_viridis(discrete=T, option='D')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(x=expression('f'[dm]))

NAL_fdm<-read.table("NA_L_100_200.txt",header=T)
NAL_fdm$Chromosome <- chrom$Chromosome[ match(NAL_fdm$chr, chrom$Contig)]
NAL_fdm <- na.omit(NAL_fdm)

nal_fdm_plot<-ggplot(NAL_fdm, aes(x=f_dM, y=Chromosome, fill=Chromosome))+
  geom_density_ridges2(alpha=0.7, scale=1.7)+
  theme_bw()+geom_vline(xintercept = 0.0, linetype='dashed')+
  coord_cartesian(clip = "off")+
  scale_y_discrete(expand = expansion(add = c(0.2, 2)))+
  scale_fill_viridis(discrete=T, option='D')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')+labs(x=expression('f'[dm]))

png('distributions_fdm_chromosomes_bfm.png', width=3, height=3, units='in', res=320)
BFM_fdm_plot
dev.off()


mid<-rowMeans(BFM_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
BFM_fdm <- cbind(BFM_fdm, mid)
mid<-rowMeans(MLF_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
MLF_fdm <- cbind(MLF_fdm, mid)
mid<-rowMeans(KEL_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
KEL_fdm <- cbind(KEL_fdm, mid)
mid<-rowMeans(NAL_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
NAL_fdm <- cbind(NAL_fdm, mid)


# ggplot(data = BFM_fdm, aes(x =Chromosome, y = f_dM, colour=Chromosome), group=mid)+
#   geom_jitter(size=3,alpha=0.7)+geom_hline(yintercept = 0.0, 
#                                            linetype='dashed')+
#   theme_bw()+scale_colour_viridis(option='D', discrete = T)
# 
# 
# ggplot(data = KEL_fdm, aes(x =Chromosome, y = f_dM, colour=Chromosome), group=mid)+
#   geom_jitter(size=3,alpha=0.7)+geom_hline(yintercept = 0.0, 
#                                            linetype='dashed')+
#   theme_bw()+scale_colour_viridis(option='D', discrete = T)




#now filter the fdm dataset.
BFM_fdm <- filter(BFM_fdm, !f_dM==0)
BFM_fdm$Introgression <- ifelse(BFM_fdm$f_dM >=0, "FM", "BM") ##assign labels to direction based on fdm.



rm(list = ls())
####FOR BFM
#get mid positions
mid<-rowMeans(BFM_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
BFM_fdm <- cbind(BFM_fdm, mid)
mid<-rowMeans(MLF_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
MLF_fdm <- cbind(MLF_fdm, mid)
mid<-rowMeans(KEL_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
KEL_fdm <- cbind(KEL_fdm, mid)
mid<-rowMeans(NAL_fdm[c('windowStart', 'windowEnd')], na.rm=TRUE)
NAL_fdm <- cbind(NAL_fdm, mid)
#get scaffold length
chrom<-chrom[order(chrom$Contig),]
chrom$add<-c(0,cumsum(chrom$End)[1:72])
#now add scaffold length to the relevant dataset(s)
BFM_fdm$Add <- chrom$add[match(BFM_fdm$chr, chrom$Contig)]
MLF_fdm$Add <- chrom$add[match(MLF_fdm$chr, chrom$Contig)]
KEL_fdm$Add <- chrom$add[match(KEL_fdm$chr, chrom$Contig)]
NAL_fdm$Add <- chrom$add[match(NAL_fdm$chr, chrom$Contig)]
#Next step is then add chromosome lengths to mid to get a relative genomic position.
BFM_fdm$added <- BFM_fdm$mid+BFM_fdm$Add
MLF_fdm$added <- MLF_fdm$mid+MLF_fdm$Add
KEL_fdm$added <- KEL_fdm$mid+KEL_fdm$Add
NAL_fdm$added <- NAL_fdm$mid+NAL_fdm$Add


aggregate(MLF_fdm[, 4], list(MLF_fdm$Chromosome), mean)
aggregate(BFM_fdm[, 4], list(BFM_fdm$Chromosome), mean)
aggregate(KEL_fdm[, 4], list(KEL_fdm$Chromosome), mean)
aggregate(NAL_fdm[, 4], list(NAL_fdm$Chromosome), mean)

ml_fdm_aov<-aov(f_dM ~ Chromosome, data = MLF_fdm)
bf_fdm_aov<-aov(f_dM ~ Chromosome, data = BFM_fdm)
ke_fdm_aov<-aov(f_dM ~ Chromosome, data = KEL_fdm)
na_fdm_aov<-aov(f_dM ~ Chromosome, data = NAL_fdm)

anova(ml_fdm_aov)
anova(bf_fdm_aov)
anova(ke_fdm_aov)
anova(na_fdm_aov)



TukeyHSD(ml_fd_aov)




##### plotting genomewide fdm ######

png('BFM_fdm.png', width=6, height=2, units='in', res=320)
ggplot(BFM_fdm, aes(x=chr, y=f_dM, colour=f_dM))+
  geom_jitter(size=1.75,alpha=0.7)+
  geom_hline(yintercept = 0.0)+scale_colour_viridis(option='D', discrete=F)+
  labs(y=expression('f'[dm] (B-F-M)))+
  theme_minimal_hgrid()+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()

png('MLF_fdm.png', width=6, height=2, units='in', res=320)
ggplot(MLF_fdm, aes(x=chr, y=f_dM, colour=f_dM))+
  geom_jitter(size=1.75,alpha=0.7)+
  geom_hline(yintercept = 0.0)+scale_colour_viridis(option='D', discrete=F)+
  labs(y=expression('f'[dm] (M-L-F)))+
  theme_minimal_hgrid()+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()


png('KEL_fdm.png', width=6, height=2, units='in', res=320)
ggplot(KEL_fdm, aes(x=chr, y=f_dM, colour=f_dM))+
  geom_jitter(size=1.75,alpha=0.7)+
  geom_hline(yintercept = 0.0)+scale_colour_viridis(option='D', discrete=F)+
  labs(y=expression('f'[dm] (K-E-L)))+
  theme_minimal_hgrid()+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()

png('NAL_fdm.png', width=6, height=2, units='in', res=320)
ggplot(NAL_fdm, aes(x=chr, y=f_dM, colour=f_dM))+
  geom_jitter(size=1.75,alpha=0.7)+
  geom_hline(yintercept = 0.0)+scale_colour_viridis(option='D', discrete=F)+
  labs(y=expression('f'[dm] (N-A-L)))+
  theme_minimal_hgrid()+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))
dev.off()


#Now plot BFM
png('BFM_fdm.png', width=6, height=2, units='in', res=320)
ggplot(BFM_fdm, aes(x=added/1000000, y=f_dM, colour=f_dM))+
  geom_jitter(size=2)+
  scale_colour_viridis(discrete = F, option='D', alpha = 1)+
  xlab('Genomic Position (MB)')+labs(y=expression('f'[dm]))+
  theme_bw()+geom_hline(yintercept = 0.0, linetype='dashed')+
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=11,face="bold"))+
  labs(color='')
dev.off()



png('MLF_fdm.png', width=6, height=2, units='in', res=320)
ggplot(MLF_fdm, aes(x=added/1000000, y=f_dM, colour=f_dM))+
  geom_jitter(size=2)+
  scale_colour_viridis(discrete = F, option='D', alpha = 1)+
  xlab('Genomic Position (MB)')+labs(y=expression('f'[dm]))+
  theme_bw()+geom_hline(yintercept = 0.0, linetype='dashed')+
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=11,face="bold"))+
  labs(color='')
dev.off()

png('KEL_fdm.png', width=6, height=2, units='in', res=320)
ggplot(KEL_fdm, aes(x=added/1000000, y=f_dM, colour=f_dM))+
  geom_jitter(size=2)+
  scale_colour_viridis(discrete = F, option='D', alpha = 1)+
  xlab('Genomic Position (MB)')+labs(y=expression('f'[dm]))+
  theme_bw()+geom_hline(yintercept = 0.0, linetype='dashed')+
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=11,face="bold"))+
  labs(color='')
dev.off()

png('NAL_fdm.png', width=6, height=2, units='in', res=320)
ggplot(NAL_fdm, aes(x=added/1000000, y=f_dM, colour=f_dM))+
  geom_jitter(size=2)+
  scale_colour_viridis(discrete = F, option='D', alpha = 1)+
  xlab('Genomic Position (MB)')+labs(y=expression('f'[dm]))+
  theme_bw()+geom_hline(yintercept = 0.0, linetype='dashed')+
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=11,face="bold"))+
  labs(color='')
dev.off()

KEL_fdm$Introgression <- ifelse(KEL_fdm$f_dM >=0, "EL", "KL")
MLF_fdm$Introgression <- ifelse(MLF_fdm$f_dM >=0, "LF", "MF")
BFM_fdm$Introgression <- ifelse(BFM_fdm$f_dM >=0, "FM", "BM")
NAL_fdm$Introgression <- ifelse(NAL_fdm$f_dM >=0, "AL", "NL")

pie_colours=c('#FF6347','#4169E1')
NAL_pie<-count(NAL_fdm, Introgression)
KEL_pie<-count(KEL_fdm, Introgression)
MLF_pie<-count(MLF_fdm, Introgression)
BFM_pie<-count(BFM_fdm, Introgression)

png('BFM_pie.png', width=2, height=2, units='in', res=320)
ggplot(BFM_pie, aes(x="", y=n, fill=Introgression))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+
  blank_theme+theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid  = element_blank())+
  scale_fill_manual(values=pie_colours)
dev.off()

png('MLF_pie.png', width=2, height=2, units='in', res=320)
ggplot(MLF_pie, aes(x="", y=n, fill=Introgression))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+
  blank_theme+theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid  = element_blank())+
  scale_fill_manual(values=pie_colours)
dev.off()

png('KEL_pie.png', width=2, height=2, units='in', res=320)
ggplot(KEL_pie, aes(x="", y=n, fill=Introgression))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+
  blank_theme+theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid  = element_blank())+
  scale_fill_manual(values=pie_colours)
dev.off()

png('NAL_pie.png', width=2, height=2, units='in', res=320)
ggplot(NAL_pie, aes(x="", y=n, fill=Introgression))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+
  blank_theme+theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid  = element_blank())+
  scale_fill_manual(values=pie_colours)
dev.off()


ML_scaf00003<-Montana_Lacicola_genic %>%
  filter(sequence_id == 'monCan3F9.00003')
MLF_Fdm_scaf00003<-MLF_fdm %>%
  filter(chr == 'monCan3F9.00003')
BF_scaf00003<-Borealis_Flavomontana_genic %>%
  filter(sequence_id == 'monCan3F9.00003')
EK_scaf00003<-Kanekoi_Ezoana_genic %>%
  filter(sequence_id == 'monCan3F9.00003')
NA_scaf00003<-Novamexicana_Americana_genic %>%
  filter(sequence_id == 'monCan3F9.00003')


ggplot(ML_scaf00003, aes(x=rel_pos/1000000, y=dxy))+
  geom_point(size=4, alpha=0.7)+theme_bw()+
  xlab('Genomic Position (Mbp)')+labs(y=expression('d'[XY] (M-L)))

png('nocte_introgress.png', width=3, height=2, units='in', res=320)
ggplot(MLF_Fdm_scaf00003, aes(x=mid/1000000, y=f_dM))+
  geom_point(size=4.5, alpha=0.9)+theme_minimal()+
  geom_hline(yintercept = 0.0,
                                            colour='red', linetype='dashed')+
  xlab('Genomic Position (Mb)')+labs(y=expression('f'[dm]))+
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=12.5,face="bold"))
dev.off()

ggplot(ML_scaf00003, aes(x=centre/1000000))+
  geom_line(data=ML_scaf00003, aes(y=pi_Montana), colour='red', size=3, alpha=0.7)+
  geom_line(data=ML_scaf00003, aes(y=pi_Lacicola), colour='blue', size=3, alpha=0.7)+
  theme_bw()+ylab('pi')

ggplot(BF_scaf00003, aes(x=centre/1000000))+
  geom_line(data=BF_scaf00003, aes(y=pi_Borealis), colour='red', size=3, alpha=0.7)+
  geom_line(data=BF_scaf00003, aes(y=pi_Flavomontana), colour='blue', size=3, alpha=0.7)+
  theme_bw()+ylab('pi')





#DOES THIS WORK?
FST <- read.table("monCan3F.posselectreference.fstsliding.txt", header=T)
#filter dataset to only include reliable regions
FST <- filter(FST, min_cov > 25)
FST <- filter(FST, snps > 5)


#FST$fdm <- ifelse(sapply(FST$pos, function(p) 
 # any(BFM_fdm$windowStart <= p & BFM_fdm$windowEnd >= p)), BFM_fdm$f_dM, 0)
#FST$Chromosome <- chrom$Chromosome[ match(FST$Scaffold, chrom$Contig)]
#FST <- na.omit(FST)
#FST$Introgression <- ifelse(FST$fdm >=0, "FM", "BM")

fst_aggregate<-aggregate(FST[, 6:9], list(FST$Scaffold), mean)

fst_agg_melt<- melt(fst_aggregate, id.vars=c('Group.1', 'fdm'), measure.vars=c('Col_Van', 
                                                                               'Van_Oul', 'Col_Oul'))

fst_agg_melt$Chromosome <- chrom$Chromosome[ match(fst_agg_melt$Group.1, chrom$Contig)]
fst_agg_melt <- na.omit(fst_agg_melt)

ggplot(fst_agg_melt, aes(x=fdm,y=value, colour=variable, shape=variable))+
  geom_jitter(alpha=0.7, size=2)+
  scale_colour_viridis(option='D', discrete=T)+
  theme_bw()+labs(x=expression('f'[dm]))+ylab('FST')+
  geom_vline(xintercept = 0.0, linetype='dashed')



########## for top 20 regions in fdm 

NAL_fdm_top_10 <- NAL_fdm %>%
  select(chr,windowStart,windowEnd,f_dM) %>%
  top_n(10) 

BFM_topandbottom_introgressed<-merge(BFM_fdm_top_10,BFM_fdm_bottom_10,
      by=c('chr','windowStart','windowEnd','f_dM'), all=T)

MLF_topandbottom_introgressed<-merge(MLF_fdm_top_10,MLF_fdm_bottom_10,
                                     by=c('chr','windowStart','windowEnd','f_dM'), all=T)
KEL_topandbottom_introgressed<-merge(KEL_fdm_top_10,KEL_fdm_bottom_10,
                                     by=c('chr','windowStart','windowEnd','f_dM'), all=T)
NAL_topandbottom_introgressed<-merge(NAL_fdm_top_10,NAL_fdm_bottom_10,
                                     by=c('chr','windowStart','windowEnd','f_dM'), all=T)


write.table(BFM_topandbottom_introgressed, file = "BFM_topandbottom_introgressed.txt", append = FALSE, quote = TRUE, sep = "\t")
write.table(MLF_topandbottom_introgressed, file = "MLF_topandbottom_introgressed.txt", append = FALSE, quote = TRUE, sep = "\t")
write.table(KEL_topandbottom_introgressed, file = "KEL_topandbottom_introgressed.txt", append = FALSE, quote = TRUE, sep = "\t")
write.table(NAL_topandbottom_introgressed, file = "NAL_topandbottom_introgressed.txt", append = FALSE, quote = TRUE, sep = "\t")


##############################
##KEL

MLF_GO_significant$trio <- "MLF"
BFM_GO_significant$trio <- "BFM"
KEL_GO_significant$trio <- "KEL"
NAL_GO_significant$trio <- "NAL"



Altogether_significant_GO<-Reduce(function(...) merge(..., all=TRUE), list(MLF_GO_significant, BFM_GO_significant, 
                                                KEL_GO_significant, NAL_GO_significant))
introgressed_newterms<-sapply(strsplit(Altogether_significant_GO$Term, "_GO"),head, 1)
Altogether_significant_GO <- cbind(Altogether_significant_GO, introgressed_newterms)
Altogether_significant_GO$introgressed_newterms<-gsub( "_", " ", as.character(Altogether_significant_GO$introgressed_newterms))
Altogether_significant_GO$trio <- factor(Altogether_significant_GO$trio, levels = Altogether_significant_GO$trio)



png('Altogether.png', width=8.1, height=8, units='in', res=320)
ggplot(Altogether_significant_GO, aes(x=Combined_Score, y=trio, group=introgressed_newterms,
                               label=introgressed_newterms, fill=Adjusted_P.value))+
  geom_bar(stat = 'identity', position="dodge", width=0.9, size=0.5, alpha=0.8, colour='black')+
  scale_y_discrete(expand = c(0, 0))+scale_fill_viridis(option='D')+
  labs(fill = "FDR")+theme(legend.text = element_text(face='bold'))+theme_classic()+
  scale_x_continuous(limits = c(0,250), expand = c(0, 0))+
  ylab('GO Biological Terms')+xlab('Combined Score')+
  geom_text(size=3.5, angle=0, hjust=-0.06 ,colour='black', fontface='bold', position = position_dodge(width = 0.91))+
  theme(axis.text.y = element_blank())+theme(axis.ticks.y = element_blank())+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"))
dev.off()


#################################





KEL_GO_introgressed <- read.table('Introgression/GO_BIOLOGICAL_KEL_introgressed.txt', header=T, sep='', fill=T)
KEL_GO_introgressed$Term<-as.character(KEL_GO_introgressed$Term)
introgressed_KEL_newterms<-sapply(strsplit(KEL_GO_introgressed$Term, "_GO"),head, 1)
KEL_GO_introgressed <- cbind(KEL_GO_introgressed, introgressed_KEL_newterms)
KEL_GO_introgressed$introgressed_KEL_newterms<-gsub( "_", " ", as.character(KEL_GO_introgressed$introgressed_KEL_newterms))
KEL_GO_significant <- KEL_GO_introgressed %>%
  top_n(10, Combined_Score)


png('KEL_introgressed_GO.png', width=8.5, height=3, units='in', res=320)
ggplot(KEL_GO_significant, aes(x=Combined_Score, y=introgressed_KEL_newterms,
                               label=introgressed_KEL_newterms))+
  geom_bar(stat = 'identity', width=0.9, size=0.6)+
  scale_y_discrete(expand = c(0, 0))+scale_fill_continuous()+
  theme(legend.position = 'none')+theme_classic()+
  scale_x_continuous(limits = c(0,50), expand = c(0, 0))+
  ylab('GO Biological Terms')+xlab('Combined Score')+
  geom_text(size=3.5, angle=0, hjust=-0.05,vjust='bottom' ,colour='black', fontface='bold', position = position_dodge(width = 1))+
  theme(axis.text.y = element_blank())+theme(axis.ticks.y = element_blank())+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"))
dev.off()
  

ggplot(KEL_GO_significant, aes(x=Combined_Score, y=reorder(introgressed_newterms, Combined_Score), colour=Adjusted_P.value))+
  geom_point(size=KEL_GO_significant$Combined_Score)+scale_colour_viridis(option='D', discrete = F)+
  theme_minimal()+xlab('Combined Score')+ylab('Gene ontology (Biological process)')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))

##BFM
BFM_GO_introgressed <- read.table('Introgression/GO_Biological_Process_BFM.txt', header=T, sep='', fill=T)
BFM_GO_introgressed$Term<-as.character(BFM_GO_introgressed$Term)
introgressed_BFM_newterms<-sapply(strsplit(BFM_GO_introgressed$Term, "_GO"),head, 1)
BFM_GO_introgressed <- cbind(BFM_GO_introgressed, introgressed_BFM_newterms)
BFM_GO_introgressed$introgressed_BFM_newterms<-gsub( "_", " ", as.character(BFM_GO_introgressed$introgressed_BFM_newterms))
BFM_GO_significant <- BFM_GO_introgressed %>%
  top_n(10, Combined_Score)



png('BFM_introgressed_GO.png', width=8.5, height=3, units='in', res=320)
ggplot(BFM_GO_significant, aes(x=Combined_Score, y=introgressed_BFM_newterms,
                               label=introgressed_BFM_newterms))+
  geom_bar(stat = 'identity', width=0.9, size=0.6)+
  scale_y_discrete(expand = c(0, 0))+scale_fill_continuous()+
  theme(legend.position = 'none')+theme_classic()+
  scale_x_continuous(limits = c(0,200), expand = c(0, 0))+
  ylab('GO Biological Terms')+xlab('Combined Score')+
  geom_text(size=3.5, angle=0, hjust=-0.05,vjust='bottom' ,colour='black', fontface='bold', position = position_dodge(width = 1))+
  theme(axis.text.y = element_blank())+theme(axis.ticks.y = element_blank())+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"))
dev.off()
           

##BFM
MLF_GO_introgressed <- read.table('Introgression/GO_biological_terms_flyenrichr_MLF_introgresed.txt', header=T, sep='', fill=T)
MLF_GO_introgressed$Term<-as.character(MLF_GO_introgressed$Term)
introgressed_MLF_newterms<-sapply(strsplit(MLF_GO_introgressed$Term, "_GO"),head, 1)
MLF_GO_introgressed <- cbind(MLF_GO_introgressed, introgressed_MLF_newterms)
MLF_GO_introgressed$introgressed_MLF_newterms<-gsub( "_", " ", as.character(MLF_GO_introgressed$introgressed_MLF_newterms))
MLF_GO_significant <- MLF_GO_introgressed %>%
  top_n(10, Combined_Score)


png('MLF_introgressed_GO.png', width=8.5, height=3, units='in', res=320)
ggplot(MLF_GO_significant, aes(x=Combined_Score, y=introgressed_MLF_newterms,
                               label=introgressed_MLF_newterms))+
  geom_bar(stat = 'identity', width=0.9, size=0.6)+
  scale_y_discrete(expand = c(0, 0))+scale_fill_continuous()+
  theme(legend.position = 'none')+theme_classic()+
  scale_x_continuous(limits = c(0,50), expand = c(0, 0))+
  ylab('GO Biological Terms')+xlab('Combined Score')+
  geom_text(size=4.5, angle=0, hjust=-0.05,vjust='bottom' ,colour='black', fontface='bold', position = position_dodge(width = 1))+
  theme(axis.text.y = element_blank())+theme(axis.ticks.y = element_blank())+
  theme(axis.text=element_text(size=14, face='bold'), 
        axis.title=element_text(size=14,face="bold"))
dev.off()



##NAL
NAL_GO_introgressed <- read.table('Introgression/GO_Biological_Process_NAL_introgressed.txt', header=T, sep='', fill=T)
NAL_GO_introgressed$Term<-as.character(NAL_GO_introgressed$Term)
introgressed_NAL_newterms<-sapply(strsplit(NAL_GO_introgressed$Term, "_GO"),head, 1)
NAL_GO_introgressed <- cbind(NAL_GO_introgressed, introgressed_NAL_newterms)
NAL_GO_introgressed$introgressed_NAL_newterms<-gsub( "_", " ", as.character(NAL_GO_introgressed$introgressed_NAL_newterms))
NAL_GO_significant <- NAL_GO_introgressed %>%
  top_n(10, Combined_Score)


png('NAL_introgressed_GO.png', width=8.5, height=3, units='in', res=320)
ggplot(NAL_GO_significant, aes(x=Combined_Score, y=introgressed_NAL_newterms,
                               label=introgressed_NAL_newterms))+
  geom_bar(stat = 'identity', width=0.9, size=0.6)+
  scale_y_discrete(expand = c(0, 0))+scale_fill_continuous()+
  theme(legend.position = 'none')+theme_classic()+
  scale_x_continuous(limits = c(0,50), expand = c(0, 0))+
  ylab('GO Biological Terms')+xlab('Combined Score')+
  geom_text(size=3.5, angle=0, hjust=-0.05,vjust='bottom' ,colour='black', fontface='bold', position = position_dodge(width = 1))+
  theme(axis.text.y = element_blank())+theme(axis.ticks.y = element_blank())+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"))
dev.off()


ggplot(BFM_fdm,aes(x=Chromosome, y=f_dM, fill=Chromosome))+
  geom_boxplot()+theme_bw()


################# fd vs dxy ################
MLF_fd<-read.table("fd_MLF_100_200.txt",header=T)
BFM_fd<-read.table("fd_BFM_100_200.txt",header=T)
KEL_fd<-read.table("fd_KEL_100_200.txt",header=T)
NAL_fd<-read.table("fd_NAL_100_200.txt",header=T)

MLF_fd$Chromosome <- chrom$Chromosome[ match(MLF_fd$chr, chrom$Contig)]
BFM_fd$Chromosome <- chrom$Chromosome[ match(BFM_fd$chr, chrom$Contig)]
KEL_fd$Chromosome <- chrom$Chromosome[ match(KEL_fd$chr, chrom$Contig)]
NAL_fd$Chromosome <- chrom$Chromosome[ match(NAL_fd$chr, chrom$Contig)]

#omit na
MLF_fd[(MLF_fd$f_d<0 | MLF_fd$f_d>1) & !is.na(MLF_fd$f_d),"f_d"]<-NA
MLF_fd <- na.omit(MLF_fd)
BFM_fd[(BFM_fd$f_d<0 | BFM_fd$f_d>1) & !is.na(BFM_fd$f_d),"f_d"]<-NA
BFM_fd <- na.omit(BFM_fd)
KEL_fd[(KEL_fd$f_d<0 | KEL_fd$f_d>1) & !is.na(KEL_fd$f_d),"f_d"]<-NA
KEL_fd <- na.omit(KEL_fd)
NAL_fd[(NAL_fd$f_d<0 | NAL_fd$f_d>1) & !is.na(NAL_fd$f_d),"f_d"]<-NA
NAL_fd <- na.omit(NAL_fd)

##################################################################
########### Find overlapping outliers for introgressed regions #############

#### MLF 
quants <- quantile(MLF_fd$f_d, c(0.95, 0.99))
MLF_fd$quant  <- with(MLF_fd, factor(ifelse(f_d < quants[1], 0, 
                                                                            ifelse(f_d< quants[2], 1, 2))))


MLF_fd_outliers <- filter(MLF_fd, quant != '0')

#### BFM 

quants <- quantile(BFM_fd$f_d, c(0.95, 0.99))
BFM_fd$quant  <- with(BFM_fd, factor(ifelse(f_d < quants[1], 0, 
                                            ifelse(f_d< quants[2], 1, 2))))

BFM_fd_outliers <- filter(BFM_fd, quant != '0')
##### KEL 

quants <- quantile(KEL_fd$f_d, c(0.95, 0.99))
KEL_fd$quant  <- with(KEL_fd, factor(ifelse(f_d < quants[1], 0, 
                                            ifelse(f_d< quants[2], 1, 2))))

KEL_fd_outliers <- filter(KEL_fd, quant != '0')

##### NAL

quants <- quantile(NAL_fd$f_d, c(0.95, 0.99))
NAL_fd$quant  <- with(NAL_fd, factor(ifelse(f_d < quants[1], 0, 
                                            ifelse(f_d< quants[2], 1, 2))))

NAL_fd_outliers <- filter(NAL_fd, quant != '0')


ggplot(MLF_fd, aes(x=Chromosome, y=f_d, fill=Chromosome))+
theme_pubr()+
  geom_boxplot(width=.2, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('f'[d] (La-F)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(BFM_fd, aes(x=Chromosome, y=f_d, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('f'[d] (F-M)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(KEL_fd, aes(x=Chromosome, y=f_d, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('f'[d] (E-Li)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

ggplot(NAL_fd, aes(x=Chromosome, y=f_d, fill=Chromosome, colour=Chromosome))+
  geom_boxplot(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_quasirandom(dodge.width=1, alpha=0.5) +scale_fill_viridis(option='D', discrete=T)+scale_colour_viridis(option='D', discrete=T)+
  labs(y=expression('f'[d] (A-Lu)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')


mlf_fd_lm<-lm(f_d ~ Chromosome, data = MLF_fd)
mlf_fd_aov<-aov(f_d ~ Chromosome, data = MLF_fd)

sd(MLF_fdm$f_dM)
mean(MLF_fd$f_d)
sd(BFM_fdm$f_dM)
mean(BFM_fd$f_d)
sd(KEL_fdm$f_dM)
mean(KEL_fd$f_d)
sd(NAL_fdm$f_dM)
mean(NAL_fd$f_d)
  
BFM_fd%>%
  group_by(Chromosome) %>%
  summarise(meanfd = mean(f_d))

MLF_fd%>%
  group_by(Chromosome) %>%
  summarise(meanfd = mean(f_d))

KEL_fd%>%
  group_by(Chromosome) %>%
  summarise(meanfd = mean(f_d))

NAL_fd%>%
  group_by(Chromosome) %>%
  summarise(meanfd = mean(f_d))


anova(mlf_fd_lm)
TukeyHSD(mlf_fd_aov)

bfm_fd_lm<-lm(f_d ~ Chromosome, data = BFM_fd)
bfm_fd_aov<-aov(f_d ~ Chromosome, data = BFM_fd)

anova(bfm_fd_lm)
TukeyHSD(bfm_fd_aov)

kel_fd_lm<-lm(f_d ~ Chromosome, data = KEL_fd)
kel_fd_aov<-aov(f_d ~ Chromosome, data = KEL_fd)

anova(kel_fd_lm)
TukeyHSD(kel_fd_aov)


nal_fd_lm<-lm(f_d ~ Chromosome, data = NAL_fd)
nal_fd_aov<-aov(f_d ~ Chromosome, data = NAL_fd)

anova(nal_fd_lm)
TukeyHSD(nal_fd_aov)



################################################

#get means for each scaffold and rename columns in fd
MLF_fd_mean<-aggregate(MLF_fd[, 4], list(MLF_fd$chr), mean)
names(MLF_fd_mean) <- c("scaffolds", "fd")
BFM_fd_mean<-aggregate(BFM_fd[, 4], list(BFM_fd$chr), mean)
names(BFM_fd_mean) <- c("scaffolds", "fd")
KEL_fd_mean<-aggregate(KEL_fd[, 4], list(KEL_fd$chr), mean)
names(KEL_fd_mean) <- c("scaffolds", "fd")
NAL_fd_mean<-aggregate(NAL_fd[, 4], list(NAL_fd$chr), mean)
names(NAL_fd_mean) <- c("scaffolds", "fd")

#get means for each scaffold and rename columns in dxy
MLF_dxy<-aggregate(Montana_Lacicola_genic[, 9], list(Montana_Lacicola_genic$sequence_id), mean)
names(MLF_dxy) <- c("scaffolds", "dxy")
BFM_dxy<-aggregate(Borealis_Flavomontana_genic[, 9], list(Borealis_Flavomontana_genic$sequence_id), mean)
names(BFM_dxy) <- c("scaffolds", "dxy")
KEL_dxy<-aggregate(Kanekoi_Ezoana_genic[, 9], list(Kanekoi_Ezoana_genic$sequence_id), mean)
names(KEL_dxy) <- c("scaffolds", "dxy")
NAL_dxy<-aggregate(Novamexicana_Americana_genic[, 9], list(Novamexicana_Americana_genic$sequence_id), mean)
names(NAL_dxy) <- c("scaffolds", "dxy")


MLF_inter<-aggregate(Montana_Lacicola_intergenic[, 9], list(Montana_Lacicola_intergenic$sequence_id), mean)
names(MLF_inter) <- c("scaffolds", "dxy")
BFM_inter<-aggregate(Borealis_Flavomontana_intergenic[, 9], list(Borealis_Flavomontana_intergenic$sequence_id), mean)
names(BFM_inter) <- c("scaffolds", "dxy")
KEL_inter<-aggregate(Kanekoi_Ezoana_intergenic[, 9], list(Kanekoi_Ezoana_intergenic$sequence_id), mean)
names(KEL_inter) <- c("scaffolds", "dxy")
NAL_inter<-aggregate(Novamexicana_Americana_intergenic[, 9], list(Novamexicana_Americana_intergenic$sequence_id), mean)
names(NAL_inter) <- c("scaffolds", "dxy")


MLF_pi<-aggregate(Montana_Lacicola_genic[, mean(7:8)], list(Montana_Lacicola_genic$sequence_id), mean)
names(MLF_pi) <- c("scaffolds", "pi")
BFM_pi<-aggregate(Borealis_Flavomontana_genic[, mean(7:8)], list(Borealis_Flavomontana_genic$sequence_id), mean)
names(BFM_pi) <- c("scaffolds", "pi")
KEL_pi<-aggregate(Kanekoi_Ezoana_genic[, mean(7:8)], list(Kanekoi_Ezoana_genic$sequence_id), mean)
names(KEL_pi) <- c("scaffolds", "pi")
NAL_pi<-aggregate(Novamexicana_Americana_genic[, mean(7:8)], list(Novamexicana_Americana_genic$sequence_id), mean)
names(NAL_pi) <- c("scaffolds", "pi")


MLF_dxy_fd_inter<-merge(MLF_inter,MLF_fd_mean,
                  by=c('scaffolds'), all=T)
BFM_dxy_fd_inter<-merge(BFM_inter,BFM_fd_mean,
                  by=c('scaffolds'), all=T)
KEL_dxy_fd_inter<-merge(KEL_inter,KEL_fd_mean,
                  by=c('scaffolds'), all=T)
NAL_dxy_fd_inter<-merge(NAL_inter,NAL_fd_mean,
                  by=c('scaffolds'), all=T)


MLF_dxy_coding_noncoding<-merge(MLF_inter,MLF_dxy,
                        by=c('scaffolds'), all=T)
BFM_dxy_coding_noncoding<-merge(BFM_inter,BFM_dxy,
                        by=c('scaffolds'), all=T)
KEL_dxy_coding_noncoding<-merge(KEL_inter,KEL_dxy,
                        by=c('scaffolds'), all=T)
NAL_dxy_coding_noncoding<-merge(NAL_inter,NAL_dxy,
                        by=c('scaffolds'), all=T)



MLF_dxy_fd<-merge(MLF_dxy,MLF_fd_mean,
                                     by=c('scaffolds'), all=T)
BFM_dxy_fd<-merge(BFM_dxy,BFM_fd_mean,
                  by=c('scaffolds'), all=T)
KEL_dxy_fd<-merge(KEL_dxy,KEL_fd_mean,
                  by=c('scaffolds'), all=T)
NAL_dxy_fd<-merge(NAL_dxy,NAL_fd_mean,
                  by=c('scaffolds'), all=T)


MLF_pi_fd<-merge(MLF_pi,MLF_fd_mean,
                  by=c('scaffolds'), all=T)
BFM_pi_fd<-merge(BFM_pi,BFM_fd_mean,
                  by=c('scaffolds'), all=T)
KEL_pi_fd<-merge(KEL_pi,KEL_fd_mean,
                  by=c('scaffolds'), all=T)
NAL_pi_fd<-merge(NAL_pi,NAL_fd_mean,
                  by=c('scaffolds'), all=T)




MLF_dxy_fd <- na.omit(MLF_dxy_fd)
BFM_dxy_fd <- na.omit(BFM_dxy_fd)
KEL_dxy_fd <- na.omit(KEL_dxy_fd)
NAL_dxy_fd <- na.omit(NAL_dxy_fd)

MLF_dxy_fd_inter <- na.omit(MLF_dxy_fd_inter)
BFM_dxy_fd_inter <- na.omit(BFM_dxy_fd_inter)
KEL_dxy_fd_inter <- na.omit(KEL_dxy_fd_inter)
NAL_dxy_fd_inter <- na.omit(NAL_dxy_fd_inter)

MLF_dxy_coding_noncoding <- na.omit(MLF_dxy_coding_noncoding)
BFM_dxy_coding_noncoding <- na.omit(BFM_dxy_coding_noncoding)
KEL_dxy_coding_noncoding <- na.omit(KEL_dxy_coding_noncoding)
NAL_dxy_coding_noncoding <- na.omit(NAL_dxy_coding_noncoding)


MLF_pi_fd <- na.omit(MLF_pi_fd)
BFM_pi_fd <- na.omit(BFM_pi_fd)
KEL_pi_fd <- na.omit(KEL_pi_fd)
NAL_pi_fd <- na.omit(NAL_pi_fd)


names(MLF_dxy_coding_noncoding)[2] <- "noncoding_dxy"
names(MLF_dxy_coding_noncoding)[3] <- "coding_dxy"

names(BFM_dxy_coding_noncoding)[2] <- "noncoding_dxy"
names(BFM_dxy_coding_noncoding)[3] <- "coding_dxy"

names(KEL_dxy_coding_noncoding)[2] <- "noncoding_dxy"
names(KEL_dxy_coding_noncoding)[3] <- "coding_dxy"

names(NAL_dxy_coding_noncoding)[2] <- "noncoding_dxy"
names(NAL_dxy_coding_noncoding)[3] <- "coding_dxy"


ggplot(MLF_dxy_coding_noncoding, aes(x=noncoding_dxy,y=coding_dxy))+
  geom_point(size=3, colour='black', pch=21)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (M-L)))+
  labs(y=expression('d'[XY] (M-L)))

ggplot(BFM_dxy_coding_noncoding, aes(x=noncoding_dxy,y=coding_dxy))+
  geom_point(size=3, colour='black', pch=21)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (B-F)))+
  labs(y=expression('d'[XY] (B-F)))

ggplot(KEL_dxy_coding_noncoding, aes(x=noncoding_dxy,y=coding_dxy))+
  geom_point(size=3, colour='black', pch=21)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (K-E)))+
  labs(y=expression('d'[XY] (K-E)))

ggplot(NAL_dxy_coding_noncoding, aes(x=noncoding_dxy,y=coding_dxy))+
  geom_point(size=3, colour='black', pch=17, alpha=0.7)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (N-A)))+
  labs(y=expression('d'[XY] (N-A)))





genic_cor_dxy_fd_mlf<-ggplot(MLF_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (M-L)))+
  labs(y=expression('f'[d] (M-L-F)))



genic_cor_dxy_fd_bfm<-ggplot(BFM_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson", label.x = 0.012, label.y=0.075)+
  labs(x=expression('d'[XY] (B-F)))+
  labs(y=expression('f'[d] (B-F-M)))



genic_cor_dxy_fd_kel<-ggplot(KEL_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (K-E)))+
         labs(y=expression('f'[d] (K-E-L)))



genic_cor_dxy_fd_nal<-ggplot(NAL_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (N-A)))+
  labs(y=expression('f'[d] (N-A-L)))

png('genic_correlations_introgression_dxy.png', width=6.25, height=4, units='in', res=320)
ggarrange(genic_cor_dxy_fd_mlf, genic_cor_dxy_fd_bfm, genic_cor_dxy_fd_kel,
          genic_cor_dxy_fd_nal, ncol=2, nrow=2)
dev.off()



genic_cor_pi_fd_ml<-ggplot(MLF_pi_fd, aes(x=pi,y=fd))+
  geom_point(size=3, colour='black', pch=21)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('pi M-L'))+
  labs(y=expression('f'[d] (M-L-F)))


genic_cor_pi_fd_bf<-ggplot(BFM_pi_fd, aes(x=pi,y=fd))+
  geom_point(size=3, colour='black', pch=21)+
  theme_bw()+stat_cor(method = "pearson",label.y=0.079)+
  labs(x=expression('pi B-F'))+
  labs(y=expression('f'[d] (B-F-M)))


genic_cor_pi_fd_ke<-ggplot(KEL_pi_fd, aes(x=pi,y=fd))+
  geom_point(size=3, colour='black', pch=21)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('pi K-E'))+
  labs(y=expression('f'[d] (M-L-F)))

genic_cor_pi_fd_na<-ggplot(NAL_pi_fd, aes(x=pi,y=fd))+
  geom_point(size=3, colour='black', pch=21)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('pi N-A'))+
  labs(y=expression('f'[d] (M-L-F)))

png('genic_correlations_introgression_pi.png', width=6.25, height=4, units='in', res=320)
ggarrange(genic_cor_dxy_fd_ml, genic_cor_dxy_fd_bf, genic_cor_dxy_fd_ke,
          genic_cor_dxy_fd_na, ncol=2, nrow=2)
dev.off()


inter_cor_dxy_fd_ml<-ggplot(MLF_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (M-L)))+
  labs(y=expression('f'[d] (M-L-F)))

inter_cor_dxy_fd_bf<-ggplot(BFM_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson",label.y=0.078)+
  labs(x=expression('d'[XY] (B-F)))+
  labs(y=expression('f'[d] (B-F-M)))

inter_cor_dxy_fd_ke<-ggplot(KEL_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (K-E)))+
  labs(y=expression('f'[d] (K-E-L)))

inter_cor_dxy_fd_na<-ggplot(NAL_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.8)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (N-A)))+
  labs(y=expression('f'[d] (N-A-L)))



png('intergenic_correlations_introgression_dxy.png', width=6.25, height=4, units='in', res=320)
ggarrange(inter_cor_dxy_fd_ml, inter_cor_dxy_fd_bf, inter_cor_dxy_fd_ke,
          inter_cor_dxy_fd_na, ncol=2, nrow=2)
dev.off()



######## Dmin analyses (replot) ########
options(scipen = 999)
#dmin_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")



Dmin <- read.table('SETS_all_wit_tree_Dmin.txt', header=T)
Dmin<-within(Dmin, Interaction <- paste(P2, P3, sep='_'))
fdr_pvalues<-p.adjust(Dmin$p.value, method = 'bonferroni', n = 120)
Dmin <- cbind(Dmin, fdr_pvalues)
Dmin_ordered <- Dmin[order(Dmin$Interaction, -abs(Dmin$Dstatistic) ), ]
Dmin_rmdups<-Dmin_ordered[ !duplicated(Dmin_ordered$Interaction), ]  

mean(Dmin$Dstatistic)


pal <- wes_palette("FantasticFox1",  52, type = "continuous")


Dmin <- filter(Dmin_rmdups, fdr_pvalues < 0.05)
Dmin <- filter(Dmin, f4.ratio > 0.01)
#Dmin <- filter(Dmin_rmdups, Dstatistic > 0.1)
Dmin<-within(Dmin, Interaction <- paste(P1,P2, P3, sep='_'))

ggplot(Dmin, aes(x=Interaction, y=Dstatistic, colour=p.value))+
  geom_point(size=3)+theme_minimal()+
  scale_colour_gradientn(colours = pal)+
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size=6))



P3_order <- factor(Dmin$P3, level = c('Mon', 'Lac', 'Bor','Flavo','Lum', 'Amer', 'Nova', 'Litt','Ez', 'Kan'))
P2_order <- factor(Dmin$P2, level = c('Mon', 'Lac', 'Bor','Flavo', 'Amer', 'Nova', 'Litt','Ez', 'Kan'))

png('Dstatistic.png', width=5, height=5, units='in', res=320)
ggballoonplot(Dmin_rmdups, x = P2_order, y = level_order,
             size = "Dstatistic", fill = "Dstatistic") +
  scale_fill_gradientn(colors = pal) +
  guides(size = FALSE)+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))
dev.off()

mean(Dmin$Dstatistic)
mean(Dmin$f4.ratio)

##### Dstatistic tree ##########################
Dtree <- read.table('SETS_all_with_tree_tree.txt', header=T)
Dtree<-within(Dtree, Interaction <- paste(P2, P3, sep='_'))
fdr_pvalues<-p.adjust(Dtree$p.value, method = 'bonferroni', n = 120)
Dtree <- cbind(Dtree, fdr_pvalues)
Dtree_ordered <- Dtree[order(Dtree$Interaction, -abs(Dtree$Dstatistic) ), ]
Dtree_rmdups<-Dtree_ordered[ !duplicated(Dtree_ordered$Interaction), ]  



pal <- wes_palette("FantasticFox1",  52, type = "continuous")


Dtree <- filter(Dtree_rmdups, fdr_pvalues < 0.05)
Dtree <- filter(Dtree, f_G > 0.01)
#Dmin <- filter(Dmin_rmdups, Dstatistic > 0.1)
Dtree<-within(Dtree, Interaction <- paste(P1,P2, P3, sep='_'))

#################################################
### June 2021 

f4.ratio_Dmin<-ggplot(data = Dmin, aes(x=P2, y=P3, fill=f4.ratio)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "grey3", 
                       midpoint = 0.008, limit = c(0.0001,0.04), space = "Lab", 
                       name="f4-ratio\n") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(size=12))+
  coord_fixed()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=8),
    legend.direction = "vertical")


Dminimum+theme(legend.direction = "vertical")
png('Dstatistic.png', width=10, height=10, units='in', res=300)
ggarrange(Dminimum+theme(legend.direction = "vertical",legend.justification = c(0.5, 0.5)), f4.ratio_Dmin,nrow = 1)
dev.off()
####################### Parker FST ####################################
#### First retrieve my set of genes
my_gene_set <- read.table('UNIPROT_2_FLYBASE_posselect.txt')

###now retrieve the Parker gene set 
Parker2018 <- read.csv('FST_each_gene_mel_Parker2018.csv')

###Introgressed variation upload!?
Introgressed_gene_set <- read.table('Introgression/FlyBase_Introgressed_altogether_IDs.txt')

#Now find and match for positive selec genes
Fst_my_gene_set<-Parker2018[(Parker2018$FB_vir_gene_id%in%my_gene_set$V1),]
names(Parker2018) <- c("mel", "Colorado_Oulanka", "Colorado_Vancouver", "Oulanka_Vancouver", "vir")
names(Fst_my_gene_set) <- c("mel", "Colorado_Oulanka", "Colorado_Vancouver", "Oulanka_Vancouver", "vir")
#Now find and match for introgressed genes
Fst_introgressed_gene_set <- Parker2018[(Parker2018$FB_vir_gene_id%in%Introgressed_gene_set$V1),]


#melt Fst for my geneset to be able to plot it
Melted_Fst_gene_set<-melt(Fst_my_gene_set, id.vars=c('mel', 'vir'), measure.vars=c('Colorado_Oulanka', 
                                                                'Oulanka_Vancouver', 'Colorado_Vancouver'))
Melted_Fst_background<-melt(Parker2018, id.vars=c('mel', 'vir'), measure.vars=c('Colorado_Oulanka', 
                                                                                               'Oulanka_Vancouver', 'Colorado_Vancouver'))
Melted_Fst_background<-na.omit(Melted_Fst_background)
Background.summary <- aggregate(. ~ variable, mean, data=Melted_Fst_background)
Foreground.summary <- aggregate(. ~ variable, mean, data=Melted_Fst_gene_set)

png('FST_pos_vs_background.png', width=7, height=2, units='in', res=320)
ggplot(Melted_Fst_gene_set, aes(x=variable, y=value,colour=variable))+
  geom_jitter(size=2, alpha=0.6)+
  geom_crossbar(data=Background.summary, aes(ymin = value, ymax = value),
                size=0.5,col="red", width = .5)+
  theme_bw()+xlab('')+labs(y=expression('F'[ST]))+coord_flip()+
  scale_colour_manual(values=wes_palette("IsleofDogs1"))+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))
dev.off()



BF_mean_pi<-Borealis_Flavomontana_genic[, mean(7:8)]
ML_mean_pi<-Montana_Lacicola_genic[, mean(7:8)]
KE_mean_pi<-Kanekoi_Ezoana_genic[, mean(7:8)]



Borealis_Flavomontana_genic <- cbind(Borealis_Flavomontana_genic, BF_mean_pi)
Montana_Lacicola_genic <- cbind(Montana_Lacicola_genic, ML_mean_pi)
Kanekoi_Ezoana_genic <- cbind(Kanekoi_Ezoana_genic, KE_mean_pi)

ggplot(Borealis_Flavomontana_genic, aes(x=dxy, y=BF_mean_pi, label=sequence_id))+
  geom_point()+
  geom_text(aes(label=ifelse(dxy>0.03,as.character(sequence_id),'')),hjust=0,vjust=0)

ggplot(Kanekoi_Ezoana_genic, aes(x=dxy, y=KE_mean_pi, label=sequence_id))+
  geom_point()+
  geom_text(aes(label=ifelse(dxy>0.03,as.character(sequence_id),'')),hjust=0,vjust=0)






########### Concordance factors ##########
tab2=read.table('concord.cf.stat_tree',header=TRUE, na.strings = "NA")
tab3=read.table('concord_scf.cf.stat_loci',header=TRUE, na.strings = "NA")
gene_concord<-na.omit(tab2)
site_concord<-na.omit(tab3)
gene_con_melt <- melt(gene_concord, id.vars=c('ID', 'TreeID'), measure.vars=c('gC', 
                                                                                 'gD1', 'gD2'))
site_con_melt <- melt(site_concord, id.vars=c('ID', 'PartID'), measure.vars=c('sC', 
                                                                              'sD1', 'sD2'))

ggplot(gene_con_melt, aes(x=variable,y=value, fill=variable))+
  geom_bar(stat='identity', alpha=0.9)+facet_wrap(~ID)+scale_fill_viridis(option='D', discrete=T)+
  theme_minimal()+labs(x='Number of decisive gene trees', y='Topology',colour = "Topology")
  

### Where are the genes that are most discordant found?
locations_concord=read.table('Coordinates_BUSCO_geneconcord.txt')
gene_con_melt$scaffold <- locations_concord$V2[ match(gene_con_melt$TreeID, locations_concord$V1)]
site_con_melt$scaffold <- locations_concord$V2[ match(site_con_melt$PartID, locations_concord$V1)]
gene_con_melt$Chromosome <- chrom$Chromosome[ match(gene_con_melt$scaffold, chrom$Contig)]
site_con_melt$Chromosome <- chrom$Chromosome[ match(site_con_melt$scaffold, chrom$Contig)]


gene_con_melt<-na.omit(gene_con_melt)
site_con_melt<-na.omit(site_con_melt)

png('gene_concordance_factors_bychromosome.png', width=7, height=5, units='in', res=320)
ggplot(gene_con_melt, aes(x=Chromosome,y=value, fill=variable))+
  geom_bar(stat='identity')+facet_wrap(~ID)+scale_fill_viridis(option='D', discrete=T)+
  theme_pubr()+labs(y='count', fill='topologies')
dev.off()

png('site_concordance_factors_bychromosome.png', width=7, height=5, units='in', res=320)
ggplot(site_con_melt, aes(x=Chromosome,y=value, fill=variable))+
  geom_bar(stat='identity')+facet_wrap(~ID)+scale_fill_viridis(option='D', discrete=T)+
  theme_pubr()+labs(y='count', fill='topologies')
dev.off()

kruskal.test(value ~ variable, data=gene_con_melt)
kruskal.test(value ~ variable, data=site_con_melt)

gene_con_lm<-lm(value ~ Chromosome, data = gene_con_melt)
gene_con_aov<-aov(value ~ Chromosome, data = gene_con_melt)

anova(gene_con_lm)
TukeyHSD(gene_con_aov)


#########  BIG TABLE COMPARISON


Big_table_all_info <- read.csv('dxy_isolation.csv')
Big_table_all_info$Dstatistic <- Dmin_rmdups$Dstatistic[ match(Big_table_all_info$Comparisons, Dmin_rmdups$Interaction)]

Big_table_all_info_ave



ggplot(Big_table_all_info, aes(x=dxy, y=Dstatistic, colour=Crossability))+
  geom_point(size=2)+facet_wrap(~Chr)+theme_bw()
  



#####################################
#### MAP
world <- map_data("world")
species_sampling_info<-read.csv('18genomes_v2_Leeban.csv')

png('worldmap_virilissamples.png', width=7, height=5, units='in', res=320)
ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill='darkgrey',) + 
  coord_fixed(1.3)+coord_sf(ylim = c(20, 100))+theme_pubclean()+scale_colour_viridis(option='D', discrete=T)+
  geom_jitter(data = species_sampling_info, mapping = aes(x = longitude, y = latitude, color = species), alpha=0.8, size=3,width = 2)
dev.off()
