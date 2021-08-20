#Reanalysis of dxy using simon martin scripts, for species pairs and more.

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
library(gt)

setwd("~/Desktop/PhD_Labbook/Comparative_analysis/virilis_tree")
#Get all the files read in
Genic_SimonMartin_sumstats <- read.csv(gzfile("genic_summarystats_output_20920.csv.gz"))
Intergenic_SimonMartin_sumstats <- read.csv(gzfile("intergenic_summarystats_output_20920.csv.gz"))
chrom<-read.table("monCan3F9_chromosomes.txt",header=T)

#Match chrom to scaffolds
Genic_SimonMartin_sumstats$Chromosome <- chrom$Chromosome[ match(Genic_SimonMartin_sumstats$scaffold, chrom$Contig)]
Intergenic_SimonMartin_sumstats$Chromosome <- chrom$Chromosome[ match(Intergenic_SimonMartin_sumstats$scaffold, chrom$Contig)]


#Get rid of NA in rows
Genic_SimonMartin_sumstats<-na.omit(Genic_SimonMartin_sumstats)
Intergenic_SimonMartin_sumstats<-na.omit(Intergenic_SimonMartin_sumstats)

###Genic

BOX_ML<-ggplot(Genic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Lacicola_Montana, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (M-L)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_NA<-ggplot(Genic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Americana_Novamexicana, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (N-A)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_KE<-ggplot(Genic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Ezoana_Kanekoi, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (K-E)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_BF<-ggplot(Genic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Borealis_Flavomontana, fill=Chromosome))+
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

#Statistics for divergence -- determining significance 
ml_dxy_lm<-lm(dxy_Lacicola_Montana ~ Chromosome, data = Genic_SimonMartin_sumstats)
ml_dxy_aov<-aov(dxy_Lacicola_Montana ~ Chromosome, data = Genic_SimonMartin_sumstats)

anova(ml_dxy_lm)
TukeyHSD(ml_dxy_aov)

#Borealis_Flavomontana
bf_dxy_lm<-lm(dxy_Borealis_Flavomontana ~ Chromosome, data = Genic_SimonMartin_sumstats)
anova(bf_dxy_lm)
bf_dxy_aov<-aov(dxy_Borealis_Flavomontana ~ Chromosome, data = Genic_SimonMartin_sumstats)
TukeyHSD(bf_dxy_aov)


#Kanekoi_Ezoana
ke_dxy_lm<-lm(dxy_Ezoana_Kanekoi ~ Chromosome, data = Genic_SimonMartin_sumstats)

anova(ke_dxy_lm)
ke_dxy_aov<-aov(dxy_Ezoana_Kanekoi ~ Chromosome, data = Genic_SimonMartin_sumstats)
TukeyHSD(ke_dxy_aov)

#Novamexicana_americana
na_dxy_lm<-lm(dxy_Americana_Novamexicana ~ Chromosome, data = Genic_SimonMartin_sumstats)
anova(na_dxy_lm)
na_dxy_aov<-aov(dxy_Americana_Novamexicana ~ Chromosome, data = Genic_SimonMartin_sumstats)
TukeyHSD(na_dxy_aov)



##Intergenic

BOX_ML_int<-ggplot(Intergenic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Lacicola_Montana, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (M-L)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_NA_int<-ggplot(Intergenic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Americana_Novamexicana, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (N-A)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_KE_int<-ggplot(Intergenic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Ezoana_Kanekoi, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (K-E)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

BOX_BF_int<-ggplot(Intergenic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Borealis_Flavomontana, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (B-F)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

png('dxy_chromosomes_ml_bf_int.png', width=5, height=2.5, units='in', res=320)
ggarrange(BOX_ML_int, BOX_BF_int, ncol = 2)
dev.off()
png('dxy_chromosomes_ke_na_int.png', width=5, height=2.5, units='in', res=320)
ggarrange(BOX_KE_int, BOX_NA_int, ncol = 2)
dev.off()

#######Now is that the case for intergenic regions?
#Montana-Lacicola
ml_int_dxy_lm<-lm(dxy_Lacicola_Montana ~ Chromosome, data = Intergenic_SimonMartin_sumstats)
ml_inter_dxy_aov<-aov(dxy_Lacicola_Montana ~ Chromosome, data = Intergenic_SimonMartin_sumstats)

anova(ml_int_dxy_lm)
TukeyHSD(ml_inter_dxy_aov)

#Borealis_Flavomontana
bf_dxy_lm_inter<-lm(dxy_Borealis_Flavomontana ~ Chromosome, data = Intergenic_SimonMartin_sumstats)
anova(bf_dxy_lm_inter)
bf_dxy_aov_inter<-aov(dxy_Borealis_Flavomontana ~ Chromosome, data = Intergenic_SimonMartin_sumstats)
TukeyHSD(bf_dxy_aov_inter)


#Kanekoi_Ezoana
ke_dxy_lm<-lm(dxy_Ezoana_Kanekoi ~ Chromosome, data = Intergenic_SimonMartin_sumstats)
anova(ke_dxy_lm)
ke_dxy_aov<-aov(dxy_Ezoana_Kanekoi ~ Chromosome, data = Intergenic_SimonMartin_sumstats)
TukeyHSD(ke_dxy_aov)

#Novamexicana_americana
na_dxy_lm<-lm(dxy_Americana_Novamexicana ~ Chromosome, data = Intergenic_SimonMartin_sumstats)
anova(na_dxy_lm)
na_dxy_aov<-aov(dxy_Americana_Novamexicana ~ Chromosome, data = Intergenic_SimonMartin_sumstats)
TukeyHSD(na_dxy_aov)



#### Now look for dxy genomewide
##################DXY CODING GENOME-WIDE #################

## MONTANA LACICOLA

mean(Genic_SimonMartin_sumstats$dxy_Lacicola_Montana)
mean(Intergenic_SimonMartin_sumstats$dxy_Lacicola_Montana)

meanData_ML = Genic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Lacicola_Montana)) 

meanData_ML

quants <- quantile(Genic_SimonMartin_sumstats$dxy_Lacicola_Montana, c(0.95, 0.99))
Genic_SimonMartin_sumstats$quant  <- with(Genic_SimonMartin_sumstats, factor(ifelse(dxy_Lacicola_Montana < quants[1], 0, 
                                                                            ifelse(dxy_Lacicola_Montana < quants[2], 1, 2))))


Outliers_spp_pairs_genic_ML<-Genic_SimonMartin_sumstats %>%
  filter(quant %in% c('1','2'),
         !quant_bf %in% c('1','2'),
         !quant_ke %in% c('1','2'),
         !quant_na %in% c('1','2'))

Outliers_spp_pairs_genic_ML<-Outliers_spp_pairs_genic_ML %>%
  select(scaffold, start, end)

write.table(Outliers_spp_pairs_genic_ML, file='ML_dxy_outliers.txt', sep='\t')

dxy_colours=c('gray10', 'orange1', 'orangered2')


png('genomewide_dxy_ml.png', width=6, height=2, units='in', res=320)
ggplot(Genic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Lacicola_Montana))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant))+
  geom_hline(data=meanData_ML, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (M-L)))+
  theme_minimal_hgrid()+ylim(0.0, 0.09
                             )+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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


Genic_species_pairs<-Genic_SimonMartin_sumstats %>% 
  select(scaffold, start, end, mid, sites,dxy_Ezoana_Kanekoi,dxy_Borealis_Flavomontana, dxy_Lacicola_Montana,
         dxy_Americana_Novamexicana, Chromosome)

#### BOREALIS flavomontana

meanData_bf = Genic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Borealis_Flavomontana)) 

meanData_bf

quants_bf <- quantile(Genic_SimonMartin_sumstats$dxy_Borealis_Flavomontana, c(0.95, 0.99))
Genic_SimonMartin_sumstats$quant_bf  <- with(Genic_SimonMartin_sumstats, factor(ifelse(dxy_Borealis_Flavomontana < quants_bf[1], 0, 
                                                                                         ifelse(dxy_Borealis_Flavomontana < quants_bf[2], 1, 2))))

png('genomewide_dxy_bf.png', width=6, height=2, units='in', res=320)
ggplot(Genic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Borealis_Flavomontana))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_bf))+
  geom_hline(data=meanData_bf, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (B-F)))+
  theme_minimal_hgrid()+ylim(0.0, 0.1)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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

mean(Genic_SimonMartin_sumstats$dxy_Ezoana_Kanekoi)
mean(Intergenic_SimonMartin_sumstats$dxy_Ezoana_Kanekoi)

meanData_ke = Genic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Ezoana_Kanekoi)) 

meanData_ke

quants_ke <- quantile(Genic_SimonMartin_sumstats$dxy_Ezoana_Kanekoi, c(0.95, 0.99))
Genic_SimonMartin_sumstats$quant_ke  <- with(Genic_SimonMartin_sumstats, factor(ifelse(dxy_Ezoana_Kanekoi < quants_ke[1], 0, 
                                                                           ifelse(dxy_Ezoana_Kanekoi < quants_ke[2], 1, 2))))

png('genomewide_dxy_ke.png', width=6, height=2, units='in', res=320)
ggplot(Genic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Ezoana_Kanekoi))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_ke))+
  geom_hline(data=meanData_ke, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (K-E)))+
  theme_minimal_hgrid()+ylim(0.02, 0.11)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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

meanData_na = Genic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Americana_Novamexicana)) 

meanData_na

quants_na <- quantile(Genic_SimonMartin_sumstats$dxy_Americana_Novamexicana, c(0.95, 0.99))
Genic_SimonMartin_sumstats$quant_na <- with(Genic_SimonMartin_sumstats, factor(ifelse(dxy_Americana_Novamexicana< quants_na[1], 0, 
                                                                                          ifelse(dxy_Americana_Novamexicana< quants_na[2], 1, 2))))

png('genomewide_dxy_na.png', width=6, height=2, units='in', res=320)
ggplot(Genic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Americana_Novamexicana))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_na))+
  geom_hline(data=meanData_na, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (N-A)))+
  theme_minimal_hgrid()+ylim(0.0, 0.08)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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



##################DXY CODING GENOME-WIDE #################

## MONTANA LACICOLA

meanData_ML_int = Intergenic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Lacicola_Montana)) 

meanData_ML_int

quants <- quantile(Intergenic_SimonMartin_sumstats$dxy_Lacicola_Montana, c(0.95, 0.99))
Intergenic_SimonMartin_sumstats$quant  <- with(Intergenic_SimonMartin_sumstats, factor(ifelse(dxy_Lacicola_Montana < quants[1], 0, 
                                                                                    ifelse(dxy_Lacicola_Montana < quants[2], 1, 2))))



png('genomewide_dxy_ml_int.png', width=6, height=2, units='in', res=320)
ggplot(Intergenic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Lacicola_Montana))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant))+
  geom_hline(data=meanData_ML, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (M-L)))+
  theme_minimal_hgrid()+ylim(0.0, 0.1
  )+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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

Genic_species_pairs<-Genic_SimonMartin_sumstats %>% 
  select(scaffold, start, end, mid, sites,dxy_Ezoana_Kanekoi,dxy_Borealis_Flavomontana, dxy_Lacicola_Montana,
         dxy_Americana_Novamexicana, Chromosome)

#### BOREALIS flavomontana

meanData_bf_int = Intergenic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Borealis_Flavomontana)) 

meanData_bf_int

quants_bf <- quantile(Intergenic_SimonMartin_sumstats$dxy_Borealis_Flavomontana, c(0.95, 0.99))
Intergenic_SimonMartin_sumstats$quant_bf  <- with(Intergenic_SimonMartin_sumstats, factor(ifelse(dxy_Borealis_Flavomontana < quants_bf[1], 0, 
                                                                                       ifelse(dxy_Borealis_Flavomontana < quants_bf[2], 1, 2))))

png('genomewide_dxy_bf_int.png', width=6, height=2, units='in', res=320)
ggplot(Intergenic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Borealis_Flavomontana))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_bf))+
  geom_hline(data=meanData_bf, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (B-F)))+
  theme_minimal_hgrid()+ylim(0.0, 0.1)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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

meanData_ke_int = Intergenic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Ezoana_Kanekoi)) 

meanData_ke_int

quants_ke <- quantile(Intergenic_SimonMartin_sumstats$dxy_Ezoana_Kanekoi, c(0.95, 0.99))
Intergenic_SimonMartin_sumstats$quant_ke  <- with(Intergenic_SimonMartin_sumstats, factor(ifelse(dxy_Ezoana_Kanekoi < quants_ke[1], 0, 
                                                                                       ifelse(dxy_Ezoana_Kanekoi < quants_ke[2], 1, 2))))

png('genomewide_dxy_ke_int.png', width=6, height=2, units='in', res=320)
ggplot(Intergenic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Ezoana_Kanekoi))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_ke))+
  geom_hline(data=meanData_ke, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (K-E)))+
  theme_minimal_hgrid()+ylim(0.02, 0.13)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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

meanData_na = Intergenic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Americana_Novamexicana)) 

meanData_na

quants_na <- quantile(Intergenic_SimonMartin_sumstats$dxy_Americana_Novamexicana, c(0.95, 0.99))
Intergenic_SimonMartin_sumstats$quant_na <- with(Intergenic_SimonMartin_sumstats, factor(ifelse(dxy_Americana_Novamexicana< quants_na[1], 0, 
                                                                                      ifelse(dxy_Americana_Novamexicana< quants_na[2], 1, 2))))

png('genomewide_dxy_na_int.png', width=6, height=2, units='in', res=320)
ggplot(Intergenic_SimonMartin_sumstats, aes(x=scaffold, y=dxy_Americana_Novamexicana))+
  geom_jitter(size=1.75,alpha=0.7, aes(colour = quant_na))+
  geom_hline(data=meanData_na, aes(yintercept=meanCDR), colour='red', linetype='dashed')+
  scale_colour_manual(values=dxy_colours)+labs(y=expression('d'[XY] (N-A)))+
  theme_minimal_hgrid()+ylim(0.0, 0.08)+theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
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

#########################################################
################ dxy across species ###########

Genic_SimonMartin_sumstats_dxymean<-Genic_SimonMartin_sumstats %>%
  select(starts_with("dxy_")) %>%
  summarise(across(everything(), list(mean)))

Intergenic_SimonMartin_sumstats_dxymean<-Intergenic_SimonMartin_sumstats %>%
  select(starts_with("dxy_")) %>%
  summarise(across(everything(), list(mean)))


Genic_SimonMartin_sumstats_dxymean<-melt(Genic_SimonMartin_sumstats_dxymean)
Intergenic_SimonMartin_sumstats_dxymean<-melt(Intergenic_SimonMartin_sumstats_dxymean)
#Genic_SimonMartin_sumstats_dxymean$variable<-gsub("dxy_","",Genic_SimonMartin_sumstats_dxymean$variable)
Big_table_with_averaged_dxy_across_chroms$Dxy_coding <- Genic_SimonMartin_sumstats_dxymean$value[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, Genic_SimonMartin_sumstats_dxymean$variable)]
Big_table_with_averaged_dxy_across_chroms$Dxy_noncoding <- Intergenic_SimonMartin_sumstats_dxymean$value[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, Intergenic_SimonMartin_sumstats_dxymean$variable)]


Big_table_all_info[Big_table_all_info==""] <- NA

Plotting_Biogeog<-filter(Big_table_all_info, Chr== '5')

ggplot(Big_table_all_info, aes(x=Node_depth, y=New_dxy))+
  geom_point()+theme_bw()

ggplot(Big_table_with_averaged_dxy_across_chroms, aes(x=PrematingIso, y=Dstatistic, colour=PrematingIso))+
  geom_point()+theme_bw()

ggplot(Plotting_Biogeog, aes(x=New_dxy, y=Dstatistic))+
  geom_point(size=2, alpha=0.7)+scale_colour_viridis(option='D', discrete=T)+theme_bw()

  
############# More simple summaries - X vs A
# Questions 
# X-A difference?

AutSex_Comparison<-Genic_SimonMartin_sumstats %>%
  select(Chromosome, starts_with('dxy_'))%>%
  mutate(SexAut=case_when(
    Chromosome=='2' ~ 'Autosome',
    Chromosome=='3' ~ 'Autosome',
    Chromosome=='4' ~ 'Autosome',
    Chromosome=='5' ~ 'Autosome',
    Chromosome=='X' ~ 'X'))

glimpse(AutSex_Comparison)



AutSex_Comparison<-pivot_longer(AutSex_Comparison,
             cols=dxy_Americana_Borealis:dxy_Novamexicana_Virilis,
             names_to="Comparisons",
             values_to="dxy")

AutSex_Comparison$Comparisons<-gsub("dxy_","",AutSex_Comparison$Comparisons)
AutSex_Comparison$Biogeography <- big_table$Biogeography[match(big_table$Comparisons, AutSex_Comparison$Comparisons)]

AutSex_Comparison<-AutSex_Comparison %>%
  group_by(SexAut,Comparisons) %>%
  mutate(Average_dxy=mean(dxy)) %>%
  select(Comparisons, Average_dxy, SexAut) %>%
  distinct(Average_dxy)
  



ggplot(AutSex_Comparison, aes(x=SexAut, y=Average_dxy, fill=SexAut))+
  geom_boxplot()+theme_bw()+labs(x='Chromosome', y='Dxy', fill='')
  

ggplot(AutSex_Comparison, aes(Average_dxy, fill=SexAut))+
  geom_area(stat='bin')

wilcox.test(Average_dxy ~ SexAut, AutSex_Comparison)


ggplot(AutSex_Comparison, aes(x = SexAut, y = Average_dxy, group=Comparisons, colour=SexAut)) +
  geom_point(size=2, alpha=0.7) + geom_line(colour='black') +
  theme_bw()+
  scale_color_hue(direction = -1)+
  labs(y='Dxy', colour='', x='Chromosome')

AutSex_Comparison$Comparisons<-gsub("dxy_","",AutSex_Comparison$Comparisons)
big_table$X_A_ratio<- X_A_Ratio$ratio[ match(big_table$Comparisons, AutSex_Comparison$Comparisons)]
AutSex_Comparison$Comparisons<-gsub("-"," vs. ",AutSex_Comparison$Comparisons)

### Plot showing difference between autosomes and X chromosome

png('dxycoding_difference_sex_aut.png', width=7, height=8, units='in', res=320)
ggplot(AutSex_Comparison, aes(x=Comparisons,y=Average_dxy, colour=SexAut))+
  geom_line(colour='black')+geom_point(alpha=1, size=2)+theme_minimal()+
  coord_flip()+labs(y='Dxy', colour='', x='')+
  scale_color_hue(direction = -1)+
  theme(text = element_text(size=10))+
  theme(text=element_text(face = "bold"))
dev.off()

X_A_Ratio<-AutSex_Comparison %>%
  group_by(Comparisons) %>%
  summarise(ratio = Average_dxy[SexAut=="X"]/Average_dxy[SexAut=="Autosome"])

X_A_Ratio$Comparisons<-gsub("-","_",X_A_Ratio$Comparisons)
Big_table_with_averaged_dxy_across_chroms$X_A_Ratio <- X_A_Ratio$ratio[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, X_A_Ratio$Comparisons)]


###Intergenic

AutSex_Comparison_Int<-Intergenic_SimonMartin_sumstats %>%
  select(Chromosome, starts_with('dxy_'))%>%
  mutate(SexAut=case_when(
    Chromosome=='2' ~ 'Autosome',
    Chromosome=='3' ~ 'Autosome',
    Chromosome=='4' ~ 'Autosome',
    Chromosome=='5' ~ 'Autosome',
    Chromosome=='X' ~ 'X'))

glimpse(AutSex_Comparison_Int)

AutSex_Comparison_Int<-pivot_longer(AutSex_Comparison_Int,
                                cols=dxy_Americana_Borealis:dxy_Novamexicana_Virilis,
                                names_to="Comparisons",
                                values_to="dxy")

AutSex_Comparison_Int<-AutSex_Comparison_Int %>%
  group_by(SexAut,Comparisons) %>%
  mutate(Average_dxy=mean(dxy)) %>%
  select(Comparisons, Average_dxy, SexAut) %>%
  distinct(Average_dxy)

wilcox.test(Average_dxy ~ SexAut, AutSex_Comparison_Int)

ggplot(AutSex_Comparison, aes(x=SexAut, y=Average_dxy, fill=SexAut))+
  geom_boxplot()+theme_bw()+labs(x='Chromosome', y='Dxy', fill='')

ggplot(AutSex_Comparison_Int, aes(x=Comparisons,y=Average_dxy, colour=SexAut))+
  geom_point(alpha=0.8, size=2)+theme_bw()+geom_line(colour='black')+
  coord_flip()+labs(y='Dxy', colour='')+
  scale_color_hue(direction = -1)+
  theme(text = element_text(size=10))



#################### Making the big table #################

Big_table_with_averaged_dxy_across_chroms <- read.csv('DXY_genomewide.csv')
Big_table_all_info$New_dxy <- Genic_SimonMartin_sumstats_dxymean$value[ match(Big_table_all_info$Comparisons, Genic_SimonMartin_sumstats_dxymean$variable)]


#Get just the dxy means
Genic_SimonMartin_sumstats_dxymean<-Genic_SimonMartin_sumstats %>%
  select(starts_with("dxy_")) %>%
  summarise(across(everything(), list(mean)))

Intergenic_SimonMartin_sumstats_dxymean<-Intergenic_SimonMartin_sumstats %>%
  select(starts_with("dxy_")) %>%
  summarise(across(everything(), list(mean)))

#melt these tables
Genic_SimonMartin_sumstats_dxymean<-melt(Genic_SimonMartin_sumstats_dxymean)
Intergenic_SimonMartin_sumstats_dxymean<-melt(Intergenic_SimonMartin_sumstats_dxymean)

#remove the annoying _1 that's been added to the variable header
Genic_SimonMartin_sumstats_dxymean$variable<-gsub("_1","",Genic_SimonMartin_sumstats_dxymean$variable)
Intergenic_SimonMartin_sumstats_dxymean$variable<-gsub("_1","",Intergenic_SimonMartin_sumstats_dxymean$variable)

#match column headers and add dxy_means to the big table.
Big_table_with_averaged_dxy_across_chroms$Dxy_coding <- Genic_SimonMartin_sumstats_dxymean$value[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, Genic_SimonMartin_sumstats_dxymean$variable)]
Big_table_with_averaged_dxy_across_chroms$Dxy_noncoding <- Intergenic_SimonMartin_sumstats_dxymean$value[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, Intergenic_SimonMartin_sumstats_dxymean$variable)]

#######Now do FST

Genic_SimonMartin_sumstats_fstmean<-Genic_SimonMartin_sumstats %>%
  select(starts_with("Fst_")) %>%
  summarise(across(everything(), list(mean)))

Intergenic_SimonMartin_sumstats_fstmean<-Intergenic_SimonMartin_sumstats %>%
  select(starts_with("Fst_")) %>%
  summarise(across(everything(), list(mean)))





#####################
#melt these tables
Genic_SimonMartin_sumstats_fstmean<-melt(Genic_SimonMartin_sumstats_fstmean)
Intergenic_SimonMartin_sumstats_fstmean<-melt(Intergenic_SimonMartin_sumstats_fstmean)

#remove the annoying _1 that's been added to the variable header
Genic_SimonMartin_sumstats_fstmean$variable<-gsub("_1","",Genic_SimonMartin_sumstats_fstmean$variable)
Intergenic_SimonMartin_sumstats_fstmean$variable<-gsub("_1","",Intergenic_SimonMartin_sumstats_fstmean$variable)

Big_table_with_averaged_dxy_across_chroms$Comparisons<-gsub("dxy_","Fst_",Big_table_with_averaged_dxy_across_chroms$Comparisons)
#match column headers and add dxy_means to the big table.
Big_table_with_averaged_dxy_across_chroms$Fst_coding <- Genic_SimonMartin_sumstats_fstmean$value[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, Genic_SimonMartin_sumstats_fstmean$variable)]
Big_table_with_averaged_dxy_across_chroms$Fst_noncoding <- Intergenic_SimonMartin_sumstats_fstmean$value[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, Intergenic_SimonMartin_sumstats_fstmean$variable)]


Big_table_with_averaged_dxy_across_chroms$Comparisons<-gsub("Fst_","",Big_table_with_averaged_dxy_across_chroms$Comparisons)
##Now add Dstatistic 
#change back to dxy_ to match for dstats
Big_table_with_averaged_dxy_across_chroms$Comparisons<-gsub("dxy_","",Big_table_with_averaged_dxy_across_chroms$Comparisons)
Dmin_rmdups$Interaction<-gsub("Mon","Montana",Dmin_rmdups$Interaction)
Dmin_rmdups
Big_table_with_averaged_dxy_across_chroms$Dstatistic <- Dmin_rmdups$Dstatistic[ match(Big_table_with_averaged_dxy_across_chroms$Comparisons, Dmin_rmdups$Interaction)]
Big_table_with_averaged_dxy_across_chroms$PrematingIso[is.na(Big_table_with_averaged_dxy_across_chroms$PrematingIso)] <- "Not estimated"


Big_table<-Big_table_with_averaged_dxy_across_chroms %>%
  select(Comparisons, Crossability, PrematingIso, Biogeography, Dxy_coding, Dxy_noncoding, Fst_coding
         , Fst_noncoding, Dstatistic)



write.csv(Big_table, file='Big_table_virilis_group_with_divergence_geneflow.csv')



###GET fd outlier positional information 
write.table(BFM_fd_outliers, file='BFM_fd_outliers.txt', sep='\t')
write.table(MLF_fd_outliers, file='MLF_fd_outliers.txt', sep='\t')
write.table(KEL_fd_outliers, file='KEL_fd_outliers.txt', sep='\t')
write.table(NAL_fd_outliers, file='NAL_fd_outliers.txt', sep='\t')








################# Biogeography ####################

big_table <- read.csv('Big_table_virilis_group_with_divergence_geneflow.csv', na = c("Not estimated", "", "Not tested"))
X_A_$Comparisons<-gsub("-"," vs. ",X_A_Ratio$Comparisons)
big_table$X_A <- X_A_Ratio$ratio[ match(X_A_Ratio$Comparisons, big_table$Comparisons)]
big_table<-big_table[!is.na(big_table$Biogeography), ]


ggplot(big_table, aes(x=big_table$Biogeography..Yukilevich.2014..Throckmorton.1982, y=big_table$TMRCA))+
  geom_boxplot()

ggplot(big_table, aes(x=X_A, y=reorder(Comparisons, TMRCA), colour=big_table$Biogeography..Yukilevich.2014..Throckmorton.1982.))+
  geom_point(size=3, alpha=0.75)+theme_minimal()+
  labs(y='Comparisons',
      x='X Chromosome/Autosome Dxy ratio', 
      colour='Biogeography')+
  xlim(0.7,1.4)+geom_vline(xintercept = 0.75, linetype='dashed')

#
png('dxycoding_sympatry_allopatry.png', width=5, height=2.5, units='in', res=320)
ggplot(big_table, aes(x=Biogeography, y=big_table$dXY..coding., fill=Biogeography, colour=Biogeography))+
  geom_boxplot(alpha=0.6)+geom_jitter(size=2, alpha=0.6)+theme_bw()+
  labs(x='Biogeography',
       y='Dxy (coding)')+ylim(0,0.08)+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))
  
dev.off()

wilcox.test(DXY ~ Geog, big_table)

ggplot(big_table, aes(x=Comparisons, y=X_A_Ratio, colour=Biogeography))+


##################### Species pairs X_A ratio
spp_pairs_info<-read.csv('Species_pairs_allinfo.csv')
spp_pairs_info$Species_pairs<-gsub("_"," vs. ",spp_pairs_info$Species_pairs)



png('dxyxa_divergencetime.png', width=7, height=5, units='in', res=320)
ggplot(spp_pairs_info, aes(x=Species_pairs, y=Divergence_time, colour=X_A)) +
  geom_pointrange(aes(ymin = Div_time_lower, ymax = Div_time_upper), size=0.8)+coord_flip()+
  theme_bw()+labs(x='Species pairs',
                               y= 'Divergence time (MYA)',
                               colour='X / Autosome dxy (coding)')+scale_colour_gradient(low = "grey", high = "black")+
  theme_minimal()+theme(axis.text=element_text(size=10, face='bold'), axis.title=element_text(size=12,face="bold"))+ theme(legend.position = c(0.85, 0.2))
dev.off()


ggplot(big_table, aes(x=big_table$Biogeography..Yukilevich.2014..Throckmorton.1982., y=X_A, fill=big_table$Biogeography..Yukilevich.2014..Throckmorton.1982., colour=big_table$Biogeography..Yukilevich.2014..Throckmorton.1982.))+
  geom_boxplot(alpha=0.6)+geom_jitter(size=2, alpha=0.6)+theme_bw()+
  labs(x='Biogeography',
       y='X/A dxy ratio)')+
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))


PreIso <-big_table$Premating.Isolation..Yukilevich.2014.
Geog <- big_table$Biogeography..Yukilevich.2014..Throckmorton.1982.
DXY <- big_table$dXY..coding.

sym_allo_model<-lm(PreIso ~ Geog + DXY + Geog:DXY, data = big_table)

png('qqplots1.png', width=4, height=3, units='in', res=320)
plot(sym_allo_model, add.smooth = FALSE, which = 1)
dev.off()
png('qqplots2.png', width=4, height=3, units='in', res=320)
plot(sym_allo_model, which = 2)
dev.off()
png('qqplots3.png', width=4, height=3, units='in', res=320)
plot(sym_allo_model, add.smooth = FALSE, which = 3)
dev.off
png('hist_symallo.png', width=4, height=3, units='in', res=320)
hist(sym_allo_model$residuals)
dev.off()

png('qqplots.png', width=9, height=5, units='in', res=320)
qq_1
ggarrange(qq_1, qq_2, qq_3)
dev.off()

anova(sym_allo_model)
cor.test(PreIso, DXY, method='pearson')

png('Preisolating_DXY_Biogeog.png', width=7, height=5, units='in', res=320)
ggplot(big_table, aes(x=DXY, y=PreIso, colour=Geog))+
  geom_point(size=3)+theme_bw()+geom_smooth(method='lm')+
  labs(x='Dxy',
       y='Premating Isolation',
       colour='Biogeography')+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=14,face="bold"))
dev.off()




#################### Showing reticulated evolution ##############

Outliers_ML_dxy <- Genic_SimonMartin_sumstats %>%
  filter(quant %in% c('1','2')) %>%
  select(scaffold, start, end)

Outlier_fdm_MLF <- MLF_fdm %>%
  filter(f_dM > 0.1) %>%
  select(chr, windowStart, windowEnd)

Outlier_fdm_MLF_negative <- MLF_fdm %>%
  filter(f_dM < -0.1) %>%
  select(chr, windowStart, windowEnd)
  
write.table(Outliers_ML_dxy, file="Outlier_ML_dxy.txt", sep="\t", col.names = F, row.names = F)
write.table(Outlier_fdm_MLF, file="Outlier_MLF_fdm.txt", sep="\t", col.names = F, row.names = F)
write.table(Outlier_fdm_MLF_negative, file="Outlier_MLF_fdm_negative.txt", sep="\t", col.names = F, row.names = F)



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
which( colnames(Genic_SimonMartin_sumstats)=="dxy_Americana_Novamexicana" )
MLF_dxy<-aggregate(Genic_SimonMartin_sumstats[, 59], list(Genic_SimonMartin_sumstats$scaffold), mean)
names(MLF_dxy) <- c("scaffolds", "dxy")
BFM_dxy<-aggregate(Genic_SimonMartin_sumstats[, 28], list(Genic_SimonMartin_sumstats$scaffold), mean)
names(BFM_dxy) <- c("scaffolds", "dxy")
KEL_dxy<-aggregate(Genic_SimonMartin_sumstats[, 37], list(Genic_SimonMartin_sumstats$scaffold), mean)
names(KEL_dxy) <- c("scaffolds", "dxy")
NAL_dxy<-aggregate(Genic_SimonMartin_sumstats[, 25], list(Genic_SimonMartin_sumstats$scaffold), mean)
names(NAL_dxy) <- c("scaffolds", "dxy")


MLF_dxy_inter<-aggregate(Intergenic_SimonMartin_sumstats[, 59], list(Intergenic_SimonMartin_sumstats$scaffold), mean)
names(MLF_dxy_inter) <- c("scaffolds", "dxy")
BFM_dxy_inter<-aggregate(Intergenic_SimonMartin_sumstats[, 28], list(Intergenic_SimonMartin_sumstats$scaffold), mean)
names(BFM_dxy_inter) <- c("scaffolds", "dxy")
KEL_dxy_inter<-aggregate(Intergenic_SimonMartin_sumstats[, 37], list(Intergenic_SimonMartin_sumstats$scaffold), mean)
names(KEL_dxy_inter) <- c("scaffolds", "dxy")
NAL_dxy_inter<-aggregate(Intergenic_SimonMartin_sumstats[, 25], list(Intergenic_SimonMartin_sumstats$scaffold), mean)
names(NAL_dxy_inter) <- c("scaffolds", "dxy")


MLF_dxy_fd_inter<-merge(MLF_dxy_inter,MLF_fd_mean,
                        by=c('scaffolds'), all=T)
BFM_dxy_fd_inter<-merge(BFM_dxy_inter,BFM_fd_mean,
                        by=c('scaffolds'), all=T)
KEL_dxy_fd_inter<-merge(KEL_dxy_inter,KEL_fd_mean,
                        by=c('scaffolds'), all=T)
NAL_dxy_fd_inter<-merge(NAL_dxy_inter,NAL_fd_mean,
                        by=c('scaffolds'), all=T)


MLF_dxy_fd<-merge(MLF_dxy,MLF_fd_mean,
                  by=c('scaffolds'), all=T)
BFM_dxy_fd<-merge(BFM_dxy,BFM_fd_mean,
                  by=c('scaffolds'), all=T)
KEL_dxy_fd<-merge(KEL_dxy,KEL_fd_mean,
                  by=c('scaffolds'), all=T)
NAL_dxy_fd<-merge(NAL_dxy,NAL_fd_mean,
                  by=c('scaffolds'), all=T)



MLF_dxy_fd <- na.omit(MLF_dxy_fd$fd)
BFM_dxy_fd <- na.omit(BFM_dxy_fd)
KEL_dxy_fd <- na.omit(KEL_dxy_fd)
NAL_dxy_fd <- na.omit(NAL_dxy_fd)

MLF_dxy_fd_inter <- na.omit(MLF_dxy_fd_inter)
BFM_dxy_fd_inter <- na.omit(BFM_dxy_fd_inter)
KEL_dxy_fd_inter <- na.omit(KEL_dxy_fd_inter)
NAL_dxy_fd_inter <- na.omit(NAL_dxy_fd_inter)


genic_cor_dxy_fd_mlf<-ggplot(MLF_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (M-L)))+
  labs(y=expression('f'[d] (M-L-F)))



genic_cor_dxy_fd_bfm<-ggplot(BFM_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (B-F)))+
  labs(y=expression('f'[d] (B-F-M)))



genic_cor_dxy_fd_kel<-ggplot(KEL_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (K-E)))+
  labs(y=expression('f'[d] (K-E-L)))



genic_cor_dxy_fd_nal<-ggplot(NAL_dxy_fd, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='blue', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (N-A)))+
  labs(y=expression('f'[d] (N-A-L)))

png('genic_correlations_introgression_dxy.png', width=6.25, height=4, units='in', res=320)
ggarrange(genic_cor_dxy_fd_mlf, genic_cor_dxy_fd_bfm, genic_cor_dxy_fd_kel,
          genic_cor_dxy_fd_nal, ncol=2, nrow=2)
dev.off()



inter_cor_dxy_fd_ml<-ggplot(MLF_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (M-L)))+
  labs(y=expression('f'[d] (M-L-F)))

inter_cor_dxy_fd_bf<-ggplot(BFM_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (B-F)))+
  labs(y=expression('f'[d] (B-F-M)))

inter_cor_dxy_fd_ke<-ggplot(KEL_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (K-E)))+
  labs(y=expression('f'[d] (K-E-L)))

inter_cor_dxy_fd_na<-ggplot(NAL_dxy_fd_inter, aes(x=dxy,y=fd))+
  geom_point(size=3, colour='black', pch=21, fill='red', alpha=0.6)+
  theme_bw()+stat_cor(method = "pearson")+
  labs(x=expression('d'[XY] (N-A)))+
  labs(y=expression('f'[d] (N-A-L)))



png('intergenic_correlations_introgression_dxy.png', width=6.25, height=4, units='in', res=320)
ggarrange(inter_cor_dxy_fd_ml, inter_cor_dxy_fd_bf, inter_cor_dxy_fd_ke,
          inter_cor_dxy_fd_na, ncol=2, nrow=2)
dev.off()





############ TESTING HYPOTHESIS THAT INTROGRESSION INTO BOREALIS AND FLAVOMONTANA cause increase in divergence
ggplot(Intergenic_SimonMartin_sumstats, aes(x=Chromosome, y=dxy_Flavomontana_Lacicola, fill=Chromosome))+
  geom_violin(position='dodge',size=0.75, width=0.5, alpha=0.75)+theme_pubr()+
  geom_boxplot(width=.1, position='dodge', coef = 6)+scale_fill_viridis(option='D', discrete=T)+
  labs(y=expression('d'[XY] (M-L)))+theme(axis.text=element_text(size=11), axis.title=element_text(size=11,face="bold"))+
  theme(legend.position='none')

Intergenic_SimonMartin_sumstats %>% group_by(Chromosome) %>%
  summarise(meanCDR = mean(dxy_Flavomontana_Lacicola)) 





chrom_f_d_bfm<-BFM_fd %>% group_by(Chromosome) %>%
  summarise('mean(BFM)' = mean(f_d),
            'SD(BFM)'= sd(f_d)) 

chrom_f_d_mlf<-MLF_fd %>% group_by(Chromosome) %>%
  summarise('mean(MLF)' = mean(f_d), 
            'SD(MLF)'= sd(f_d)) 

chrom_f_d_kel<-KEL_fd %>% group_by(Chromosome) %>%
  summarise('mean(KEL)' = mean(f_d),
            'SD(KEL)'= sd(f_d))

chrom_f_d_nal<-NAL_fd %>% group_by(Chromosome) %>%
  summarise('mean(NAL)' = mean(f_d),
            'SD(NAL)'= sd(f_d)) 


Altogether_chrom_fd_info<-cbind(chrom_f_d_bfm, chrom_f_d_mlf, chrom_f_d_kel, chrom_f_d_nal)
Altogether_chrom_fd_info[[4]] <- NULL
Altogether_chrom_fd_info[[6]] <- NULL
Altogether_chrom_fd_info[[8]] <- NULL
Altogether_chrom_fd_info<-na.omit(Altogether_chrom_fd_info)

sd_mean_chrom_fd_gt<-gt(Altogether_chrom_fd_info)

gtsave(sd_mean_chrom_fd_gt, filename = 'fd_chromosome_mean_fd_gt_table.png')


