library(ggplot2)
library(Hmisc)
library(ggrepel)
library(data.table)
library(tidyverse)
library(fgsea)
library(edgeR)



#FIGURE 2A
#drawing in illustrator



#FIGURE 2B

ggplot(Fig2B_summary_evo_CNAS_wreps_absolutem_merge,aes(x=(factor(id, level = c('pop-a','pop-b',  'ctrl1','ctrlm-a','ctrlm-b','ctrlm-c','bq-a','bq-b','eh-a','eh-b','ae','af-a','af-b','af-c','df','dc-a','dc-b','da','ef-a','ef-b','ee-a','ee-b','ec-a','ec-b','ed-a','ed-b','ed-c','eg-a','eg-b','ea-a','ea-b','eb-a','eb-c','FS','FG-a','FG-b','FG-c','FA-a','FA-b','FA-c','FC-a','FC-b','FC-c','FE-a','FE-b','FE-c','FM-a','FM-b','FM-c','FV-a','FV-b','FV-c','FQ-a','FQ-b','FQ-c','FR-a','FR-b','FR-c','FY-a','FY-b','FY-c','FX-a','FX-b','FX-c','FN-a','FN-b','FN-c','FI','FD','CQ-L-a','CQ-L-b','CQ-J','CQ-I-a','CQ-I-b','CQ-I-c','CQ-G-a','CQ-G-b','CQ-G-c','CQ-E-a','CQ-E-b','CQ-E-c','CQ-F','CQ-M','CQ-C-a','CQ-C-b','CQ-C-c','CQ-H-a','CQ-H-b','CQ-H-c'))),y=factor(variable, level=rev(c('X1q','X8q', 'X16p', 'X20q',   'X5p', 'X7p', 'X3q', 'X20p', 'X17q', 'X7q',  'X10p', 'X6p','X19q',  'X12p','X2p', 'X21q', 'Xp', 'X12q', 'Xq',  'X19p', 'X2q', 'X11p', 'X5q', 'X14q', 'X18p', 'X9q', 'X9p',   'X1p',       'X4q', 'X15q', 'X10q', 'X6q',  'X3p', 'X18q', 'X4p',         'X13q', 'X11q', 'X22q', 'X8p',   'X17p' , 'X16q'))), fill=as.numeric(original_state) )) + geom_tile(colour="black",size=0.25) + scale_fill_gradientn(colours=c("blue","lightblue", "grey98","brown1", "red", "firebrick"), breaks=c(-2,-1,0,1,2,4), na.value = "grey98",limits = c(-2, 4)) + theme_bw() + geom_point(data=Fig2B_summary_evo_CNAS_wreps_absolutem_merge, aes(x=id, y=variable, shape=as.factor(evo_direction), color=as.numeric(evo_final_state)), size=2) + scale_color_gradientn(colours=c("blue","lightblue", "grey98","brown1", "red", "firebrick"), breaks=c(-2,-1,0,1,2,4), na.value = "grey98",limits = c(-2, 4)) + scale_shape_manual(values=c(25, 25,24,24,24)) + theme(axis.text.x = element_text(angle = 45, hjust=1, size=8))



#FIGURE 2C

ggplot(all_by_all_corr_tissues_HMEC_TrueArms_xxx2, aes(x= HMEC_amp_minus_del_arms , y = BRCA_amp_minus_del_arms)) +geom_point(alpha=0.7)  + geom_smooth(method="lm",se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 3)) +theme_bw() +geom_text_repel(data=subset(all_by_all_corr_tissues_HMEC_TrueArms_xxx2,  Chromosome_arm == "8q" | Chromosome_arm == "20q" | Chromosome_arm == "1q" | Chromosome_arm == "16p"), aes(label = Chromosome_arm), color = "red", size=4) +geom_text_repel(data=subset(all_by_all_corr_tissues_HMEC_TrueArms_xxx2,  Chromosome_arm == "8p" | Chromosome_arm == "22q" | Chromosome_arm == "17p" | Chromosome_arm == "16q"), aes(label = Chromosome_arm), color = "blue", size=4) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)



#FIGURE 2D

my_corp <- rcorr(as.matrix(all_by_all_corr_tissues_HMEC_TrueArms_xxx2[,-c(1)]),type="pearson")
my_corp <- my_corp$P
my_corp[is.na(my_corp)] <- 1
my_corp <- -log10(my_corp)
my_corp <- data.frame(my_corp)
HMEC_vs_Tissues_TrueArms_corp <- my_corp$HMEC_amp_minus_del_arms
HMEC_vs_Tissues_TrueArms_corp <- data.frame(HMEC_vs_Tissues_TrueArms_corp)
HMEC_vs_Tissues_TrueArms_corp$Tissue <- rownames(my_corp)
ggplot(HMEC_vs_Tissues_TrueArms_corp, aes(x=Tissue, y=1, width = 1)) + geom_tile(aes(fill = HMEC_vs_Tissues_TrueArms_corp )) + scale_fill_gradient2(low = "white", mid = "khaki", high = "red",  midpoint = 3, limits = c(0,6.1)) +theme_bw()



#FIGURE 2E


BRCA_CNA_summaryNEW_0.5 <- ChrArm_CNA_freq_TCGA_puritycorrection(BRCA_CNA, BRCA_purity, "BRCA_CNAs_0.5_-0.41_0.32", 0.5, -0.41, 0.32, chromdata_peri_cent2)
my.list16pfwa_BRCANEW0.5_CNorm0.7 <- Arm_Level_CNorm_diffexp2(BRCA_RSEM, BRCA_purity, BRCA_CNA, BRCA_CNA_summaryNEW_0.5, 10, 0.7, "16p", "gain", 5, 30)

ranks <- my.list16pfwa_BRCANEW0.5_CNorm0.7[[1]]$neglogP
names(ranks) <- my.list16pfwa_BRCANEW0.5_CNorm0.7[[1]]$genes
fgseamy.list16pfwa_BRCANEW <- fgsea(RP_immune_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

my_breaks2 <- c(-4, -3.5,-3, -2.5, 0, 2.5, 3, 3.5, 4 )
my_palette3 <- myPalette(low = "blue4", high = "red3", mid="white", k =8)
ggplot(BRCA_CNorm0.7_fgsea_RPimmunesets_TopChrArms, aes(x=ChrArm, y=gene_set)) + geom_tile(aes(fill = neglogPadj)) + scale_fill_gradientn(colours = my_palette3,breaks = my_breaks2, limits = c(-3.3,3.3) )



