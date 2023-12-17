library(ggplot2)
library(fgsea)
library(ggpmisc)
library(ggpubr)

#Ext Data Figure 10a

ExtData_Fig10a_data$id <- factor(ExtData_Fig10a_data$id,levels = c("bq", "bn", "FY_ev2", "FC_ev2", "FM_ev2", "FN_ev3", "FG_ev1", "FR_ev3", "CQ_ev_I_ev3", "FS", "FD_ev_A", "CQ_ev_L_ev1", "FX_ev1_B"))
ExtData_Fig10a_data$variable <- factor(ExtData_Fig10a_data$variable,levels = c("bq", "bn", "FY_ev2", "FC_ev2", "FM_ev2", "FN_ev3", "FG_ev1", "FR_ev3", "CQ_ev_I_ev3", "FS", "FD_ev_A", "CQ_ev_L_ev1", "FX_ev1_B"))
ExtData_Fig10a_data$value <- ifelse(ExtData_Fig10a_data$value > 6.5, 6,5, ExtData_Fig10a_data$value)
ExtData_Fig10a_data$value <- ifelse(ExtData_Fig10a_data$value < -6, -6, ExtData_Fig10a_data$value)
ggplot(ExtData_Fig10a_data, aes(id, variable, fill=value)) +
geom_tile() +
scale_fill_viridis(discrete=FALSE, direction=1)

#Ext Data Figure 10b
#copy number models from Ext Data Figure 5
#please see BioProject: PRJNA634423 for gDNA-seq data corresponding to these samples
#see separate repository https://github.com/emmavwatson/CNAplot for all raw gDNA sequencing processing, copy number calling, and plotting code 



#Ext Data Figure 10c
my_comparisons2 <- list( c("1q-dip", "wt-dip"), c("1q-wt", "wt-wt"), c("wt-1q_orig", "1q-1q_orig"), c("1q-1q", "wt-1q") )
ggplot(ExtData_Fig10c_data, aes(x=pair_type3, y=diff_growth_APC2, color = pair_type3)) + geom_hline(yintercept = 0, alpha=0.3) + geom_boxplot()  + stat_compare_means(comparisons = my_comparisons2, na.rm=TRUE, method="t.test", aes(label = ..p.signif..))  +theme_bw() + scale_color_manual(values = c("blue","red","blue","red","blue","red","blue","red","blue","red","blue","red","blue","red","blue","red")) + facet_grid(cols = vars(condition))


#Ext Data Figure 10d

ranks <- ExtData_Fig10d_data$neglogP_RNAi
names(ranks) <- ExtData_Fig10d_data$genes
plotEnrichment(Notch_curated_sets[["proNotch"]], ranks)


ranks <- ExtData_Fig10d_data$neglogP_CRISPR
names(ranks) <- ExtData_Fig10d_data$genes
plotEnrichment(Notch_curated_sets[["proNotch"]], ranks)


#Ext Data Figure 10e
#conceptual illustration of predictions from Notch lateral inhibition model with respect to clonal takeover of +1q subclones