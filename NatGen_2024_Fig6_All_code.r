library(ggplot2)
library(edgeR)
library(fgsea)

#FIGURE 6A
#drawing in illustrator



#FIGURE 6B
#drawing in illustrator



#FIGURE 6C


ggplot(Sim_1q1.0_vs_wt, aes(x=group, y=mean, fill = Notch_status)) +
    geom_bar(stat="identity", color = "white", position=position_dodge()) +
    
    scale_fill_manual(values=c("blue4", "yellow2")) 


ggplot(Sim_1q2.2_vs_wt, aes(x=group, y=mean, fill = Notch_status)) +
    geom_bar(stat="identity", color = "white", position=position_dodge()) +
    
    scale_fill_manual(values=c("blue4", "yellow2")) 


ggplot(Sim_1q1.0_vs_1q1.0, aes(x=group, y=mean, fill = Notch_status)) +
    geom_bar(stat="identity", color = "white", position=position_dodge()) +
    
    scale_fill_manual(values=c("blue4", "yellow2")) 





#FIGURE 6D

ggplot(Sim_combo, aes(x=group, y=mean, fill = Notch_status)) +
    geom_bar(stat="identity", color = "white") +
    facet_wrap(~poising_factor, nrow = 1) +
    scale_fill_manual(values=c("blue4", "yellow2")) +ylim(0,1)



#FIGURE 6E

ggplot(subset(Sim_combo_proportions, experiment != 3 & experiment != 2.8), aes(x=group, y=mean, fill = Notch_status)) +
    geom_bar(stat="identity", color = "black") +
    facet_wrap(~fraction_1q, nrow = 1) +
    scale_fill_manual(values=c("blue4", "yellow2")) +ylim(0,1) +theme_bw()


#FIGURE 6G
#coculture RNAseq experiment; processing and filtering raw reads based on RPKM

BRCA_RSEM_genes <- BRCA_RSEM[,c(1)]
BRCA_RSEM_genes <- data.frame(BRCA_RSEM_genes)
BRCA_RSEM_genes <- BRCA_RSEM_genes[-1,]
BRCA_RSEM_genes <- data.frame(BRCA_RSEM_genes)
BRCA_RSEM_genes <- separate(data = BRCA_RSEM_genes, col = BRCA_RSEM_genes, into = c('Gene', 'x'), sep = "\\|")
BRCA_RSEM_genes <- BRCA_RSEM_genes[-c(1:29),]

RNAseq_041121_counts <- RNAseq_041121
RNAseq_041121_counts$Length <- as.numeric(as.character(RNAseq_041121_counts$Length))
RNAseq_041121_counts$Length <- RNAseq_041121_counts$Length/1000
rownames(RNAseq_041121_counts) <- RNAseq_041121_counts$Geneid
RNAseq_041121_counts_geninfo <- RNAseq_041121_counts[,c(1:6)]

for(i in 7:ncol(RNAseq_041121_counts)){
RNAseq_041121_counts[,i] <- as.numeric(as.character(RNAseq_041121_counts[,i]))
}

RNAseq_041121_counts2 <- subset(RNAseq_041121_counts, Geneid %in% BRCA_RSEM_genes$Gene)

sample_counts_tot <- colSums(RNAseq_041121_counts2[,c(7:ncol(RNAseq_041121_counts2))])
sample_counts_tot <- sample_counts_tot/1000000
RNAseq_041121_counts2[,c(7:ncol(RNAseq_041121_counts2))] <- RNAseq_041121_counts2[,c(7:ncol(RNAseq_041121_counts2))]/sample_counts_tot
for(i in 7:ncol(RNAseq_041121_counts2)) {
RNAseq_041121_counts2[,i] <- RNAseq_041121_counts2[,i]/RNAseq_041121_counts2$Length
}
RNAseq_041121_counts2_lcrm <- RNAseq_041121_counts2
RNAseq_041121_counts2_lcrm$AVG <- rowMeans(RNAseq_041121_counts2_lcrm[,c(7:ncol(RNAseq_041121_counts2_lcrm))])
RNAseq_041121_counts2_lcrm <- subset(RNAseq_041121_counts2_lcrm, RNAseq_041121_counts2_lcrm$AVG > 2)

#edgeR 

RNAseq_2020_sva_codes4_temp <- subset(coculture_experiments_samples, (color_collapse == "1q_mono" | color_collapse == "1q_co_dip") )
RNAseq_2020_sva2_temp2_dblNORM_temp <- RNAseq_041121_counts2_lcrm[,names(RNAseq_041121_counts2_lcrm) %in% RNAseq_2020_sva_codes4_temp$id]
yEndtoEnd2 <- DGEList(counts=RNAseq_2020_sva2_temp2_dblNORM_temp, genes=rownames(RNAseq_2020_sva2_temp2_dblNORM_temp))
labels <- data.frame("Name"=colnames(RNAseq_2020_sva2_temp2_dblNORM_temp), "TimePoint"= RNAseq_2020_sva_codes4_temp$color_collapse2)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_1q_codip_vs_1q_mono_RPKM2 <- lrtres
limma_1q_codip_vs_1q_mono_RPKM2[limma_1q_codip_vs_1q_mono_RPKM2 == "NaN"] <- 0
ranks <- limma_1q_codip_vs_1q_mono_RPKM2$neglogP
names(ranks) <- limma_1q_codip_vs_1q_mono_RPKM2$genes
plotEnrichment(pathwaysHx[["proNotch"]], ranks)
fgsealimma_1q_codip_vs_1q_mono_RPKM2 <- fgsea(Notch_curated_sets31_NOGS3_18_MSprotBreastGEX4_new12b.gmx, ranks, minSize=1, maxSize = 500, nperm=10000)



#FIGURE 6H

RNAseq_2020_sva_codes4_temp <- subset(coculture_experiments_samples, (color_collapse == "1q_mono" | color_collapse == "1q_co_wt") )
RNAseq_2020_sva2_temp2_dblNORM_temp <- RNAseq_041121_counts2_lcrm[,names(RNAseq_041121_counts2_lcrm) %in% RNAseq_2020_sva_codes4_temp$id]
yEndtoEnd2 <- DGEList(counts=RNAseq_2020_sva2_temp2_dblNORM_temp, genes=rownames(RNAseq_2020_sva2_temp2_dblNORM_temp))
labels <- data.frame("Name"=colnames(RNAseq_2020_sva2_temp2_dblNORM_temp), "TimePoint"= RNAseq_2020_sva_codes4_temp$color_collapse2)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_1q_cowt_vs_1q_mono_RPKM2 <- lrtres
limma_1q_cowt_vs_1q_mono_RPKM2[limma_1q_cowt_vs_1q_mono_RPKM2 == "NaN"] <- 0
ranks <- limma_1q_cowt_vs_1q_mono_RPKM2$neglogP
names(ranks) <- limma_1q_cowt_vs_1q_mono_RPKM2$genes
plotEnrichment(pathwaysHx[["proNotch"]], ranks)
fgseallimma_1q_cowt_vs_1q_mono_RPKM2 <- fgsea(Notch_curated_sets31_NOGS3_18_MSprotBreastGEX4_new12b.gmx, ranks, minSize=1, maxSize = 500, nperm=10000)


#plots color bars
my_palette3 <- myPalette(low = "blue", high = "red", mid="white", k =200)
limma_1q_cowt_vs_1q_mono_RPKM2 <- limma_1q_cowt_vs_1q_mono_RPKM2[order(-limma_1q_cowt_vs_1q_mono_RPKM2$neglogP),]
limma_1q_cowt_vs_1q_mono_RPKM2$genes2 <- factor(limma_1q_cowt_vs_1q_mono_RPKM2$genes, levels = limma_1q_cowt_vs_1q_mono_RPKM2$genes)
ggplot(limma_1q_cowt_vs_1q_mono_RPKM2, aes(x=genes2, y=1)) + geom_tile(aes(fill = neglogP)) + scale_fill_gradientn(colours = my_palette3,limits=c(-3.2,3.2)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

limma_1q_codip_vs_1q_mono_RPKM2 <- limma_1q_codip_vs_1q_mono_RPKM2[order(-limma_1q_codip_vs_1q_mono_RPKM2$neglogP),]
limma_1q_codip_vs_1q_mono_RPKM2$genes2 <- factor(limma_1q_codip_vs_1q_mono_RPKM2$genes, levels = limma_1q_codip_vs_1q_mono_RPKM2$genes)
ggplot(limma_1q_codip_vs_1q_mono_RPKM2, aes(x=genes2, y=1)) + geom_tile(aes(fill = neglogP)) + scale_fill_gradientn(colours = my_palette3,limits=c(-3.2,3.2)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())


#FIGURE 6I

library(ggbeeswarm)
my_comparisons3 <- list( c("0wt-wt", "1q-wt"))
ggplot(subset(competition_assay_summary2_and_3, pair_type3 != "1q-1q0" & pair_type3 !="wt-wt0"), aes(x=pair_type3, y=diff_avg2)) + geom_boxplot() + geom_beeswarm(alpha=0.2) + stat_compare_means(comparisons = my_comparisons3)  +theme_bw()



#FIGURE 6J

ggplot(DepMap_RNAi_CRISP_NotchPathway_SynLet_ARM_ampdels_gt10perc2, aes(x = pro_Notch_NegLogP_RNAi_perm, y = pro_Notch_NegLogP_CRISP_perm, color = ampdel.x)) +geom_hline(yintercept=0) +geom_vline(xintercept=0) + geom_point() + geom_text_repel(data= DepMap_RNAi_CRISP_NotchPathway_SynLet_ARM_ampdels_gt10perc2, aes(label=ChrArm, color = ampdel.x)) +theme_bw() + geom_vline(xintercept=-2) + geom_hline(yintercept=-2) +ylim(-4,4) +xlim(-4,4)
