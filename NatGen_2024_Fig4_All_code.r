library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(fgsea)
library(edgeR)


#FIGURE 4A

ggplot(subset(dfHMEC, strain !="diploid" ), aes(x=frequency_score, y=growth_rate)) +
    geom_smooth(method="lm", se=F, color = "black") +
    theme_bw() +
    stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat((r.squared)), after_stat(p.value), sep="\n" ),size = 2)) +
    geom_point(data = subset(dfHMEC, strain !="diploid" ), size = 3, alpha=0.5 )+
    ylim(0,1.3)+
    xlim(0,0.22)+
    geom_hline(yintercept = 1.18)+
    geom_hline(yintercept = 1.18 + 0.0765)+
    geom_hline(yintercept = 1.18 - 0.0765)



#FIGURE 4B

ggplot(subset(dfRPTEC, strain !="diploid" ), aes(x=frequency_score, y=growth_rate)) +
    geom_smooth(method="lm", se=F, color = "black") +
    theme_bw() +
    stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2))  +
    geom_point(data = subset(dfRPTEC, strain !="diploid" ), size = 3, alpha=0.5)+
    ylim(0,1.3)+
    xlim(0,0.3)+
    geom_hline(yintercept = 0.833)+
    geom_hline(yintercept = 0.833 + 0.0264) +
    geom_hline(yintercept = 0.833 - 0.0264)




#FIGURE 4C

ggplot(growth_rate_combined4, aes(x = status1, y=growth_rate, fill=status1)) +     geom_boxplot() + geom_hline(yintercept = 1.28) +  theme_bw() + scale_fill_manual(values = c("goldenrod2", "maroon")) + facet_grid(. ~ strain5) +
    stat_compare_means(comparisons = list(c("evo", "orig"))) +geom_point(alpha=0.4) +ylim(0,1.5)




#FIGURE 4D

ggplot(time_to_clonal_takeover, aes(x=orig_growth_rate, y=log2FC_growth)) + geom_point(aes(color = time_to_subclone), size = 3) +  geom_smooth(method="lm", se=F, color = "black", size = 0.5) + theme_bw() + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) + scale_colour_gradientn(colours=c("yellow","orange", "red3"))




#FIGURE 4E

RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2020_sva_codes7, evo == "evo" | class == "diploid")
RNAseq_redo_010222d_temp <- RNAseq_redo_010222a2_RPKM_lcrm[,names(RNAseq_redo_010222a2_RPKM_lcrm) %in% RNAseq_2020_sva_codes2_temp$unq_id]
yEndtoEnd2 <- DGEList(counts=RNAseq_redo_010222d_temp, genes=rownames(RNAseq_redo_010222d_temp))
labels <- data.frame("Name"=colnames(RNAseq_redo_010222d_temp), "TimePoint"=RNAseq_2020_sva_codes2_temp$class)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(1,-1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_HMEC_Aneup_evo_vs_diploid_RPKM3 <- lrtres

colnames(fgsealimma_HMEC_Aneup_orig_vs_diploid_RPKM3) <- paste(colnames(fgsealimma_HMEC_Aneup_orig_vs_diploid_RPKM3), "evo", sep = '_')
RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2020_sva_codes7, evo == "orig" | class == "diploid")
RNAseq_redo_010222d_temp <- RNAseq_redo_010222a2_RPKM_lcrm[,names(RNAseq_redo_010222a2_RPKM_lcrm) %in% RNAseq_2020_sva_codes2_temp$unq_id]
yEndtoEnd2 <- DGEList(counts=RNAseq_redo_010222d_temp, genes=rownames(RNAseq_redo_010222d_temp))
labels <- data.frame("Name"=colnames(RNAseq_redo_010222d_temp), "TimePoint"=RNAseq_2020_sva_codes2_temp$class)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(1,-1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_HMEC_Aneup_orig_vs_diploid_RPKM3 <- lrtres

limma_HMEC_Aneup_evo_vs_diploid_RPKM3[limma_HMEC_Aneup_evo_vs_diploid_RPKM3 == "NaN"] <- 0
limma_HMEC_Aneup_orig_vs_diploid_RPKM3[limma_HMEC_Aneup_orig_vs_diploid_RPKM3 == "NaN"] <- 0
limma_HMEC_Aneup_evo_vs_diploid_RPKM3[limma_HMEC_Aneup_evo_vs_diploid_RPKM3 == "NaN"] <- 0

ranks <- limma_HMEC_Aneup_evo_vs_diploid_RPKM3$neglogP
names(ranks) <- limma_HMEC_Aneup_evo_vs_diploid_RPKM3$genes
fgsealimma_HMEC_Aneup_evo_vs_diploid_RPKM3 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)
limma_HMEC_Aneup_orig_vs_diploid_RPKM3[limma_HMEC_Aneup_orig_vs_diploid_RPKM3 == "NaN"] <- 0
ranks <- limma_HMEC_Aneup_orig_vs_diploid_RPKM3$neglogP
names(ranks) <- limma_HMEC_Aneup_orig_vs_diploid_RPKM3$genes
fgsealimma_HMEC_Aneup_orig_vs_diploid_RPKM3 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)
fgsealimma_HMEC_Aneup_orig_vs_diploid_RPKM3$evo <- "orig"
fgsealimma_HMEC_Aneup_evo_vs_diploid_RPKM3$evo <- "evo"
fgsealimma_HMEC_Aneup_evo_orig_vs_diploid_RPKM3 <- rbind(fgsealimma_HMEC_Aneup_orig_vs_diploid_RPKM3, fgsealimma_HMEC_Aneup_evo_vs_diploid_RPKM3)


ggplot(fgsealimma_HMEC_Aneup_evo_orig_vs_diploid_RPKM3, aes(x=pathway,y=neglogP, fill=evo)) +geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))