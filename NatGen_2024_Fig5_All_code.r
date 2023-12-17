library(edgeR)
library(fgsea)
library(pheatmap)
library(ggplot2)

#FIGURE 5A

#TCGA BRCA analysis; 8q and 1q

BRCA_CNA_summaryNEW_0.5_nopur <- ChrArm_CNA_freq_TCGA_nopuritycorrection(BRCA_CNA, "BRCA_CNAs_nopur_0.5_-0.41_0.32", 0.5, -0.41, 0.32, chromdata_peri_cent2)
my.list8qfwa_BRCANEW0.5_noCNorm_nopur <- Arm_Level_CNorm_diffexp2_nopurity(BRCA_RSEM, BRCA_CNA, BRCA_CNA_summaryNEW_0.5_nopur, 10, 1, "8q", "gain", 5, 30)
ranks <- my.list8qfwa_BRCANEW0.5_noCNorm_nopur[[1]]$neglogP
names(ranks) <- my.list8qfwa_BRCANEW0.5_noCNorm_nopur[[1]]$genes
fgsea_my.list8qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)
fgsea_my.list8qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks <- fgsea_my.list8qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks[,-8]
colnames(fgsea_my.list8qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks) <- c("pathway", "pval_BRCA_8q", "padj_BRCA_8q", "ES_BRCA_8q", "NES_BRCA_8q", "nMoreExtreme_BRCA_8q", "size_BRCA_8q" )

my.list1qfwa_BRCANEW0.5_noCNorm_nopur <- Arm_Level_CNorm_diffexp2_nopurity(BRCA_RSEM, BRCA_CNA, BRCA_CNA_summaryNEW_0.5_nopur, 10, 1, "1q", "gain", 5, 30)
ranks <- my.list1qfwa_BRCANEW0.5_noCNorm_nopur[[1]]$neglogP
names(ranks) <- my.list1qfwa_BRCANEW0.5_noCNorm_nopur[[1]]$genes
fgsea_my.list1qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)
fgsea_my.list1qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks <- fgsea_my.list1qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks[,-8]
colnames(fgsea_my.list1qfwa_BRCANEW0.5_noCNorm_nopur_Hallmarks) <- c("pathway", "pval_BRCA_1q", "padj_BRCA_1q", "ES_BRCA_1q", "NES_BRCA_1q", "nMoreExtreme_BRCA_1q", "size_BRCA_1q" )

#HMEC RNAseq analysis; +1q vs WT 1q and +8q vs WT 8q

RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2020_sva_codes5, class != "diploid" & class != "RPTEC" )
RNAseq_redo_010222d_temp <- RNAseq_redo_010222a2_RPKM_lcrm[,names(RNAseq_redo_010222a2_RPKM_lcrm) %in% RNAseq_2020_sva_codes2_temp$unq_id]
yEndtoEnd2 <- DGEList(counts=RNAseq_redo_010222d_temp, genes=rownames(RNAseq_redo_010222d_temp))
labels <- data.frame("Name"=colnames(RNAseq_redo_010222d_temp), "TimePoint"=RNAseq_2020_sva_codes2_temp$chr1q)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_HMEC_1q_new_RPKM3 <- lrtres
ranks <- limma_HMEC_1q_new_RPKM3$neglogP
names(ranks) <- limma_HMEC_1q_new_RPKM3$genes
fgsealimma_HMEC_1q_new_RPKM3 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2020_sva_codes5, class != "diploid" & class != "RPTEC" & chr8q_discrep != 1 )
RNAseq_redo_010222d_temp <- RNAseq_redo_010222a2_RPKM_lcrm[,names(RNAseq_redo_010222a2_RPKM_lcrm) %in% RNAseq_2020_sva_codes2_temp$unq_id]
yEndtoEnd2 <- DGEList(counts=RNAseq_redo_010222d_temp, genes=rownames(RNAseq_redo_010222d_temp))
labels <- data.frame("Name"=colnames(RNAseq_redo_010222d_temp), "TimePoint"=RNAseq_2020_sva_codes2_temp$chr8q)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_HMEC_8q_new_RPKM3 <- lrtres
ranks <- limma_HMEC_8q_new_RPKM3$neglogP
names(ranks) <- limma_HMEC_8q_new_RPKM3$genes
fgsealimma_HMEC_8q_new_RPKM3 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

#HMEC RNAseq analysis; +1q vs pre-evo 1q and +8q vs pre-evo 8q
RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2020_sva_codes5, class != "diploid" & class != "RPTEC" & id %like% "G23_C1_27")
RNAseq_redo_010222d_temp <- RNAseq_redo_010222a2_RPKM_lcrm[,names(RNAseq_redo_010222a2_RPKM_lcrm) %in% RNAseq_2020_sva_codes2_temp$unq_id]
yEndtoEnd2 <- DGEList(counts=RNAseq_redo_010222d_temp, genes=rownames(RNAseq_redo_010222d_temp))
labels <- data.frame("Name"=colnames(RNAseq_redo_010222d_temp), "TimePoint"=RNAseq_2020_sva_codes2_temp$chr1q)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27 <- lrtres

limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27 == "Inf"] <- 200
limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27 == "-Inf"] <- -200
limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27 == "NaN"] <- 0
ranks <- limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27$neglogP
names(ranks) <- limma_HMEC_1q_new_lineage_RPKM3_G23_C1_27$genes
fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_27 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22 == "Inf"] <- 200
limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22 == "-Inf"] <- -200
limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22 == "NaN"] <- 0
ranks <- limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22$neglogP
names(ranks) <- limma_HMEC_1q_new_lineage_RPKM3_G23_C1_22$genes
fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_22 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)


limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19 == "Inf"] <- 200
limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19 == "-Inf"] <- -200
limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19[limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19 == "NaN"] <- 0
ranks <- limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19$neglogP
names(ranks) <- limma_HMEC_1q_new_lineage_RPKM3_G23_C1_19$genes
fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_19 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)


fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG <- fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_19
fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG$pval <- (fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_19$pval + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_22$pval + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_27$pval )/3
fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG$padj <- (fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_19$padj + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_22$padj + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_27$padj )/3
fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG$ES <- (fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_19$ES + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_22$ES + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_27$ES )/3
fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG$NES <- (fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_19$NES + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_22$NES + fgsealimma_HMEC_1q_new_lineage_RPKM3_G23_C1_27$NES )/3

RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2020_sva_codes5, id == "102B" | id == "102B_evo3")
RNAseq_redo_010222d_temp <- RNAseq_redo_010222a2_RPKM_lcrm[,names(RNAseq_redo_010222a2_RPKM_lcrm) %in% RNAseq_2020_sva_codes2_temp$unq_id]
yEndtoEnd2 <- DGEList(counts=RNAseq_redo_010222d_temp, genes=rownames(RNAseq_redo_010222d_temp))
labels <- data.frame("Name"=colnames(RNAseq_redo_010222d_temp), "TimePoint"=RNAseq_2020_sva_codes2_temp$chr8q)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_HMEC_8q_new_lineage_RPKM3_102B <- lrtres #repeat for all lineages, then combine using averages


limma_HMEC_8q_new_lineage_RPKM3_102B[limma_HMEC_8q_new_lineage_RPKM3_102B == "Inf"] <- 200
limma_HMEC_8q_new_lineage_RPKM3_102B[limma_HMEC_8q_new_lineage_RPKM3_102B == "-Inf"] <- -200
limma_HMEC_8q_new_lineage_RPKM3_102B[limma_HMEC_8q_new_lineage_RPKM3_102B == "NaN"] <- 0
ranks <- limma_HMEC_8q_new_lineage_RPKM3_102B$neglogP
names(ranks) <- limma_HMEC_8q_new_lineage_RPKM3_102B$genes
fgsealimma_HMEC_8q_new_lineage_RPKM3_102B <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

limma_HMEC_8q_new_lineage_RPKM3_HMNRC4[limma_HMEC_8q_new_lineage_RPKM3_HMNRC4 == "Inf"] <- 200
limma_HMEC_8q_new_lineage_RPKM3_HMNRC4[limma_HMEC_8q_new_lineage_RPKM3_HMNRC4 == "-Inf"] <- -200
limma_HMEC_8q_new_lineage_RPKM3_HMNRC4[limma_HMEC_8q_new_lineage_RPKM3_HMNRC4 == "NaN"] <- 0
ranks <- limma_HMEC_8q_new_lineage_RPKM3_HMNRC4$neglogP
names(ranks) <- limma_HMEC_8q_new_lineage_RPKM3_HMNRC4$genes
fgsealimma_HMEC_8q_new_lineage_RPKM3_HMNRC4 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

limma_HMEC_8q_new_lineage_RPKM3_C05[limma_HMEC_8q_new_lineage_RPKM3_C05 == "Inf"] <- 200
limma_HMEC_8q_new_lineage_RPKM3_C05[limma_HMEC_8q_new_lineage_RPKM3_C05 == "-Inf"] <- -200
limma_HMEC_8q_new_lineage_RPKM3_C05[limma_HMEC_8q_new_lineage_RPKM3_C05 == "NaN"] <- 0
ranks <- limma_HMEC_8q_new_lineage_RPKM3_C05$neglogP
names(ranks) <- limma_HMEC_8q_new_lineage_RPKM3_C05$genes
fgsealimma_HMEC_8q_new_lineage_RPKM3_C05 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

limma_HMEC_8q_new_lineage_RPKM3_K16[limma_HMEC_8q_new_lineage_RPKM3_K16 == "Inf"] <- 200
limma_HMEC_8q_new_lineage_RPKM3_K16[limma_HMEC_8q_new_lineage_RPKM3_K16 == "-Inf"] <- -200
limma_HMEC_8q_new_lineage_RPKM3_K16[limma_HMEC_8q_new_lineage_RPKM3_K16 == "NaN"] <- 0
ranks <- limma_HMEC_8q_new_lineage_RPKM3_K16$neglogP
names(ranks) <- limma_HMEC_8q_new_lineage_RPKM3_K16$genes
fgsealimma_HMEC_8q_new_lineage_RPKM3_K16 <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

limma_HMEC_8q_new_lineage_RPKM3_0.4_1A[limma_HMEC_8q_new_lineage_RPKM3_0.4_1A == "Inf"] <- 200
limma_HMEC_8q_new_lineage_RPKM3_0.4_1A[limma_HMEC_8q_new_lineage_RPKM3_0.4_1A == "-Inf"] <- -200
limma_HMEC_8q_new_lineage_RPKM3_0.4_1A[limma_HMEC_8q_new_lineage_RPKM3_0.4_1A == "NaN"] <- 0
ranks <- limma_HMEC_8q_new_lineage_RPKM3_0.4_1A$neglogP
names(ranks) <- limma_HMEC_8q_new_lineage_RPKM3_0.4_1A$genes
fgsealimma_HMEC_8q_new_lineage_RPKM3_0.4_1A <- fgsea(Hallmarks_gene_sets, ranks, minSize=1, maxSize = 500, nperm=10000)

fgsealimma_HMEC_8q_new_lineage_RPKM3_AVG <- fgsealimma_HMEC_8q_new_lineage_RPKM3_102B
fgsealimma_HMEC_8q_new_lineage_RPKM3_AVG$pval <- (fgsealimma_HMEC_8q_new_lineage_RPKM3_102B$pval + fgsealimma_HMEC_8q_new_lineage_RPKM3_0.4_1A$pval + fgsealimma_HMEC_8q_new_lineage_RPKM3_C05$pval + fgsealimma_HMEC_8q_new_lineage_RPKM3_K16$pval + fgsealimma_HMEC_8q_new_lineage_RPKM3_HMNRC4$pval)/5
fgsealimma_HMEC_8q_new_lineage_RPKM3_AVG$padj <- (fgsealimma_HMEC_8q_new_lineage_RPKM3_102B$padj + fgsealimma_HMEC_8q_new_lineage_RPKM3_0.4_1A$padj + fgsealimma_HMEC_8q_new_lineage_RPKM3_C05$padj + fgsealimma_HMEC_8q_new_lineage_RPKM3_K16$padj + fgsealimma_HMEC_8q_new_lineage_RPKM3_HMNRC4$padj)/5
fgsealimma_HMEC_8q_new_lineage_RPKM3_AVG$ES <- (fgsealimma_HMEC_8q_new_lineage_RPKM3_102B$ES + fgsealimma_HMEC_8q_new_lineage_RPKM3_0.4_1A$ES + fgsealimma_HMEC_8q_new_lineage_RPKM3_C05$ES + fgsealimma_HMEC_8q_new_lineage_RPKM3_K16$ES + fgsealimma_HMEC_8q_new_lineage_RPKM3_HMNRC4$ES)/5
fgsealimma_HMEC_8q_new_lineage_RPKM3_AVG$NES <- (fgsealimma_HMEC_8q_new_lineage_RPKM3_102B$NES + fgsealimma_HMEC_8q_new_lineage_RPKM3_0.4_1A$NES + fgsealimma_HMEC_8q_new_lineage_RPKM3_C05$NES + fgsealimma_HMEC_8q_new_lineage_RPKM3_K16$NES + fgsealimma_HMEC_8q_new_lineage_RPKM3_HMNRC4$NES)/5


fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW <- merge(fgsealimma_HMEC_1q_new_RPKM2, fgsealimma_HMEC_1q_new_lineage_RPKM_AVG, by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW, fgsealimma_HMEC_8q_new_RPKM, by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW, fgsealimma_HMEC_8q_new_lineage_RPKM_AVG, by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW, fgsea_limma_BRCANEW0.5_1q_noCNorm_nopur_Hallmarks[,c(1:7)], by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW, fgsea_limma_BRCANEW0.5_8q_noCNorm_nopur_Hallmarks[,c(1:7)], by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW <- data.frame(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW)

rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW) <- fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW$pathway
rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW) <- tolower(rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW))
rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW) <- gsub("hallmark_","",as.character(rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW)))
rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW) <- gsub("_"," ",as.character(rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW)))


fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST <- fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW
colnames(fgsealimma_HMEC_8q_new_lineage_RPKM3_AVG)[4] <- "ES_HMEC_8q_lineage2"
colnames(fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG)[4] <- "ES_HMEC_1q_lineage2"
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST, fgsealimma_HMEC_8q_new_lineage_RPKM3_AVG[,c(1,4)], by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST, fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG[,c(1,4)], by = "pathway")

colnames(fgsealimma_HMEC_8q_new_RPKM3)[4] <- "ES_HMEC_8q_2"
colnames(fgsealimma_HMEC_1q_new_RPKM3)[4] <- "ES_HMEC_1q_2"
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST, fgsealimma_HMEC_8q_new_RPKM3[,c(1,4)], by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST, fgsealimma_HMEC_1q_new_RPKM3[,c(1,4)], by = "pathway")
colnames(fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG)[4] <- "ES_HMEC_1q_lineage3"
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST <- merge(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST, fgsealimma_HMEC_1q_new_lineage_RPKM3_AVG[,c(1,4)], by = "pathway")
fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST <- fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST[order(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST$pathway),]
rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST) <- rownames(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEW)

my_palette3 <- myPalette(low = "darkblue", high = "red3", mid="white", k =10)
my_breaks <- c(-1, -0.7, -0.6, -0.5, -0.3, 0, 0.3, 0.5, 0.6, 0.7, 1 )


fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST_forsubmission <- fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST[,c(28,34,54, 53, 52, 50)]

pheatmap(fgsea_limma_HMEC_RPKM_vs_BRCA_Hallmarks_NEWEST_forsubmission, cluster_rows=TRUE, cluster_cols=TRUE, color = my_palette3, show_rownames=TRUE,  show_colnames=TRUE, border_col=NA, fontsize=6, breaks=my_breaks)


#FIGURE 5B
#drawing in illustrator


#FIGURE 5C

RNA_seq_060921counts2 <- cbind(RNA_seq_060921counts,RNA_seq_060921_2counts[,-c(1:6)])
rownames(RNA_seq_060921counts2) <- RNA_seq_060921counts2$Geneid
RNA_seq_060921counts2_geninfo <- RNA_seq_060921counts2[,c(1:6)]
RNA_seq_060921counts3 <- RNA_seq_060921counts2[,-c(1:6)]
for(i in 7:ncol(RNA_seq_060921counts2)){
RNA_seq_060921counts2[,i] <- as.numeric(as.character(RNA_seq_060921counts2[,i]))
}
RNA_seq_060921counts2$Length <- as.numeric(as.character(RNA_seq_060921counts2$Length))
RNA_seq_060921counts2$Length <- RNA_seq_060921counts2$Length/1000
RNA_seq_060921counts2 <- subset(RNA_seq_060921counts2, Geneid %in% BRCA_RSEM_genes$Gene)
nrow(RNA_seq_060921counts2)
sample_counts_tot <- colSums(RNA_seq_060921counts2[,c(7:ncol(RNA_seq_060921counts2))])
sample_counts_tot <- sample_counts_tot/1000000

RNA_seq_060921counts2[,c(7:ncol(RNA_seq_060921counts2))] <- RNA_seq_060921counts2[,c(7:ncol(RNA_seq_060921counts2))]/sample_counts_tot

for(i in 7:ncol(RNA_seq_060921counts2)) {
RNA_seq_060921counts2[,i] <- RNA_seq_060921counts2[,i]/RNA_seq_060921counts2$Length
}
RNA_seq_060921counts2_lcrm <- RNA_seq_060921counts2
RNA_seq_060921counts2_lcrm$AVG <- rowMeans(RNA_seq_060921counts2_lcrm[,c(7:ncol(RNA_seq_060921counts2_lcrm))])
RNA_seq_060921counts2_lcrm <- subset(RNA_seq_060921counts2_lcrm, RNA_seq_060921counts2_lcrm$AVG > 2)


RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2021_sva_codes5, id %like% "K16" & id != "K16_GSI" & id != "K16_DLL_GSI" )
RNAseq_2020_sva2_temp2_dblNORM_temp <- RNA_seq_060921counts2_lcrm[,names(RNA_seq_060921counts2_lcrm) %in% RNAseq_2020_sva_codes2_temp$unq_id]
yEndtoEnd2 <- DGEList(counts=RNAseq_2020_sva2_temp2_dblNORM_temp, genes=rownames(RNAseq_2020_sva2_temp2_dblNORM_temp))
labels <- data.frame("Name"=colnames(RNAseq_2020_sva2_temp2_dblNORM_temp), "TimePoint"=RNAseq_2020_sva_codes2_temp$id)
timepoint<-factor(labels$TimePoint)
design<-model.matrix(~0+timepoint)
xglm2 <-  estimateDisp(yEndtoEnd2, design)
fit <- glmFit(xglm2, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
lrtres <- data.frame(lrt$genes,lrt$table)
lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
limma_K16_vs_K16_DLL_test2 <- lrtres
limma_K16_vs_K16_DLL_test2[limma_K16_vs_K16_DLL_test2 == "Inf"] <- 200
limma_K16_vs_K16_DLL_test2[limma_K16_vs_K16_DLL_test2 == "-Inf"] <- -200
ranks <- limma_K16_vs_K16_DLL_test2$neglogP
names(ranks) <-limma_K16_vs_K16_DLL_test2$genes
fgsealimma_K16_vs_K16_DLL_test2_notch <- fgsea(Notch_curated_sets, ranks, minSize=1, maxSize = 500, nperm=10000)
fgsealimma_K16_vs_K16_DLL_test2_notch

plotEnrichment(Notch_curated_sets[["antiNotch"]], ranks)

plotEnrichment(Notch_curated_sets[["proNotch"]], ranks)


#FIGURE 5D

#CCLE data analysis
testrun_1qAVG_CNA_vs_GEX_cervix_test <- TCGA_slopeNorm_func_tissue3(CCLE_expression_genes_t2, "c1q", "Cervical Cancer")
ranks <- testrun_1qAVG_CNA_vs_GEX_cervix_test$neglogP
names(ranks) <- testrun_1qAVG_CNA_vs_GEX_cervix_test$genes
fgsea_testrun_1qAVG_CNA_vs_GEX_cervix_test_Notch <- fgsea(pathwaysHx, ranks, minSize=1, maxSize = 10000, nperm=10000)

#TCGA data for 1q NotchUP and NotchDN analysis
BLCA_CNA_summaryNEW_0.5_nopur <- ChrArm_CNA_freq_TCGA_nopuritycorrection(BLCA_CNA, "BLCA_CNAs_nopur_0.5_-0.41_0.32", 0.5, -0.41, 0.32, chromdata_peri_cent2)
my.list1qfwa_BLCANEW0.5_noCNorm_nopur <- Arm_Level_CNorm_diffexp2_nopurity(BLCA_RSEM, BLCA_CNA, BLCA_CNA_summaryNEW_0.5_nopur, 10, 1, "1q", "gain", 5, 30)
ranks <- my.list1qfwa_BLCANEW0.5_noCNorm_nopur[[1]]$neglogP
names(ranks) <- my.list1qfwa_BLCANEW0.5_noCNorm_nopur[[1]]$genes
fgsea_limma_BLCANEW0.5_noCNorm_nopur_Notch <- fgsea(pathwaysHx, ranks, minSize=1, maxSize = 10000, nperm=10000)

my_breaks2 <- c(-2.5,  -1.3, -0.5, 0, 0.5, 1.3, 2.5 )

ggplot(subset(fgsea_summary_TCGA_nopur_Notch_UP2_DN2_forsubmission, class == "TCGA" & gene_set == "Notch_UP"), aes(x=tissue, y=1, width = 1)) + geom_tile(aes(fill = neglogP)) + scale_fill_gradientn(colours = my_palette3, breaks = my_breaks2, limits=c(-2.5,2.5))
ggplot(subset(fgsea_summary_TCGA_nopur_Notch_UP2_DN2_forsubmission, class == "CCLE" & gene_set == "Notch_UP"), aes(x=tissue, y=1, width = 1)) + geom_tile(aes(fill = neglogP)) + scale_fill_gradientn(colours = my_palette3, breaks = my_breaks2, limits=c(-2.5,2.5))
ggplot(subset(fgsea_summary_TCGA_nopur_Notch_UP2_DN2_forsubmission, class == "TCGA" & gene_set == "Notch_DN"), aes(x=tissue, y=1, width = 1)) + geom_tile(aes(fill = neglogP)) + scale_fill_gradientn(colours = my_palette3, breaks = my_breaks2, limits=c(-2.5,2.5))
ggplot(subset(fgsea_summary_TCGA_nopur_Notch_UP2_DN2_forsubmission, class == "CCLE" & gene_set == "Notch_DN"), aes(x=tissue, y=1, width = 1)) + geom_tile(aes(fill = neglogP)) + scale_fill_gradientn(colours = my_palette3, breaks = my_breaks2, limits=c(-2.5,2.5))

#FIGURE 5E
# EGTA notch cleavage experiment; western blot



#FIGURE 5F
# Quantification of western blots


ggplot(Notch_Cleavage_Exp_Summary, aes(x = status_1q, y=N1ICD_normalized, color=status_1q)) +     geom_boxplot() + geom_hline(yintercept = 0) +  theme_bw() + scale_color_manual(values = c("cornflowerblue", "red")) + facet_grid(. ~ Exp) + geom_point() + stat_compare_means(comparisons = list(c("WT", "gain 1q"))) +geom_point(alpha=0.4)



#FIGURE 5G
# Drawing in illustrator



#FIGURE 5H
#mRNA::DNA correlations, matched tumor/normal data


ggplot(mRNA_DNA_correlations_chr1, aes(x=Rank, y = neglogP)) +geom_point(aes(color = neglogP)) + geom_hline(yintercept = c(0,10.5), size = 0.25) + scale_colour_gradientn(colours=c("white","pink", "red3")) +geom_point(data = subset(mRNA_DNA_correlations_chr1, genes %in% c("APH1A", "NCSTN", "AKT3", "MDM4", "MCL1", "KDM5B", "PSEN2"), color = "black"), pch = 1) + geom_text_repel(data=subset(mRNA_DNA_correlations_chr1, genes %in% c("APH1A", "NCSTN", "AKT3", "MDM4", "MCL1", "KDM5B", "PSEN2")), aes(label =  genes), color = "red", size=3.2, max.overlaps=100) + theme_bw()


#FIGURE 5I
#Notch_Activation_Assay_NCSTN_KO_110622.csv

Notch_Activation_Assay_NCSTN_KO_110622n <- Notch_Activation_Assay_NCSTN_KO_110622
Notch_Activation_Assay_NCSTN_KO_110622n$NCSTN_vs_GAPDH <- (Notch_Activation_Assay_NCSTN_KO_110622n$NCSTN_vs_GAPDH/1.850608)
Notch_Activation_Assay_NCSTN_KO_110622n$N1ICD_vs_GAPDH <- (Notch_Activation_Assay_NCSTN_KO_110622n$N1ICD_vs_GAPDH/0.8222574)
my_comparisons <- list( c("dip_ctrl", "K26_ctrl"), c("dip_ctrl", "dip_sgNCSTN_1"), c("dip_ctrl", "dip_sgNCSTN_2"),  c("dip_ctrl", "K26_sgNCSTN_1"),  c("K26_ctrl", "K26_sgNCSTN_1"), c("K26_ctrl", "K26_sgNCSTN_2"), c("dip_sgNCSTN_1", "K26_sgNCSTN_2"))
dev.off()
pdf(file = "Notch_Activation_Assay_NCSTN_KO_110622_NCSTN.pdf", width = 4, height = 6)
ggboxplot(Notch_Activation_Assay_NCSTN_KO_110622n, x = "guide", y = "NCSTN_vs_GAPDH", color = "guide", palette = c("grey60", "yellow2", "orange2", "grey60", "yellow2", "orange2")) + geom_point() + stat_compare_means(comparisons = my_comparisons, method = "t.test") + stat_compare_means(label.y = 4) + geom_hline(yintercept = c(1, 0.5)) + ylim(0,5)
dev.off()



#FIGURE 5J

dev.off()
pdf(file = "Notch_Activation_Assay_NCSTN_KO_110622_N1ICD2.pdf", width = 4, height = 6)
ggboxplot(Notch_Activation_Assay_NCSTN_KO_110622n, x = "guide", y = "N1ICD_vs_GAPDH", color = "guide", palette = c("grey60", "yellow2", "orange2", "grey60", "yellow2", "orange2")) + geom_point() + stat_compare_means(comparisons = my_comparisons, method = "t.test") + stat_compare_means(label.y = 4) + geom_hline(yintercept = c(1, 0.5)) + ylim(0,5)
dev.off()




#FIGURE 5K


dev.off()
pdf(file = "Notch_Activation_Assay_NCSTN_KO_110622_N1ICD_vs_NCSTN.pdf", width = 6, height = 6)
ggplot(data=Notch_Activation_Assay_NCSTN_KO_110622n, aes(x=NCSTN_vs_GAPDH, y=N1ICD_vs_GAPDH)) +geom_point(aes(color=guide,shape=guide ), size = 3) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 3)) +theme_light() + geom_text_repel(data= Notch_Activation_Assay_NCSTN_KO_110622n, aes(label =  line), color = "black", size=4, max.overlaps=100)
dev.off()




#functions

TCGA_slopeNorm_func_tissue3 <- function(RNAi_combined_gene_dep_scores2_t, gene_target, tissue) {
    RNAi_combined_gene_dep_scores2_t <- subset(RNAi_combined_gene_dep_scores2_t, primary_disease == tissue)
    limma_gene_DNA_RNA_fidelity_summaryTCGA3 <- matrix(nrow=(ncol(RNAi_combined_gene_dep_scores2_t)), ncol = 3)
    for (j in 54:(ncol(RNAi_combined_gene_dep_scores2_t)) ) {
        testdf <- matrix(nrow = nrow(RNAi_combined_gene_dep_scores2_t), ncol = 2)
        testdf[,1] <- RNAi_combined_gene_dep_scores2_t[, gene_target]
        testdf[,2] <- RNAi_combined_gene_dep_scores2_t[,(j)]
        testdf <- data.frame(testdf)
        M.lm=lm(X1~X2,data=testdf)
        if (nrow(summary(M.lm)$coefficients) < 2)
            next
        limma_gene_DNA_RNA_fidelity_summaryTCGA3[j,1] <- as.character(colnames(RNAi_combined_gene_dep_scores2_t[j]))
        limma_gene_DNA_RNA_fidelity_summaryTCGA3[j,2] <- as.numeric(as.character(summary(M.lm)$coefficients[2,1]))
        limma_gene_DNA_RNA_fidelity_summaryTCGA3[j,3] <- as.numeric(as.character(summary(M.lm)$coefficients[2,4]))
    }
    limma_gene_DNA_RNA_fidelity_summaryTCGA3 <- data.frame(limma_gene_DNA_RNA_fidelity_summaryTCGA3)
    colnames( limma_gene_DNA_RNA_fidelity_summaryTCGA3) <- c("genes", "slope", "pval")
    limma_gene_DNA_RNA_fidelity_summaryTCGA3$slope <- as.numeric(as.character(limma_gene_DNA_RNA_fidelity_summaryTCGA3$slope))
    limma_gene_DNA_RNA_fidelity_summaryTCGA3$pval <- as.numeric(as.character(limma_gene_DNA_RNA_fidelity_summaryTCGA3$pval))
    limma_gene_DNA_RNA_fidelity_summaryTCGA3 <- limma_gene_DNA_RNA_fidelity_summaryTCGA3[-c(1:53),]
    limma_gene_DNA_RNA_fidelity_summaryTCGA3$neglogP <- (limma_gene_DNA_RNA_fidelity_summaryTCGA3$slope/abs(limma_gene_DNA_RNA_fidelity_summaryTCGA3$slope))*-log10(limma_gene_DNA_RNA_fidelity_summaryTCGA3$pval)
    limma_gene_DNA_RNA_fidelity_summaryTCGA3$genes <- gsub("[.]y", "", limma_gene_DNA_RNA_fidelity_summaryTCGA3$genes)
    return(limma_gene_DNA_RNA_fidelity_summaryTCGA3)
}


ChrArm_CNA_freq_TCGA_nopuritycorrection <- function(TCGA_CNA, y, z, w, g, chromdata) { 
    #prep files
    print("... preparing data...")
    
    TCGA_CNA <- separate(data = TCGA_CNA, col = Sample, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_CNA$id <- paste(TCGA_CNA$TCGA, TCGA_CNA$TSS, TCGA_CNA$Participant, TCGA_CNA$SampleVial, sep='-')
    TCGA_CNA <- subset(TCGA_CNA, TCGA_CNA$SampleVial == "01A" | TCGA_CNA$SampleVial == "01B" )
    TCGA_CNA <- TCGA_CNA[,-c(1:7)]
    print("success")
    TCGA_CNA_purity <- TCGA_CNA
    TCGA_CNA_purity$Segment_Mean[is.na(TCGA_CNA_purity$Segment_Mean)] <- 0
    #TCGA_CNA_purity <- TCGA_CNA_purity[,-7]
    print("success")
    #assign segments to arms, deal with centromere-spanning segments
    print("... assigning segments to arms ...")
    TCGA_CNA_purity <- merge(TCGA_CNA_purity, chromdata, by = "Chromosome")
    TCGA_CNA_purity$length <- TCGA_CNA_purity$End - TCGA_CNA_purity$Start
    TCGA_CNA_purity$length <- as.numeric(as.character(TCGA_CNA_purity$length))
    TCGA_CNA_purity$Mid_Segment <- TCGA_CNA_purity$End - (TCGA_CNA_purity$End - TCGA_CNA_purity$Start)/2
    TCGA_CNA_purity$Mid_Segment <- as.numeric(as.character(TCGA_CNA_purity$Mid_Segment))
    TCGA_CNA_purity$Segment_Mean <- as.numeric(as.character(TCGA_CNA_purity$Segment_Mean))
    TCGA_CNA_purity$centromere_start <- as.numeric(as.character(TCGA_CNA_purity$centromere_start))
    TCGA_CNA_purity$Plength <- as.numeric(as.character(TCGA_CNA_purity$Plength))
    TCGA_CNA_purity$Qlength <- as.numeric(as.character(TCGA_CNA_purity$Qlength))
    TCGA_CNA_purity$WC <- ifelse(TCGA_CNA_purity$Start < TCGA_CNA_purity$centromere_start & TCGA_CNA_purity$End > TCGA_CNA_purity$centromere_end, 1, 0)
    TCGA_CNA_purity <- rbind(TCGA_CNA_purity, TCGA_CNA_purity %>% filter(WC == 1) %>% mutate(WC = 2))
    TCGA_CNA_purity$ChromArm <- ifelse(TCGA_CNA_purity$WC == 2, paste(TCGA_CNA_purity$Chromosome, "q", sep ="", collapse=NULL), ifelse(TCGA_CNA_purity$WC == 1, paste(TCGA_CNA_purity$Chromosome, "p", sep ="", collapse=NULL), ifelse(TCGA_CNA_purity$Mid_Segment < TCGA_CNA_purity$centromere_start, paste(TCGA_CNA_purity$Chromosome, "p", sep ="", collapse=NULL), paste(TCGA_CNA_purity$Chromosome, "q", sep ="", collapse=NULL))))
    TCGA_CNA_purity$ChromArm_length <- ifelse(TCGA_CNA_purity$ChromArm %like% "p", TCGA_CNA_purity$Plength, TCGA_CNA_purity$Qlength)
    TCGA_CNA_purity$ChromArm_length <- as.numeric(as.character(TCGA_CNA_purity$ChromArm_length))
    TCGA_CNA_purity$length <- ifelse(TCGA_CNA_purity$WC == 0, TCGA_CNA_purity$length, ifelse(TCGA_CNA_purity$WC == 2, (TCGA_CNA_purity$End - TCGA_CNA_purity$centromere_end), (TCGA_CNA_purity$centromere_start - TCGA_CNA_purity$Start) ))
    TCGA_CNA_purity$length <- as.numeric(as.character(TCGA_CNA_purity$length))
    TCGA_CNA_purity <- TCGA_CNA_purity[!(TCGA_CNA_purity$ChromArm=="13p" | TCGA_CNA_purity$ChromArm=="14p" | TCGA_CNA_purity$ChromArm=="15p" | TCGA_CNA_purity$ChromArm=="21p" | TCGA_CNA_purity$ChromArm=="22p"),]
    ids <- split(TCGA_CNA_purity, TCGA_CNA_purity$id)
    print("success")
    #arm-level loss summary
    print("... preparing chromosome loss summary ...")
    {  
        TCGA_CNA_purity_losses <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean < w)
        TCGA_CNA_purity_losses_chrX <- split(TCGA_CNA_purity_losses, TCGA_CNA_purity_losses$ChromArm)
        TCGA_CNA_purity_losses_chrX <- TCGA_CNA_purity_losses_chrX[sapply(TCGA_CNA_purity_losses_chrX, function(x) dim(x)[1]) > 0]
        TCGA_deletion_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_deletion_frequencies_ALL <- data.frame(TCGA_deletion_frequencies_ALL)
        colnames(TCGA_deletion_frequencies_ALL) <- c('Chromosome', 'del_frequencies', 'avg_segMean_deleted')
        Deletion_summary <- matrix(ncol=4, nrow=1)
        colnames(Deletion_summary) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'ChrArm')
        for (j in 1:length(TCGA_CNA_purity_losses_chrX)) {
            TCGA_deletion_frequencies <- matrix(ncol = 3, nrow = 1)
            colnames(TCGA_deletion_frequencies) <- c('Chromosome', 'del_frequencies', 'avg_segMean_deleted')
            TCGA_CNA_purity_losses_chrJ <- TCGA_CNA_purity_losses_chrX[[j]]
            TCGA_CNA_purity_losses_chrJ_sample <- split(TCGA_CNA_purity_losses_chrJ, TCGA_CNA_purity_losses_chrJ$id)
            xout <- matrix(ncol = 4, nrow = length(TCGA_CNA_purity_losses_chrJ_sample))
            for (i in 1:length(TCGA_CNA_purity_losses_chrJ_sample)) {
                TCGA_CNA_purity_losses_chrJ_sample_patX <- TCGA_CNA_purity_losses_chrJ_sample[[i]]
                TCGA_CNA_purity_losses_chrJ_sample_patX$weightedSegMean <- TCGA_CNA_purity_losses_chrJ_sample_patX$length*TCGA_CNA_purity_losses_chrJ_sample_patX$Segment_Mean
                xout[i, 1] <- TCGA_CNA_purity_losses_chrJ_sample_patX$id[1]
                xout[i, 2] <- sum(TCGA_CNA_purity_losses_chrJ_sample_patX$weightedSegMean)/sum(TCGA_CNA_purity_losses_chrJ_sample_patX$length)
                xout[i, 3] <- sum(TCGA_CNA_purity_losses_chrJ_sample_patX$length)/TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm_length[1]
                xout <- data.frame(xout)
                xout$X1 <- as.character(xout$X1)
                xout$X2 <- as.numeric(as.character(xout$X2))
                xout$X3 <- as.numeric(as.character(xout$X3))
            }
            xout_sig <- subset(xout, xout$X2 < w & xout$X3 > z)
            xout_sig$X4 <- as.character(xout_sig$X4)
            if(nrow(xout_sig) == 0) xout_sig[1,] <- c(0,0,0,TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1])
            xout_sig$X4 <- TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1]
            colnames(xout_sig) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'ChrArm')
            write.table(xout_sig, file = paste(y, "Chromosome_", TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1], "_deletion_samples", ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
            Deletion_summary <- rbind(Deletion_summary, xout_sig)
            TCGA_deletion_frequencies[1, 1] <- TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1]
            TCGA_deletion_frequencies[1, 2] <- nrow(subset(xout_sig, id != 0))/length(ids)
            TCGA_deletion_frequencies[1, 3] <- mean(xout_sig$avg_segMean_deleted)
            TCGA_deletion_frequencies <- data.frame(TCGA_deletion_frequencies)
            TCGA_deletion_frequencies_ALL <- rbind(TCGA_deletion_frequencies_ALL, TCGA_deletion_frequencies)
        }
        Deletion_summary <- Deletion_summary[-1,]
        TCGA_deletion_frequencies_ALL <- TCGA_deletion_frequencies_ALL[-1,]
        write.table(TCGA_deletion_frequencies_ALL, file = paste(y, "_deletion_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        write.table(Deletion_summary, file = paste(y, "_Deletion_summary", ".txt",sep=""), sep="\t", row.names=FALSE,  quote = FALSE)
        
    }
    print("success")
    #arm-level gain summary
    print("... preparing chromosome gain summary ...")
    {
        TCGA_CNA_purity_gains <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean > g)
        TCGA_CNA_purity_gains_chrX <- split(TCGA_CNA_purity_gains, TCGA_CNA_purity_gains$ChromArm)
        TCGA_CNA_purity_gains_chrX <- TCGA_CNA_purity_gains_chrX[sapply(TCGA_CNA_purity_gains_chrX, function(x) dim(x)[1]) > 0]
        
        TCGA_amplification_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_amplification_frequencies_ALL <- data.frame(TCGA_amplification_frequencies_ALL)
        colnames(TCGA_amplification_frequencies_ALL) <- c('Chromosome', 'amp_frequencies', 'avg_segMean_amplified')
        Amplification_summary <- matrix(ncol = 4, nrow = 1)
        colnames(Amplification_summary) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'ChrArm')
        for (j in 1:length(TCGA_CNA_purity_gains_chrX)) {
            TCGA_amplification_frequencies <- matrix(ncol = 3, nrow = 1)
            colnames(TCGA_amplification_frequencies) <- c('Chromosome', 'amp_frequencies', 'avg_segMean_amplified')
            TCGA_CNA_purity_gains_chrJ <- TCGA_CNA_purity_gains_chrX[[j]]
            TCGA_CNA_purity_gains_chrJ_sample <- split(TCGA_CNA_purity_gains_chrJ, TCGA_CNA_purity_gains_chrJ$id)
            xout2 <- matrix(ncol = 4, nrow = length(TCGA_CNA_purity_gains_chrJ_sample))
            
            for (i in 1:length(TCGA_CNA_purity_gains_chrJ_sample)) {
                TCGA_CNA_purity_gains_chrJ_sample_patX <- TCGA_CNA_purity_gains_chrJ_sample[[i]]
                TCGA_CNA_purity_gains_chrJ_sample_patX$weightedSegMean <- TCGA_CNA_purity_gains_chrJ_sample_patX$length*TCGA_CNA_purity_gains_chrJ_sample_patX$Segment_Mean
                xout2[i, 1] <- TCGA_CNA_purity_gains_chrJ_sample_patX$id[1]
                xout2[i, 2] <- sum(TCGA_CNA_purity_gains_chrJ_sample_patX$weightedSegMean)/sum(TCGA_CNA_purity_gains_chrJ_sample_patX$length)
                xout2[i, 3] <- sum(TCGA_CNA_purity_gains_chrJ_sample_patX$length)/TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm_length[1]
                xout2 <- data.frame(xout2)
                xout2$X1 <- as.character(xout2$X1)
                xout2$X2 <- as.numeric(as.character(xout2$X2))
                xout2$X3 <- as.numeric(as.character(xout2$X3))
            }
            
            xout_sig2 <- subset(xout2, xout2$X2 > g & xout2$X3 > z)
            xout_sig2$X4 <- as.character(xout_sig2$X4)
            if(nrow(xout_sig2) == 0) xout_sig2[1,] <- c(0,0,0,TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1])
            xout_sig2$X4 <- TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1]
            colnames(xout_sig2) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'ChrArm')
            write.table(xout_sig2, file = paste(y, "Chromosome_", TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1], "_amplification_samples", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
            Amplification_summary <- rbind(Amplification_summary, xout_sig2)
            TCGA_amplification_frequencies[1, 1] <- TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1]
            TCGA_amplification_frequencies[1, 2] <- nrow(subset(xout_sig2, id != 0))/length(ids)
            TCGA_amplification_frequencies[1, 3] <- mean(xout_sig2$avg_segMean_amplified)
            TCGA_amplification_frequencies <- data.frame(TCGA_amplification_frequencies)
            TCGA_amplification_frequencies_ALL <- rbind(TCGA_amplification_frequencies_ALL, TCGA_amplification_frequencies)
        }
        Amplification_summary <- Amplification_summary[-1,]
        TCGA_amplification_frequencies_ALL <- TCGA_amplification_frequencies_ALL[-1,]
        write.table(TCGA_amplification_frequencies_ALL, file = paste(y, "_amplification_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        write.table(Amplification_summary, file = paste(y, "_Amplification_summary", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        
    }
    final_output <- list(Deletion_summary, TCGA_deletion_frequencies_ALL, Amplification_summary, TCGA_amplification_frequencies_ALL)
    final_output[[1]] <- final_output[[1]][final_output[[1]]$id != 0, ]
    final_output[[2]][is.na(final_output[[2]])] <- 0
    final_output[[3]] <- final_output[[3]][final_output[[3]]$id != 0, ]
    final_output[[4]][is.na(final_output[[4]])] <- 0
    return(final_output)
    print("success")
}


Arm_Level_CNorm_diffexp2_nopurity <- function(TCGA_RSEM,  TCGA_CNA, ChrArm_CNA_summary, RSEM_min_count, cohort_subsample_minFrac, arm, gain_loss, acceptable_difference_percent, iter) {
    #define function
    CNA_distributions_plot_func <- function(x) {
        CNA_distributions_ChrArm2 <- x[,c(1,4)]
        CNA_distributions_ChrArm3 <- x[,c(1,5)]
        CNA_distributions_ChrArm2$status <- "yes"
        CNA_distributions_ChrArm3$status <- "no"
        colnames(CNA_distributions_ChrArm3)[2] <- "yes"
        CNA_distributions_ChrArm4 <- rbind(CNA_distributions_ChrArm2,CNA_distributions_ChrArm3)
        return(CNA_distributions_ChrArm4) }
    #preparing files
    print("... preparing RSEM file ... ")
    TCGA_RSEM <- TCGA_RSEM[-1,]
    TCGA_RSEM <- separate(data = TCGA_RSEM, col = Hybridization.REF, into = c('Gene', 'x'), sep = "\\|")
    TCGA_RSEM$Gene <- ifelse(TCGA_RSEM$Gene =="?", TCGA_RSEM$x, TCGA_RSEM$Gene)
    TCGA_RSEM <- TCGA_RSEM[,-2]
    TCGA_RSEM <- data.frame(TCGA_RSEM)
    rownames(TCGA_RSEM) <- make.names(TCGA_RSEM$Gene, unique=TRUE)
    TCGA_RSEM <- TCGA_RSEM[,-1]

    for(i in 1:ncol(TCGA_RSEM)){
        TCGA_RSEM[,i] <- as.numeric(TCGA_RSEM[,i])
    }
    TCGA_RSEM_ids <- colnames(TCGA_RSEM)
    TCGA_RSEM_ids <- data.frame(TCGA_RSEM_ids)
    colnames(TCGA_RSEM_ids)[1] <- c("id")
    TCGA_RSEM_ids <- separate(data = TCGA_RSEM_ids, col = id, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\.")
    TCGA_RSEM_ids$id <- paste(TCGA_RSEM_ids$TCGA, TCGA_RSEM_ids$TSS, TCGA_RSEM_ids$Participant, TCGA_RSEM_ids$SampleVial, sep='-')
    TCGA_RSEM_ids$SampleVial <- gsub('.{1}$', '', TCGA_RSEM_ids$SampleVial)
    TCGA_RSEM_ids$id <- gsub('.{1}$', '', TCGA_RSEM_ids$id)
    TCGA_RSEM_ids$Sample <- colnames(TCGA_RSEM)
    TCGA_RSEM_ids$SampleVial <- as.numeric(as.character(TCGA_RSEM_ids$SampleVial))
    TCGA_RSEM_ids <- subset(TCGA_RSEM_ids, SampleVial < 6)
    TCGA_RSEM2 <- TCGA_RSEM[, names(TCGA_RSEM) %in% TCGA_RSEM_ids$Sample]
    print("success")
    #remove samples with no CNA or purity data
    print("... preparing merge of CNA, RSEM, purity data ... ")
    TCGA_CNA <- TCGA_CNA[,c(1,2)]
    TCGA_CNA <- separate(data = TCGA_CNA, col = Sample, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_CNA$id <- paste(TCGA_CNA$TCGA, TCGA_CNA$TSS, TCGA_CNA$Participant, TCGA_CNA$SampleVial, sep='-')
    TCGA_CNA$id <- gsub('.{1}$', '', TCGA_CNA$id)
    TCGA_CNA$SampleVial <- gsub('.{1}$', '', TCGA_CNA$SampleVial)
    TCGA_CNA <- subset(TCGA_CNA, SampleVial < 6)
    TCGA_CNA_samples <- unique(TCGA_CNA$id)
    TCGA_CNA_samples <- data.frame(TCGA_CNA_samples)
    colnames(TCGA_CNA_samples) <- c("id")
    TCGA_CNA_samples$CNAcheck <- 1
    TCGA_RSEM_ids <- merge(TCGA_RSEM_ids, TCGA_CNA_samples, by="id")
    TCGA_RSEM_ids <- distinct(TCGA_RSEM_ids,Sample, .keep_all= TRUE)
    TCGA_RSEM2 <- TCGA_RSEM2[, names(TCGA_RSEM2) %in% TCGA_RSEM_ids$Sample]
    print("success")
    w <- ncol(TCGA_RSEM2)/2
    #remove gene rows with low counts
    TCGA_RSEM3 <- TCGA_RSEM2[rowSums(TCGA_RSEM2 < RSEM_min_count) <=  w, , drop = FALSE]
    #grouping based on CNA status
    print("... grouping based on CNA status ...")
    `%ni%` <- Negate(`%in%`)
    ChrArm_summary_gain <- ChrArm_CNA_summary[[3]]
    ChrArm_summary_loss <- ChrArm_CNA_summary[[1]]
    ChrArm_summary_gain$id <- gsub('.{1}$', '', ChrArm_summary_gain$id)
    ChrArm_summary_loss$id <- gsub('.{1}$', '', ChrArm_summary_loss$id)
    ChrArm_summary_gain$status <- "gain"
    ChrArm_summary_loss$status <- "loss"
    ChrArm_summary_gain <- ChrArm_summary_gain[,c(1,4:5)]
    ChrArm_summary_loss <- ChrArm_summary_loss[,c(1,4:5)]
    ChrArm_summary <- rbind(ChrArm_summary_gain, ChrArm_summary_loss)
    ChrArm_summary$ChrArm <- as.character(ChrArm_summary$ChrArm)
    ChrArm_summary$ChrArm2 <- paste(ChrArm_summary$ChrArm, ChrArm_summary$status, sep="_")
    ChrArm_summary <- subset(ChrArm_summary, id %in% TCGA_RSEM_ids$id)
    ChrArmY <- paste(arm, gain_loss, sep="_")
    ChrArm_summary2 <- ChrArm_summary
    ids_to_grab1 <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)$id
    ChrArm_summary_subset1 <- subset(ChrArm_summary2, id %in% ids_to_grab1)
    ChrArm_summary_subset_ELSE1 <- subset(ChrArm_summary2, id %ni% ids_to_grab1)
    min3subset <- round(length(unique(ChrArm_summary_subset1$id))*cohort_subsample_minFrac)
    min3ELSE <- round(length(unique(ChrArm_summary_subset_ELSE1$id))*cohort_subsample_minFrac)
    print("success")
    ChrArm_list <- unique(ChrArm_summary$ChrArm2)
    #my.list2 <- list()
    print("... performing CNorm to generate custom cohort ...")
    for(i in 1:iter) {
        ids_to_grab <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)$id
        ChrArm_summary_subset <- subset(ChrArm_summary2, id %in% ids_to_grab)
        ChrArm_summary_subset_ELSE <- subset(ChrArm_summary2, id %ni% ids_to_grab)
        CNA_distributions_ChrArm <- matrix(ncol=5, nrow=length(unique(ChrArm_summary$ChrArm2)))
        colnames(CNA_distributions_ChrArm) <- c("ChrArm_status","yes_count", "no_count", "yes", "no")
        for(i in 1:length(ChrArm_list)){
            ChrArm_X <- ChrArm_list[i]
            CNA_distributions_ChrArm[i,1] <- ChrArm_X
            if (ChrArm_X == ChrArmY)
                next
            test_ids <- subset(ChrArm_summary_subset_ELSE, ChrArm2 == ChrArm_X)
            test_ids2 <- subset(ChrArm_summary_subset, ChrArm2 == ChrArm_X)
            CNA_distributions_ChrArm[i,2] <- nrow(test_ids2)
            CNA_distributions_ChrArm[i,3] <- nrow(test_ids)
            CNA_distributions_ChrArm[i,4] <- 100*nrow(test_ids2)/length(unique(ChrArm_summary_subset$id))
            CNA_distributions_ChrArm[i,5] <- 100*nrow(test_ids)/length(unique(ChrArm_summary_subset_ELSE$id))
        }
        CNA_distributions_ChrArm <- data.frame(CNA_distributions_ChrArm)
        CNA_distributions_ChrArm$yes <- as.numeric(as.character(CNA_distributions_ChrArm$yes))
        CNA_distributions_ChrArm$no <- as.numeric(as.character(CNA_distributions_ChrArm$no))
        CNA_distributions_ChrArm$yes_count <- as.numeric(as.character(CNA_distributions_ChrArm$yes_count))
        CNA_distributions_ChrArm$no_count <- as.numeric(as.character(CNA_distributions_ChrArm$no_count))
        CNA_distributions_ChrArm$diff <- CNA_distributions_ChrArm$yes - CNA_distributions_ChrArm$no
        CNA_distributions_ChrArm$yes_to_remove <- CNA_distributions_ChrArm$yes_count - (CNA_distributions_ChrArm$no/100)*(length(unique(ChrArm_summary_subset$id)))
        CNA_distributions_ChrArm$no_to_remove <- CNA_distributions_ChrArm$no_count - (CNA_distributions_ChrArm$yes/100)*(length(unique(ChrArm_summary_subset_ELSE$id)))
        CNA_distributions_ChrArm$yes_to_remove <- round(CNA_distributions_ChrArm$yes_to_remove)
        CNA_distributions_ChrArm$no_to_remove <- round(CNA_distributions_ChrArm$no_to_remove)
        CNA_distributions_ChrArm2 <- subset(CNA_distributions_ChrArm, abs(diff) > acceptable_difference_percent)
        CNA_distributions_ChrArm2 <- subset(CNA_distributions_ChrArm2, yes_count != 0 & no_count !=0)
        if (cohort_subsample_minFrac == 1)
            break
        if (nrow(CNA_distributions_ChrArm2) == 0)
            break
        if (length(unique(ChrArm_summary_subset$id)) < min3subset)
            break
        if (length(unique(ChrArm_summary_subset_ELSE$id)) < min3ELSE)
            break
        CNA_distributions_ChrArm2$same_or_else <- ifelse(CNA_distributions_ChrArm2$no_to_remove > 0, "same", "else")
        CNA_distributions_ChrArm2$ChrArm_status <- as.character(CNA_distributions_ChrArm2$ChrArm_status)
        CNA_distributions_ChrArm2 <- CNA_distributions_ChrArm2[sample(nrow(CNA_distributions_ChrArm2)),]
        ChrArm_summary_subset_ELSE <- ChrArm_summary_subset_ELSE[sample(nrow(ChrArm_summary_subset_ELSE)),]
        ChrArm_summary_subset <- ChrArm_summary_subset[sample(nrow(ChrArm_summary_subset)),]
        ELSE_ids_removable <- sort(table(subset(ChrArm_summary_subset_ELSE, ChrArm2 %ni% subset(CNA_distributions_ChrArm2, same_or_else == "else")$ChrArm_status )$id), decreasing = TRUE)
        ELSE_ids_removable <- data.frame(ELSE_ids_removable)
        ELSE_ids_removable$Var1 <- as.character(ELSE_ids_removable$Var1)
        ELSE_ids_removable <- ELSE_ids_removable$Var1
        ELSE_ids_removable2 <- sort(table(subset(ChrArm_summary_subset_ELSE, ChrArm2 %in% subset(CNA_distributions_ChrArm2, same_or_else == "same")$ChrArm_status )$id), decreasing = TRUE)
        ELSE_ids_removable2 <- data.frame(ELSE_ids_removable2)
        ELSE_ids_removable2$Var1 <- as.character(ELSE_ids_removable2$Var1)
        ELSE_ids_removable2 <- ELSE_ids_removable2$Var1
        subset_ids_removable <- sort(table(subset(ChrArm_summary_subset, ChrArm2 %in% subset(CNA_distributions_ChrArm2, yes_to_remove > 0)$ChrArm_status )$id), decreasing = TRUE)
        subset_ids_removable <- data.frame(subset_ids_removable)
        subset_ids_removable$Var1 <- as.character(subset_ids_removable$Var1)
        subset_ids_removable <- subset_ids_removable$Var1
        if (CNA_distributions_ChrArm2$same_or_else[1] == "same") {
            ids_to_remove <- ELSE_ids_removable2[1:round(abs(CNA_distributions_ChrArm2$no_to_remove[1])/3)]
            ChrArm_summary2 <- subset(ChrArm_summary2, id %ni%  ids_to_remove)
        } else {
            ids_to_remove <- subset_ids_removable[1:round(abs(CNA_distributions_ChrArm2$yes_to_remove[1])/3)]
            ChrArm_summary2 <- subset(ChrArm_summary2, id %ni%  ids_to_remove)
        }
    }
    print("success")
    print("... performing edgeR analysis ...")
    CNA_distributions_ChrArm5 <-  CNA_distributions_plot_func(CNA_distributions_ChrArm)
    CNA_distributions_ChrArmX <- CNA_distributions_ChrArm
    after_plot <- ggpar(ggbarplot(CNA_distributions_ChrArm5, x="ChrArm_status", y ="yes", fill="status", size = 0, palette = c("#00AFBB", "#FC4E07"), position = position_dodge(0.9)), font.tickslab = c(6), xtickslab.rt = 45)
    ChrArm_summary_subset <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)
    CNA_ids_list <- unique(ChrArm_summary2$id)
    TCGA_RSEM_ids$CNA_status <- ifelse(TCGA_RSEM_ids$id %in% ChrArm_summary_subset$id, 2, 1)
    TCGA_RSEM_ids$remove_status <- ifelse(TCGA_RSEM_ids$id %in% ChrArm_summary2$id, 1, 0)
    to.remove <- subset(TCGA_RSEM_ids, remove_status==0)$Sample
    TCGA_RSEM3 <- subset(TCGA_RSEM3,select = names(TCGA_RSEM3) %ni% to.remove)
    TCGA_RSEM_ids <- subset(TCGA_RSEM_ids, remove_status ==1 )
    yEndtoEnd2 <- DGEList(counts=TCGA_RSEM3, genes=rownames(TCGA_RSEM3))
    labels <- data.frame("Name"=colnames(yEndtoEnd2), "TimePoint"=TCGA_RSEM_ids$CNA_status)
    timepoint<-factor(labels$TimePoint)
    design<-model.matrix(~0+timepoint)
    xglm2 <-  estimateDisp(yEndtoEnd2, design)
    fit <- glmFit(xglm2, design)
    lrt <- glmLRT(fit, contrast=c(-1,1))
    lrtres <- data.frame(lrt$genes,lrt$table)
    lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
    print("success")
    print("... preparing plot ...")
    ChrArm_summary2 <- ChrArm_summary
    ChrArm_list <- unique(ChrArm_summary$ChrArm2)
    ids_to_grab <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)$id
    ChrArm_summary_subset <- subset(ChrArm_summary2, id %in% ids_to_grab)
    ChrArm_summary_subset_ELSE <- subset(ChrArm_summary2, id %ni% ids_to_grab)
    CNA_distributions_ChrArm <- matrix(ncol=5, nrow=length(unique(ChrArm_summary$ChrArm2)))
    colnames(CNA_distributions_ChrArm) <- c("ChrArm_status","yes_count", "no_count", "yes", "no")
    for(i in 1:length(ChrArm_list)){
        ChrArm_X <- ChrArm_list[i]
        CNA_distributions_ChrArm[i,1] <- ChrArm_X
        if (ChrArm_X == ChrArmY)
            next
        test_ids <- subset(ChrArm_summary_subset_ELSE, ChrArm2 == ChrArm_X)
        test_ids2 <- subset(ChrArm_summary_subset, ChrArm2 == ChrArm_X)
        CNA_distributions_ChrArm[i,2] <- nrow(test_ids2)
        CNA_distributions_ChrArm[i,3] <- nrow(test_ids)
        CNA_distributions_ChrArm[i,4] <- 100*nrow(test_ids2)/length(unique(ChrArm_summary_subset$id))
        CNA_distributions_ChrArm[i,5] <- 100*nrow(test_ids)/length(unique(ChrArm_summary_subset_ELSE$id))
    }
    CNA_distributions_ChrArm <- data.frame(CNA_distributions_ChrArm)
    CNA_distributions_ChrArm$yes <- as.numeric(as.character(CNA_distributions_ChrArm$yes))
    CNA_distributions_ChrArm$no <- as.numeric(as.character(CNA_distributions_ChrArm$no))
    CNA_distributions_ChrArm$yes_count <- as.numeric(as.character(CNA_distributions_ChrArm$yes_count))
    CNA_distributions_ChrArm$no_count <- as.numeric(as.character(CNA_distributions_ChrArm$no_count))
    CNA_distributions_ChrArm$diff <- CNA_distributions_ChrArm$yes - CNA_distributions_ChrArm$no
    CNA_distributions_ChrArm$yes_to_remove <- CNA_distributions_ChrArm$yes_count - (CNA_distributions_ChrArm$no/100)*(length(unique(ChrArm_summary_subset$id)))
    CNA_distributions_ChrArm$no_to_remove <- CNA_distributions_ChrArm$no_count - (CNA_distributions_ChrArm$yes/100)*(length(unique(ChrArm_summary_subset_ELSE$id)))
    CNA_distributions_ChrArm$yes_to_remove <- round(CNA_distributions_ChrArm$yes_to_remove)
    CNA_distributions_ChrArm$no_to_remove <- round(CNA_distributions_ChrArm$no_to_remove)
    CNA_distributions_ChrArm4 <-  CNA_distributions_plot_func(CNA_distributions_ChrArm)
    before_plot <- ggpar(ggbarplot(CNA_distributions_ChrArm4, x="ChrArm_status", y ="yes", fill="status", size = 0, palette = c("#00AFBB", "#FC4E07"), position = position_dodge(0.9)), font.tickslab = c(6), xtickslab.rt = 45)
    theplot <- ggarrange(before_plot, after_plot,
                         labels = c("before normalization", "after normalization"),
                         ncol = 1, nrow = 2)
    print(theplot)
    my.list <- list(lrtres, TCGA_RSEM_ids, CNA_distributions_ChrArmX, CNA_distributions_ChrArm5, CNA_distributions_ChrArm, CNA_distributions_ChrArm4)
    return(my.list)
}



