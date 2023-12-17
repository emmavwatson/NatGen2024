library(edgeR)
library(ggplot2)


#limma_HMEC_108B_evo_vs_diploid_RPKM3, evo2, FR-ev2
RNAseq_2020_sva_codes2_temp <- subset(RNAseq_2020_sva_codes7, (id == "108B_evo" | class == "diploid") & run == 2)
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
limma_HMEC_108B_evo_vs_diploid_RPKM3 <- lrtres
limma_HMEC_108B_evo_vs_diploid_RPKM3[limma_HMEC_108B_evo_vs_diploid_RPKM3 == "NaN"] <- 0
limma_HMEC_108B_evo_vs_diploid_RPKM3 <- merge(limma_HMEC_108B_evo_vs_diploid_RPKM3, Gene_Chrom_pos_hg19_3, by = "genes")
limma_HMEC_108B_evo_vs_diploid_RPKM3 <- merge(limma_HMEC_108B_evo_vs_diploid_RPKM3, chromdata3_grch37, by = 'Chromosome')
limma_HMEC_108B_evo_vs_diploid_RPKM3$genpos <- limma_HMEC_108B_evo_vs_diploid_RPKM3$MID + limma_HMEC_108B_evo_vs_diploid_RPKM3$add
limma_HMEC_108B_evo_vs_diploid_RPKM3 <-  limma_HMEC_108B_evo_vs_diploid_RPKM3[!duplicated(limma_HMEC_108B_evo_vs_diploid_RPKM3$genes),]

# sum(segments108B_evo_bwag37.sorted.bam$copy.numberNew2*segments108B_evo_bwag37.sorted.bam$width)/sum(segments108B_evo_bwag37.sorted.bam$width)
# avg ploidy: 4.1742

ggplot(data = limma_HMEC_108B_evo_vs_diploid_RPKM3, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(limma_HMEC_108B_evo_vs_diploid_RPKM3, Chromosome != "7" & Chromosome != "12" & Chromosome != "19" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(limma_HMEC_108B_evo_vs_diploid_RPKM3, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red3") +
    geom_point(data = subset(limma_HMEC_108B_evo_vs_diploid_RPKM3, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(limma_HMEC_108B_evo_vs_diploid_RPKM3, Chromosome == "19" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(limma_HMEC_108B_evo_vs_diploid_RPKM3, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red3") +
    theme_classic() +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.1742), log2(4/4.1742), log2(5/4.1742),  log2(6/4.1742)))


#repeat edgeR analysis for each line relative to diploid and avg ploidy calculation and plot


#Ext Data Figure 8a

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ctrl13_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC)
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "20"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red")  +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.021372), log2(2/2.021372), log2(3/2.021372)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ctrl13_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "20", 3, 2)
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.021372), log2(2/2.021372))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ae_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC)
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "20" & Chromosome != "10" & Chromosome != "15"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red")  +  geom_point(data = subset(Fig_S8a_data, Chromosome == "15" ), size = 1.2, alpha= 0.3, color = "red")    + geom_point(data = subset(Fig_S8a_data, Chromosome == "10" & Karyotype_band %like% "p" ), size = 1.2, alpha= 0.3, color = "blue") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "10" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "grey60")  +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.032952), log2(2/2.032952), log2(3/2.032952)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ae_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "20" | Fig_S8a_data$Chromosome == "15", 3, ifelse(Fig_S8a_data$Chromosome == "10" & Fig_S8a_data$Karyotype_band %like% "p", "1", 2))
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.032952), log2(2/2.032952), log2(1/2.032952))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "af_ev_b_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC)
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "8" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" ), size = 1.2, alpha= 0.3, color = "red")   +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.049693), log2(2/2.049693), log2(3/2.049693)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "af_ev_b_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "8" , 3, 2)
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.049693), log2(2/2.049693))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ae_ev_d_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC)
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "20" & Chromosome != "10" & Chromosome != "8" & Chromosome != "15"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" ), size = 1.2, alpha= 0.3, color = "red")  +  geom_point(data = subset(Fig_S8a_data, Chromosome == "15" ), size = 1.2, alpha= 0.3, color = "red")    + geom_point(data = subset(Fig_S8a_data, Chromosome == "10" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "blue") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "10" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "grey60")  +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.08514),log2(2/2.08514), log2(3/2.08514)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ae_ev_d_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "20" | Fig_S8a_data$Chromosome == "15" | Fig_S8a_data$Chromosome == "8", 3, ifelse(Fig_S8a_data$Chromosome == "10" & Fig_S8a_data$Karyotype_band %like% "p", "1", 2))
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.08514), log2(2/2.08514), log2(1/2.08514))) + xlim(-2,2) +theme_bw()


Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "da_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC)
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "5" & Chromosome != "17" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "5" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "17" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red3")  +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.132),log2(2/2.132), log2(3/2.132), log2(4/2.132)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "da_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "5" | Fig_S8a_data$Chromosome == "17", 3, ifelse(Fig_S8a_data$Chromosome == "20", 4, 2))
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(4/2.132), log2(3/2.132), log2(2/2.132))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "dc_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC) 
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "9" & Chromosome != "12" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "9" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red")   +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.09),log2(2/2.09), log2(3/2.09)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "dc_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "9" | Fig_S8a_data$Chromosome == "12", 3, 2)
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.09), log2(2/2.09))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "dc_ev1_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC) 
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "8" & Chromosome != "20" & Chromosome != "9" & Chromosome != "12" ), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "9" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    theme_classic() +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red3")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60")    +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.162643), log2(2/2.162643), log2(3/2.162643),log2(4/2.162643)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "dc_ev1_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "20" | Fig_S8a_data$Chromosome == "9" | Fig_S8a_data$Chromosome == "12" , 3, ifelse(Fig_S8a_data$Chromosome == "8" & Fig_S8a_data$Karyotype_band %like% "q", 4, 2 ) )
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(4/2.162643), log2(3/2.162643), log2(2/2.162643))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "df_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC) 
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "2" & Chromosome != "8" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "2" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" ), size = 1.2, alpha= 0.3, color = "red")   +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.131),log2(2/2.131), log2(3/2.131)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "df_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "2" | Fig_S8a_data$Chromosome == "8", 3, 2)
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.131), log2(2/2.131))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "bq_ev_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC) 
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "8" & Chromosome != "20"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" & Karyotype_band %like% "p" ), size = 1.2, alpha= 0.3, color = "grey")  +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    theme_classic() +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.036393), log2(2/2.036393), log2(3/2.036393)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "bq_ev_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "20" | (Fig_S8a_data$Chromosome == "8" & Fig_S8a_data$Karyotype_band %like% "q"), 3, 2)
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.036393), log2(2/2.036393))) + xlim(-2,2) +theme_bw()



Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ae_ev_e_vs_diploid")
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC > 2, 2, Fig_S8a_data$logFC)
Fig_S8a_data$logFC <- ifelse(Fig_S8a_data$logFC < -2, -2, Fig_S8a_data$logFC)
ggplot(data = Fig_S8a_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8a_data, Chromosome != "20" & Chromosome != "10" & Chromosome != "8" & Chromosome != "15"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red")  +  geom_point(data = subset(Fig_S8a_data, Chromosome == "15" ), size = 1.2, alpha= 0.3, color = "red")    + geom_point(data = subset(Fig_S8a_data, Chromosome == "10" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "blue") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "10" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "grey60") + geom_point(data = subset(Fig_S8a_data, Chromosome == "8" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8a_data, Chromosome == "8" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red")  +
    theme_classic()  +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(1/2.08514),log2(2/2.08514), log2(3/2.08514)))

Fig_S8a_data  <- subset(Fig_S8a_data_all, sample == "ae_ev_d_vs_diploid")
Fig_S8a_data$copy.number <- ifelse(Fig_S8a_data$Chromosome == "20" | Fig_S8a_data$Chromosome == "15" | (Fig_S8a_data$Chromosome == "8" &  Fig_S8a_data$Karyotype_band %like% "q"), 3, ifelse(Fig_S8a_data$Chromosome == "10" & Fig_S8a_data$Karyotype_band %like% "p", 1, 2))
ggplot(Fig_S8a_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/2.08514), log2(2/2.08514), log2(1/2.08514))) + xlim(-2,2) +theme_bw()



#Ext Data Figure 8b


Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FG-ev2_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC)
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "2" & Chromosome != "3" & Chromosome != "7" & Chromosome != "12" & Chromosome != "16" & Chromosome != "20" & Chromosome != "22"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "2" & genpos > 380846708), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "2" & genpos < 380846708), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "3" ), size = 1.2, alpha= 0.3, color = "blue2") +
    theme_classic() +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "16" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "22" ), size = 1.2, alpha= 0.3, color = "blue2") +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.949671),log2(4/3.949671), log2(5/3.949671)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FG-ev2_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( Fig_S8b_data$Chromosome == "7" | Fig_S8b_data$Chromosome == "12"  | Fig_S8b_data$Chromosome == "20" , 5, ifelse( (Fig_S8b_data$Chromosome == "2" & Fig_S8b_data$genpos > 380846708 ) | Fig_S8b_data$Chromosome == "16" | Fig_S8b_data$Chromosome == "22" , 3, 4))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.949671), log2(4/3.949671), log2(5/3.949671))) + xlim(-2,2) +theme_bw()


Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FQ_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC)
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "6" & Chromosome != "13" & Chromosome != "17" & Chromosome != "18" & Chromosome != "3"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2") +
    theme_classic() +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "13" & genpos > 2149807307 ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "13" & genpos < 2149807307 ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "17" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "18" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "3" ), size = 1.2, alpha= 0.3, color = "blue2")   +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.863242), log2(4/3.863242), log2(5/3.863242)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FQ_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( Fig_S8b_data$Chromosome == "17" | Fig_S8b_data$Chromosome == "18"  | Fig_S8b_data$Chromosome == "3" |  Fig_S8b_data$Chromosome == "6", 3, ifelse( (Fig_S8b_data$Chromosome == "13" & Fig_S8b_data$genpos > 2149807307) , 5, 4))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.863242), log2(4/3.863242), log2(5/3.863242))) + xlim(-2,2) +theme_bw()


Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FQ_ev2_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC)
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "6" & Chromosome != "8" & Chromosome != "13" & Chromosome != "17" & Chromosome != "18" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "8" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "8" &  Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red") +
    theme_classic() +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "13" & genpos > 2149807307), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "13" & genpos < 2149807307), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "17" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "18" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "blue2") +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.931876),log2(4/3.931876), log2(5/3.931876)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FQ_ev2_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( Fig_S8b_data$Chromosome == "6" | (Fig_S8b_data$Chromosome == "8" & Fig_S8b_data$Karyotype_band %like% "p")  | Fig_S8b_data$Chromosome == "17" |  Fig_S8b_data$Chromosome == "18" | (Fig_S8b_data$Chromosome == "20" & Fig_S8b_data$Karyotype_band %like% "p"), 3, ifelse( (Fig_S8b_data$Chromosome == "8" & Fig_S8b_data$Karyotype_band %like% "q") | (Fig_S8b_data$Chromosome == "13" & Fig_S8b_data$genpos > 2149807307) | (Fig_S8b_data$Chromosome == "20" & Fig_S8b_data$Karyotype_band %like% "q"), 5, 4))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.931876), log2(4/3.931876), log2(5/3.931876))) + xlim(-2,2) +theme_bw()


Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FM-ev3_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC)
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "2" & Chromosome != "4" & Chromosome != "6" & Chromosome != "8" & Chromosome != "9" & Chromosome != "12" & Chromosome != "20"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "2" & Karyotype_band %like% "q" ), size = 1.2, alpha= 0.3, color = "blue2")+
    geom_point(data = subset(Fig_S8b_data, Chromosome == "2" & Karyotype_band %like% "p"  ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "4" & Karyotype_band %like% "q" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "4" & Karyotype_band %like% "p" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2")+
    geom_point(data = subset(Fig_S8b_data, Chromosome == "9" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "8" & genpos > 1469807307 ), size = 1.2, alpha= 0.3, color = "red4") + 
    geom_point(data = subset(Fig_S8b_data, Chromosome == "8" & genpos < 1469807307 ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    theme_classic() +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.95022), log2(4/3.95022), log2(5/3.95022),  log2(6/3.95022), log2(8/3.95022)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FM-ev3_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( (Fig_S8b_data$Chromosome == "2" & Fig_S8b_data$Karyotype_band %like% "q") | (Fig_S8b_data$Chromosome == "4" & Fig_S8b_data$Karyotype_band %like% "q")  | Fig_S8b_data$Chromosome == "6" |  Fig_S8b_data$Chromosome == "9" | (Fig_S8b_data$Chromosome == "8" & Fig_S8b_data$genpos < 1469807307), 3, ifelse(  Fig_S8b_data$Chromosome == "12"  | Fig_S8b_data$Chromosome == "20" , 5, ifelse((Fig_S8b_data$Chromosome == "8" & Fig_S8b_data$genpos > 1469807307), 8,4)))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.95022), log2(4/3.95022), log2(5/3.95022), log2(8/3.95022))) + xlim(-2,2) +theme_bw()

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FM_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC)
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data,  Chromosome != "4" & Chromosome != "6" & Chromosome != "9" & Chromosome != "12" & Chromosome != "20"), size = 1.2, alpha= 0.3, color = "grey60")   +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "4" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "4" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "9" ), size = 1.2, alpha= 0.3, color = "blue2") +
    theme_classic() +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.918948), log2(4/3.918948), log2(5/3.918948)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FM_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( (Fig_S8b_data$Chromosome == "2" & Fig_S8b_data$Karyotype_band %like% "q") | (Fig_S8b_data$Chromosome == "4" & Fig_S8b_data$Karyotype_band %like% "q")  | Fig_S8b_data$Chromosome == "6" |  Fig_S8b_data$Chromosome == "9" | (Fig_S8b_data$Chromosome == "8" & Fig_S8b_data$Karyotype_band %like% "p"), 3, ifelse( (Fig_S8b_data$Chromosome == "8" & Fig_S8b_data$Karyotype_band %like% "q") | Fig_S8b_data$Chromosome == "12"  | Fig_S8b_data$Chromosome == "20" , 5, 4))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.918948), log2(4/3.918948), log2(5/3.918948))) + xlim(-2,2) +theme_bw()



Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FN-ev2_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC)
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "3" & Chromosome != "6" & Chromosome != "8" & Chromosome != "11" & Chromosome != "15" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2") +
    theme_classic() +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "3" & Karyotype_band %like% "q" ), size = 1.2, alpha= 0.3, color = "blue2")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "8" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "11" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "15" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "3" & Karyotype_band %like% "p" ), size = 1.2, alpha= 0.3, color = "grey60")   +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.996421), log2(4/3.996421), log2(5/3.996421)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FN-ev2_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( (Fig_S8b_data$Chromosome == "3" & Fig_S8b_data$Karyotype_band %like% "q") |  Fig_S8b_data$Chromosome == "6" |  Fig_S8b_data$Chromosome == "15"  , 3, ifelse(  Fig_S8b_data$Chromosome == "8"  | Fig_S8b_data$Chromosome == "20" , 5, 4))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.996421), log2(4/3.996421), log2(5/3.996421))) + xlim(-2,2) +theme_bw()



Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FD-ev-A_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC) 
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "1" & Chromosome != "2" &  Chromosome != "3" & Chromosome != "5" & Chromosome != "6" & Chromosome != "7" & Chromosome != "8" & Chromosome != "10" & Chromosome != "15" & Chromosome != "16" & Chromosome != "17" & Chromosome != "22"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "1"  & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "1"  & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "blue2")  +
    theme_classic() +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "2"  ), size = 1.2, alpha= 0.3, color = "blue2")   +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "3" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "3" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "5" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "blue2")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "8" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "10" ), size = 1.2, alpha= 0.3, color = "blue2")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "15" ), size = 1.2, alpha= 0.3, color = "blue2")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "16" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "17" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "22" ), size = 1.2, alpha= 0.3, color = "blue2") +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.580342), log2(4/3.580342), log2(5/3.580342), log2(6/3.580342)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FD-ev-A_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( (Fig_S8b_data$Chromosome == "3" & Fig_S8b_data$Karyotype_band %like% "q") | (Fig_S8b_data$Chromosome == "1" & Fig_S8b_data$Karyotype_band %like% "p") |  Fig_S8b_data$Chromosome == "6" |  Fig_S8b_data$Chromosome == "5" |  Fig_S8b_data$Chromosome == "7" |  Fig_S8b_data$Chromosome == "10" |  Fig_S8b_data$Chromosome == "15" |  Fig_S8b_data$Chromosome == "16" |  Fig_S8b_data$Chromosome == "17" |  Fig_S8b_data$Chromosome == "22" , 3, ifelse(  (Fig_S8b_data$Chromosome == "1" & Fig_S8b_data$Karyotype_band %like% "q")  | Fig_S8b_data$Chromosome == "8" , 5, 4))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.580342), log2(4/3.580342), log2(5/3.580342))) + xlim(-2,2) +theme_bw()


Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FM-ev2_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC) 
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "2" & Chromosome != "4" & Chromosome != "6"  & Chromosome != "9" & Chromosome != "16" & Chromosome != "12" & Chromosome != "20"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "2" & Karyotype_band %like% "q" ), size = 1.2, alpha= 0.3, color = "blue2")+
    geom_point(data = subset(Fig_S8b_data, Chromosome == "2" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "4" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "4" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2")+
    geom_point(data = subset(Fig_S8b_data, Chromosome == "9" ), size = 1.2, alpha= 0.3, color = "blue2")  +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red")+
    geom_point(data = subset(Fig_S8b_data, Chromosome == "16" ), size = 1.2, alpha= 0.3, color = "blue2") +
    theme_classic() +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/3.868495), log2(4/3.868495), log2(5/3.868495)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FM-ev2_vs_diploid")
Fig_S8b_data$copy.number <- ifelse( (Fig_S8b_data$Chromosome == "2" & Fig_S8b_data$Karyotype_band %like% "q") | (Fig_S8b_data$Chromosome == "4" & Fig_S8b_data$Karyotype_band %like% "q") |  Fig_S8b_data$Chromosome == "6" |  Fig_S8b_data$Chromosome == "9" |  Fig_S8b_data$Chromosome == "16" , 3, ifelse(  (Fig_S8b_data$Chromosome == "2" & Fig_S8b_data$Karyotype_band %like% "p")  | Fig_S8b_data$Chromosome == "12" | Fig_S8b_data$Chromosome == "20" , 5, 4))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/3.868495), log2(4/3.868495), log2(5/3.868495))) + xlim(-2,2) +theme_bw()



Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FX-ev2-B_vs_diploid")
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC > 2, 2, Fig_S8b_data$logFC)
Fig_S8b_data$logFC <- ifelse(Fig_S8b_data$logFC < -2, -2, Fig_S8b_data$logFC)
ggplot(data = Fig_S8b_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8b_data, Chromosome != "3" & Chromosome != "5" &  Chromosome != "6" & Chromosome != "7" & Chromosome != "8" & Chromosome != "9" & Chromosome != "11" & Chromosome != "12" & Chromosome != "13" & Chromosome != "15" & Chromosome != "18" & Chromosome != "19" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "3" ), size = 1.2, alpha= 0.3, color = "blue2")  +
    theme_classic() +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "5"  ), size = 1.2, alpha= 0.3, color = "red")   +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "6" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "7"  ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "8" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "9" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "11" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "11"  & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red")    +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "13" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "15" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "18" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "19" ), size = 1.2, alpha= 0.3, color = "blue2") +
    geom_point(data = subset(Fig_S8b_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red3") +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.089825), log2(4/4.089825), log2(5/4.089825), log2(6/4.089825)))

Fig_S8b_data  <- subset(Fig_S8b_data_all, sample == "FX-ev2-B_vs_diploid")
Fig_S8b_data$copy.number <- ifelse(    Fig_S8b_data$Chromosome == "6" |  Fig_S8b_data$Chromosome == "3" |  Fig_S8b_data$Chromosome == "13" |  Fig_S8b_data$Chromosome == "15" |  Fig_S8b_data$Chromosome == "18" |  Fig_S8b_data$Chromosome == "19" , 3, ifelse(  (Fig_S8b_data$Chromosome == "11" & Fig_S8b_data$Karyotype_band %like% "p")  | Fig_S8b_data$Chromosome == "8" | Fig_S8b_data$Chromosome == "7" | Fig_S8b_data$Chromosome == "9" | Fig_S8b_data$Chromosome == "12"  , 5, ifelse(  Fig_S8b_data$Chromosome == "20" & Fig_S8b_data$genpos < 2769807307, 6, ifelse( Fig_S8b_data$Chromosome == "20" & Fig_S8b_data$genpos > 2769807307 ,7,4) ) ))
ggplot(Fig_S8b_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/4.089825), log2(4/4.089825), log2(5/4.089825), log2(6/4.089825), log2(7/4.089825))) + xlim(-2,2) +theme_bw()



#Ext Data Figure 8c

Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-L_vs_diploid")
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC > 2, 2, Fig_S8c_data$logFC)
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC < -2, -2, Fig_S8c_data$logFC)
ggplot(data = Fig_S8c_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8c_data, Chromosome != "7" & Chromosome != "8" & Chromosome != "11" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red3") +
    theme_classic() +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "8"  ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "11" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red")   +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.225384), log2(4/4.225384), log2(5/4.225384), log2(6/4.225384)))

Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-L_vs_diploid")
Fig_S8c_data$copy.number <- ifelse(  Fig_S8c_data$Chromosome == "8" |  Fig_S8c_data$Chromosome == "11" |  Fig_S8c_data$Chromosome == "20" , 5, ifelse(  Fig_S8c_data$Chromosome == "7" , 6, 4))
ggplot(Fig_S8c_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(6/4.225384), log2(4/4.225384), log2(5/4.225384))) + xlim(-2,2) +theme_bw()


Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-L-ev1_vs_diploid")
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC > 2, 2, Fig_S8c_data$logFC)
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC < -2, -2, Fig_S8c_data$logFC)
ggplot(data = Fig_S8c_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8c_data, Chromosome != "7" & Chromosome != "8" &  Chromosome != "20" &  Chromosome != "1" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "1" & Karyotype_band %like% "q" ), size = 1.2, alpha= 0.3, color = "red3") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "1" & Karyotype_band %like% "p" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red3")  +
    theme_classic() +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "8"  ), size = 1.2, alpha= 0.3, color = "red")   +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red3")   +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.297832), log2(4/4.297832), log2(5/4.297832), log2(6/4.297832)))

Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-L-ev1_vs_diploid")
Fig_S8c_data$copy.number <- ifelse(   Fig_S8c_data$Chromosome == "8" , 5, ifelse(  Fig_S8c_data$Chromosome == "7" |  Fig_S8c_data$Chromosome == "20" | (Fig_S8c_data$Chromosome == "1" & Fig_S8c_data$Karyotype_band %like% "q" ) , 6, 4))
ggplot(Fig_S8c_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(6/4.297832), log2(4/4.297832), log2(5/4.297832))) + xlim(-2,2) +theme_bw()


Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-G_vs_diploid")
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC > 2, 2, Fig_S8c_data$logFC)
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC < -2, -2, Fig_S8c_data$logFC)
ggplot(data = Fig_S8c_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8c_data, Chromosome != "7" & Chromosome != "8" &  Chromosome != "20" & Chromosome != "11" & Chromosome != "12" & Chromosome != "16" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red3")  +
    theme_classic() +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "8"  ), size = 1.2, alpha= 0.3, color = "red3")   +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "11" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "16" ), size = 1.2, alpha= 0.3, color = "blue2")    +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.28993), log2(4/4.28993), log2(5/4.28993), log2(6/4.28993)))

Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-G_vs_diploid")
Fig_S8c_data$copy.number <- ifelse(   Fig_S8c_data$Chromosome == "12" |  Fig_S8c_data$Chromosome == "11" |  Fig_S8c_data$Chromosome == "20" , 5, ifelse(  Fig_S8c_data$Chromosome == "8" | Fig_S8c_data$Chromosome == "7" , 6, ifelse( Fig_S8c_data$Chromosome == "16",3, 4) )  )
ggplot(Fig_S8c_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/4.28993), log2(6/4.28993), log2(4/4.28993), log2(5/4.28993))) + xlim(-2,2) +theme_bw()



Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-G-ev1_vs_diploid")
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC > 2, 2, Fig_S8c_data$logFC)
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC < -2, -2, Fig_S8c_data$logFC)
ggplot(data = Fig_S8c_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8c_data, Chromosome != "7" & Chromosome != "8" &  Chromosome != "20" & Chromosome != "12" & Chromosome != "16" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red3")  +
    theme_classic() +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "8"  ), size = 1.2, alpha= 0.3, color = "red3")   +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "1" & Karyotype_band %like% "q"  ), size = 1.2, alpha= 0.3, color = "red4") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "1" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "16" ), size = 1.2, alpha= 0.3, color = "blue2")    +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.407001), log2(4/4.407001), log2(5/4.407001), log2(6/4.407001), log2(8/4.407001)))

Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-G-ev1_vs_diploid")
Fig_S8c_data$copy.number <- ifelse(   Fig_S8c_data$Chromosome == "12" |  Fig_S8c_data$Chromosome == "20" , 5, ifelse(  Fig_S8c_data$Chromosome == "8" | Fig_S8c_data$Chromosome == "7" , 6, ifelse( Fig_S8c_data$Chromosome == "16",3, ifelse( Fig_S8c_data$Chromosome == "1" & Fig_S8c_data$Karyotype_band %like% "q" ,8,4) ) )  )
ggplot(Fig_S8c_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/4.407001),log2(8/4.407001), log2(6/4.407001), log2(4/4.407001), log2(5/4.407001))) + xlim(-2,2) +theme_bw()



Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-I_vs_diploid")
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC > 2, 2, Fig_S8c_data$logFC)
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC < -2, -2, Fig_S8c_data$logFC)
ggplot(data = Fig_S8c_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8c_data, Chromosome != "7" & Chromosome != "8" &  Chromosome != "20" &  Chromosome != "11" & Chromosome != "12" & Chromosome != "16" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red3")  +
    theme_classic() +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "8"  ), size = 1.2, alpha= 0.3, color = "red3")   +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "11" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "16" ), size = 1.2, alpha= 0.3, color = "blue2")    +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.289887), log2(4/4.289887), log2(5/4.289887), log2(6/4.289887)))

Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-I_vs_diploid")
Fig_S8c_data$copy.number <- ifelse(   Fig_S8c_data$Chromosome == "12" |  Fig_S8c_data$Chromosome == "11" |  Fig_S8c_data$Chromosome == "20" , 5, ifelse(  Fig_S8c_data$Chromosome == "8" | Fig_S8c_data$Chromosome == "7" , 6, ifelse( Fig_S8c_data$Chromosome == "16",3, 4) )  )
ggplot(Fig_S8c_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/4.28993), log2(6/4.28993), log2(4/4.28993), log2(5/4.28993))) + xlim(-2,2) +theme_bw()



Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-I-ev1_vs_diploid")
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC > 2, 2, Fig_S8c_data$logFC)
Fig_S8c_data$logFC <- ifelse(Fig_S8c_data$logFC < -2, -2, Fig_S8c_data$logFC)
ggplot(data = Fig_S8c_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(Fig_S8c_data, Chromosome != "1" & Chromosome != "7" & Chromosome != "8" &  Chromosome != "20"  & Chromosome != "12" & Chromosome != "16" & Chromosome != "20" ), size = 1.2, alpha= 0.3, color = "grey60") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "7" ), size = 1.2, alpha= 0.3, color = "red3")  +
    theme_classic() +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "8"  ), size = 1.2, alpha= 0.3, color = "red3")   +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "20" ), size = 1.2, alpha= 0.3, color = "red")  +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "1" & Karyotype_band %like% "q" ), size = 1.2, alpha= 0.3, color = "red3") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "1" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60")+
    geom_point(data = subset(Fig_S8c_data, Chromosome == "12" ), size = 1.2, alpha= 0.3, color = "red") +
    geom_point(data = subset(Fig_S8c_data, Chromosome == "16" ), size = 1.2, alpha= 0.3, color = "blue2")    +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome, expand = c(0,0)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("log2FC") +
    xlab("Chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(3/4.315126), log2(4/4.315126), log2(5/4.315126), log2(6/4.315126)))

Fig_S8c_data  <- subset(Fig_S8c_data_all, sample == "CQ-ev-I_vs_diploid")
Fig_S8c_data$copy.number <- ifelse(   Fig_S8c_data$Chromosome == "12" |   Fig_S8c_data$Chromosome == "20" , 5, ifelse(  (Fig_S8c_data$Chromosome == "1" & Fig_S8c_data$Karyotype_band %like% "q") | Fig_S8c_data$Chromosome == "8" | Fig_S8c_data$Chromosome == "7" , 6, ifelse( Fig_S8c_data$Chromosome == "16",3, 4) )  )
ggplot(Fig_S8c_data, aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/4.315126), log2(6/4.315126), log2(4/4.315126), log2(5/4.315126))) + xlim(-2,2) +theme_bw()


#Ext Data Figure 8d

ggplot(subset(All_logFC_2N_aneuploids, copy.number != "NA"), aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(1/2), log2(2/2), log2(3/2),  log2(4/2))) + xlim(-2,2) +theme_bw()


#Ext Data Figure 8e

ggplot(subset(All_logFC_4N_aneuploids, copy.number != "NA"), aes(x=logFC, color=as.character(copy.number))) +
    geom_density() +geom_vline(xintercept = c(log2(3/4), log2(2/2), log2(5/4),  log2(6/4) , log2(7/4) ,  log2(8/4))) + xlim(-2,2) +theme_bw()



