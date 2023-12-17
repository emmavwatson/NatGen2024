library(ggplot2)
library(ggcyto)

#EXT DATA FIGURE 1A
ggplot(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm, aes(x=logFC, color=tissue_specificity)) + geom_density() +  theme_bw()
ggplot(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm, aes(x=logFC_invitro, color=tissue_specificity)) + geom_density() +  theme_bw()

#EXT DATA FIGURE 1B
ggplot(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm, aes(x=logFC, y = logFC_invitro )) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(data=subset(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm, tissue_specificity == "nonspecific"), color="grey20", size=1, alpha = 0.3) + geom_point(data=subset(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm, tissue_specificity == "Breast" | genes == "ESR1" ), color="red",  size=1) + geom_point(data=subset(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm, tissue_specificity == "Kidney"), color="blue",  size=1) + theme_bw() + geom_text_repel(data=subset(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm, (tissue_specificity == "Breast" & (logFC > 0 &  logFC_invitro > 0 )) | genes == "ESR1" | genes == "KRT5" | genes == "KRT17" | genes %like% "KRT6" | genes %like% "CLCA2") , aes(label =  genes), color = "red", size=3, max.overlaps=100) + geom_text_repel(data=subset(limma_RPTEC_vs_HMEC_RPKM_vs_humanTumor_lcrm,  (tissue_specificity == "Kidney" & (logFC < -4 & logFC_invitro < -2 )) | (tissue_specificity == "Kidney" & (logFC < -1.5 & logFC_invitro < -7))) , aes(label =  genes), color = "blue", size=3, max.overlaps=100) + geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = after_stat(p.value), sep="\n" ))

#EXT DATA FIGURE 1C
#see separate repository https://github.com/emmavwatson/CNAplot for all raw gDNA sequencing processing, copy number calling, and plotting code 

#EXT DATA FIGURE 1D
#HMEC single cells

#file uploading and pre-processing

#fwalist_081622scp4 <- list.files("~/Documents/gDNA_HMEC_scWGA4_081622p4_vwr")
#fwalist_081622scp2 <- list.files("~/Documents/gDNA_HMEC_scWGA4_081622p2_vwr")
#fwalist_081622scp3 <- list.files("~/Documents/gDNA_HMEC_scWGA4_081622p3_vwr")
#fwalist_081622scp1 <- list.files("~/Documents/gDNA_HMEC_scWGA4_081622p1_vwr")
#fwalist_081622scp1 <- data.frame(fwalist_081622scp1)
#fwalist_081622scp2 <- data.frame(fwalist_081622scp2)
#fwalist_081622scp3 <- data.frame(fwalist_081622scp3)
#fwalist_081622scp4 <- data.frame(fwalist_081622scp4)
#colnames(fwalist_081622scp4)[1] <- "fwa"
#colnames(fwalist_081622scp3)[1] <- "fwa"
#colnames(fwalist_081622scp2)[1] <- "fwa"
#colnames(fwalist_081622scp1)[1] <- "fwa"
#fwalist_081622scp_all <- rbind(fwalist_081622scp1, fwalist_081622scp2, fwalist_081622scp3, fwalist_081622scp4)
#fwalist_081622scp_all$fwa <- substring(fwalist_081622scp_all$fwa, 28)
#fwalist_081622scp_all$fwa <- substr(fwalist_081622scp_all$fwa,1,nchar(fwalist_081622scp_all$fwa)-12)
#gDNAscWGA4_081622scp1_vwr_min20_wFilter <- Merging_CNA_bin_files("~/Documents/gDNA_HMEC_scWGA4_081622_ALL_segments")
#gDNAscWGA4_081622scp1_vwr_min20_wFilter$sample<- substring(gDNAscWGA4_081622scp1_vwr_min20_wFilter$sample, 20)
#gDNAscWGA4_081622scp1_vwr_min20_wFilter$sample <- substr(gDNAscWGA4_081622scp1_vwr_min20_wFilter$sample,1,nchar(gDNAscWGA4_081622scp1_vwr_min20_wFilter$sample)-4)

#gDNAscWGA4_081622scp1_vwr_min20_wFilter <- subset(gDNAscWGA4_081622scp1_vwr_min20_wFilter, sample %in% fwalist_081622scp_all$fwa)
#gDNAscWGA4_081622scp1_vwr_min20_wFilter$copy.number2 <- gDNAscWGA4_081622scp1_vwr_min20_wFilter$copy.number
#colnames(gDNAscWGA4_081622scp1_vwr_min20_wFilter)[11] <- "Sample_ID"
#colnames(gDNAscWGA4_081622scp1_vwr_min20_wFilter)[1] <- "Chromosome"
#colnames(gDNAscWGA4_081622scp1_vwr_min20_wFilter)[2] <- "Start"
#colnames(gDNAscWGA4_081622scp1_vwr_min20_wFilter)[3] <- "End"
#colnames(gDNAscWGA4_081622scp1_vwr_min20_wFilter)[4] <- "length"

#gDNAscWGA4_081622scp1_vwr_min20_wFilter2 <- gDNAscWGA4_081622scp1_vwr_min20_wFilter
#gDNAscWGA4_081622scp1_vwr_min20_wFilter2 <- gDNAscWGA4_081622scp1_vwr_min20_wFilter2[,c(1,11,12)]
#gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast <- dcast(gDNAscWGA4_081622scp1_vwr_min20_wFilter2, Sample_ID~Chromosome, fun.aggregate = mean)
#rownames(gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast) <- gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast$Sample_ID
#gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast <- gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast[,-1]
#gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast[is.na(gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast)] <- 0
#gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast_hc <- hclust(d=dist(gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast))
#fwa <- gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast_hc$labels[gDNAscWGA4_081622scp1_vwr_min20_wFilter2_dcast_hc$order]
#fwa <- data.frame(fwa)
#fwa$order <- rownames(fwa)
#colnames(fwa)[1] <- "Sample_ID"
#gDNAscWGA4_081622scp1_vwr_min20_wFilter <- merge(gDNAscWGA4_081622scp1_vwr_min20_wFilter, fwa, by = "Sample_ID", all= TRUE)
#gDNAscWGA4_081622scp1_vwr_min20_wFilter$order <- as.numeric(as.character(gDNAscWGA4_081622scp1_vwr_min20_wFilter$order))
#gDNAscWGA4_081622scp1_vwr_min20_wFilter <- gDNAscWGA4_081622scp1_vwr_min20_wFilter[order(gDNAscWGA4_081622scp1_vwr_min20_wFilter$order),]


dev.off()
pdf(file = "HMEC_single_cell_profiles_hclust.pdf", width = 8, height = 8)
ggplot(data = gDNAscWGA4_081622scp1_vwr_min20_wFilter, aes(x = (start.genome+end.genome)/2, y = copy.number2)) + geom_point(size = 1.2, alpha= 0.3) + theme_classic() + geom_segment(data=gDNAscWGA4_081622scp1_vwr_min20_wFilter, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4, color='copy.number2'), size=10) + scale_x_continuous(breaks=chromdata_hg19$chromlength/2+df.chroms[-24,1], labels=chromdata_hg19$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=1) + ylim(4,4) + geom_segment(data=whiteout, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4), color="white", size=10) + theme(legend.position="none") + theme(plot.margin=grid::unit(c(5,0,5,0), "cm")) + ylab("") + xlab("") + scale_colour_gradient2(limits = c(0.96, 3.05), low="blue", mid="gray90", high="red", na.value="firebrick", midpoint=2) + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.margin = unit(0, "lines")) +facet_wrap(~order, nrow=109) + theme(strip.text.x = element_blank())
dev.off()



#EXT DATA FIGURE 1E
#RPTEC single cells

#file uploading and pre-processing

#fwalist_080322scp2 <- list.files("~/Documents/gDNA_RPTEC_scWGA4_080322p2_vwr")
#fwalist_080322scp3 <- list.files("~/Documents/gDNA_RPTEC_scWGA4_080322p3_vwr")

#fwalist_080322scp2 <- data.frame(fwalist_080322scp2)
#fwalist_080322scp3 <- data.frame(fwalist_080322scp3)
#colnames(fwalist_080322scp2)[1] <- "fwa"
#colnames(fwalist_080322scp3)[1] <- "fwa"

#fwalist_080322scp_all <- rbind(fwalist_080322scp2, fwalist_080322scp3)
#fwalist_080322scp_all$fwa <- substring(fwalist_080322scp_all$fwa, 28)
#fwalist_080322scp_all$fwa <- substr(fwalist_080322scp_all$fwa,1,nchar(fwalist_080322scp_all$fwa)-12)

#gDNAscWGA4_080322scp1_vwr_min20_wFilter <- Merging_CNA_bin_files("~/Documents/gDNA_RPTEC_scWGA4_080322_ALL_segments")
#gDNAscWGA4_080322scp1_vwr_min20_wFilter$sample<- substring(gDNAscWGA4_080322scp1_vwr_min20_wFilter$sample, 20)
#gDNAscWGA4_080322scp1_vwr_min20_wFilter$sample <- substr(gDNAscWGA4_080322scp1_vwr_min20_wFilter$sample,1,nchar(gDNAscWGA4_080322scp1_vwr_min20_wFilter$sample)-4)
#gDNAscWGA4_080322scp1_vwr_min20_wFilter <- subset(gDNAscWGA4_080322scp1_vwr_min20_wFilter, sample %in% fwalist_080322scp_all$fwa)

#gDNAscWGA4_080322scp1_vwr_min20_wFilter$copy.number2 <- gDNAscWGA4_080322scp1_vwr_min20_wFilter$copy.number
#colnames(gDNAscWGA4_080322scp1_vwr_min20_wFilter)[11] <- "Sample_ID"
#colnames(gDNAscWGA4_080322scp1_vwr_min20_wFilter)[1] <- "Chromosome"
#colnames(gDNAscWGA4_080322scp1_vwr_min20_wFilter)[2] <- "Start"
#colnames(gDNAscWGA4_080322scp1_vwr_min20_wFilter)[3] <- "End"
#colnames(gDNAscWGA4_080322scp1_vwr_min20_wFilter)[4] <- "length"

#gDNAscWGA4_080322scp1_vwr_min20_wFilter2 <- gDNAscWGA4_080322scp1_vwr_min20_wFilter
#gDNAscWGA4_080322scp1_vwr_min20_wFilter2 <- gDNAscWGA4_080322scp1_vwr_min20_wFilter2[,c(1,11,12)]
#gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast <- dcast(gDNAscWGA4_080322scp1_vwr_min20_wFilter2, Sample_ID~Chromosome, fun.aggregate = mean)
#rownames(gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast) <- gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast$Sample_ID
#gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast <- gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast[,-1]
#gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast[is.na(gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast)] <- 0
#gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast_hc <- hclust(d=dist(gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast))
#fwa <- gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast_hc$labels[gDNAscWGA4_080322scp1_vwr_min20_wFilter2_dcast_hc$order]
#fwa <- data.frame(fwa)
#fwa$order <- rownames(fwa)
#colnames(fwa)[1] <- "Sample_ID"
#gDNAscWGA4_080322scp1_vwr_min20_wFilter <- merge(gDNAscWGA4_080322scp1_vwr_min20_wFilter, fwa, by = "Sample_ID", all= TRUE)
#gDNAscWGA4_080322scp1_vwr_min20_wFilter$order <- as.numeric(as.character(gDNAscWGA4_080322scp1_vwr_min20_wFilter$order))
#gDNAscWGA4_080322scp1_vwr_min20_wFilter <- gDNAscWGA4_080322scp1_vwr_min20_wFilter[order(gDNAscWGA4_080322scp1_vwr_min20_wFilter$order),]


dev.off()
pdf(file = "RPTEC_single_cell_profiles_hclust.pdf", width = 8, height = 8)
ggplot(data = gDNAscWGA4_080322scp1_vwr_min20_wFilter, aes(x = (start.genome+end.genome)/2, y = copy.number2)) + geom_point(size = 1.2, alpha= 0.3) + theme_classic() + geom_segment(data=gDNAscWGA4_080322scp1_vwr_min20_wFilter, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4, color='copy.number2'), size=10) + scale_x_continuous(breaks=chromdata_hg19$chromlength/2+df.chroms[-24,1], labels=chromdata_hg19$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=1) + ylim(4,4) + geom_segment(data=whiteout, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4), color="white", size=10) + theme(legend.position="none") + theme(plot.margin=grid::unit(c(5,0,5,0), "cm")) + ylab("") + xlab("") + scale_colour_gradient2(limits = c(0.96, 3.05), low="blue", mid="gray90", high="red", na.value="firebrick", midpoint=2) + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.margin = unit(0, "lines")) +facet_wrap(~order, nrow=109) + theme(strip.text.x = element_blank())
dev.off()

#EXT DATA FIGURE 1F
#microscopy images
#PI stain flow data

#file upload and pre-processing

#fs <- read.flowSet(path = "~/Downloads/PI_test", pattern = ".fcs", alter.names = T)
#for (i in 1:length(fs)) {
#    fs[[i]]  <-  transform(fs[[i]],
#                           `SSC.A`=log10(`SSC.A`),
#    )
#}

g.singlets <- polygonGate(filterId = "singlets","FSC.A"=c(40, 400, 550, 115),"PerCP.Cy5.5.A"=c(125, 125, 450, 450))
gs <- GatingSet(fs)
add(gs,g.singlets)
recompute(gs)
ggcyto(gs[[2]],aes(x=PerCP.Cy5.5.A ,y=FSC.A),subset="singlets")+geom_hex(bins = 200) + ggcyto_par_set(limits = list(x = c(0,1000), y = c(0, 1000))) +theme_bw()
ggcyto(gs[[2]], aes(x = PerCP.Cy5.5.A), subset = "singlets") + geom_density(aes(y = ..density..)) +theme_bw() + ggcyto_par_set(limits = list(x = c(0,1000), y = c(0, 0.016)))


g.singlets4N <- polygonGate(filterId = "singlets4N","FSC.A"=c(80, 1020, 1020, 230),"PerCP.Cy5.5.A"=c(250, 250, 800, 800))
gs <- GatingSet(fs)
add(gs,g.singlets4N)
recompute(gs)
ggcyto(gs[[3]],aes(x=PerCP.Cy5.5.A ,y=FSC.A),subset="singlets4N")+geom_hex(bins = 200) + ggcyto_par_set(limits = list(x = c(0,1000), y = c(0, 1000))) +theme_bw()
ggcyto(gs[[3]], aes(x = PerCP.Cy5.5.A), subset = "singlets4N") + geom_density(aes(y = ..density..)) +theme_bw() + ggcyto_par_set(limits = list(x = c(0,1000), y = c(0, 0.016)))


#EXT DATA FIGURE 1G

g.singlets <- polygonGate(filterId = "singlets","FSC.A"=c(40, 400, 550, 115),"PerCP.Cy5.5.A"=c(125, 125, 450, 450))
gs <- GatingSet(fs)
add(gs,g.singlets)
recompute(gs)
ggcyto(gs[[2]], aes(x = PerCP.Cy5.5.A), subset = "singlets") + geom_density(aes(y = ..density..)) +theme_bw() + ggcyto_par_set(limits = list(x = c(0,1000), y = c(0, 0.016)))


g.singlets4N <- polygonGate(filterId = "singlets4N","FSC.A"=c(80, 1020, 1020, 230),"PerCP.Cy5.5.A"=c(250, 250, 800, 800))
gs <- GatingSet(fs)
add(gs,g.singlets4N)
recompute(gs)
ggcyto(gs[[3]], aes(x = PerCP.Cy5.5.A), subset = "singlets4N") + geom_density(aes(y = ..density..)) +theme_bw() + ggcyto_par_set(limits = list(x = c(0,1000), y = c(0, 0.016)))


#EXT DATA FIGURE 1H
#cell size analysis from imaging in ImageJ


#EXT DATA FIGURE 1I
ggplot(X2N_4N_PI_vs_diploid_ctrl, aes(x=FC_FSC2N, y = FC_perc2N)) + geom_point(aes(color = ctrls)) +theme_bw()

#EXT DATA FIGURE 1J TOP 
#screen #1 HMEC tets segment files

plots_ALL_FIG1_HMECtets <- heatmapGenomewide9b(returned_models_102517_tets3, files_102517_tets3, namessegments_102317_plate1_new4, 10000000, 4, 0.2, 0.2, 1.98, 6.05, "euclidean", "complete", 20)


#EXT DATA FIGURE 1J BOTTOM 
#screen #2 HMEC tets segment files

#file upload and pre-processing

#fwalist_121222_tets <- list.files("~/Documents/gDNA_120122_vwr_tet_models")
#fwalist_121222_tets <- data.frame(fwalist_121222_tets)
#colnames(fwalist_121222_tets)[1] <- "fwa"
#fwalist_121222_tets$fwa <- substring(fwalist_121222_tets$fwa, 28)
#fwalist_121222_tets$fwa <- substr(fwalist_121222_tets$fwa,1,nchar(fwalist_121222_tets$fwa)-12)
#gDNAscreen2_tets_vwr_min20_wFilter <- Merging_CNA_bin_files("~/Documents/gDNA_HMEC_tets_screen2_ALL_segments")
#gDNAscreen2_tets_vwr_min20_wFilter$sample<- substring(gDNAscreen2_tets_vwr_min20_wFilter$sample, 20)
#gDNAscreen2_tets_vwr_min20_wFilter$sample <- substr(gDNAscreen2_tets_vwr_min20_wFilter$sample,1,nchar(gDNAscreen2_tets_vwr_min20_wFilter$sample)-4)
#gDNAscreen2_tets_vwr_min20_wFilter <- subset(gDNAscreen2_tets_vwr_min20_wFilter, sample %in% fwalist_121222_tets$fwa)
#gDNAscreen2_tets_vwr_min20_wFilter$copy.number2 <- gDNAscreen2_tets_vwr_min20_wFilter$copy.number
#colnames(gDNAscreen2_tets_vwr_min20_wFilter)[11] <- "Sample_ID"
#colnames(gDNAscreen2_tets_vwr_min20_wFilter)[1] <- "Chromosome"
#colnames(gDNAscreen2_tets_vwr_min20_wFilter)[2] <- "Start"
#colnames(gDNAscreen2_tets_vwr_min20_wFilter)[3] <- "End"
#colnames(gDNAscreen2_tets_vwr_min20_wFilter)[4] <- "length"
#gDNAscreen2_tets_vwr_min20_wFilter2 <- gDNAscreen2_tets_vwr_min20_wFilter
#gDNAscreen2_tets_vwr_min20_wFilter2 <- gDNAscreen2_tets_vwr_min20_wFilter2[,c(1,11,12)]
#gDNAscreen2_tets_vwr_min20_wFilter2_dcast <- dcast(gDNAscreen2_tets_vwr_min20_wFilter2, Sample_ID~Chromosome, fun.aggregate = mean)
#rownames(gDNAscreen2_tets_vwr_min20_wFilter2_dcast) <- gDNAscreen2_tets_vwr_min20_wFilter2_dcast$Sample_ID
#gDNAscreen2_tets_vwr_min20_wFilter2_dcast <- gDNAscreen2_tets_vwr_min20_wFilter2_dcast[,-1]
#gDNAscreen2_tets_vwr_min20_wFilter2_dcast[is.na(gDNAscreen2_tets_vwr_min20_wFilter2_dcast)] <- 0
#gDNAscreen2_tets_vwr_min20_wFilter2_dcast_hc <- hclust(d=dist(gDNAscreen2_tets_vwr_min20_wFilter2_dcast))
#fwa <- gDNAscreen2_tets_vwr_min20_wFilter2_dcast_hc$labels[gDNAscreen2_tets_vwr_min20_wFilter2_dcast_hc$order]
#fwa <- data.frame(fwa)
#fwa$order <- rownames(fwa)
#colnames(fwa)[1] <- "Sample_ID"
#gDNAscreen2_tets_vwr_min20_wFilter <- merge(gDNAscreen2_tets_vwr_min20_wFilter, fwa, by = "Sample_ID", all= TRUE)
#gDNAscreen2_tets_vwr_min20_wFilter$order <- as.numeric(as.character(gDNAscreen2_tets_vwr_min20_wFilter$order))
#gDNAscreen2_tets_vwr_min20_wFilter <- gDNAscreen2_tets_vwr_min20_wFilter[order(gDNAscreen2_tets_vwr_min20_wFilter$order),]


dev.off()
pdf(file = "HMEC_screen2_cell_profiles_hclust_tets.pdf", width = 8, height = 8)
ggplot(data = gDNAscreen2_tets_vwr_min20_wFilter, aes(x = (start.genome+end.genome)/2, y = copy.number)) + geom_point(size = 1.2, alpha= 0.3) + theme_classic() + geom_segment(data=gDNAscreen2_tets_vwr_min20_wFilter, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4, color='copy.number'), size=10) + scale_x_continuous(breaks=chromdata_hg19$chromlength/2+df.chroms[-24,1], labels=chromdata_hg19$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=1) + ylim(4,4) + geom_segment(data=whiteout, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4), color="white", size=10) + theme(legend.position="none") + theme(plot.margin=grid::unit(c(6,0,6,0), "cm")) + ylab("") + xlab("") + scale_colour_gradient2(limits = c(1.96, 6.05), low="blue", mid="gray90", high="red", na.value="firebrick", midpoint=4) + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.margin = unit(0, "lines")) +facet_wrap(~order, nrow=38) + theme(strip.text.x = element_blank())
dev.off()


#EXT DATA FIGURE 1L
#RPTEC screen replicate # 2 

#file upload and pre-processing

#fwalist_112822_RPTECscreen2 <- list.files("~/Documents/gDNA_RPTEC_screen2a_2b_vwr")
#fwalist_112822_RPTECscreen2 <- data.frame(fwalist_112822_RPTECscreen2)
#colnames(fwalist_112822_RPTECscreen2)[1] <- "fwa"
#fwalist_112822_RPTECscreen2$fwa <- substring(fwalist_112822_RPTECscreen2$fwa, 28)
#fwalist_112822_RPTECscreen2$fwa <- substr(fwalist_112822_RPTECscreen2$fwa,1,nchar(fwalist_112822_RPTECscreen2$fwa)-12)
#gDNAscreen2_RPTEC_vwr_min20_wFilter <- Merging_CNA_bin_files("~/Documents/gDNA_RPTEC_screen2a_2b_vwr_segments")
#gDNAscreen2_RPTEC_vwr_min20_wFilter$sample<- substring(gDNAscreen2_RPTEC_vwr_min20_wFilter$sample, 20)
#gDNAscreen2_RPTEC_vwr_min20_wFilter$sample <- substr(gDNAscreen2_RPTEC_vwr_min20_wFilter$sample,1,nchar(gDNAscreen2_RPTEC_vwr_min20_wFilter$sample)-4)
#gDNAscreen2_RPTEC_vwr_min20_wFilter <- subset(gDNAscreen2_RPTEC_vwr_min20_wFilter, sample %in% fwalist_112822_RPTECscreen2$fwa)
#gDNAscreen2_RPTEC_vwr_min20_wFilter$copy.number2 <- gDNAscreen2_RPTEC_vwr_min20_wFilter$copy.number
#colnames(gDNAscreen2_RPTEC_vwr_min20_wFilter)[11] <- "Sample_ID"
#colnames(gDNAscreen2_RPTEC_vwr_min20_wFilter)[1] <- "Chromosome"
#colnames(gDNAscreen2_RPTEC_vwr_min20_wFilter)[2] <- "Start"
#colnames(gDNAscreen2_RPTEC_vwr_min20_wFilter)[3] <- "End"
#colnames(gDNAscreen2_RPTEC_vwr_min20_wFilter)[4] <- "length"
#gDNAscreen2_RPTEC_vwr_min20_wFilter2 <- gDNAscreen2_RPTEC_vwr_min20_wFilter
#gDNAscreen2_RPTEC_vwr_min20_wFilter2 <- gDNAscreen2_RPTEC_vwr_min20_wFilter2[,c(1,11,12)]
#gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast <- dcast(gDNAscreen2_RPTEC_vwr_min20_wFilter2, Sample_ID~Chromosome, fun.aggregate = mean)
#rownames(gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast) <- gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast$Sample_ID
#gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast <- gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast[,-1]
#gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast[is.na(gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast)] <- 0
#gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast_hc <- hclust(d=dist(gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast))
#fwa <- gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast_hc$labels[gDNAscreen2_RPTEC_vwr_min20_wFilter2_dcast_hc$order]
#fwa <- data.frame(fwa)
#fwa$order <- rownames(fwa)
#colnames(fwa)[1] <- "Sample_ID"
#gDNAscreen2_RPTEC_vwr_min20_wFilter <- merge(gDNAscreen2_RPTEC_vwr_min20_wFilter, fwa, by = "Sample_ID", all= TRUE)
#gDNAscreen2_RPTEC_vwr_min20_wFilter$order <- as.numeric(as.character(gDNAscreen2_RPTEC_vwr_min20_wFilter$order))
#gDNAscreen2_RPTEC_vwr_min20_wFilter <- gDNAscreen2_RPTEC_vwr_min20_wFilter[order(gDNAscreen2_RPTEC_vwr_min20_wFilter$order),]

dev.off()
pdf(file = "RPTEC_screen2_cell_profiles_hclust.pdf", width = 8, height = 8)
ggplot(data = gDNAscreen2_RPTEC_vwr_min20_wFilter, aes(x = (start.genome+end.genome)/2, y = copy.number)) + geom_point(size = 1.2, alpha= 0.3) + theme_classic() + geom_segment(data=gDNAscreen2_RPTEC_vwr_min20_wFilter, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4, color='copy.number'), size=10) + scale_x_continuous(breaks=chromdata_hg19$chromlength/2+df.chroms[-24,1], labels=chromdata_hg19$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=1) + ylim(4,4) + geom_segment(data=whiteout, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4), color="white", size=10) + theme(legend.position="none") + theme(plot.margin=grid::unit(c(4.7807,0,4.7807,0), "cm")) + ylab("") + xlab("") + scale_colour_gradient2(limits = c(0.96, 3.05), low="blue", mid="gray90", high="red", na.value="firebrick", midpoint=2) + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.margin = unit(0, "lines")) +facet_wrap(~order, nrow=60+54) + theme(strip.text.x = element_blank())
dev.off()



#EXT DATA FIGURE 1K
#HMEC screen 2 dips

#file upload and pre-processing

#fwalist_HMEC_screen2_dips <- list.files("~/Documents/gDNA_HMEC_screen2_ALL_dips")
#fwalist_HMEC_screen2_dips <- data.frame(fwalist_HMEC_screen2_dips)
#colnames(fwalist_HMEC_screen2_dips)[1] <- "fwa"
#fwalist_HMEC_screen2_dips$fwa <- substring(fwalist_HMEC_screen2_dips$fwa, 28)
#fwalist_HMEC_screen2_dips$fwa <- substr(fwalist_HMEC_screen2_dips$fwa,1,nchar(fwalist_HMEC_screen2_dips$fwa)-12)
#gDNAscreen2_combo_dips_vwr_min20_wFilter <- Merging_CNA_bin_files("~/Documents/gDNA_HMEC_screen2_ALL_dips_good_segments")
#gDNAscreen2_combo_dips_vwr_min20_wFilter$sample<- substring(gDNAscreen2_combo_dips_vwr_min20_wFilter$sample, 20)
#gDNAscreen2_combo_dips_vwr_min20_wFilter$sample <- substr(gDNAscreen2_combo_dips_vwr_min20_wFilter$sample,1,nchar(gDNAscreen2_combo_dips_vwr_min20_wFilter$sample)-4)
#gDNAscreen2_combo_dips_vwr_min20_wFilter <- subset(gDNAscreen2_combo_dips_vwr_min20_wFilter, sample %in% fwalist_HMEC_screen2_dips$fwa)
#gDNAscreen2_combo_dips_vwr_min20_wFilter$copy.number2 <- gDNAscreen2_combo_dips_vwr_min20_wFilter$copy.number
#colnames(gDNAscreen2_combo_dips_vwr_min20_wFilter)[11] <- "Sample_ID"
#colnames(gDNAscreen2_combo_dips_vwr_min20_wFilter)[1] <- "Chromosome"
#colnames(gDNAscreen2_combo_dips_vwr_min20_wFilter)[2] <- "Start"
#colnames(gDNAscreen2_combo_dips_vwr_min20_wFilter)[3] <- "End"
#colnames(gDNAscreen2_combo_dips_vwr_min20_wFilter)[4] <- "length"
#gDNAscreen2_combo_dips_vwr_min20_wFilter2 <- gDNAscreen2_combo_dips_vwr_min20_wFilter
#gDNAscreen2_combo_dips_vwr_min20_wFilter2 <- gDNAscreen2_combo_dips_vwr_min20_wFilter2[,c(1,11,12)]
#gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast <- dcast(gDNAscreen2_combo_dips_vwr_min20_wFilter2, Sample_ID~Chromosome, fun.aggregate = mean)
#rownames(gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast) <- gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast$Sample_ID
#gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast <- gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast[,-1]
#gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast[is.na(gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast)] <- 0
#gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast_hc <- hclust(d=dist(gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast))
#fwa <- gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast_hc$labels[gDNAscreen2_combo_dips_vwr_min20_wFilter2_dcast_hc$order]
#fwa <- data.frame(fwa)
#fwa$order <- rownames(fwa)
#colnames(fwa)[1] <- "Sample_ID"
#gDNAscreen2_combo_dips_vwr_min20_wFilter <- merge(gDNAscreen2_combo_dips_vwr_min20_wFilter, fwa, by = "Sample_ID", all= TRUE)
#gDNAscreen2_combo_dips_vwr_min20_wFilter$order <- as.numeric(as.character(gDNAscreen2_combo_dips_vwr_min20_wFilter$order))
#gDNAscreen2_combo_dips_vwr_min20_wFilter <- gDNAscreen2_combo_dips_vwr_min20_wFilter[order(gDNAscreen2_combo_dips_vwr_min20_wFilter$order),]

dev.off()
pdf(file = "HMEC_screen2_combo_cell_profiles_hclust_dips.pdf", width = 8, height = 8)
ggplot(data = gDNAscreen2_combo_dips_vwr_min20_wFilter, aes(x = (start.genome+end.genome)/2, y = copy.number)) + geom_point(size = 1.2, alpha= 0.3) + theme_classic() + geom_segment(data=gDNAscreen2_combo_dips_vwr_min20_wFilter, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4, color='copy.number'), size=10) + scale_x_continuous(breaks=chromdata_hg19$chromlength/2+df.chroms[-24,1], labels=chromdata_hg19$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=1) + ylim(4,4) + geom_segment(data=whiteout, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4), color="white", size=10) + theme(legend.position="none") + theme(plot.margin=grid::unit(c(6,0,6,0), "cm")) + ylab("") + xlab("") + scale_colour_gradient2(limits = c(0.96, 3.05), low="blue", mid="gray90", high="red", na.value="firebrick", midpoint=2) + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.margin = unit(0, "lines")) +facet_wrap(~order, nrow=66) + theme(strip.text.x = element_blank())
dev.off()



