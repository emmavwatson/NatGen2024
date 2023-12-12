# FIGURE 1A
#drawing in Illustrator




#FIGURE 1B
#HMEC screen #1 dips and tets


plots_ALL_FIG1_HMECdips <- heatmapGenomewide9b(ALL_FIG1_HMEC_dips2, ALL_FIG1_HMEC_dips_files2, namessegments_080717_plate1_new5, 10000000, 2, 0.2, 0.2, 0.98, 3.05, "euclidean", "complete", 22.5)

plots_ALL_FIG1_HMECtets <- heatmapGenomewide9b(returned_models_ALL_FIG1_HMEC_tets_1e5, ALL_FIG1_HMEC_tets_files, namessegments_080717_plate1_new5, 10000000, 4, 0.2, 0.2, 1.98, 6.05, "euclidean", "complete", 9)



#FIGURE 1C
#RPTEC screen #1 dips

plots_ALL_FIG1_RPTEC <- heatmapGenomewide9b(models_ALL_RPTEC3, files_ALL_RPTEC,  DNA_shredding_021020_RPTEC_clones, 10000000, 2, 0.2, 0.2, 0.98, 3.05, "euclidean", "complete", 27)


# FIGURE 1D
#WC_CNA_freq_TCGA_puritycorrection_OV used for the Ovarian TCGA sample set which has no purity data available, all other datasets use WC_CNA_freq_TCGA_puritycorrection

#example run for LUAD; repeat for all tumor types except OV
# LUAD_WC_CNA_summaryNEW_0.75 <- WC_CNA_freq_TCGA_puritycorrection(LUAD_CNA, LUAD_purity, "LUAD_CNAs_0.5_-0.41_0.32", 0.75, -0.41, 0.32, chromdata_peri_cent2)
#example run for OV
# OV_WC_CNA_summaryNEW_0.75 <- WC_CNA_freq_TCGA_puritycorrection_OV(OV_CNA, "OV_CNAs_0.5_-0.41_0.32", 0.75, -0.41, 0.32, chromdata_peri_cent2)

#merging all TCGA data
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- BRCA_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)]
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[2] <- "BRCA_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, OV_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[3] <- "OV_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, DLBC_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[4] <- "DLBC_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, READ_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[5] <- "READ_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, COAD_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[6] <- "COAD_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, UCEC_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[7] <- "UCEC_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, HNSC_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[8] <- "HNSC_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, LUSC_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[9] <- "LUSC_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, ESCA_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[10] <- "ESCA_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, BLCA_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[11] <- "BLCA_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, KIRP_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[12] <- "KIRP_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, KIRC_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[13] <- "KIRC_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, PRAD_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[14] <- "PRAD_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, GBM_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[15] <- "GBM_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, LGG_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[16] <- "LGG_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, LUAD_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[17] <- "LUAD_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, LIHC_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[18] <- "LIHC_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, PAAD_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[19] <- "PAAD_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, SKCM_WC_CNA_summaryNEW_0.75[[4]][,c(1,2)], by = "Chromosome")
# colnames(all_by_all_corr_tissues_HMEC_RPTEC_WC)[20] <- "SKCM_amp_freq"
# all_by_all_corr_tissues_HMEC_RPTEC_WC$Chromosome[23] <- "X"

# all_by_all_corr_tissues_HMEC_RPTEC_WC_FINAL <- merge(all_by_all_corr_tissues_HMEC_RPTEC_WC, Combined_HMEC_RPTEC_screen_data, by = "Chromosome")

#"all_by_all_corr_tissues_HMEC_RPTEC_WC_FINAL.csv"

library(ggrepel)
library(ggpmisc)
library(ggpubr)

ggplot(all_by_all_corr_tissues_HMEC_RPTEC_WC_FINAL, aes(x= BRCA_amp_freq, y=  HMEC_amp_freq_s1_s2_avg)) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(sqrt(r.squared)), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=all_by_all_corr_tissues_HMEC_RPTEC_WC_FINAL, aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)


# FIGURE 1E


library(ggrepel)
library(ggpmisc)
library(ggpubr)

ggplot(all_by_all_corr_tissues_HMEC_RPTEC_WC_FINAL, aes(x= KIRC_amp_freq, y=  RPTEC_amp_freq_s1_s2_avg)) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(sqrt(r.squared)), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=all_by_all_corr_tissues_HMEC_RPTEC_WC_FINAL, aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)




# FIGURE 1F
#correlations all-by-all and plotting
#"my_cor_WC.csv"

library(gplots)

pal <- colorpanel(1000, "white", "khaki", "red")
my_cor_WC <- cor(all_by_all_corr_tissues_HMEC_RPTEC_WC_FINAL[,c(2:22)])
my_cor_WC[my_cor_WC <0] <- 0
heatmap.2(my_cor_WC, scale='none', col =  pal, trace = "none")





# FIGURE 1G
# WGD calls from PCAWG, from Jake
# 'Tissue_specific_rates_Tet.csv'

ggplot(Tissue_specific_rates_Tet, aes(x = Tissue, y=percent_cohort, fill=Ploidy)) +
    geom_bar(position="stack", stat="identity") +  theme_bw() + scale_fill_manual(values = c("maroon", "goldenrod2")) + facet_grid(. ~ exp)



#functions required 

filterSegments <- function(segments, min.seg.width) {

if (is.null(segments)) {
return(NULL)
}
if (min.seg.width<=0) {
return(segments)
}

replace.index <- which(width(segments) < min.seg.width)
repl.segments <- segments[replace.index]
keep.index <- which(width(segments) >= min.seg.width)
keep.segments <- segments[keep.index]
nearest.index <- nearest(repl.segments, keep.segments)
na.mask <- is.na(nearest.index)
nearest.index <- nearest.index[!na.mask]
replace.index <- replace.index[!na.mask]
if (length(nearest.index)>0) {
nearest.keep.segments <- keep.segments[nearest.index]
segments$state[replace.index] <- nearest.keep.segments$state
segments$mstate[replace.index] <- nearest.keep.segments$mstate
segments$pstate[replace.index] <- nearest.keep.segments$pstate
}
segments.df <- as.data.frame(segments)
segments.df <- collapseBins(segments.df, column2collapseBy='state', columns2drop=c('width'))
segments.filtered <- as(segments.df, 'GRanges')
seqlevels(segments.filtered) <- seqlevels(segments) 
seqlengths(segments.filtered) <- seqlengths(segments)
return(segments.filtered)
}


heatmapGenomewide9b <- function (x, files, keepers, msw, w, maxdif, maxvar, low, high,methoddist,methodhclust, height1) {
    rm(dataset1)
    rm(temp_dataset1)
    hmms <- x
    listfiles <- data.frame()
    for(i1 in 1:length(hmms)){
        listfiles[i1, 1] <- basename(files[i1])
        listfiles[i1, 2] <- i1
        hmms[[i1]]$ID <- basename(files[i1])
    }
    colnames(listfiles) <- c("file_name", "order")
    listfiles <- merge(listfiles, keepers, by ="file_name")

    hmms <- list.filter(hmms, ID %in% listfiles$file_name)
    hmms1 <- hmms
    for(i in 1:length(hmms1)) {
        if (!exists("dataset1")){
            dataset1 <- data.frame(hmms1[[i]]$bins)
            dataset1$IRanges <- paste(dataset1$seqnames, dataset1$start, dataset1$end, sep="_")
            dataset1$ID <- i
        }
        if (exists("dataset1")){
            temp_dataset1 <- data.frame(hmms1[[i]]$bins)
            temp_dataset1$IRanges <- paste(temp_dataset1$seqnames, temp_dataset1$start, temp_dataset1$end, sep="_")
            temp_dataset1$ID <- i
            dataset1 <-rbind(dataset1, temp_dataset1)
            rm(temp_dataset1)
        }
    }
    hmm_combo_bins <- dataset1
    hmm_combo_bins <- hmm_combo_bins  %>% distinct()
    df.wide <- reshape2::dcast(hmm_combo_bins, ID ~ IRanges, value.var = "copy.number", factorsAsStrings = FALSE, fun.aggregate=sum, drop=FALSE)
    df.wide <- df.wide[ , colSums(is.na(df.wide)) == 0]
    hc <- stats::hclust(stats::dist(data.matrix(df.wide[-1]), method=methoddist), method=methodhclust)
    hmms <- hmms[hc$order]

    listfiles <- listfiles[(hc$order),]
    rm(df1)
    for(i1 in 1:length(hmms)) { 
        hmms[[i1]]$segments <- filterSegments(hmms[[i1]]$segments, msw)
        dfplot.seg <- as.data.frame(transCoord(hmms[[i1]]$segments))

        dfplot.seg$copy.number <- substr(dfplot.seg$state, 1, 1)
        dfplot.seg$copy.number <- as.numeric(dfplot.seg$copy.number)
        dfplot.seg$counts.CNV <- hmms[[i1]]$distributions[as.character(dfplot.seg$state),'mu']
        dfplot.seg$width <- dfplot.seg$end.genome - dfplot.seg$start.genome
        binsmodelgrTC <- transCoord(hmms[[i1]]$bins)
        binsmodelgrTC <- as.data.frame(binsmodelgrTC)
        counts <- as.data.frame(hmms[[i1]]$bincounts)
        binsmodelgrTC <- cbind(binsmodelgrTC, counts)
        colnames(binsmodelgrTC) <- make.unique(names(binsmodelgrTC))  
        dfplot.segX <- subcloneCorrection_after_filterSeg3(dfplot.seg, w, binsmodelgrTC, maxdif)
        colnames(dfplot.segX) <- make.unique(names(dfplot.segX))
        rm(dataset1)
        rm(temp_dataset1)
        rm(chr)
        chr <- split(dfplot.segX, dfplot.segX$seqnames)
        for (i in 1:length(chr)) {
            if (!exists("dataset1")){
                dataset1 <- data.frame(chr[[i]])
                if (nrow(dataset1) <3) {next}
                if (nrow(dataset1) >6) {
                    dataset1$copy.number <- mean(dataset1$copy.number) 
                    next}
                if(var(dataset1$copy.number) > maxvar) {next}
                dataset1$copy.number <- mean(dataset1$copy.number)
            }
            
            if (exists("dataset1")){
                temp_dataset1 <- data.frame(chr[[i]])
                if (nrow(temp_dataset1) <3) {
                    dataset1 <-rbind(dataset1, temp_dataset1)
                    next}
                if (nrow(temp_dataset1) >6) {
                    temp_dataset1$copy.number <- mean(temp_dataset1$copy.number)
                    dataset1 <-rbind(dataset1, temp_dataset1)
                    next}
                if(var(temp_dataset1$copy.number) > maxvar) {
                    dataset1 <-rbind(dataset1, temp_dataset1)
                    next}
                temp_dataset1$copy.number <- mean(temp_dataset1$copy.number)
                dataset1 <-rbind(dataset1, temp_dataset1)
                rm(temp_dataset1)
            }
        }
        dfplot.segX <- dataset1
        dfplot.segX$ID <- i1
        dfplot.segX$ID2 <- hmms[[i1]]$ID
        dfplot.segX$group <- subset(listfiles, file_name==hmms[[i1]]$ID)[1,ncol(listfiles)]
        if (!exists("df1")) {
            df1 <- dfplot.segX
        }
        if (exists("df1")) {
            temp_df1 <- dfplot.segX
            df1 <-rbind(df1, temp_df1)
            rm(temp_df1)
        }
    }

    pltlist <- list()
    ggplt <- ggplot(df1) + geom_linerange(aes_string(ymin = "start.genome", 
                                                     ymax = "end.genome", x = "ID", col = "copy.number"), size = 5) + scale_y_continuous(breaks = label.pos, 
                                                                                                                                         labels = names(label.pos), expand=c(0,0)) + scale_x_continuous(name = "ID", 
                                                                                                                                                                                                        breaks = 1:length(unique(df1$ID2)), labels = unique(df1$ID2), expand=c(0.01,0.01))
    ggplt <- ggplt + scale_colour_gradient2(limits = c(low, high), low="blue", mid="gray90", high="red", na.value="firebrick", midpoint=w)
    ggplt <- ggplt + theme(panel.background = element_blank(), 
                           axis.ticks.x = element_blank(), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), 
                           axis.line = element_blank(), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor.y = element_line(colour = "black"))
    ggplt <- ggplt + geom_segment(aes_string(x = 0.7, xend = length(hmms)+0.3, 
                                             y = "x", yend = "x"), data = df.chroms, col = "black")
    ggplt <- ggplt + coord_flip()
    width.heatmap <- 50
    height <- height1
    pltlist[["heatmap"]] <- ggplt
    widths["heatmap"] <- width.heatmap

    p2 <- ggplot(df1, aes(x=ID,y=1,fill=group))+geom_tile() + scale_x_continuous(name = "ID", 
                                                                                 breaks = 1:length(unique(df1$ID2)), labels = unique(df1$ID2), expand=c(0.0065,0.0065))
    p2 <-  p2 + theme(panel.background = element_blank(), axis.ticks = element_blank(), 
                      axis.text = element_blank(), 
                      axis.line = element_blank(), axis.title = element_blank())
    p2 <- p2 + coord_flip()
    pltlist[["dendro1"]] <- p2
    width.dendro1 <- 1.3
    widths["dendro1"] <- width.dendro1

    dhc <- stats::as.dendrogram(hc)
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    ggdndr <- ggplot(ddata$segments) + geom_segment(aes_string(x = "x", 
                                                               xend = "xend", y = "y", yend = "yend")) + scale_y_reverse() + scale_x_continuous(expand=c(0.015,0.015))
    ggdndr <- ggdndr + coord_flip()
    ggdndr <- ggdndr + theme(panel.background = element_blank(), 
                             axis.ticks = element_blank(), axis.text = element_blank(), 
                             axis.line = element_blank(), axis.title = element_blank(), plot.margin = unit(c(1, 0, 1, 0), "cm"))
    width.dendro2 <- 5
    pltlist[["dendro2"]] <- ggdndr
    widths["dendro2"] <- width.dendro2
    cowplt <- cowplot::plot_grid(plotlist = rev(pltlist), align = "h", 
                                 ncol = length(pltlist), rel_widths = rev(widths))
    ggsave("plottest.pdf", cowplt, width = (width.heatmap+width.dendro1+width.dendro2), height = height, units = "cm", limitsize = FALSE)
    return(df1)
}


subcloneCorrection_after_filterSeg3 <- function(fwa,w,binsmodelgrTC,maxdif) {  
    dfplot.segX <- fwa
    w <- w
    for (j in 1:nrow(dfplot.segX)) {
        dfplot.segX$avg.countspbin[j] <- sum(subset(binsmodelgrTC, start.genome >= dfplot.segX$start.genome[j] & end.genome <= dfplot.segX$end.genome[j])$counts)/length(subset(binsmodelgrTC, start.genome >= dfplot.segX$start.genome[j] & end.genome <= dfplot.segX$end.genome[j])$counts)
    }
    dfplot.segX$copy.number <- ifelse( ((dfplot.segX$avg.countspbin > (1+maxdif)*dfplot.segX$counts.CNV) | (dfplot.segX$avg.countspbin < (1-maxdif)*dfplot.segX$counts.CNV)), dfplot.segX$copy.number*(dfplot.segX$avg.countspbin/dfplot.segX$counts.CNV), dfplot.segX$copy.number )
    return(dfplot.segX)
}


transCoord <- function(gr) {
     cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
     cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
     names(cum.seqlengths.0) <- seqlevels(gr)
     gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
     gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
     return(gr)
 }


WC_CNA_freq_TCGA_puritycorrection <- function(TCGA_CNA, TCGA_purity, y, z, w, g, chromdata) { 
    #prep files
    print("... preparing data...")
    TCGA_purity <- separate(data = TCGA_purity, col = SampleName, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_purity$id <- paste(TCGA_purity$TCGA, TCGA_purity$TSS, TCGA_purity$Participant, TCGA_purity$SampleVial, sep='-')
    TCGA_purity <- TCGA_purity[,-c(1:7)]
    TCGA_purity = TCGA_purity[!duplicated(TCGA_purity$id),]
    TCGA_CNA <- separate(data = TCGA_CNA, col = Sample, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_CNA$id <- paste(TCGA_CNA$TCGA, TCGA_CNA$TSS, TCGA_CNA$Participant, TCGA_CNA$SampleVial, sep='-')
    TCGA_CNA <- subset(TCGA_CNA, TCGA_CNA$SampleVial == "01A" | TCGA_CNA$SampleVial == "01B" )
    TCGA_CNA <- TCGA_CNA[,-c(1:7)]
    print("success")
    #correct CNAs based on purity
    print("... correcting CNAs based on purity ...")
    TCGA_CNA_purity <- merge(TCGA_CNA, TCGA_purity, by = "id")
    TCGA_CNA_purity$Segment_Mean <- log2(((2^(TCGA_CNA_purity$Segment_Mean) )/TCGA_CNA_purity$Purity_InfiniumPurify) - ((1-TCGA_CNA_purity$Purity_InfiniumPurify)/(TCGA_CNA_purity$Purity_InfiniumPurify)))
    TCGA_CNA_purity$Segment_Mean[is.na(TCGA_CNA_purity$Segment_Mean)] <- 0
    TCGA_CNA_purity <- TCGA_CNA_purity[,-7]
    print("success")
    #assign segments to arms, deal with centromere-spanning segments
    print("... assigning segments to arms ...")
    TCGA_CNA_purity <- merge(TCGA_CNA_purity, chromdata, by = "Chromosome")
    TCGA_CNA_purity$length <- TCGA_CNA_purity$End - TCGA_CNA_purity$Start
    TCGA_CNA_purity$length <- as.numeric(as.character(TCGA_CNA_purity$length))
    #TCGA_CNA_purity$chromlength <- TCGA_CNA_purity$chromlength - TCGA_CNA_purity$cent_length
    TCGA_CNA_purity$chromlength <- ifelse( (TCGA_CNA_purity$Chromosome == "13" | TCGA_CNA_purity$Chromosome == "14" | TCGA_CNA_purity$Chromosome == "15" | TCGA_CNA_purity$Chromosome == "21" | TCGA_CNA_purity$Chromosome == "22"), TCGA_CNA_purity$Qlength, TCGA_CNA_purity$chromlength)
    
    ids <- split(TCGA_CNA_purity, TCGA_CNA_purity$id)
    print("success")
    #arm-level loss summary
    print("... preparing chromosome loss summary ...")
    {  
        TCGA_CNA_purity_losses <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean < w)
        TCGA_CNA_purity_losses_chrX <- split(TCGA_CNA_purity_losses, TCGA_CNA_purity_losses$Chromosome)
        TCGA_CNA_purity_losses_chrX <- TCGA_CNA_purity_losses_chrX[sapply(TCGA_CNA_purity_losses_chrX, function(x) dim(x)[1]) > 0]
        TCGA_deletion_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_deletion_frequencies_ALL <- data.frame(TCGA_deletion_frequencies_ALL)
        colnames(TCGA_deletion_frequencies_ALL) <- c('Chromosome', 'del_frequencies', 'avg_segMean_deleted')
        Deletion_summary <- matrix(ncol=4, nrow=1)
        colnames(Deletion_summary) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'Chrom')
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
                xout[i, 3] <- sum(TCGA_CNA_purity_losses_chrJ_sample_patX$length)/TCGA_CNA_purity_losses_chrJ_sample_patX$chromlength[1]
                xout <- data.frame(xout)
                xout$X1 <- as.character(xout$X1)
                xout$X2 <- as.numeric(as.character(xout$X2))
                xout$X3 <- as.numeric(as.character(xout$X3))
            }
            xout_sig <- subset(xout, xout$X2 < w & xout$X3 > z)
            xout_sig$X4 <- as.character(xout_sig$X4)
            if(nrow(xout_sig) == 0) xout_sig[1,] <- c(0,0,0,TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1])
            xout_sig$X4 <- TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1]
            colnames(xout_sig) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'Chrom')
            #write.table(xout_sig, file = paste(y, "Chromosome_", TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1], "_deletion_samples", ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
            Deletion_summary <- rbind(Deletion_summary, xout_sig)
            TCGA_deletion_frequencies[1, 1] <- TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1]
            TCGA_deletion_frequencies[1, 2] <- nrow(subset(xout_sig, id != 0))/length(ids)
            TCGA_deletion_frequencies[1, 3] <- mean(xout_sig$avg_segMean_deleted)
            TCGA_deletion_frequencies <- data.frame(TCGA_deletion_frequencies)
            TCGA_deletion_frequencies_ALL <- rbind(TCGA_deletion_frequencies_ALL, TCGA_deletion_frequencies)
        }
        Deletion_summary <- Deletion_summary[-1,]
        TCGA_deletion_frequencies_ALL <- TCGA_deletion_frequencies_ALL[-1,]
        #write.table(TCGA_deletion_frequencies_ALL, file = paste(y, "_deletion_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        #write.table(Deletion_summary, file = paste(y, "_Deletion_summary", ".txt",sep=""), sep="\t", row.names=FALSE,  quote = FALSE)
        
    }
    print("success")
    #arm-level gain summary
    print("... preparing chromosome gain summary ...")
    {
        TCGA_CNA_purity_gains <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean > g)
        TCGA_CNA_purity_gains_chrX <- split(TCGA_CNA_purity_gains, TCGA_CNA_purity_gains$Chromosome)
        TCGA_CNA_purity_gains_chrX <- TCGA_CNA_purity_gains_chrX[sapply(TCGA_CNA_purity_gains_chrX, function(x) dim(x)[1]) > 0]
        
        TCGA_amplification_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_amplification_frequencies_ALL <- data.frame(TCGA_amplification_frequencies_ALL)
        colnames(TCGA_amplification_frequencies_ALL) <- c('Chromosome', 'amp_frequencies', 'avg_segMean_amplified')
        Amplification_summary <- matrix(ncol = 4, nrow = 1)
        colnames(Amplification_summary) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'Chrom')
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
                xout2[i, 3] <- sum(TCGA_CNA_purity_gains_chrJ_sample_patX$length)/TCGA_CNA_purity_gains_chrJ_sample_patX$chromlength[1]
                xout2 <- data.frame(xout2)
                xout2$X1 <- as.character(xout2$X1)
                xout2$X2 <- as.numeric(as.character(xout2$X2))
                xout2$X3 <- as.numeric(as.character(xout2$X3))
            }
            
            xout_sig2 <- subset(xout2, xout2$X2 > g & xout2$X3 > z)
            xout_sig2$X4 <- as.character(xout_sig2$X4)
            if(nrow(xout_sig2) == 0) xout_sig2[1,] <- c(0,0,0,TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1])
            xout_sig2$X4 <- TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1]
            colnames(xout_sig2) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'Chrom')
            #write.table(xout_sig2, file = paste(y, "Chromosome_", TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1], "_amplification_samples", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
            Amplification_summary <- rbind(Amplification_summary, xout_sig2)
            TCGA_amplification_frequencies[1, 1] <- TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1]
            TCGA_amplification_frequencies[1, 2] <- nrow(subset(xout_sig2, id != 0))/length(ids)
            TCGA_amplification_frequencies[1, 3] <- mean(xout_sig2$avg_segMean_amplified)
            TCGA_amplification_frequencies <- data.frame(TCGA_amplification_frequencies)
            TCGA_amplification_frequencies_ALL <- rbind(TCGA_amplification_frequencies_ALL, TCGA_amplification_frequencies)
        }
        Amplification_summary <- Amplification_summary[-1,]
        TCGA_amplification_frequencies_ALL <- TCGA_amplification_frequencies_ALL[-1,]
        #write.table(TCGA_amplification_frequencies_ALL, file = paste(y, "_amplification_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        #write.table(Amplification_summary, file = paste(y, "_Amplification_summary", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        
    }
    final_output <- list(Deletion_summary, TCGA_deletion_frequencies_ALL, Amplification_summary, TCGA_amplification_frequencies_ALL)
    final_output[[1]] <- final_output[[1]][final_output[[1]]$id != 0, ]
    final_output[[2]][is.na(final_output[[2]])] <- 0
    final_output[[3]] <- final_output[[3]][final_output[[3]]$id != 0, ]
    final_output[[4]][is.na(final_output[[4]])] <- 0
    return(final_output)
    print("success")
}


WC_CNA_freq_TCGA_puritycorrection_OV <- function(TCGA_CNA, y, z, w, g, chromdata) { 
    #prep files
    print("... preparing data...")
    
    TCGA_CNA <- separate(data = TCGA_CNA, col = Sample, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_CNA$id <- paste(TCGA_CNA$TCGA, TCGA_CNA$TSS, TCGA_CNA$Participant, TCGA_CNA$SampleVial, sep='-')
    TCGA_CNA <- subset(TCGA_CNA, TCGA_CNA$SampleVial == "01A" | TCGA_CNA$SampleVial == "01B" )
    TCGA_CNA <- TCGA_CNA[,-c(1:7)]
    print("success")
    #correct CNAs based on purity
    print("... correcting CNAs based on purity ...")
    TCGA_CNA_purity <- TCGA_CNA
    #TCGA_CNA_purity$Segment_Mean <- log2(((2^(TCGA_CNA_purity$Segment_Mean) )/TCGA_CNA_purity$Purity_InfiniumPurify) - ((1-TCGA_CNA_purity$Purity_InfiniumPurify)/(TCGA_CNA_purity$Purity_InfiniumPurify)))
    TCGA_CNA_purity$Segment_Mean[is.na(TCGA_CNA_purity$Segment_Mean)] <- 0
    TCGA_CNA_purity <- TCGA_CNA_purity[,-7]
    print("success")
    #assign segments to arms, deal with centromere-spanning segments
    print("... assigning segments to arms ...")
    TCGA_CNA_purity <- merge(TCGA_CNA_purity, chromdata, by = "Chromosome")
    TCGA_CNA_purity$length <- TCGA_CNA_purity$End - TCGA_CNA_purity$Start
    TCGA_CNA_purity$length <- as.numeric(as.character(TCGA_CNA_purity$length))
    #TCGA_CNA_purity$chromlength <- TCGA_CNA_purity$chromlength - TCGA_CNA_purity$cent_length
    TCGA_CNA_purity$chromlength <- ifelse( (TCGA_CNA_purity$Chromosome == "13" | TCGA_CNA_purity$Chromosome == "14" | TCGA_CNA_purity$Chromosome == "15" | TCGA_CNA_purity$Chromosome == "21" | TCGA_CNA_purity$Chromosome == "22"), TCGA_CNA_purity$Qlength, TCGA_CNA_purity$chromlength)
    
    ids <- split(TCGA_CNA_purity, TCGA_CNA_purity$id)
    print("success")
    #arm-level loss summary
    print("... preparing chromosome loss summary ...")
    {  
        TCGA_CNA_purity_losses <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean < w)
        TCGA_CNA_purity_losses_chrX <- split(TCGA_CNA_purity_losses, TCGA_CNA_purity_losses$Chromosome)
        TCGA_CNA_purity_losses_chrX <- TCGA_CNA_purity_losses_chrX[sapply(TCGA_CNA_purity_losses_chrX, function(x) dim(x)[1]) > 0]
        TCGA_deletion_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_deletion_frequencies_ALL <- data.frame(TCGA_deletion_frequencies_ALL)
        colnames(TCGA_deletion_frequencies_ALL) <- c('Chromosome', 'del_frequencies', 'avg_segMean_deleted')
        Deletion_summary <- matrix(ncol=4, nrow=1)
        colnames(Deletion_summary) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'Chrom')
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
                xout[i, 3] <- sum(TCGA_CNA_purity_losses_chrJ_sample_patX$length)/TCGA_CNA_purity_losses_chrJ_sample_patX$chromlength[1]
                xout <- data.frame(xout)
                xout$X1 <- as.character(xout$X1)
                xout$X2 <- as.numeric(as.character(xout$X2))
                xout$X3 <- as.numeric(as.character(xout$X3))
            }
            xout_sig <- subset(xout, xout$X2 < w & xout$X3 > z)
            xout_sig$X4 <- as.character(xout_sig$X4)
            if(nrow(xout_sig) == 0) xout_sig[1,] <- c(0,0,0,TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1])
            xout_sig$X4 <- TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1]
            colnames(xout_sig) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'Chrom')
            #write.table(xout_sig, file = paste(y, "Chromosome_", TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1], "_deletion_samples", ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
            Deletion_summary <- rbind(Deletion_summary, xout_sig)
            TCGA_deletion_frequencies[1, 1] <- TCGA_CNA_purity_losses_chrJ_sample_patX$Chromosome[1]
            TCGA_deletion_frequencies[1, 2] <- nrow(subset(xout_sig, id != 0))/length(ids)
            TCGA_deletion_frequencies[1, 3] <- mean(xout_sig$avg_segMean_deleted)
            TCGA_deletion_frequencies <- data.frame(TCGA_deletion_frequencies)
            TCGA_deletion_frequencies_ALL <- rbind(TCGA_deletion_frequencies_ALL, TCGA_deletion_frequencies)
        }
        Deletion_summary <- Deletion_summary[-1,]
        TCGA_deletion_frequencies_ALL <- TCGA_deletion_frequencies_ALL[-1,]
        #write.table(TCGA_deletion_frequencies_ALL, file = paste(y, "_deletion_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        #write.table(Deletion_summary, file = paste(y, "_Deletion_summary", ".txt",sep=""), sep="\t", row.names=FALSE,  quote = FALSE)
        
    }
    print("success")
    #arm-level gain summary
    print("... preparing chromosome gain summary ...")
    {
        TCGA_CNA_purity_gains <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean > g)
        TCGA_CNA_purity_gains_chrX <- split(TCGA_CNA_purity_gains, TCGA_CNA_purity_gains$Chromosome)
        TCGA_CNA_purity_gains_chrX <- TCGA_CNA_purity_gains_chrX[sapply(TCGA_CNA_purity_gains_chrX, function(x) dim(x)[1]) > 0]
        
        TCGA_amplification_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_amplification_frequencies_ALL <- data.frame(TCGA_amplification_frequencies_ALL)
        colnames(TCGA_amplification_frequencies_ALL) <- c('Chromosome', 'amp_frequencies', 'avg_segMean_amplified')
        Amplification_summary <- matrix(ncol = 4, nrow = 1)
        colnames(Amplification_summary) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'Chrom')
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
                xout2[i, 3] <- sum(TCGA_CNA_purity_gains_chrJ_sample_patX$length)/TCGA_CNA_purity_gains_chrJ_sample_patX$chromlength[1]
                xout2 <- data.frame(xout2)
                xout2$X1 <- as.character(xout2$X1)
                xout2$X2 <- as.numeric(as.character(xout2$X2))
                xout2$X3 <- as.numeric(as.character(xout2$X3))
            }
            
            xout_sig2 <- subset(xout2, xout2$X2 > g & xout2$X3 > z)
            xout_sig2$X4 <- as.character(xout_sig2$X4)
            if(nrow(xout_sig2) == 0) xout_sig2[1,] <- c(0,0,0,TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1])
            xout_sig2$X4 <- TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1]
            colnames(xout_sig2) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'Chrom')
            #write.table(xout_sig2, file = paste(y, "Chromosome_", TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1], "_amplification_samples", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
            Amplification_summary <- rbind(Amplification_summary, xout_sig2)
            TCGA_amplification_frequencies[1, 1] <- TCGA_CNA_purity_gains_chrJ_sample_patX$Chromosome[1]
            TCGA_amplification_frequencies[1, 2] <- nrow(subset(xout_sig2, id != 0))/length(ids)
            TCGA_amplification_frequencies[1, 3] <- mean(xout_sig2$avg_segMean_amplified)
            TCGA_amplification_frequencies <- data.frame(TCGA_amplification_frequencies)
            TCGA_amplification_frequencies_ALL <- rbind(TCGA_amplification_frequencies_ALL, TCGA_amplification_frequencies)
        }
        Amplification_summary <- Amplification_summary[-1,]
        TCGA_amplification_frequencies_ALL <- TCGA_amplification_frequencies_ALL[-1,]
        #write.table(TCGA_amplification_frequencies_ALL, file = paste(y, "_amplification_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        #write.table(Amplification_summary, file = paste(y, "_Amplification_summary", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        
    }
    final_output <- list(Deletion_summary, TCGA_deletion_frequencies_ALL, Amplification_summary, TCGA_amplification_frequencies_ALL)
    final_output[[1]] <- final_output[[1]][final_output[[1]]$id != 0, ]
    final_output[[2]][is.na(final_output[[2]])] <- 0
    final_output[[3]] <- final_output[[3]][final_output[[3]]$id != 0, ]
    final_output[[4]][is.na(final_output[[4]])] <- 0
    return(final_output)
    print("success")
}
