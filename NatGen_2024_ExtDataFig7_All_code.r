library(reshape2)
library(reshape)
library(ggplot2)


#Ext Data Figure 7a
#please see BioProject: PRJNA634423 for deep gDNA-seq data corresponding to this sample

#Ext Data Figure 7b
#please see BioProject: PRJNA634423 for HiC data corresponding to these samples
#see separate repository https://github.com/emmavwatson/sparseHiC for all raw HiC data processing and plotting code, also please see Fig 3 code for additional guidance

ExtDataFigure7bb <- melt(ExtDataFigure7b, na.rm = FALSE, value.name = "name", id = 'Chromosome')
ExtDataFigure7bb$variable <- gsub("^.{0,1}", "", ExtDataFigure7bb$variable)
colnames(ExtDataFigure7bb) <- c("chromosome1", "chromosome2", "ES")
ExtDataFigure7bb$chromosome1 <- as.character(ExtDataFigure7bb$chromosome1)
ExtDataFigure7bb$chromosome2 <- as.character(ExtDataFigure7bb$chromosome2)
ExtDataFigure7bb$chromosome1 <- factor(ExtDataFigure7bb$chromosome1,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"))
ExtDataFigure7bb$chromosome2 <- factor(ExtDataFigure7bb$chromosome2,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"))
ggplot(ExtDataFigure7bb, aes(chromosome1, chromosome2, fill= ES)) + geom_tile() + scale_fill_gradientn(colours=brewer.pal(6,"YlGnBu"))

#Ext Data Figure 7c
#please see BioProject: PRJNA634423 for HiC data corresponding to these samples
#see separate repository https://github.com/emmavwatson/sparseHiC for all raw HiC data processing and plotting code, also please see Fig 3 code for additional guidance

ExtDataFigure7cc <- melt(ExtDataFigure7c, na.rm = FALSE, value.name = "name", id = 'Chromosome')
ExtDataFigure7cc$variable <- gsub("^.{0,1}", "", ExtDataFigure7cc$variable)
colnames(ExtDataFigure7cc) <- c("chromosome1", "chromosome2", "ES")
ExtDataFigure7cc$chromosome1 <- as.character(ExtDataFigure7cc$chromosome1)
ExtDataFigure7cc$chromosome2 <- as.character(ExtDataFigure7cc$chromosome2)
ExtDataFigure7cc$chromosome1 <- factor(ExtDataFigure7cc$chromosome1,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"))
ExtDataFigure7cc$chromosome2 <- factor(ExtDataFigure7cc$chromosome2,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"))
ggplot(ExtDataFigure7cc, aes(chromosome1, chromosome2, fill= ES)) + geom_tile() + scale_fill_gradientn(colours=brewer.pal(6,"YlGnBu"))

#Ext Data Figure 7d
#please see BioProject: PRJNA634423 for HiC data corresponding to these samples
#see separate repository https://github.com/emmavwatson/sparseHiC for all raw HiC data processing and plotting code, also please see Fig 3 code for additional guidance

ExtDataFigure7dd <- melt(ExtDataFigure7d, na.rm = FALSE, value.name = "name", id = 'Chromosome')
ExtDataFigure7dd$variable <- gsub("^.{0,1}", "", ExtDataFigure7dd$variable)
colnames(ExtDataFigure7dd) <- c("chromosome1", "chromosome2", "ES")
ExtDataFigure7dd$chromosome1 <- as.character(ExtDataFigure7dd$chromosome1)
ExtDataFigure7dd$chromosome2 <- as.character(ExtDataFigure7dd$chromosome2)
ExtDataFigure7dd$chromosome1 <- factor(ExtDataFigure7dd$chromosome1,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"))
ExtDataFigure7dd$chromosome2 <- factor(ExtDataFigure7dd$chromosome2,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"))
ggplot(ExtDataFigure7dd, aes(chromosome1, chromosome2, fill= ES)) + geom_tile() + scale_fill_gradientn(colours=brewer.pal(6,"YlGnBu"))

#Ext Data Figure 7e
#please see BioProject: PRJNA634423 for deep gDNA-seq data corresponding to this sample


#Ext Data Figure 7f-h
#microscopy images