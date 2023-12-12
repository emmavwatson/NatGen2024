library(dplyr)
library(stringr)

# Pedigree plot (using source data for fig 3c as input)
df <- read.csv("[filename]", header=T, as.is=T, sep="\t")
catalog7 <- c("#CD8500", "#104E8B", "#6495ED", "#DC143C", "#000000", "#FF82AB", "#E3E3E3")

pdf("stackedbar_signatures_pedigree.pdf", height=9, width=9)
par(mar=c(5.1,12.1,4.1,2.1)) # bottom, left, top, right
yy <- barplot(as.matrix(t(df[,c(2,4,6,3,5,7,9)])), border = NA, col = catalog7, names = df$cl, horiz = T, las=1, space=1, xlim = c(0,2700))
dev.off()

# Box plot (using source data for fig 3e as input)
df <- read.csv("[filename]", header=T, as.is=T, sep="\t")
twocolors <- c("#1874CD", "#FF9912")

pdf("boxplot_HMEC_mutations.pdf", height=6, width=6)
par(mfrow = c(2, 3))
boxplot((df$num_snv/(df$num_div*df$ploidy)) ~ as.factor(df$group), frame = FALSE, ylim = c(0,12), outline = F, ylab = "Haploid Mutation Rate", las = 1, col = twocolors, names=c("diploid", "WGD"))
stripchart((df$num_snv/(df$num_div*df$ploidy)) ~ as.factor(df$group), pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_indel/(df$num_div*df$ploidy)) ~ as.factor(df$group), frame = FALSE, ylim = c(0,1.1), outline = F, ylab = "Haploid Mutation Rate", las = 1, col = twocolors, names=c("diploid", "WGD"))
stripchart((df$num_indel/(df$num_div*df$ploidy)) ~ as.factor(df$group), pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_sv/(df$num_div*df$ploidy)) ~ as.factor(df$group), frame = FALSE, ylim = c(0, 0.12), outline = F, ylab = "Haploid Mutation Rate", las = 1, col = twocolors, names=c("diploid", "WGD"))
stripchart((df$num_sv/(df$num_div*df$ploidy)) ~ as.factor(df$group), pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_snv/df$num_div) ~ as.factor(df$group), frame = FALSE, ylim = c(0,55), outline = F, ylab = "Unadjusted Mutation Rate", las = 1, col = twocolors, names=c("diploid", "WGD"))
stripchart((df$num_snv/df$num_div) ~ as.factor(df$group), pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_indel/df$num_div) ~ as.factor(df$group), frame = FALSE, ylim = c(0,4.5), outline = F, ylab = "Unadjusted Mutation Rate", las = 1, col = twocolors, names=c("diploid", "WGD"))
stripchart((df$num_indel/df$num_div) ~ as.factor(df$group), pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_sv/df$num_div) ~ as.factor(df$group), frame = FALSE, ylim = c(0, 0.5), outline = F, ylab = "Unadjusted Mutation Rate", las = 1, col = twocolors, names=c("diploid", "WGD"))
stripchart((df$num_sv/df$num_div) ~ as.factor(df$group), pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
dev.off()

t.test((df$num_snv/(df$num_div*df$ploidy)) ~ as.factor(df$group))
t.test((df$num_indel/(df$num_div*df$ploidy)) ~ as.factor(df$group))
t.test((df$num_sv/(df$num_div*df$ploidy)) ~ as.factor(df$group))
t.test((df$num_snv/(df$num_div)) ~ as.factor(df$group))
t.test((df$num_indel/(df$num_div)) ~ as.factor(df$group))
t.test((df$num_sv/(df$num_div)) ~ as.factor(df$group))

# Box plot (using source data for fig 3f as input)
df <- read.csv("[filename]", header=T, as.is=T, sep="\t")
twocolors <- c("#1874CD", "#FF9912")

pdf("boxplot_PCAWG_mutplots.pdf", height=6, width=6)
par(mfrow = c(2, 3))
boxplot((df$num_snv) ~ df$wgd, frame = FALSE, ylim = c(0, max(df$num_snv)), outline = F, ylab = "Number of substitutions", las = 1, col = twocolors, names = c("no WGD\n(n=105)", "WGD\n(n=103)"))
stripchart((df$num_snv) ~ df$wgd, pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_indel) ~ df$wgd, frame = FALSE, ylim = c(0, max(df$num_indel)), outline = F, ylab = "Number of indels", las = 1, col = twocolors, names = c("no WGD\n(n=105)", "WGD\n(n=103)"))
stripchart((df$num_indel) ~ df$wgd, pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_sv) ~ df$wgd, frame = FALSE, ylim = c(0, max(df$num_sv)), outline = F, ylab = "Number of rearrangements", las = 1, col = twocolors, names = c("no WGD\n(n=105)", "WGD\n(n=103)"))
stripchart((df$num_sv) ~ df$wgd, pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_snv/df$ploidy) ~ df$wgd, frame = FALSE, ylim = c(0, max(df$num_snv/df$ploidy)), outline = F, ylab = "Number of substitutions in haploid genome", las = 1, col = twocolors, names = c("no WGD\n(n=105)", "WGD\n(n=103)"))
stripchart((df$num_snv/df$ploidy) ~ df$wgd, pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_indel/df$ploidy) ~ df$wgd, frame = FALSE, ylim = c(0, max(df$num_indel/df$ploidy)), outline = F, ylab = "Number of indels in haploid genome", las = 1, col = twocolors, names = c("no WGD\n(n=105)", "WGD\n(n=103)"))
stripchart((df$num_indel/df$ploidy) ~ df$wgd, pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
boxplot((df$num_sv/df$ploidy) ~ df$wgd, frame = FALSE, ylim = c(0, max(df$num_sv/df$ploidy)), outline = F, ylab = "Number of rearrangements in haploid genome", las = 1, col = twocolors, names = c("no WGD\n(n=105)", "WGD\n(n=103)"))
stripchart((df$num_sv/df$ploidy) ~ df$wgd, pch = 19, col = rgb(0,0,0,.3), vertical=T, method="jitter", add = T)
dev.off()

t.test((df$num_snv) ~ df$wgd)
t.test((df$num_indel) ~ df$wgd)
t.test((df$num_sv) ~ df$wgd)
t.test((df$num_snv/df$ploidy) ~ df$wgd)
t.test((df$num_indel/df$ploidy) ~ df$wgd)
t.test((df$num_sv/df$ploidy) ~ df$wgd)

# Heatmap plot (using source data for extended data fig 6c)
library(ComplexHeatmap)
library(viridis)

## for chromosome 20
snp20 <- read.csv("[filename]", header=T, as.is=T, sep="\t")
rownames(snp20) = snp20[,1]
snp20 <- snp20[,c(2:ncol(snp20))]

pdf("Chromosome20_heatmap.pdf", height=4.5, width=5)
Heatmap(as.matrix(snp20), col=rev(viridis(25)), column_order = c("D02", "D04", "D03", "D05", "A01", "A02", "A03", "A04", "A05", "T01", "T02", "T04", "T03", "T05", "T06", "T07", "B01", "B02", "B03", "B04", "B05", "B06", "B07"), column_title = "Chromosome 20")
dev.off()

## for chromosome 1q
snp1q <- read.csv("[filename", header=T, as.is=T, sep="\t")
rownames(snp1q) = snp1q[,1]
snp1q <- snp1q[,c(2:ncol(snp1q))]

pdf("Chromosome1q_heatmap.pdf", height=4.5, width=5)
Heatmap(as.matrix(snp1q), col=rev(viridis(25)), column_order = c("D02", "D04", "D03", "D05", "A01", "A02", "A03", "A04", "A05", "T01", "T02", "T04", "T03", "T05", "T06", "T07", "B01", "B02", "B03", "B04", "B05", "B06", "B07"), column_title = "Chromosome 1q")
dev.off()

# Mutational spectra plot (using source data for extended data fig 6defg)
df <- read.csv("[filename]", header=T, as.is=T, sep="\t")

ord1 <- c("AC>AA","AC>AC","AC>AG","AC>AT","CC>AA","CC>AC","CC>AG","CC>AT","GC>AA","GC>AC","GC>AG","GC>AT","TC>AA","TC>AC","TC>AG","TC>AT",
          "AC>GA","AC>GC","AC>GG","AC>GT","CC>GA","CC>GC","CC>GG","CC>GT","GC>GA","GC>GC","GC>GG","GC>GT","TC>GA","TC>GC","TC>GG","TC>GT",
          "AC>TA","AC>TC","AC>TG","AC>TT","CC>TA","CC>TC","CC>TG","CC>TT","GC>TA","GC>TC","GC>TG","GC>TT","TC>TA","TC>TC","TC>TG","TC>TT",
          "AT>AA","AT>AC","AT>AG","AT>AT","CT>AA","CT>AC","CT>AG","CT>AT","GT>AA","GT>AC","GT>AG","GT>AT","TT>AA","TT>AC","TT>AG","TT>AT",
          "AT>CA","AT>CC","AT>CG","AT>CT","CT>CA","CT>CC","CT>CG","CT>CT","GT>CA","GT>CC","GT>CG","GT>CT","TT>CA","TT>CC","TT>CG","TT>CT",
          "AT>GA","AT>GC","AT>GG","AT>GT","CT>GA","CT>GC","CT>GG","CT>GT","GT>GA","GT>GC","GT>GG","GT>GT","TT>GA","TT>GC","TT>GG","TT>GT")

ord2 <- c("TG>TT","GG>TT","CG>TT","AG>TT","TG>TG","GG>TG","CG>TG","AG>TG","TG>TC","GG>TC","CG>TC","AG>TC","TG>TA","GG>TA","CG>TA","AG>TA",
          "TG>CT","GG>CT","CG>CT","AG>CT","TG>CG","GG>CG","CG>CG","AG>CG","TG>CC","GG>CC","CG>CC","AG>CC","TG>CA","GG>CA","CG>CA","AG>CA",
          "TG>AT","GG>AT","CG>AT","AG>AT","TG>AG","GG>AG","CG>AG","AG>AG","TG>AC","GG>AC","CG>AC","AG>AC","TG>AA","GG>AA","CG>AA","AG>AA",
          "TA>TT","GA>TT","CA>TT","AA>TT","TA>TG","GA>TG","CA>TG","AA>TG","TA>TC","GA>TC","CA>TC","AA>TC","TA>TA","GA>TA","CA>TA","AA>TA",
          "TA>GT","GA>GT","CA>GT","AA>GT","TA>GG","GA>GG","CA>GG","AA>GG","TA>GC","GA>GC","CA>GC","AA>GC","TA>GA","GA>GA","CA>GA","AA>GA",
          "TA>CT","GA>CT","CA>CT","AA>CT","TA>CG","GA>CG","CA>CG","AA>CG","TA>CC","GA>CC","CA>CC","AA>CC","TA>CA","GA>CA","CA>CA","AA>CA")

mutspec <- data.frame(ord1, ord2)
mutspec$ord1 <- as.character(mutspec$ord1)
mutspec$ord2 <- as.character(mutspec$ord2)

subsig <- df[df$group %in% c("A", "D"),]
pyrdf <- data.frame(table(subsig$hg19_spectra[subsig$cat_refbase == "pyrimidine"]))
colnames(pyrdf)[1] = "ord1"
pyrdf$ord1 <- as.character(pyrdf$ord1)
purdf <- data.frame(table(subsig$hg19_spectra[subsig$cat_refbase == "purine"]))
colnames(purdf)[1] = "ord2"
purdf$ord2 <- as.character(purdf$ord2)
tempdf <- left_join(mutspec, pyrdf, by = "ord1")
tempdf <- left_join(tempdf, purdf, by = "ord2")
tempdf$Freq.x[is.na(tempdf$Freq.x)] = 0
tempdf$Freq.y[is.na(tempdf$Freq.y)] = 0
tempdf$count = tempdf$Freq.x + tempdf$Freq.y
pdf(paste0("mutational_spectra_group_AD.pdf"), height=5, width=10)
plottemp <- barplot(tempdf$count[1:96], space=1, border=NA, col=sigcolors, las=2, xlab = "Trinucleotide Context", ylab= "Number of Mutations", main=paste0("Substitutions: ", sum(tempdf$count)))
axis(1, plottemp, labels=tempdf$ord1[1:96], cex.axis=0.5, las=2)
dev.off()

subsig <- df[df$group %in% c("B", "T"),]
pyrdf <- data.frame(table(subsig$hg19_spectra[subsig$cat_refbase == "pyrimidine"]))
colnames(pyrdf)[1] = "ord1"
pyrdf$ord1 <- as.character(pyrdf$ord1)
purdf <- data.frame(table(subsig$hg19_spectra[subsig$cat_refbase == "purine"]))
colnames(purdf)[1] = "ord2"
purdf$ord2 <- as.character(purdf$ord2)
tempdf <- left_join(mutspec, pyrdf, by = "ord1")
tempdf <- left_join(tempdf, purdf, by = "ord2")
tempdf$Freq.x[is.na(tempdf$Freq.x)] = 0
tempdf$Freq.y[is.na(tempdf$Freq.y)] = 0
tempdf$count = tempdf$Freq.x + tempdf$Freq.y
pdf(paste0("mutational_spectra_group_BT.pdf"), height=5, width=10)
plottemp <- barplot(tempdf$count[1:96], space=1, border=NA, col=sigcolors, las=2, xlab = "Trinucleotide Context", ylab= "Number of Mutations", main=paste0("Substitutions: ", sum(tempdf$count)))
axis(1, plottemp, labels=tempdf$ord1[1:96], cex.axis=0.5, las=2)
dev.off()

