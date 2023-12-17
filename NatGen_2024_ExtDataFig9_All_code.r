library(ggplot2)
library(pheatmap)
library(ggpmisc)
library(tidyverse)
library(ggpubr)
library(rstatix)

#Ext Data Figure 9a

ExtData_Fig9a_Top_data <- ExtData_Fig9a_Top_data[!duplicated(ExtData_Fig9a_Top_data$genes), ]
ggplot(data = ExtData_Fig9a_Top_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(ExtData_Fig9a_Top_data, Chromosome != "1"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(ExtData_Fig9a_Top_data, Chromosome == "1" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red3") +
    geom_point(data = subset(ExtData_Fig9a_Top_data, Chromosome == "1" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    theme_classic() +
    scale_color_manual(values = chromColors) +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("copy number") +
    xlab("chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(0.5),0, log2(1.5)))

ExtData_Fig9a_Top_data$copy <- ifelse(ExtData_Fig9a_Top_data$Chromosome == "1" & ExtData_Fig9a_Top_data$Karyotype_band %like% "q", "3", "2")
ggplot(ExtData_Fig9a_Top_data, aes(x=logFC, color=copy)) +
    geom_density() +geom_vline(xintercept = c(log2(2/2.05),log2(3/2.05))) + xlim(-2,2) +theme_bw()


ExtData_Fig9a_Bottom_data <- ExtData_Fig9a_Bottom_data[!duplicated(ExtData_Fig9a_Bottom_data$genes), ]
ggplot(data = ExtData_Fig9a_Bottom_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(ExtData_Fig9a_Bottom_data, Chromosome != "8"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(ExtData_Fig9a_Bottom_data, Chromosome == "8" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red3") +
    geom_point(data = subset(ExtData_Fig9a_Bottom_data, Chromosome == "8" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    theme_classic() +
    scale_color_manual(values = chromColors) +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("copy number") +
    xlab("chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(0.5),0, log2(1.5)))

ExtData_Fig9a_Bottom_data$copy <- ifelse(ExtData_Fig9a_Bottom_data$Chromosome == "8" & ExtData_Fig9a_Bottom_data$Karyotype_band %like% "q", "3", "2")
ggplot(ExtData_Fig9a_Bottom_data, aes(x=logFC, color=copy)) +
    geom_density() +geom_vline(xintercept = c(log2(2/2.05),log2(3.3/2.05))) + xlim(-2,2) +theme_bw()



#Ext Data Figure 9b

ExtData_Fig9b_Top_data <- ExtData_Fig9b_Top_data[!duplicated(ExtData_Fig9b_Top_data$genes), ]
ggplot(data = ExtData_Fig9b_Top_data, aes(x = gen_pos, y = logFC)) +
    geom_point(data = subset(ExtData_Fig9b_Top_data, Chromosome != "1"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(ExtData_Fig9b_Top_data, Chromosome == "1" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red3") +
    geom_point(data = subset(ExtData_Fig9b_Top_data, Chromosome == "1" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    theme_classic() +
    scale_color_manual(values = chromColors) +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("copy number") +
    xlab("chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(0.5),0, log2(1.5)))

ExtData_Fig9b_Top_data$copy <- ifelse(ExtData_Fig9b_Top_data$Chromosome == "1" & ExtData_Fig9b_Top_data$Karyotype_band %like% "q", "3", "2")
ggplot(ExtData_Fig9b_Top_data, aes(x=logFC, color=copy)) +
    geom_density() +geom_vline(xintercept = c(log2(2/2.05),log2(3/2.05))) + xlim(-2,2) +theme_bw()

ExtData_Fig9b_Bottom_data <- ExtData_Fig9b_Bottom_data[!duplicated(ExtData_Fig9b_Bottom_data$genes), ]
ggplot(data = ExtData_Fig9b_Bottom_data, aes(x = genpos, y = logFC)) +
    geom_point(data = subset(ExtData_Fig9b_Bottom_data, Chromosome != "8"), size = 1.2, alpha= 0.3, color = "grey60")  +
    geom_point(data = subset(ExtData_Fig9b_Bottom_data, Chromosome == "8" & Karyotype_band %like% "q"), size = 1.2, alpha= 0.3, color = "red3") +
    geom_point(data = subset(ExtData_Fig9b_Bottom_data, Chromosome == "8" & Karyotype_band %like% "p"), size = 1.2, alpha= 0.3, color = "grey60") +
    theme_classic() +
    scale_color_manual(values = chromColors) +
    scale_x_continuous(breaks=(chromdata3_grch37_3$add + (chromdata3_grch37_3$chromlength/2)), labels=chromdata3_grch37_3$Chromosome) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(-2, 2) + ylab("copy number") +
    xlab("chromosome") +
    theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) +
    theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=12), axis.title = element_text(size=13)) + geom_hline(yintercept = c(log2(0.5),0, log2(1.5)))

ExtData_Fig9b_Bottom_data$copy <- ifelse(ExtData_Fig9b_Bottom_data$Chromosome == "8" & ExtData_Fig9b_Bottom_data$Karyotype_band %like% "q", "3", "2")
ggplot(ExtData_Fig9b_Bottom_data, aes(x=logFC, color=copy)) +
    geom_density() +geom_vline(xintercept = c(log2(2/2.05),log2(3.3/2.05))) + xlim(-2,2) +theme_bw()

#Ext Data Figure 9c

my_palette3 <- myPalette(low = "darkblue", high = "red3", mid="white", k =12)
my_breaks <- c(-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5,3 )

pheatmap(ExtData_Fig9c_data[,-1], cluster_rows=TRUE, cluster_cols=TRUE, color = my_palette3, show_rownames=TRUE,  show_colnames=TRUE, border_col=NA, fontsize=6, breaks=my_breaks)

#Ext Data Figure 9d


`%!like%` = Negate(`%like%`)

stat.test <- subset(fwaDF, id_rep %!like% "K16" & exp != "GSI") %>%
    group_by(exp) %>%
    t_test( logFC_NotchUP~ chr1q) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test <- stat.test %>%
    add_xy_position(x = "exp", dodge = 0.8)

ExtData_Fig9d_data$chr1q <- factor(ExtData_Fig9d_data$chr1q,levels = c("WT", "gain_1q"))
ggplot(subset(ExtData_Fig9d_data, id_rep %!like% "K16" & exp != "GSI"), aes(x=exp,y=logFC_NotchUP, color=chr1q)) + geom_boxplot(width=.5) +theme_bw() + geom_point(position=position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1)) + stat_pvalue_manual(
    stat.test,  label = "p", tip.length = 0)


#Ext Data Figure 9e


ggplot(subset(ExtData_Fig9e_data, genes == "NCSTN"), aes(x= DNA_logFC, y= mRNA_logFC )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() +ylim(-2.2,2.2) + xlim(-2.2,2.2) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
ggplot(subset(ExtData_Fig9e_data, genes == "APH1A"), aes(x= DNA_logFC, y= mRNA_logFC )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() +ylim(-2.2,2.2) + xlim(-2.2,2.2) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
ggplot(subset(ExtData_Fig9e_data, genes == "PSEN2"), aes(x= DNA_logFC, y= mRNA_logFC )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() +ylim(-2.2,2.2) + xlim(-2.2,2.2) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
   
#Ext Data Figure 9f

ggplot(ExtData_Fig9f_data, aes(x=as.character(gain_1q), y=APH1A_RPKM)) + geom_beeswarm(alpha=0.2)  +theme_bw() +stat_compare_means(method = "t.test")
ggplot(ExtData_Fig9f_data, aes(x=as.character(gain_1q), y=NCSTN_RPKM)) + geom_beeswarm(alpha=0.2)  +theme_bw() +stat_compare_means(method = "t.test")
ggplot(ExtData_Fig9f_data, aes(x=as.character(gain_1q), y=PSEN2_RPKM)) + geom_beeswarm(alpha=0.2)  +theme_bw() +stat_compare_means(method = "t.test")

