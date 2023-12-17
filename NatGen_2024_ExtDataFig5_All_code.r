library(ggplot2)
library(ggpmisc)
library(ggpubr)


#FIGURE S5A-C

#please see BioProject: PRJNA634423 for gDNA-seq data corresponding to these samples
#see separate repository https://github.com/emmavwatson/CNAplot for all raw gDNA sequencing processing, copy number calling, and plotting code 


#FIGURE S5D

ggplot(Fig_S5d_data, aes(x=BRCA_amp_minus_delNEW0.5, y = HMEC_evo_amp_minus_del_freq)) +geom_point(alpha=0.7)  + geom_smooth(method="lm",se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 3)) +theme_bw() +geom_text_repel(data=subset(Fig_S5d_data, HMEC_evo_amp_minus_del_freq > 0.1 | HMEC_evo_amp_minus_del_freq < -0.07 | Chromosome == "8p" | Chromosome == "17p"), aes(label = Chromosome), color = "red", size=4) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)