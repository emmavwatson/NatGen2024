library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggrepel)

#all HMEC comparisons to tissues, individual plots

myplots2 <- vector('list', 15)
for (i in 2:16) {
    message(i)
    myplots2[[i]] <- local({
        i <- i
        p1 <- ggplot(FigS3_data, aes(x= FigS3_data[[i]], y= HMEC_amp_freq_s1_s2_avg  )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", mapping = aes(label = sprintf('r2 = %.3f \n pval = ~%.2g',after_stat(r.squared), after_stat(p.value)))) +theme_light() + geom_text_repel(data=FigS3_data, aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100) +theme(axis.title.y = element_blank(), axis.title.x = element_blank())
        print(p1)
    })
}

myplots2 <- myplots2[2:16]
dev.off()
pdf(file = "individual_lr_HMEC_vs_tumors.pdf", width = 3, height = 44)
plot_grid(myplots2[[1]], myplots2[[2]], myplots2[[3]],myplots2[[4]], myplots2[[5]],myplots2[[6]],myplots2[[7]],myplots2[[8]],myplots2[[9]],myplots2[[10]],myplots2[[11]],myplots2[[12]],myplots2[[13]],myplots2[[14]], myplots2[[15]], ncol = 1)
dev.off()


#all RPTEC comparisons to tissues, individual plots

myplots2 <- vector('list', 15)
for (i in 2:16) {
    message(i)
    myplots2[[i]] <- local({
        i <- i
        p1 <- ggplot(subset(FigS3_data, Chromosome != "15x" & Chromosome != "8x"), aes(x= FigS3_data[[i]], y= RPTEC_amp_freq_s1_s2_avg  )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", mapping = aes(label = sprintf('r2 = %.3f \n pval = ~%.2g',after_stat(r.squared), after_stat(p.value)))) +theme_light() + geom_text_repel(data=FigS3_data, aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100) +theme(axis.title.y = element_blank(), axis.title.x = element_blank())
        print(p1)
    })
}

myplots2 <- myplots2[2:16]
dev.off()
pdf(file = "individual_lr_RRPTEC_vs_tumors.pdf", width = 3, height = 44)
plot_grid(myplots2[[1]], myplots2[[2]], myplots2[[3]],myplots2[[4]], myplots2[[5]],myplots2[[6]],myplots2[[7]],myplots2[[8]],myplots2[[9]],myplots2[[10]],myplots2[[11]],myplots2[[12]],myplots2[[13]],myplots2[[14]], myplots2[[15]], ncol = 1)
dev.off()

#RPTEC vs HMEC

ggplot(FigS3_data, aes(x= HMEC_amp_freq_s1_s2_avg, y= RPTEC_amp_freq_s1_s2_avg   )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", mapping = aes(label = sprintf('r2 = %.3f \n pval = ~%.2g',after_stat(r.squared), after_stat(p.value)))) +theme_light() + geom_text_repel(data=FigS3_data, aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100) +theme(axis.title.y = element_blank(), axis.title.x = element_blank())