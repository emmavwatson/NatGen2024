library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggrepel)



#HMEC misseg vs selection
ggplot(subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(x= HMEC_misseg_freq, y= HMEC_screen_1_and_2 )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)

#HMEC screen 1 vs screen 2
ggplot(subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(x= HMEC_amp_screen1, y= HMEC_amp_screen2 )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)

#RPTEC misseg vs selection
ggplot(subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(x= RPTEC_misseg_freq, y= RPTEC_screen_1_and_2  )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)


#RPTEC screen 1 vs screen 2

ggplot(subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(x= RPTEC_amp_screen1, y= RPTEC_amp_screen2 )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)


#Kops comparison, HMECs

ggplot(subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(x= RPE1_Kops_final/100, y= HMEC_misseg_freq )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)


#Kops comparison, RPTECs

ggplot(subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(x= RPE1_Kops_final/100, y= RPTEC_misseg_freq )) + geom_point(size=1.2) +geom_smooth(method="lm", se=F) + stat_fit_glance(method = 'lm', method.args = list(formula = x~y), geom = 'text', label.x = "middle", label.y = "top", aes(label = paste(after_stat(r.squared), after_stat(p.value), sep="\n" ),size = 2)) +theme_light() + geom_text_repel(data=subset(Fig_S2a_c, Chromosome != "15x" & Chromosome != "8x"), aes(label =  Chromosome), color = "black", size=3.2, max.overlaps=100)
