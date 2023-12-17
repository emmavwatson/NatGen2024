library(ggplot2)

#Ext Data Figure 4a

cols <- c("WC_genes" = "magenta", "Centrans_genes" = "cyan", "Complex_genes" = "steelblue")
summary_cna_affected_genes_all_tissuesM2$tissue <- factor(summary_cna_affected_genes_all_tissuesM2$tissue,levels = c("LUSC", "OV", "ESCA", "LUAD", "SKCM", "BRCA", "BLCA", "STAD", "HNSC", "LIHC", "COAD", "KIRP", "GBM", "PAAD", "UCEC", "KIRC", "LGG", "PRAD"))
summary_cna_affected_genes_all_tissuesM2$variable <- factor(summary_cna_affected_genes_all_tissuesM2$variable,levels = c("WC_genes", "Centrans_genes", "Complex_genes"))
ggplot(summary_cna_affected_genes_all_tissuesM2, aes(x = tissue, y = value, fill = variable)) + geom_bar(stat = "identity") +theme_bw() +scale_fill_manual(values=cols) + theme(axis.text.x = element_text(angle = 45, hjust=1))
