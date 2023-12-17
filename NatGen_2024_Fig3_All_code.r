library(ggplot2)
library(ggpubr)
library(grid)
library(RColorBrewer)
library(circlize)


#Figure 3a

ggplot(summary_events_evolution_experiment_tally_wreps2, aes(x=ploidy, y=new_total_events_per_haploidGenome_per40PD, fill=ploidy))  +
    geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + stat_compare_means() + theme_bw() + ylim(0,2.5) + stat_summary(fun=mean, geom="point", color = "black", size=2, shape =12) 

ggplot(summary_events_evolution_experiment_tally_wreps2, aes(x=ploidy, y=new_arm_events_per_haploidGenome_per40PD, fill=ploidy))  +
    geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color="black") + stat_compare_means() + theme_bw() + ylim(0,2.5) + stat_summary(fun=mean, geom="point", color = "black", size=2, shape =12) 

ggplot(summary_events_evolution_experiment_tally_wreps2, aes(x=ploidy, y=new_WC_events_per_haploidGenome_per40PD, fill=ploidy))  +
    geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color="black") + stat_compare_means() + theme_bw() + ylim(0,2.5) + stat_summary(fun=mean, geom="point", color = "black", size=2, shape =12) 


#Figure 3b

#HMEC gains
ggplot(Fig2C_summary_final_state_CNAS_wreps_absolute_noCQ_nomulti, aes(x=arm_amps, group = wgd, fill = as.character(wgd) )) +
    geom_histogram(binwidth=2, aes(y=..density..*2), position="dodge") +
    labs(x="arm gains", y="HMEC lines", fill="wgd")  +
    theme_bw() + xlim(-1,22)

#HMEC losses
ggplot(Fig2C_summary_final_state_CNAS_wreps_absolute_noCQ_nomulti, aes(x=arm_dels, group = wgd, fill = as.character(wgd) )) +
    geom_histogram(binwidth=2, aes(y=..density..*2), position="dodge") +
    labs(x="arm losses", y="HMEC lines", fill="wgd")  +
    theme_bw() + xlim(-1,22)

#PCAWG gains
all_all <- data.frame(all_all)
colnames(all_all) <- c("id", "frac_altered", "chr","log2ratio","wgd")
all_all$log2ratio <- as.numeric(all_all$log2ratio)
all_all2_gains <- all_all
sumfwa <- matrix(nrow=length(unique(all_all2_gains$id)), ncol = 3)
for (i in 1:length(unique(all_all2_gains$id))) {
    fwasdx <- subset(all_all2_gains, id == unique(all_all2_gains$id)[i])
    sumfwa[i,1] <- fwasdx$id[1]
    sumfwa[i,2] <- ifelse (fwasdx$wgd[1] == "1", nrow(subset(fwasdx, log2ratio > log2(5/4) & frac_altered > 0.5)), nrow(subset(fwasdx, log2ratio > 0.57 & frac_altered > 0.5)))
    sumfwa[i,3] <- fwasdx$wgd[1]
}
PCAWG_gains_per_tumor <- data.frame(sumfwa)
colnames(PCAWG_gains_per_tumor) <- c("id", "gains", "wgd")
PCAWG_gains_per_tumor$gains <- as.numeric(as.character(PCAWG_gains_per_tumor$gains))
PCAWG_gains_per_tumor$gains <- ifelse(PCAWG_gains_per_tumor$gains > 20, 20, PCAWG_gains_per_tumor$gains)


ggplot(PCAWG_gains_per_tumor, aes(x=gains, fill = wgd)) +
    geom_histogram(binwidth=2, aes(y=..density..*2), position="dodge") +
    labs(x="arm gains", y="BRCA tumors", fill="group")  +
    theme_bw() +xlim(-1,22)


#PCAWG losses    
all_all <- data.frame(all_all)
colnames(all_all) <- c("id", "frac_altered", "chr","log2ratio","wgd")
all_all$log2ratio <- as.numeric(all_all$log2ratio)
all_all2_losses <- all_all
sumfwa <- matrix(nrow=length(unique(all_all2_losses$id)), ncol = 3)
for (i in 1:length(unique(all_all2_losses$id))) {
    fwasdx <- subset(all_all2_losses, id == unique(all_all2_losses$id)[i])
    sumfwa[i,1] <- fwasdx$id[1]
    sumfwa[i,2] <- ifelse (fwasdx$wgd[1] == "1", nrow(subset(fwasdx, log2ratio < log2(3/4) & frac_altered > 0.5)), nrow(subset(fwasdx, log2ratio < -0.85 & frac_altered > 0.5)))
    sumfwa[i,3] <- fwasdx$wgd[1]
}
PCAWG_losses_per_tumor <- data.frame(sumfwa)
colnames(PCAWG_losses_per_tumor) <- c("id", "losses", "wgd")
PCAWG_losses_per_tumor$losses <- as.numeric(as.character(PCAWG_losses_per_tumor$losses))
PCAWG_losses_per_tumor$losses <- ifelse(PCAWG_losses_per_tumor$losses > 20, 20, PCAWG_losses_per_tumor$losses)


ggplot(PCAWG_losses_per_tumor, aes(x=losses, fill = wgd)) +
    geom_histogram(binwidth=2, aes(y=..density..*2), position="dodge") +
    labs(x="arm losses", y="BRCA tumors", fill="group")  +
    theme_bw() +xlim(-1,22) +ylim(0,0.4)


#Figure 3c-f
#please see separate code file "plotting_code.R" for deep WGS analysis and plotting


#Figure 3g
#please see additional repository https://github.com/emmavwatson/sparseHiC for more detail on HiC data processing, analysis, and plotting pipeline

Sample_8_and_1_balanced_8_21_tri_oe_test2 <- subset(Sample_8_and_1_balanced_8_21_tri_oe_test, ( (fwa == 0 | fwa > 0) & sample == "ctrl") | ( (fwa == 0 | fwa < 0) & sample == "ae_ev2") |  ( (chrom1=="8" & chrom2 == "21") & sample == "ae_ev2") | ( (chrom1=="21" & chrom2 == "8") & sample == "ctrl") )
plotsky <- HiC_oe_plot_func_3chrom(Sample_8_and_1_balanced_8_21_tri_oe_test2, chromdata3_grch37_2,8,21,22,1)
print(plotsky,vp=viewport(angle=-45))

#making the control
test_Sample_12_all_trans_3_x_3_x_3 <- obs_vs_exp_all_trans_simple5("~/Sample_12_dedup_all_by_all_trans", Sample_12_101819_balanced_trans_exp, 3,3,3)
test_Sample_12_all_trans_3_x_3_x_3$combo <- test_Sample_12_all_trans_3_x_3_x_3$gini_start1*test_Sample_12_all_trans_3_x_3_x_3$gini_start2*test_Sample_12_all_trans_3_x_3_x_3$max_oe*test_Sample_12_all_trans_3_x_3_x_3$gini_all*test_Sample_12_all_trans_3_x_3_x_3$frac_sig_oe
#making the test sample
test_Sample_17_all_trans_3_x_3_x_3 <- obs_vs_exp_all_trans_simple5("~/Sample_17_dedup_all_by_all_trans", Sample_17_101819_balanced_trans_exp, 3,3,3)
test_Sample_17_all_trans_3_x_3_x_3$combo <- test_Sample_17_all_trans_3_x_3_x_3$gini_start1*test_Sample_17_all_trans_3_x_3_x_3$gini_start2*test_Sample_17_all_trans_3_x_3_x_3$max_oe*test_Sample_17_all_trans_3_x_3_x_3$gini_all*test_Sample_17_all_trans_3_x_3_x_3$frac_sig_oe
test_Sample_17_all_trans_3_x_3_x_3 <- merge(test_Sample_17_all_trans_3_x_3_x_3[-1,], test_Sample_12_all_trans_3_x_3_x_3[-1,], by = "compartment")
test_Sample_17_all_trans_3_x_3_x_3 <- separate(data = test_Sample_17_all_trans_3_x_3_x_3, col = compartment, into = c('x1', 'chrom1', 'x2', 'chrom2'), sep = "_")
test_Sample_17_all_trans_3_x_3_x_3 <- na.omit(test_Sample_17_all_trans_3_x_3_x_3)
#plotting
ggplot(test_Sample_17_all_trans_3_x_3_x_3, aes(chrom1, chrom2, fill= combo.x/combo.y)) + geom_tile() + scale_fill_gradientn(colours=brewer.pal(6,"YlGnBu"))


#Figure 3h

circos.initializeWithIdeogram()
set.seed(123)
bed1 <- read.delim("~/bed1.txt")
bed2 <- read.delim("~/bed2.txt")
bed1 <- bed1[,-1]
bed2 <- bed2[,-1]
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), border = NA)



#additional functions required

HiC_oe_plot_func_2chrom3 <- function(sample_x_balanced_tri_oe,chromdata,chrom1,chrom2,mincount)
{
    c1_1 <- subset(chromdata, chr ==chrom1)$add
    c1_2 <- subset(chromdata, chr ==chrom1)$chromlength
    c2_1 <- subset(chromdata, chr ==chrom2)$chromlength
    gg88 <- ggplot(subset(sample_x_balanced_tri_oe, count > mincount), aes(x=genome_pos1, y = genome_pos2, color=log2(obs_vs_exp) )) + 
        geom_point(size=0.1,shape=18) + 
        theme_bw() + 
        scale_colour_gradientn(colours=c("mediumblue","white", "red3"), na.value="white", limits = c(-3,3)) + 
        scale_x_continuous(limits = c(c1_1,(c1_1+c1_2+c2_1)), expand = c(0.001,0.001)) + 
        scale_y_continuous(limits = c(c1_1,(c1_1+c1_2+c2_1)), expand = c(0.001,0.001))+ 
        theme(axis.text.y=element_text(angle=90, hjust=0.5)) + 
        geom_vline(xintercept =c(c1_1,(c1_1+c1_2)), size=0.5,color="maroon4" )+ 
        geom_vline(xintercept =c((c1_1+c1_2),(c1_1+c1_2)), size=0.5,color="maroon4" )+ 
        geom_hline(yintercept =c((c1_1+c1_2),(c1_1+c1_2)), size=0.5,color="maroon4" ) + 
        geom_hline(yintercept =c((c1_1+c1_2),(c1_1+c1_2)), size=0.5,color="maroon4" ) +
        geom_hline(yintercept =c((c1_1+c1_2+c2_1),(c1_1+c1_2+c2_1)), size=0.5,color="maroon4" ) +
        theme_void() + 
        theme(legend.position = "none", plot.margin = unit(c(2, 2, 2, 2), "cm"))
    return(gg88)
}


HiC_oe_plot_func_3chrom <- function(sample_x_balanced_tri_oe,chromdata,chrom1,chrom2,chrom3,mincount)
{
    c1_1 <- subset(chromdata, chr ==chrom1)$add
    c1_2 <- subset(chromdata, chr ==chrom1)$chromlength
    c2_1 <- subset(chromdata, chr ==chrom2)$chromlength
    c3_1 <- subset(chromdata, chr ==chrom3)$chromlength
    gg88 <- ggplot(subset(sample_x_balanced_tri_oe, count > mincount), aes(x=genome_pos1, y = genome_pos2, color=log2_oe )) + 
        geom_point(size=0.1,shape=18) + 
        theme_bw() + 
        scale_colour_gradientn(colours=c("mediumblue","white", "red3"), na.value="white", limits = c(-3,3)) + 
        scale_x_continuous(limits = c(c1_1,(c1_1+c1_2+c2_1+c3_1)), expand = c(0.001,0.001)) + 
        scale_y_continuous(limits = c(c1_1,(c1_1+c1_2+c2_1+c3_1)), expand = c(0.001,0.001))+ 
        theme(axis.text.y=element_text(angle=90, hjust=0.5)) + 
        geom_vline(xintercept =c(c1_1,(c1_1+c1_2)), size=0.5,color="maroon4" )+ 
        geom_vline(xintercept =c((c1_1+c1_2),(c1_1+c1_2)), size=0.5,color="maroon4" )+ 
        geom_hline(yintercept =c((c1_1+c1_2), (c1_1+c1_2+c2_1+c3_1)), size=0.5,color="maroon4" )+ 
        geom_hline(yintercept =c((c1_1+c1_2),(c1_1+c1_2)), size=0.5,color="maroon4" ) +
        geom_hline(yintercept =c((c1_1+c1_2+c2_1),(c1_1+c1_2+c2_1+c3_1)), size=0.5,color="maroon4" ) +
        geom_hline(yintercept =c(c1_1,(c1_1+c1_2+c2_1+c3_1)), size=0.5,color="maroon4" ) +
        theme_void() + 
        theme(legend.position = "none", plot.margin = unit(c(2, 2, 2, 2), "cm"))
    return(gg88)
}


obs_vs_exp_all_trans_simple5 <- function(working_dir, expected_trans, min_oe, sig_oe, xx) {
    setwd(working_dir)
    rm(dataset1)
    file_list <- list.files()
    dataset1 <- matrix(ncol=6, nrow=length(file_list))
    colnames(dataset1) <- c("compartment", "gini_start1", "gini_start2", "max_oe", "gini_all", "frac_sig_oe")
    for (i in 1:length(file_list)){
        trans_a_b <- read.delim(file_list[i])
        Sample_11_101819_trans_exp1 <- expected_trans
        Sample_11_101819_trans_exp2 <- expected_trans
        Sample_11_101819_trans_exp1$diagonal <- paste("c", Sample_11_101819_trans_exp1$chrom2, "c", Sample_11_101819_trans_exp1$chrom1, sep='_')
        Sample_11_101819_trans_exp2$diagonal <- paste("c", Sample_11_101819_trans_exp2$chrom1, "c", Sample_11_101819_trans_exp2$chrom2, sep='_')
        Sample_11_101819_trans_exp <- rbind(Sample_11_101819_trans_exp1, Sample_11_101819_trans_exp2)
        Sample_11_101819_trans_exp_smol <- Sample_11_101819_trans_exp[,c(3,4,6,8)]
        out_Sample_11_101819_4_8_tri_tester1 <- trans_a_b
        out_Sample_11_101819_4_8_tri_tester1$diagonal <- paste("c", out_Sample_11_101819_4_8_tri_tester1$chrom1, "c", out_Sample_11_101819_4_8_tri_tester1$chrom2, sep ='_')
        out_Sample_11_101819_4_8_tri_tester12 <- merge(out_Sample_11_101819_4_8_tri_tester1, Sample_11_101819_trans_exp, by="diagonal")
        out_Sample_11_101819_4_8_tri_tester12$obs_vs_exp <- out_Sample_11_101819_4_8_tri_tester12$count/out_Sample_11_101819_4_8_tri_tester12$count.avg
        abc_start1 <- data.frame(table(subset(out_Sample_11_101819_4_8_tri_tester12, obs_vs_exp > min_oe)$start1))
        abc_start2 <- data.frame(table(subset(out_Sample_11_101819_4_8_tri_tester12, obs_vs_exp > min_oe)$start2))
        dataset1[i,1] <- out_Sample_11_101819_4_8_tri_tester12[1,1]
        dataset1[i,2] <- ineq(abc_start1$Freq, parameter = NULL, type = c("Gini"), na.rm = TRUE)
        dataset1[i,3] <- ineq(abc_start2$Freq, parameter = NULL, type = c("Gini"), na.rm = TRUE)
        dataset1[i,4] <- mean(tail(sort(out_Sample_11_101819_4_8_tri_tester12$obs_vs_exp),10))
        dataset1[i,5] <- ineq(subset(out_Sample_11_101819_4_8_tri_tester12, obs_vs_exp > xx)$obs_vs_exp, parameter = NULL, type = c("Gini"), na.rm = TRUE)
        dataset1[i,6] <- nrow(subset(out_Sample_11_101819_4_8_tri_tester12, obs_vs_exp > sig_oe))/nrow(subset(out_Sample_11_101819_4_8_tri_tester12, obs_vs_exp > 0))
    }
    dataset1 <- data.frame(dataset1)
    dataset1$gini_start1 <- as.numeric(as.character(dataset1$gini_start1))
    dataset1$gini_start2 <- as.numeric(as.character(dataset1$gini_start2))
    dataset1$max_oe <- as.numeric(as.character(dataset1$max_oe))
    dataset1$gini_all <- as.numeric(as.character(dataset1$gini_all))
    dataset1$frac_sig_oe <- as.numeric(as.character(dataset1$frac_sig_oe))
    dataset1$combo <- dataset1$gini_start1*dataset1$gini_start2* dataset1$max_oe*dataset1$gini_all*dataset1$frac_sig_oe
    return(dataset1)
}

