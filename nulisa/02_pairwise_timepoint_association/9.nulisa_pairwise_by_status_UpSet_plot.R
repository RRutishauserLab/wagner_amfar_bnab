#install.packages("VennDiagram")
#install.packages("UpSetR")                        # Install VennDiagram package
library(UpSetR)
library(VennDiagram) 
library(dplyr)
library(tidyr)
#library(qvalue)

plot_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/nulisa/plots"
table_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/nulisa/tables/paired_t_test/by_status/"
data_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/data"
setwd(data_dir)


bl0_1w_table <-read.csv(paste0(table_dir,"paired_t_test_targets_BL0_1w_status_effect_noOut_PID_2887.csv"))
bl0_preATI_table <-read.csv(paste0(table_dir,"paired_t_test_targets_BL0_preATI_status_effect_noOut_PID_2887.csv"))
preATI_ATI1_table <-read.csv(paste0(table_dir,"paired_t_test_targets_preATI_ATI1_status_effect_noOut_PID_2887.csv"))
preATI_preR_table <-read.csv(paste0(table_dir,"paired_t_test_targets_preATI_preR_status_effect_noOut_PID_2887.csv"))
preR_postR1_table <-read.csv(paste0(table_dir,"paired_t_test_targets_preR_postR1_status_effect_noOut_PID_2887.csv"))

#sign
bl0_1w_table_p <- bl0_1w_table %>% dplyr::filter(intercept_p<=0.05)
bl0_preATI_table_p <- bl0_preATI_table %>% dplyr::filter(intercept_p<=0.05)
preATI_ATI1_table_p <- preATI_ATI1_table %>% dplyr::filter(intercept_p<=0.05)
preATI_preR_table_p <- preATI_preR_table %>% dplyr::filter(intercept_p<=0.05)
preR_postR1_table_p <- preR_postR1_table %>% dplyr::filter(intercept_p<=0.05)


#upset plot


#UpSet plot
#combining the three lists
bl0_1w <- data.frame( Target= bl0_1w_table_p$Target, bl0_1w=1)
bl0_preATI <- data.frame( Target= bl0_preATI_table_p$Target, bl0_preATI=1)
preATI_ATI1 <- data.frame( Target= preATI_ATI1_table_p$Target, preATI_ATI1=1)
preATI_preR<- data.frame(Target= preATI_preR_table_p$Target, preATI_preR=1)
preR_postR1<- data.frame(Target= preR_postR1_table_p$Target, preR_postR1=1)


set1_set2 <- full_join(bl0_1w,unique(bl0_preATI) )
set1_set2_3 <- full_join(set1_set2,unique(preATI_ATI1) )
set1_set2_3_4 <- full_join(set1_set2_3,unique(preATI_preR) )
set1_set2_3_4_5 <- full_join(set1_set2_3_4,unique(preR_postR1) )
set1_set2_3_4_5[is.na(set1_set2_3_4_5)] <- 0

pdf(paste0(plot_dir,"/UpSet/UpSet_plot_gene_protein_accession_pairwise_timepoints_status_effect_no2887.pdf"))
print(upset(set1_set2_3_4_5, nsets = 5))
dev.off()
write.csv(set1_set2_3_4_5, file=paste0(table_dir,"UpSet_plot_gene_protein_accession_pairwise_timepoints_status_effect_no2887.csv"), quote=F, row.names=F)


bl0_1w_table <-read.csv(paste0(table_dir,"paired_t_test_targets_BL0_1w_status_effect.csv"))
bl0_preATI_table <-read.csv(paste0(table_dir,"paired_t_test_targets_BL0_preATI_status_effect.csv"))
preATI_ATI1_table <-read.csv(paste0(table_dir,"paired_t_test_targets_preATI_ATI1_status_effect.csv"))
preATI_preR_table <-read.csv(paste0(table_dir,"paired_t_test_targets_preATI_preR_status_effect.csv"))
preR_postR1_table <-read.csv(paste0(table_dir,"paired_t_test_targets_preR_postR1_status_effect.csv"))

#sign
bl0_1w_table_p <- bl0_1w_table %>% dplyr::filter(intercept_p<=0.05)
bl0_preATI_table_p <- bl0_preATI_table %>% dplyr::filter(intercept_p<=0.05)
preATI_ATI1_table_p <- preATI_ATI1_table %>% dplyr::filter(intercept_p<=0.05)
preATI_preR_table_p <- preATI_preR_table %>% dplyr::filter(intercept_p<=0.05)
preR_postR1_table_p <- preR_postR1_table %>% dplyr::filter(intercept_p<=0.05)


#upset plot


#UpSet plot
#combining the three lists
bl0_1w <- data.frame( Target= bl0_1w_table_p$Target, bl0_1w=1)
bl0_preATI <- data.frame( Target= bl0_preATI_table_p$Target, bl0_preATI=1)
preATI_ATI1 <- data.frame( Target= preATI_ATI1_table_p$Target, preATI_ATI1=1)
preATI_preR<- data.frame(Target= preATI_preR_table_p$Target, preATI_preR=1)
preR_postR1<- data.frame(Target= preR_postR1_table_p$Target, preR_postR1=1)


set1_set2 <- full_join(bl0_1w,unique(bl0_preATI) )
set1_set2_3 <- full_join(set1_set2,unique(preATI_ATI1) )
set1_set2_3_4 <- full_join(set1_set2_3,unique(preATI_preR) )
set1_set2_3_4_5 <- full_join(set1_set2_3_4,unique(preR_postR1) )
set1_set2_3_4_5[is.na(set1_set2_3_4_5)] <- 0

pdf(paste0(plot_dir,"/UpSet/UpSet_plot_gene_protein_accession_pairwise_timepoints_status_effect.pdf"))
print(upset(set1_set2_3_4_5, nsets = 5))
dev.off()
write.csv(set1_set2_3_4_5, file=paste0(table_dir,"UpSet_plot_gene_protein_accession_pairwise_timepoints_status_effect.csv"), quote=F, row.names=F)



