##Removes all objects from workspace
rm(list=ls())

#Install packages
library(tidyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(magrittr)

setwd("~/Gladstone Dropbox/Alan Xie/GB-JW-1447")
plot_dir <- "~/Gladstone Dropbox/Alan Xie/GB-JW-1447/results/proteomics/plots"
table_dir <- "~/Gladstone Dropbox/Alan Xie/GB-JW-1447/results/proteomics/tables"
data_dir <- "~/Gladstone Dropbox/Alan Xie/GB-JW-1447/data"

#uploading the proteome seq dataset
proteome <- read.csv(paste0(data_dir,"/imputed_log2_transformed_protein_intensities_v2.csv"), check.names=F)

proteome_long <- proteome %>% pivot_longer(
  cols = `3555696563_8092_ATI1`:`3555696339_8012_ATI5`, 
  names_to = "id_PID_timepoint",
  values_to = "log2_quantification") %>%
  as.data.frame()

#splitting the id_PID_timepoint in the column and create a new variable PID_timepoint
proteome_long <- proteome_long %>% 
  separate(id_PID_timepoint, c('id', 'PID', 'timepoint')) %>% 
  mutate(PID_timepoint= paste0(PID,"_", timepoint ))

unique(proteome_long$timepoint)

# Make proteome_long into data frame
# proteome_long <- data.frame(proteome_long)
# head(proteome_long)

#create a vector of timepoints of interest
t_points <- c("BL0", "BL1", "1w", "preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","ATI10",
              "ATI7",  "ATI11",    "ATI8" ,"ATI9", "preR", "postR1", "postR2",
              "postR3", "postR4")
t_points   

proteome_subset <- proteome_long %>% dplyr::filter(timepoint %in% t_points)
head(proteome_subset)

#Create a column of controllers / non controllers
non_controllers<- c("1292", "1408", "8106")
controllers <- proteome_subset$PID[!proteome_subset$PID %in% non_controllers]
no_rebound <- proteome_subset$PID[proteome_subset$PID %in% "2523"]
controllers_rebound <- controllers[ !controllers %in% no_rebound]

proteome_subset <- proteome_subset %>%
  mutate(status = case_when(PID %in% non_controllers ~ 'non_controllers',
                            PID %in% controllers ~ 'controllers')) %>%
  mutate(rebound = case_when( PID %in% no_rebound ~ 'no' , PID %in% c(controllers_rebound,non_controllers) ~ 'yes'
  ))

#target <-unique(proteome_subset$Genes)

# #Make new data frame that reorganizes as Gene + Protein Accession.
# proteome_subset_t <-  proteome_subset %>%
#   mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions )) %>%
#   dplyr::select(c(PID_timepoint, log2_quantification, Genes_protein_accessions)) %>%
#   pivot_wider(names_from = Genes_protein_accessions, values_from = log2_quantification) %>%
#   as.data.frame()
# 
# #Define variable target as unique Gene+Protein Accessions
# gene_protein_accession <- proteome_subset %>%
#   mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions )) %>%
#   dplyr::select(c(Genes_protein_accessions))
# 
# target <-unique(gene_protein_accession$Genes_protein_accessions)
# 
# #Make new data frame with columns Gene_protein_accessions2, PID_timepoint, and log2_quantification
# proteome_long2 <- proteome_subset_t %>%pivot_longer(cols = -PID_timepoint,
#   names_to = "Genes_protein_accessions", values_to = "log2_quantification")
# 
# proteome_subset_t2 <- gene_protein_accession %>% left_join(proteome_long2,
#   by = "Genes_protein_accessions")
# 
# proteome_subset_t2$timepoint <- substring(proteome_subset_t2,6)

proteome_subset <-  proteome_subset %>% 
  mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions ))

target <-unique(proteome_subset$Genes_protein_accessions)

t1="preATI"
t2="ATI4"

#paired.w.test <- data.frame()

paired.t.test <- data.frame()
linear_models <- data.frame()
paired.t.test_noOut <- data.frame()
linear_models_noOut <- data.frame()

for (i in 1:length(target)){
  proteome_subset_paired <- proteome_subset %>% dplyr::filter(timepoint %in% c(t1,t2)) %>% 
    dplyr::filter(Genes_protein_accessions %in% target[i]) 
  
  #if missing data - remove the pair
  #proteome_subset_paired %>% dplyr::filter(duplicated(PID))
  
  
  proteome_subset_paired <- proteome_subset_paired %>% dplyr::filter(!(PID %in% (((proteome_subset_paired %>%  group_by(PID) %>%
                                                                                     dplyr::filter(n() ==1))$PID)) ))
  
  t1_table <- proteome_subset_paired %>% dplyr::filter(timepoint %in% t1) %>% arrange(PID)
  t2_table <- proteome_subset_paired %>% dplyr::filter(timepoint %in% t2) %>% arrange(PID)
  
  differences <- (t1_table$log2_quantification - t2_table$log2_quantification)
  
  #perform shapiro-wilk test for normality on this vector of values
  #if(shapiro.test(differences)$p.value <=0.05){
  pvalue <- t.test(t1_table$log2_quantification, t2_table$log2_quantification, paired = TRUE)$p.value
  estimate <- t.test(t1_table$log2_quantification, t2_table$log2_quantification, paired = TRUE)$estimate
  CI_L <- t.test(t1_table$log2_quantification, t2_table$log2_quantification, paired = TRUE)$conf.int[1]
  CI_U <- t.test(t1_table$log2_quantification, t2_table$log2_quantification, paired = TRUE)$conf.int[2]
  pvalue_summary <- data.frame(Genes_protein_accessions= target[i], t1_t2_mean_diff = estimate , t1_t2_mean_diff_95CI_L=CI_L , t1_t2_mean_diff_95CI_U=CI_U ,  pvalue=pvalue, t1_t2_pair= paste0(t1,"_",t2), test="paired_t_test")
  paired.t.test <- rbind(paired.t.test, pvalue_summary)
  #}
  #if(shapiro.test(differences)$p.value >0.05){
  #pvalue <- wilcox.test(proteome_subset_paired$log2_quantification[proteome_subset_paired$timepoint == t1], proteome_subset_paired$log2_quantification[proteome_subset_paired$timepoint == t2], paired = TRUE)$p.value
  
  #pvalue_summary <- data.frame(Target= target[i], pvalue=pvalue, pair= paste0(t1,"_",t2), test="paired_wilcoxon_test")
  #paired.w.test <- rbind(paired.w.test, pvalue_summary)
  #}
  
  ###check for controller and non controller status effects using a linear model modeling the differences between the timepoints 
  #lm(diff ~ status)
  diff_table = data.frame(PID=t1_table$PID, status=t1_table$status, diff=differences)
  linear_model_target <- lm(data=diff_table, diff ~ status )
  
  linear_model_target_df <-  data.frame(Genes_protein_accessions= target[i], 
                                        t1_t2_diff= paste0(t1,"_", t2, "_diff"),
                                        intercept_beta= summary(linear_model_target)$coeff[1,1], 
                                        intercept_p= summary(linear_model_target)$coeff[1,4],
                                        non_controllers_vs_contr_beta= summary(linear_model_target)$coeff[2,1], 
                                        non_controllers_vs_contr_p= summary(linear_model_target)$coeff[2,4], 
                                        test="linear_model")
  linear_models <- rbind(linear_models, linear_model_target_df)
  
  
  ###removing sample 2887
  
  proteome_subset_paired_noOut <- proteome_subset %>% dplyr::filter(timepoint %in% c(t1,t2)) %>% 
    dplyr::filter(Genes_protein_accessions %in% target[i]) %>% dplyr::filter(!PID == "2887")
  
  #if missing data - remove the pair
  #proteome_subset_paired %>% dplyr::filter(duplicated(PID))
  
  
  proteome_subset_paired_noOut <- proteome_subset_paired_noOut %>% dplyr::filter(!(PID %in% (((proteome_subset_paired_noOut %>%  group_by(PID) %>%
                                                                                                 dplyr::filter(n() ==1))$PID)) ))
  
  t1_table_noOut <- proteome_subset_paired_noOut %>% dplyr::filter(timepoint %in% t1) %>% arrange(PID)
  t2_table_noOut <- proteome_subset_paired_noOut %>% dplyr::filter(timepoint %in% t2) %>% arrange(PID)
  
  differences_noOut <- (t1_table_noOut$log2_quantification - t2_table_noOut$log2_quantification)
  
  #perform shapiro-wilk test for normality on this vector of values
  #if(shapiro.test(differences)$p.value <=0.05){
  pvalue_noOut <- t.test(t1_table_noOut$log2_quantification, t2_table_noOut$log2_quantification, paired = TRUE)$p.value
  estimate_noOut <- t.test(t1_table_noOut$log2_quantification, t2_table_noOut$log2_quantification, paired = TRUE)$estimate
  CI_L_noOut <- t.test(t1_table_noOut$log2_quantification, t2_table_noOut$log2_quantification, paired = TRUE)$conf.int[1]
  CI_U_noOut <- t.test(t1_table_noOut$log2_quantification, t2_table_noOut$log2_quantification, paired = TRUE)$conf.int[2]
  pvalue_summary_noOut <- data.frame(Genes_protein_accessions= target[i], t1_t2_mean_diff = estimate_noOut , t1_t2_mean_diff_95CI_L=CI_L_noOut , t1_t2_mean_diff_95CI_U=CI_U_noOut ,  pvalue=pvalue_noOut, t1_t2_pair= paste0(t1,"_",t2), test="paired_t_test")
  paired.t.test_noOut <- rbind(paired.t.test_noOut, pvalue_summary_noOut)
  ###check for controller and non controller status effects using a linear model modeling the differences between the timepoints 
  #lm(diff ~ status)
  diff_table_noOut = data.frame(PID=t1_table_noOut$PID, status=t1_table_noOut$status, diff=differences_noOut)
  linear_model_target_noOut <- lm(data=diff_table_noOut, diff ~ status )
  
  linear_model_target_df_noOut <-  data.frame(Genes_protein_accessions= target[i], 
                                              t1_t2_diff= paste0(t1,"_", t2, "_diff"),
                                              intercept_beta= summary(linear_model_target_noOut)$coeff[1,1], 
                                              intercept_p= summary(linear_model_target_noOut)$coeff[1,4],
                                              non_controllers_vs_contr_beta= summary(linear_model_target_noOut)$coeff[2,1], 
                                              non_controllers_vs_contr_p= summary(linear_model_target_noOut)$coeff[2,4],
                                              test="linear_model")
  linear_models_noOut <- rbind(linear_models_noOut, linear_model_target_df_noOut)
  print(target[i])
  print(i)
}

#with outlier
#using only t test for the small sample size
combined_par_nonpar_test <- rbind(paired.t.test)
combined_par_nonpar_test_ord <- arrange(combined_par_nonpar_test, pvalue)

combined_par_nonpar_test_ord %<>% mutate(p_adj_BH = p.adjust(pvalue, method="BH"))
write.csv(combined_par_nonpar_test_ord, file=paste0(table_dir, "/paired_t_test/overall/paired_t_test_targets_",t1,"_",t2,".csv"), quote=F, row.names=F)

#linear models 
linear_models_ord <-  arrange(linear_models, intercept_p) %>% mutate(p_adj_BH = p.adjust(pvalue, method="BH"))
write.csv(linear_models_ord, file=paste0(table_dir, "/paired_t_test/by_status/paired_t_test_targets_",t1,"_",t2,"_status_effect.csv"), quote=F, row.names=F)

#without outlier

#using only t test for the small sample size
combined_par_nonpar_test_noOut <- rbind(paired.t.test_noOut )
combined_par_nonpar_test_ord_noOut <-  arrange(combined_par_nonpar_test_noOut, pvalue) %>%
  mutate(p_adj_fdr = p.adjust(pvalue, method="BH"))
write.csv(combined_par_nonpar_test_ord_noOut, file=paste0(table_dir, "/paired_t_test/overall/paired_t_test_targets_",t1,"_",t2,"_noOut_PID_2887.csv"), quote=F, row.names=F)


#linear models 
linear_models_ord_noOut <- arrange(linear_models_noOut, intercept_p) %>% mutate(p_adj_BH = p.adjust(pvalue, method="BH"))
write.csv(linear_models_ord_noOut, file=paste0(table_dir, "/paired_t_test/by_status/paired_t_test_targets_",t1,"_",t2,"_status_effect_noOut_PID_2887.csv"), quote=F, row.names=F)


##visualization all points
combined_par_nonpar_test_ord<-read.csv(paste0(table_dir, "/paired_t_test/overall/paired_t_test_targets_",t1,"_",t2,".csv"), header=T)
mean_diff_ord <- arrange(combined_par_nonpar_test_ord, Genes_protein_accessions) %>% mutate(rank = seq(1,nrow(.)))
colnames(combined_par_nonpar_test_ord)

# Make the Manhattan plot on the  results dataset
mean_diff_ord$pvalue <- mean_diff_ord$pvalue + 0.00000001

mean_diff_ord$dir <- "NA"
mean_diff_ord$dir[mean_diff_ord$t1_t2_mean_diff < 0 ] <- "NEG_mean_diff"
mean_diff_ord$dir[mean_diff_ord$t1_t2_mean_diff >= 0 ] <- "POS_mean_diff"

axis_set <- mean_diff_ord$rank

pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,".pdf"), height=10, width=40)
manhplot <- ggplot(mean_diff_ord,  aes(x = rank, y = -log10(pvalue), color = Genes_protein_accessions)) +
  geom_text(aes(label=ifelse(-log10(pvalue)>=2,as.character(Genes_protein_accessions),'')),hjust=0.5,vjust=0, angle=45) +
  geom_point(aes(size=abs(t1_t2_mean_diff), shape=dir), alpha = 0.35) + 
  scale_x_continuous(labels = mean_diff_ord$Genes_protein_accessions, 
                     breaks = seq(from = 1, by = 1, length.out = length(axis_set))) + 
  scale_y_continuous(breaks = seq(0, 8, by = 1)) +
  geom_hline(yintercept = -log10(0.01), color = "grey40", linetype = "dashed") +
  scale_shape_manual(values=c("NEG_mean_diff" = 17,"POS_mean_diff"= 19)) +
  theme(legend.position="none") +
  ggtitle(paste0(t1,"-",t2, " mean differences by target p<=0.01")) +
  labs(x = NULL, y = "-log10(P)" , size=10) +
  theme_minimal()+
  theme(legend.position="none")
# manhplot <- ggplot(mean_diff_ord,  aes(x = rank, y = -log10(pvalue), color = Genes_protein_accessions)) +
#   geom_point(aes(size=abs(t1_t2_mean_diff), shape=dir, alpha = 0.35)) +
#   geom_text(aes(label=ifelse(-log10(pvalue)>=1.3,as.character(Genes_protein_accessions),'')),hjust=0.5,vjust=0)+
#   scale_x_continuous(labels = mean_diff_ord$Genes_protein_accessions, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
#   ) +
#   scale_y_continuous(breaks = seq(0, 8, by = 1))+
#   geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") +
#   scale_shape_manual(values=c("NEG_mean_diff" = 17,"POS_mean_diff"= 19))+
#   theme(legend.position="none") +
#   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle(paste0(t1,"-",t2, " significant mean differences by target")) +
#   labs(x = NULL,
#        y = "-log10(P)" , size=10)+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), axis.title=element_text(size=20), axis.text.y = element_text(size = 20))

print(manhplot)
dev.off()


#visualization noOut
combined_par_nonpar_test_ord_noOut<-read.csv(paste0(table_dir, "/paired_t_test/overall/paired_t_test_targets_",t1,"_",t2,"_noOut_PID_2887.csv"), header=T)
mean_diff_ord_noOut <- arrange(combined_par_nonpar_test_ord_noOut, Genes_protein_accessions) %>% mutate(rank = seq(1,nrow(.)))

# Make the Manhattan plot on the  results dataset
mean_diff_ord_noOut$pvalue <- mean_diff_ord_noOut$pvalue + 0.00000001

mean_diff_ord_noOut$dir <- "NA"
mean_diff_ord_noOut$dir[mean_diff_ord_noOut$t1_t2_mean_diff < 0 ] <- "NEG_mean_diff"
mean_diff_ord_noOut$dir[mean_diff_ord_noOut$t1_t2_mean_diff >= 0 ] <- "POS_mean_diff"

axis_set <- mean_diff_ord_noOut$rank


pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,"_noOut_PID_2887.pdf"), height=10, width=40)
manhplot <- ggplot(mean_diff_ord_noOut,  aes(x = rank, y = -log10(pvalue), color = Genes_protein_accessions)) +
  geom_point(aes(size=abs(t1_t2_mean_diff), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse(-log10(pvalue)>=2,as.character(Genes_protein_accessions),'')),hjust=0.5,vjust=0, angle=45)+
  scale_x_continuous(label = mean_diff_ord_noOut$Genes_protein_accessions, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(breaks = seq(0, 2, by = 1))+
  geom_hline(yintercept = -log10(0.01), color = "grey40", linetype = "dashed") +
  scale_shape_manual(values=c("NEG_mean_diff" = 17,"POS_mean_diff"= 19))+
  
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste0(t1,"-",t2, " mean differences by target p<=0.01")) +
  labs(x = NULL,
       y = "-log10(P)" , size=10)+
  theme_minimal() +
  theme(legend.position="none") 
print(manhplot)

dev.off()


###linear models effects
linear_models_ord <- read.csv(paste0(table_dir, "/paired_t_test/by_status/paired_t_test_targets_",t1,"_",t2,"_status_effect.csv"))
non_controller_effect_ord <- arrange(linear_models_ord, Genes_protein_accessions) %>% mutate(rank = seq(1,nrow(.)))

# Make the Manhattan plot on the  results dataset
non_controller_effect_ord$non_controllers_vs_contr_p <- non_controller_effect_ord$non_controllers_vs_contr_p + 0.00000001

non_controller_effect_ord$dir <- "NA"
non_controller_effect_ord$dir[non_controller_effect_ord$non_controllers_vs_contr_beta < 0 ] <- "NEG_non_controller_diff"
non_controller_effect_ord$dir[non_controller_effect_ord$non_controllers_vs_contr_beta >= 0 ] <- "POS_non_controller_diff"

axis_set <- non_controller_effect_ord$rank


pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,"_by_status.pdf"), height=10, width=40)
manhplot <- ggplot(non_controller_effect_ord,  aes(x = rank, y = -log10(non_controllers_vs_contr_p), color = Genes_protein_accessions)) +
  geom_point(aes(size=abs(non_controllers_vs_contr_beta), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse((-log10(intercept_p)>=1.3 & -log10(non_controllers_vs_contr_p)>=1.3),as.character(Genes_protein_accessions),'')),hjust=0.5,vjust=0, angle=45)+
  scale_x_continuous(label = non_controller_effect_ord$Genes_protein_accessions, breaks = seq(from = 1, by = 1, length.out = length(non_controller_effect_ord$Genes_protein_accessions))
  ) +
  scale_y_continuous(breaks = seq(0, 2, by = 1))+
  geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") +
  scale_shape_manual(values=c("NEG_non_controller_diff" = 17,"POS_non_controller_diff"= 19))+
  theme(legend.position="none") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste0(t1,"-",t2, " mean differences by target p<=0.05 - non controllers vs controllers")) +
  labs(x = NULL,
       y = "-log10(P)" , size=10)+
  theme_minimal() +
  theme(legend.position="none")

print(manhplot)
dev.off()


###linear models effects no Out
linear_models_ord_noOut <- read.csv(paste0(table_dir, "/paired_t_test/by_status/paired_t_test_targets_",t1,"_",t2,"_status_effect_noOut_PID_2887.csv"))
non_controller_effect_ord_noOut <- arrange(linear_models_ord_noOut, Genes_protein_accessions) %>% mutate(rank = seq(1,nrow(.)))

# Make the Manhattan plot on the  results dataset
non_controller_effect_ord_noOut$non_controllers_vs_contr_p <- non_controller_effect_ord_noOut$non_controllers_vs_contr_p + 0.00000001

non_controller_effect_ord_noOut$dir <- "NA"
non_controller_effect_ord_noOut$dir[non_controller_effect_ord_noOut$non_controllers_vs_contr_beta < 0 ] <- "NEG_non_controller_diff"
non_controller_effect_ord_noOut$dir[non_controller_effect_ord_noOut$non_controllers_vs_contr_beta >= 0 ] <- "POS_non_controller_diff"

axis_set <- non_controller_effect_ord_noOut$rank


pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,"_noOut_PID_2887_by_status.pdf"), height=10, width=40)
manhplot <- ggplot(non_controller_effect_ord_noOut,  aes(x = rank, y = -log10(non_controllers_vs_contr_p), color = Genes_protein_accessions)) +
  geom_point(aes(size=abs(non_controllers_vs_contr_beta), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse((-log10(intercept_p)>=1.3 & -log10(non_controllers_vs_contr_p)>=1.3),as.character(Genes_protein_accessions),'')),hjust=0.5,vjust=0, angle=45)+
  scale_x_continuous(label = non_controller_effect_ord_noOut$Genes_protein_accessions, breaks = seq(from = 1, by = 1, length.out = length(non_controller_effect_ord_noOut$Genes_protein_accessions))
  ) +
  scale_y_continuous(breaks = seq(0, 2, by = 1))+
  geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") +
  scale_shape_manual(values=c("NEG_non_controller_diff" = 17,"POS_non_controller_diff"= 19))+
  theme(legend.position="none") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste0(t1,"-",t2, " mean differences by target p<=0.05 - non controllers vs controllers")) +
  labs(x = NULL,
       y = "-log10(P)" , size=10)+
  theme_minimal() +
  theme(legend.position="none")
print(manhplot)
dev.off()

# Filter for significant proteins (p < 0.05)
significant_proteins <- combined_par_nonpar_test_ord %>% filter(pvalue < 0.05)

# Print the names of significant proteins
print(significant_proteins$Genes_protein_accessions)


