rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(magrittr)
library(lme4)
library(lmerTest)
library(tidyr)
library(factoextra)


setwd("/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447")
plot_dir <- "/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447/results/proteomics/plots"
table_dir <- "/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447/results/proteomics/tables"
data_dir <- "/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447/data"

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

#create a vector of timepoints of interest
t_points <- c("BL0", "BL1", "1w", "preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","ATI10",
              "ATI7",  "ATI11",    "ATI8" ,"ATI9", "preR", "postR1", "postR2",
              "postR3", "postR4")

proteome_subset <- proteome_long %>% dplyr::filter(timepoint %in% t_points)
head(proteome_subset)

#adding one column
proteome_subset$PID_timepoint <- paste(proteome_subset$PID, proteome_subset$timepoint, sep="_")

#create a column of controllers / non controllers
non_controllers<- c("1292", "1408", "8106")
controllers <- proteome_subset$PID[!proteome_subset$PID %in% non_controllers]
no_rebound <- proteome_subset$PID[proteome_subset$PID %in% "2523"]
controllers_rebound <- controllers[ !controllers %in% no_rebound]

proteome_subset <- proteome_subset %>%
  mutate(status = case_when(PID %in% non_controllers ~ 'non_controllers',
                            PID %in% controllers ~ 'controllers')) %>%
  mutate(rebound = case_when(PID %in% no_rebound ~ 'no' , 
                      PID %in% c(controllers_rebound,non_controllers) ~ 'yes'))

proteome_subset <-  proteome_subset %>% 
  mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions ))

##uploading the variable of interest to be used as covariate in a univariate lnear model

confounders <- read_excel(paste0(data_dir,"/20241104_amfAR_participants_key_COMPLETE_share.xlsx"), sheet="Values")
viral_loads <- read_excel(paste0(data_dir,"/postR_viral_loads.xlsx"))

confounders <- as.data.frame(confounders)

#age_enrollment + highest_viral_load_ever (qt) + ART_start (categorical early late acute) + 
#viral_load_afterR1 (qt) +  viral_load_afterR2 (qt) + viral_load_afterR3 (qt) + 
#rebound_group(categorical , none for no rebound) + setpoint_after_ART_stopped (qt)
colnames(confounders)
viral_loads <- as.data.frame(viral_loads)

#univariate model of interest 
#lmer(data=proteome_subset, (NPQ-BL0) ~ (1|PID) + covariate1)
#get the coeff and pvalue for each covariate
#running the model for NPQ values between:
#preATI -preR
#preR-postR3
#preATI-postR3
##preATI-ATI6

colnames(viral_loads)
confounders_of_interest <- c("PID",
                             "timepoint",
                             "status",
                             "rebound",
                             "Age_enrollment",
                             "Highest_viral_load_ever",
                             "ART_start",
                             "Rebound_group" ,
                             "setpoint",
                             "postR1",
                             "postR2",
                             "postR3"
)
#extract the set of covariate of interest
confoundersForUnivariate <- confounders %>% dplyr::select(any_of(confounders_of_interest))
viralLoadsForUnivariate <- viral_loads %>% dplyr::select(any_of(confounders_of_interest))
#full set of covariates
covariatesForUnivariate <- merge(confoundersForUnivariate, viralLoadsForUnivariate, by="PID")
covariatesForUnivariate$Highest_viral_load_ever <- as.numeric(covariatesForUnivariate$Highest_viral_load_ever)
covariatesForUnivariate$ART_start <- as.factor(covariatesForUnivariate$ART_start)
covariatesForUnivariate$ART_start <- relevel(covariatesForUnivariate$ART_start, ref="EARLY")
covariatesForUnivariate$ART_start <- as.numeric(covariatesForUnivariate$ART_start)
covariatesForUnivariate$ART_start <- covariatesForUnivariate$ART_start - 1

covariatesForUnivariate$Rebound_group <- as.factor(covariatesForUnivariate$Rebound_group)
covariatesForUnivariate$Rebound_group <- relevel(covariatesForUnivariate$Rebound_group, ref="fast")
covariatesForUnivariate$Rebound_group <- as.numeric(covariatesForUnivariate$Rebound_group)
covariatesForUnivariate$Rebound_group <- covariatesForUnivariate$Rebound_group - 1

#baseline value to subtract from each NPQ

proteome_baseline0 <- proteome_subset %>% dplyr::filter(timepoint %in% "BL0") %>% arrange(PID)
head(proteome_baseline0)


t_points <- c("preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","preR", "postR1", "postR2",
              "postR3")
t_points    
#if needed subset of covariates
#covariateToInclude <- c("Age_enrollment",
# "Highest_viral_load_ever",
#"ART_start"
#)
covariateToInclude <- confounders_of_interest
proteome_all_t <- proteome_subset %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)

target_list <- unique(proteome_baseline0$Genes_protein_accessions)
pid_list <- unique(proteome_baseline0$PID)

proteome_BL0_adj_NPQ_all_proteins <- data.frame()
#for each target extract the timepoint
for (tar in 1:length(target_list)){
  proteome_target <- proteome_all_t %>%  dplyr::filter(Genes_protein_accessions %in% target_list[tar])
  proteome_baseline0_target <-  proteome_baseline0 %>% dplyr::filter(Genes_protein_accessions %in% target_list[tar])
  
  
  
  #adjust NPQ values by BL0
  proteome_BL0_adj_NPQ_target_all_pid <- data.frame()
  
  for (pid in 1:length(pid_list)){
    covariatesForUnivariate_PID <- covariatesForUnivariate[covariatesForUnivariate$PID %in% pid_list[pid],]
    covariatesForUnivariate_PID$PID <- as.character(covariatesForUnivariate_PID$PID)
    covariatesForUnivariate_PID_subset <- covariatesForUnivariate_PID %>% dplyr::select(any_of(covariateToInclude) )
    proteome_baseline0_target_NPQ_pid <- proteome_baseline0_target %>% 
      dplyr::filter(PID %in% pid_list[pid]) %>% 
      dplyr::select(log2_quantification) 
    
    
    proteome_target_NPQ_BL0_adj <-  
      proteome_target %>% 
      dplyr::filter(PID %in% pid_list[pid]) %>% 
      mutate(deltaNPQ = log2_quantification - proteome_baseline0_target_NPQ_pid$log2_quantification) %>% as.data.frame()
    
    
    
    proteome_target_BL0_adj_NPQ_df <- data.frame(target=target_list[tar], 
                                               PID=as.character(pid_list[pid]), 
                                               timepoint=proteome_target_NPQ_BL0_adj$timepoint,
                                               status=proteome_target_NPQ_BL0_adj$status,
                                               deltaNPQ=proteome_target_NPQ_BL0_adj$deltaNPQ
    )
    #adding the covariates
    proteome_target_BL0_adj_NPQ_dff_w_covar_pid <- full_join(proteome_target_BL0_adj_NPQ_df,covariatesForUnivariate_PID_subset, by="PID" )      
    #all timepoint for one given protein
    
    proteome_BL0_adj_NPQ_target_all_pid  <- rbind(proteome_BL0_adj_NPQ_target_all_pid,proteome_target_BL0_adj_NPQ_dff_w_covar_pid ) 
  }#pid
  
  proteome_BL0_adj_NPQ_all_proteins <- rbind(proteome_BL0_adj_NPQ_all_proteins, proteome_BL0_adj_NPQ_target_all_pid)
  print(paste0("Done with: ",target_list[tar],",number ", tar ))
  
}#target

write.csv(proteome_BL0_adj_NPQ_all_proteins,file=paste0(table_dir,"/input_univariate_models/proteome_all_proteins_covariates.csv"), row.names=F, quote=F)
proteome_BL0_adj_NPQ_all_proteins_2 <- read.csv(paste0(table_dir,"/input_univariate_models/proteome_all_proteins_covariates.csv"))
# proteome_BL0_adj_NPQ_all_proteins_2[1:10,1:10]
# head(proteome_BL0_adj_NPQ_all_proteins)
# tail(proteome_BL0_adj_NPQ_all_proteins)


#preATI-ATI6 - include all the individuals


t_points <- c("preATI", "ATI1", "ATI2",  "ATI3","ATI4")

proteome_preATI_ATI4 <- proteome_BL0_adj_NPQ_all_proteins_2 %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)

#transform timepoint in quantitative

proteome_preATI_ATI4_qt <- 
  proteome_preATI_ATI4 %>% mutate(timepoint = as.factor(timepoint),
                                timepoint = relevel(timepoint, ref = "preATI") ,
                                timepoint = as.numeric(timepoint),
                                timepoint = timepoint - 1)

#Add new column with PID_timepoint 
proteome_preATI_ATI4_qt <- proteome_preATI_ATI4_qt %>% mutate(PID_timepoint = paste0(PID,"_",timepoint))
                                      
proteome_preATI_ATI4_t <- proteome_preATI_ATI4_qt %>% 
  pivot_wider(names_from = target, values_from = deltaNPQ) %>%
  as.data.frame()


#variance across sample and timepoints across 250 proteins
target_var <- data.frame(Target=unique(proteome_preATI_ATI4_qt$target), var=diag(var(na.omit(proteome_preATI_ATI4_t[,-c(1:12)]))))
target_var_rank <- target_var[order(target_var$var, decreasing=T),]

median_target_var_rank <- median(target_var_rank$var)
median_target_var_rank

perc95_target_var_rank <- quantile(target_var_rank$var, 0.95)
perc95_target_var_rank

#filtering the target above the threshold of var

Target_95perc<- target_var_rank %>% 
  dplyr::filter(var >= perc95_target_var_rank) %>% 
  dplyr::select(Target)


most_var_target_95th <- dplyr::select(proteome_preATI_ATI4_t ,one_of(Target_95perc$Target))
head(most_var_target_95th)

#PCA
res.pca.norm_95th<- prcomp(most_var_target_95th,  scale=T, center = T, rank.=60)
?prcomp()
res.pca.norm_95th


PCA_95th <- data.frame(proteome_preATI_ATI4_t[c(1:12)], 
                       res.pca.norm_95th$x)

print(summary(res.pca.norm_95th)$importance[3,])
#use the first 6 PCs to perform the association with the covariates   

SummaryCovariate_preATI_ATI4  <- data.frame() 


for (cov in 4:11) {
  print(paste0(names(PCA_95th[cov])))
  
  
  PCs <- colnames(PCA_95th)
  for (pc in c(13:37)) {
    print(paste0(PCs[pc])) 
    #PC1 ~ (1|PID) + timepoint+ age
    
    lmm_univariate <- lmer(data=PCA_95th, PCA_95th[,pc] ~ (1|PID) + PCA_95th[,2] + PCA_95th[,cov])
    Sumcov <- data.frame(range="preATI-ATI4", included="all", PC=PCs[pc], covariate=colnames(PCA_95th[cov]),t(data.frame( summary(lmm_univariate)$coeff[3,])))
    rownames(Sumcov) <- ""
    SummaryCovariate_preATI_ATI4 <- rbind(SummaryCovariate_preATI_ATI4,Sumcov )
  }#PCs
}#cov

#no2523
proteome_preATI_ATI4_No2523 <- proteome_preATI_ATI4_qt  %>% dplyr::filter(PID !="2523") %>% arrange(PID)

proteome_preATI_ATI4_t <- proteome_preATI_ATI4_No2523 %>% 
  pivot_wider(names_from = target, values_from = deltaNPQ) %>%
  as.data.frame()

#variance across sample and timepoints across 250 proteins
target_var <- data.frame(Target=unique(proteome_preATI_ATI4_No2523$target), var=diag(var(na.omit(proteome_preATI_ATI4_t[,-c(1:12)]))))

#Original Code
#target_var <- data.frame(Target=unique(proteome_preATI_ATI6$target), var=diag(var(na.omit(proteome_preATI_ATI6_t[,-c(1:12)]))))

target_var_rank <- target_var[order(target_var$var, decreasing=T),]
head(target_var)
head(target_var_rank)
median_target_var_rank <- median(target_var_rank$var)
median_target_var_rank

perc95_target_var_rank <- quantile(target_var_rank$var, 0.95)
perc95_target_var_rank


#filtering the target above the threshold of var

Target_95perc<- target_var_rank %>% 
  dplyr::filter(var >= perc95_target_var_rank) %>% 
  dplyr::select(Target)


most_var_target_95th <- dplyr::select(proteome_preATI_ATI4_t ,one_of(Target_95perc$Target))
head(most_var_target_95th)

#PCA
res.pca.norm_95th<- prcomp(most_var_target_95th,  scale=T, center = T, rank.=60)
?prcomp()
res.pca.norm_95th


PCA_95th <- data.frame(proteome_preATI_ATI4_t[c(1:12)], 
                       res.pca.norm_95th$x)

print(summary(res.pca.norm_95th)$importance[3,])
#use the first 6 PCs to perform the association with the covariates   



for (cov in 4:11) {
  print(paste0(names(PCA_95th[cov])))
  
  PCs <- colnames(PCA_95th)
  for (pc in c(13:37)) {
    print(paste0(PCs[pc]))    
    lmm_univariate <- lmer(data=PCA_95th, PCA_95th[,pc] ~ (1|PID) + PCA_95th[,2] + PCA_95th[,cov])
    Sumcov <- data.frame(range="preATI-ATI4", included="no2523", PC=PCs[pc], covariate=colnames(PCA_95th[cov]),t(data.frame( summary(lmm_univariate)$coeff[3,])))
    rownames(Sumcov) <- ""
    SummaryCovariate_preATI_ATI4 <- rbind(SummaryCovariate_preATI_ATI4,Sumcov )
  }
}

# test_1 = unique(proteome_preATI_ATI6$target)
# test_2 = unique(proteome_preATI_ATI6_No2523$target)

#Remove 2523 and 2887
proteome_preATI_ATI4_No2523_No2887 <- proteome_preATI_ATI4_No2523 %>% dplyr::filter(PID !="2887")

proteome_preATI_ATI4_t <- proteome_preATI_ATI4_No2523_No2887 %>% 
  pivot_wider(names_from = target, values_from = deltaNPQ) %>%
  as.data.frame()

#variance across ssample and timepoints across 250 proteins
target_var <- data.frame(Target=unique(proteome_preATI_ATI4$target), var=diag(var(na.omit(proteome_preATI_ATI4_t[,-c(1:12)]))))
target_var_rank <- target_var[order(target_var$var, decreasing=T),]

median_target_var_rank <- median(target_var_rank$var)

perc95_target_var_rank <- quantile(target_var_rank$var, 0.95)

#filtering the target above the threshold of var

Target_95perc<- target_var_rank %>% 
  dplyr::filter(var >= perc95_target_var_rank) %>% 
  dplyr::select(Target)

most_var_target_95th <- dplyr::select(proteome_preATI_ATI4_t ,one_of(Target_95perc$Target))

#PCA
res.pca.norm_95th<- prcomp(most_var_target_95th,  scale=T, center = T, rank.=60)
res.pca.norm_95th


PCA_95th <- data.frame(proteome_preATI_ATI4_t[c(1:12)], 
                       res.pca.norm_95th$x)

print(summary(res.pca.norm_95th)$importance[3,])
#use the first 6 PCs to perform the association with the covariates   


for (cov in 4:11) {
  print(paste0(names(PCA_95th[cov])))
  
  PCs <- colnames(PCA_95th)
  for (pc in c(13:37)) {
    print(paste0(PCs[pc]))    
    lmm_univariate <- lmer(data=PCA_95th, PCA_95th[,pc] ~ (1|PID) + PCA_95th[,2] + PCA_95th[,cov])
    Sumcov <- data.frame(range="preATI-ATI4", PC=PCs[pc], included="no2523_no2887",covariate=colnames(PCA_95th[cov]),t(data.frame( summary(lmm_univariate)$coeff[3,])))
    rownames(Sumcov) <- ""
    SummaryCovariate_preATI_ATI6 <- rbind(SummaryCovariate_preATI_ATI4,Sumcov )
  }
}

head(SummaryCovariate_preATI_ATI4)

colnames(SummaryCovariate_preATI_ATI4)[6] <- "StdError" 
colnames(SummaryCovariate_preATI_ATI4)[9] <- "P"
SummaryCovariate_preATI_ATI4_ord <- SummaryCovariate_preATI_ATI4 %>% arrange(P)
head(SummaryCovariate_preATI_ATI4_ord)
write.csv(SummaryCovariate_preATI_ATI4_ord, file=paste0(table_dir,"/univariate_covariate_models/preATI_ATI4_association_between_PCs_covariates_timepoint_adjusted.csv"), row.names=F, quote=F)
