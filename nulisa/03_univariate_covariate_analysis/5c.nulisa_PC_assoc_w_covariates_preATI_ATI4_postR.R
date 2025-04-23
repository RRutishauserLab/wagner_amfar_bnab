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


setwd("~/Dropbox (Gladstone)/GB-JW-1447/data/")
plot_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/plots"
table_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/tables"
data_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/data"

#uploading the nulisa seq dataset
nulisa <- read_excel(paste0(data_dir,"/NULISA.xlsx"), sheet="NPQ Counts")
head(nulisa)
summary(nulisa)
colnames(nulisa)

unique(nulisa$timepoint)

#create a vector of timepoints of interest
#dataframe
nulisa <- data.frame(nulisa)
head(nulisa)
#basic R
#nulisa_subset <- nulisa[nulisa$timepoint %in% t_points,]
#tidyverse
t_points <- c("BL0", "BL1", "1w", "preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","ATI10",
              "ATI7",  "ATI11",    "ATI8" ,"ATI9", "preR", "postR1", "postR2",
              "postR3", "postR4")
t_points             


nulisa_subset <- nulisa %>% dplyr::filter(timepoint %in% t_points)
head(nulisa_subset)

dim(nulisa)
dim(nulisa_subset)
#adding one column
nulisa_subset$PID_timepoint <- paste(nulisa_subset$PID, nulisa_subset$timepoint, sep="_")
#nulisa_subset <- nulisa_subset %>% mutate(PID_timepoint= paste(nulisa_subset$PID, nulisa_subset$timepoint, sep="_"))
#nulisa_subset$PID_timepoint <- paste0(nulisa_subset$PID, "_",nulisa_subset$timepoint)
#create a column of controllers / non controllers
non_controllers<- c("1292", "1408", "8106")
controllers <- nulisa_subset$PID[!nulisa_subset$PID %in% non_controllers]
no_rebound <- nulisa_subset$PID[nulisa_subset$PID %in% "2523"]
controllers_rebound <- controllers[ !controllers %in% no_rebound]

nulisa_subset <- nulisa_subset %>%
   mutate(status = case_when(PID %in% non_controllers ~ 'non_controllers',
                           PID %in% controllers ~ 'controllers')) %>%
   mutate(rebound = case_when( PID %in% no_rebound ~ 'no' , PID %in% c(controllers_rebound,non_controllers) ~ 'yes'
                              ))

head(nulisa_subset)  







##uploading the variable of interest to be used as covariate in a univariate lnear model

confounders <- read_excel(paste0(data_dir,"/20241104_amfAR_participants_key_COMPLETE_share.xlsx"), sheet="Values")
viral_loads <- read_excel(paste0(data_dir,"/postR_viral_loads.xlsx"))

confounders <- as.data.frame(confounders)

                                   #age_enrollment + 
                                   #highest_viral_load_ever (qt) +
                                   #ART_start(categorical early late acute) + 
                                   #viral_load_afterR1 (qt) +  
                                   #viral_load_afterR2 (qt) + 
                                   #viral_load_afterR3 (qt) + 
                                   #rebound_group(categorical , none for no rebound) +
                                   #setpoint_after_ART_stopped (qt)
colnames(confounders)
viral_loads <- as.data.frame(viral_loads)

#univariate model of interest 
#lmer(data=nulisa_subset, (NPQ-BL0) ~ (1|PID) + covariate1)
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

#baseline value to subtract from each NQP

nulisa_baseline0 <- nulisa_subset %>% dplyr::filter(timepoint %in% "BL0") %>% arrange(PID)
head(nulisa_baseline0)


t_points <- c("preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","preR", "postR1", "postR2",
              "postR3")
t_points    
#if needed subset of covariates
#covariateToInclude <- c("Age_enrollment",
                            # "Highest_viral_load_ever",
                             #"ART_start"
                          #)
covariateToInclude <- confounders_of_interest
nulisa_all_t <- nulisa_subset %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)

target_list <- unique(nulisa_baseline0$Target)
pid_list <- unique(nulisa_baseline0$PID)

nulisa_BL0_adj_NPQ_all_proteins <- data.frame()
SummaryCovariate_preATI_postR3 <- data.frame() 
#for each target extract the timepoint
for (tar in 1:length(target_list)){
nulisa_target <- nulisa_all_t %>%  dplyr::filter(Target %in% target_list[tar])
nulisa_baseline0_target <-  nulisa_baseline0 %>% dplyr::filter(Target %in% target_list[tar])

 

  #adjust NPQ values by BL0
  nulisa_BL0_adj_NPQ_target_all_pid <- data.frame()

  for (pid in 1:length(pid_list)){
  covariatesForUnivariate_PID <- covariatesForUnivariate[covariatesForUnivariate$PID %in% pid_list[pid],]
  covariatesForUnivariate_PID$PID <- as.character(covariatesForUnivariate_PID$PID)
  covariatesForUnivariate_PID_subset <- covariatesForUnivariate_PID %>% dplyr::select(any_of(covariateToInclude) )
  nulisa_baseline0_target_NPQ_pid <- nulisa_baseline0_target %>% 
                             dplyr::filter(PID %in% pid_list[pid]) %>% 
                             dplyr::select(NPQ) 
  
  
  nulisa_target_NPQ_BL0_adj <-  
                nulisa_target %>% 
                dplyr::filter(PID %in% pid_list[pid]) %>% 
                mutate(deltaNPQ = NPQ - nulisa_baseline0_target_NPQ_pid$NPQ) %>% as.data.frame()
  

  
  nulisa_target_BL0_adj_NPQ_df <- data.frame(target=target_list[tar], 
                                                      PID=as.character(pid_list[pid]), 
                                                      timepoint=nulisa_target_NPQ_BL0_adj$timepoint,
                                                      status=nulisa_target_NPQ_BL0_adj$status,
                                                      deltaNPQ=nulisa_target_NPQ_BL0_adj$deltaNPQ
                                                       )
  #adding the covariates
  nulisa_target_BL0_adj_NPQ_dff_w_covar_pid <- full_join(nulisa_target_BL0_adj_NPQ_df,covariatesForUnivariate_PID_subset, by="PID" )      
  #all timepoint for one given protein

  nulisa_BL0_adj_NPQ_target_all_pid  <- rbind(nulisa_BL0_adj_NPQ_target_all_pid,nulisa_target_BL0_adj_NPQ_dff_w_covar_pid ) 
}#pid

  nulisa_BL0_adj_NPQ_all_proteins <- rbind(nulisa_BL0_adj_NPQ_all_proteins, nulisa_BL0_adj_NPQ_target_all_pid)


}#target


#preR-postR3 - include all the individuals

    
    t_points <- c("preATI", "ATI1", "ATI2",  "ATI3","ATI4","preR", "postR1", "postR2",
              "postR3")
    
    nulisa_preATI_postR3 <- nulisa_BL0_adj_NPQ_all_proteins %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)
    nulisa_preATI_postR3_qt <- 
    nulisa_preATI_postR3 %>% mutate(timepoint = as.factor(timepoint),
                                  timepoint = relevel(timepoint, ref = "preATI") ,
                                  timepoint = as.numeric(timepoint),
                                  timepoint = timepoint - 1)

    
    
    #
nulisa_preATI_postR3_t <- nulisa_preATI_postR3_qt %>% 
                   pivot_wider(names_from = target, values_from = deltaNPQ) %>%
                   as.data.frame()

#variance across ssample and timepoints across 250 proteins
target_var <- data.frame(Target=unique(nulisa_preATI_postR3$target), var=diag(var(na.omit(nulisa_preATI_postR3_t[,-c(1:11)]))))
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


most_var_target_95th <- dplyr::select(nulisa_preATI_postR3_t ,one_of(Target_95perc$Target))
head(most_var_target_95th)

#PCA
res.pca.norm_95th<- prcomp(most_var_target_95th,  scale=T, center = T, rank.=10)
?prcomp()
res.pca.norm_95th


PCA_95th <- data.frame(nulisa_preATI_postR3_t[c(1:11)], 
                       res.pca.norm_95th$x)
                      
print(summary(res.pca.norm_95th)$importance[3,])
#use the first 7 PCs to perform the association with the covariates    + timepoint( column 2)
    


for (cov in 4:(ncol(PCA_95th)-10)) {
    print(paste0(names(PCA_95th[cov])))
    
    
    PCs <- colnames(PCA_95th)
    #use the first 8 PCs to perform the association with the covariates   
 for (pc in c(12:18)) {
    print(paste0(PCs[pc]))    
    lmm_univariate <- lmer(data=PCA_95th, PCA_95th[,pc] ~ (1|PID) + PCA_95th[,2] + PCA_95th[,cov])
     Sumcov <- data.frame(range="preATI-postR3", included="all", PC=PCs[pc], covariate=colnames(PCA_95th[cov]),t(data.frame( summary(lmm_univariate)$coeff[3,])))
     rownames(Sumcov) <- ""
     SummaryCovariate_preATI_postR3 <- rbind(SummaryCovariate_preATI_postR3,Sumcov )
}
}
    
        #no2523
    nulisa_preATI_postR3_No2523 <- nulisa_preATI_postR3_qt  %>% dplyr::filter(PID !="2523") %>% arrange(PID)

    #
nulisa_preATI_postR3_t <- nulisa_preATI_postR3_No2523 %>% 
                   pivot_wider(names_from = target, values_from = deltaNPQ) %>%
                   as.data.frame()

#variance across ssample and timepoints across 250 proteins
target_var <- data.frame(Target=unique(nulisa_preATI_postR3$target), var=diag(var(na.omit(nulisa_preATI_postR3_t[,-c(1:11)]))))
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


most_var_target_95th <- dplyr::select(nulisa_preATI_postR3_t ,one_of(Target_95perc$Target))
head(most_var_target_95th)

#PCA
res.pca.norm_95th<- prcomp(most_var_target_95th,  scale=T, center = T, rank.=10)
?prcomp()
res.pca.norm_95th


PCA_95th <- data.frame(nulisa_preATI_postR3_t[c(1:11)], 
                       res.pca.norm_95th$x)
                      
print(summary(res.pca.norm_95th)$importance[3,])
#use the first 7 PCs to perform the association with the covariates   
    


for (cov in 4:(ncol(PCA_95th)-10)) {
    print(paste0(names(PCA_95th[cov])))
    
    
    PCs <- colnames(PCA_95th)
     for (pc in c(12:18)) {
    print(paste0(PCs[pc]))    
    lmm_univariate <- lmer(data=PCA_95th, PCA_95th[,pc] ~ (1|PID) + PCA_95th[,2] + PCA_95th[,cov])
     Sumcov <- data.frame(range="preATI-postR3", included="no2523", PC=PCs[pc], covariate=colnames(PCA_95th[cov]),t(data.frame( summary(lmm_univariate)$coeff[3,])))
     rownames(Sumcov) <- ""
     SummaryCovariate_preATI_postR3 <- rbind(SummaryCovariate_preATI_postR3,Sumcov )
}
}
    


    nulisa_preATI_postR3_No2523_No2887 <- nulisa_preATI_postR3_qt %>% dplyr::filter(PID !="2887")

    #
nulisa_preATI_postR3_t <- nulisa_preATI_postR3_No2523_No2887 %>% 
                   pivot_wider(names_from = target, values_from = deltaNPQ) %>%
                   as.data.frame()

#variance across ssample and timepoints across 250 proteins
target_var <- data.frame(Target=unique(nulisa_preATI_postR3$target), var=diag(var(na.omit(nulisa_preATI_postR3_t[,-c(1:11)]))))
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


most_var_target_95th <- dplyr::select(nulisa_preATI_postR3_t ,one_of(Target_95perc$Target))
head(most_var_target_95th)

#PCA
res.pca.norm_95th<- prcomp(most_var_target_95th,  scale=T, center = T, rank.=10)
?prcomp()
res.pca.norm_95th


PCA_95th <- data.frame(nulisa_preATI_postR3_t[c(1:11)], 
                       res.pca.norm_95th$x)
                      
print(summary(res.pca.norm_95th)$importance[3,])
#use the first 7 PCs to perform the association with the covariates   
    


for (cov in 4:(ncol(PCA_95th)-10)) {
    print(paste0(names(PCA_95th[cov])))
    
    
    PCs <- colnames(PCA_95th)
     for (pc in c(12:18)) {
    print(paste0(PCs[pc]))    
    lmm_univariate <- lmer(data=PCA_95th, PCA_95th[,pc] ~ (1|PID) + PCA_95th[,2] + PCA_95th[,cov])
     Sumcov <- data.frame(range="preATI-postR3", PC=PCs[pc], included="no2523_no2887",covariate=colnames(PCA_95th[cov]),t(data.frame( summary(lmm_univariate)$coeff[3,])))
     rownames(Sumcov) <- ""
     SummaryCovariate_preATI_postR3 <- rbind(SummaryCovariate_preATI_postR3,Sumcov )
}
}
    

colnames(SummaryCovariate_preATI_postR3)[6] <- "StdError" 
colnames(SummaryCovariate_preATI_postR3)[9] <- "P"
SummaryCovariate_preATI_postR3_ord <- SummaryCovariate_preATI_postR3 %>% arrange(P)

write.csv(SummaryCovariate_preATI_postR3_ord, file=paste0(table_dir,"/univariate_covariate_models/preATI_postR3_association_between_PCs_covariates_timepoint_adjusted.csv"), row.names=F, quote=F)



    
   