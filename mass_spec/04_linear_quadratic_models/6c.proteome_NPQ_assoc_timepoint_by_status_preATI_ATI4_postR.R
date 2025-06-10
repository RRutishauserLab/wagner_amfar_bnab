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


setwd("~/Gladstone Dropbox/Michela Traglia/GB-JW-1447")
plot_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/proteomics/plots"
table_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/proteomics/tables"
data_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/data"

#setwd("/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447")
#plot_dir <- "/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447/results/proteomics/plots/timepoint_ordinal"
#table_dir <- "/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447/results/proteomics/tables/timepoint_ordinal"
#data_dir <- "/Users/axie/Gladstone Dropbox/Alan Xie/GB-JW-1447/data"

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

#basic R
#tidyverse
t_points <- c("BL0", "BL1", "1w", "preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","ATI10",
              "ATI7",  "ATI11",    "ATI8" ,"ATI9", "preR", "postR1", "postR2",
              "postR3", "postR4")
t_points             


proteome_subset <- proteome_long %>% dplyr::filter(timepoint %in% t_points)
head(proteome_subset)

#adding one column
proteome_subset$PID_timepoint <- paste(proteome_subset$PID, proteome_subset$timepoint, sep="_")
#proteome_subset <- proteome_subset %>% mutate(PID_timepoint= paste(proteome_subset$PID, proteome_subset$timepoint, sep="_"))
#proteome_subset$PID_timepoint <- paste0(proteome_subset$PID, "_",proteome_subset$timepoint)
#create a column of controllers / non controllers
non_controllers<- c("1292", "1408", "8106")
controllers <- proteome_subset$PID[!proteome_subset$PID %in% non_controllers]
no_rebound <- proteome_subset$PID[proteome_subset$PID %in% "2523"]
controllers_rebound <- controllers[ !controllers %in% no_rebound]

proteome_subset <- proteome_subset %>%
   mutate(status = case_when(PID %in% non_controllers ~ 'non_controllers',
                           PID %in% controllers ~ 'controllers')) %>%
   mutate(rebound = case_when( PID %in% no_rebound ~ 'no' , PID %in% c(controllers_rebound,non_controllers) ~ 'yes'
                              ))

proteome_subset <-  proteome_subset %>% 
  mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions ))

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
#lmer(data=proteome_subset, (NPQ-BL0) ~ (1|PID) + covariate1)
#get the coeff and pvalue for each covariate
#running the model for NPQ values between:
#preATI -preR
#preR-ATI63
#preATI-ATI63
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
# covariatesForUnivariate$ART_start <- as.factor(covariatesForUnivariate$ART_start)
# covariatesForUnivariate$ART_start <- relevel(covariatesForUnivariate$ART_start, ref="EARLY")
# covariatesForUnivariate$ART_start <- as.numeric(covariatesForUnivariate$ART_start)
# covariatesForUnivariate$ART_start <- covariatesForUnivariate$ART_start - 1

covariatesForUnivariate$Rebound_group <- as.factor(covariatesForUnivariate$Rebound_group)
covariatesForUnivariate$Rebound_group <- relevel(covariatesForUnivariate$Rebound_group, ref="fast")
covariatesForUnivariate$Rebound_group <- as.numeric(covariatesForUnivariate$Rebound_group)
covariatesForUnivariate$Rebound_group <- covariatesForUnivariate$Rebound_group - 1

#baseline value to subtract from each NQP

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

proteome_target_all_proteins <- data.frame()
anova_summary <- data.frame()
anova_summary_no2523 <- data.frame()
anova_summary_no2523_no2887 <- data.frame()

#for each target extract the timepoint

###test witha  subset of target witht he issue of NPQ ==0
####target_list <- target_list[3300:3330]
target_list <- unique(proteome_baseline0$Genes_protein_accessions)
skipped_targets <- data.frame()

#Open bracket error
for (tar in 1:length(target_list)){
  print(paste0("target ", tar, " starting ...."))
  proteome_target <- proteome_all_t %>%  dplyr::filter(Genes_protein_accessions %in% target_list[tar])
  proteome_baseline0_target <-  proteome_baseline0 %>% dplyr::filter(Genes_protein_accessions %in% target_list[tar])


  #adjust NPQ values by BL0
  proteome_target_all_pid <- data.frame()

  for (pid in 1:length(pid_list)){
  covariatesForUnivariate_PID <- covariatesForUnivariate[covariatesForUnivariate$PID %in% pid_list[pid],]
  covariatesForUnivariate_PID$PID <- as.character(covariatesForUnivariate_PID$PID)
  covariatesForUnivariate_PID_subset <- covariatesForUnivariate_PID %>% dplyr::select(any_of(covariateToInclude) )
  proteome_baseline0_target_NPQ_pid <- proteome_baseline0_target %>% 
                             dplyr::filter(PID %in% pid_list[pid]) %>% 
                             dplyr::select(log2_quantification) 
  
  
 
  
  proteome_target_NPQ_df <- data.frame(target=target_list[tar], 
                                                      PID=as.character(pid_list[pid]), 
                                                      timepoint=proteome_target$timepoint,
                                                      status=proteome_target$status,
                                                      NPQ=proteome_target$log2_quantification
                                                       )
  #adding the covariates
  proteome_target_NPQ_dff_w_covar_pid <- full_join(proteome_target_NPQ_df,covariatesForUnivariate_PID_subset, by="PID" )      
  #all timepoint for one given protein

   proteome_target_all_pid  <- rbind(proteome_target_all_pid,proteome_target_NPQ_dff_w_covar_pid ) 



  }#pid
  proteome_target_all_proteins <- rbind(proteome_target_all_proteins, proteome_target_all_pid)

   #preATI-ATI6 - include all the individuals
    t_points <- c("preATI", "ATI1", "ATI2",  "ATI3","ATI4", "preR", "postR1", "postR2",  "postR3")

 
    
    proteome_preATI_postR3 <- proteome_target_all_pid %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)
    proteome_preATI_postR3_all_target <-   proteome_target_all_proteins %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID) #transform timepoint in quantitative
    
    proteome_preATI_postR3_qt <- 
    proteome_preATI_postR3 %>% mutate(timepoint = as.factor(timepoint),
                                  timepoint = relevel(timepoint, ref = "preATI") ,
                                  timepoint = as.numeric(timepoint),
                                  timepoint = timepoint - 1
                                  )
    if(as.vector(quantile(proteome_preATI_postR3_qt$NPQ, 0.50)) ==0)
    { 
    anova_summary <- anova_summary
    print(paste0("skipped target ", target_list[tar]))
    skipped_individual_target <- data.frame(target_w_zero_NPQ = target_list[tar])
    skipped_targets <- rbind(skipped_targets, skipped_individual_target)
    }

  if(as.vector(quantile(proteome_preATI_postR3_qt$NPQ, 0.50)) !=0)
    { 
     #all samples
    min_model_timepoint <- lmer(data=proteome_preATI_postR3_qt, NPQ ~ (1|PID) + 
                                                                        timepoint)
    min_model_timepoint2 <- lmer(data=proteome_preATI_postR3_qt, NPQ ~ (1|PID) + timepoint + I(timepoint^2))

    full_model_timepoint_status <- lmer(data=proteome_preATI_postR3_qt, NPQ ~ (1|PID) + 
                                                                                timepoint*status )
    full_model_timepoint2_status <- lmer(data=proteome_preATI_postR3_qt, NPQ ~ (1|PID) +  timepoint*status + I(timepoint^2) + 
                                                                                 I(timepoint^2):status)
    anova_timepoint  <- anova(min_model_timepoint,  full_model_timepoint_status)
    anova_timepoint2 <- anova(min_model_timepoint2,  full_model_timepoint2_status)
    
      anova_summary_df <- data.frame(range="preATI-postR3", 
                                   included="all_ind" ,
                                   variable="NPQ" ,
                                   target=target_list[tar],
                                   best_fit=ifelse(as.numeric(anova(full_model_timepoint_status, full_model_timepoint2_status)$Pr[2]) <= 0.05 , "quadratic", "linear"),
                                   estimate_status_in_linear=round(summary(full_model_timepoint_status)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_timepoint_status)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_timepoint_status)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status)$coeff[4,5], digits=5),
                                   anova_linear_p=round(anova_timepoint[2,8], digits=5),
                                   estimate_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[5,5], digits=5),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[6,1], digits=3),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[6,2], digits=3),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status)$coeff[6,5], digits=5),
                                   anova_quadratic_p=round(anova_timepoint2[2,8], digits=5)
                                   )
                                   
    anova_summary <- rbind(anova_summary, anova_summary_df)
  }
######################################################################################################
#plot quadratic
if(anova_summary_df$anova_quadratic_p <=0.05){ ###change
    fixed_effects <- fixef(full_model_timepoint2_status) ###change
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        choose_timepoints_2_predict, ##status:timepoint term
                                        choose_timepoints_2_predict^2 ##status:timepoint^2 term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "timepoint2",
                                              "status_timepoint",
                                              "status_timepoint2")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        rep(0, length(choose_timepoints_2_predict)), ##status:timepoint term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint^2 term
    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "timepoint2",
                                              "status_timepoint",
                                              "status_timepoint2")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint2_status) ###change

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/quadratic/all_preATI_postR3_associated_target_", str_split(gsub("/", "_", target_list[tar]), ";")[[1]][1],"_quadratic_lmer.pdf")) ###change
    print(ggplot(proteome_preATI_postR3_qt, aes(timepoint, NPQ, color = status)) +  ###change
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
      
    )
    
    print(ggplot(proteome_preATI_postR3_qt, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    
    dev.off() 
} 
    
######################################################################################################
######################################################################################################
#plot linear
    if(anova_summary_df$anova_linear_p <=0.05){ ###change
    fixed_effects <- fixef(full_model_timepoint_status) ###change
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict ##status:timepoint term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "status_timepoint")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint term
                                    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "status_timepoint")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint_status) ###change

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/linear/all_preATI_postR3_associated_target_", str_split(gsub("/", "_", target_list[tar]), ";")[[1]][1],"_linear_lmer.pdf")) ###change
    print(ggplot(proteome_preATI_postR3_qt, aes(timepoint, NPQ, color = status)) + ###change
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
            
    )
    
    print(ggplot(proteome_preATI_postR3_qt, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    dev.off()

    }
  
    #no2523
    proteome_preATI_postR3_qt_no2523 <- 
    proteome_preATI_postR3 %>% mutate(timepoint = as.factor(timepoint),
                                  timepoint = relevel(timepoint, ref = "preATI") ,
                                  timepoint = as.numeric(timepoint),
                                  timepoint = timepoint - 1
                                  ) %>% dplyr::filter(PID != "2523")

    if(as.vector(quantile(proteome_preATI_postR3_qt_no2523$NPQ, 0.50)) ==0)
    { 
    anova_summary_no2523 <- anova_summary_no2523
    print(paste0("skipped target ", target_list[tar]))
    skipped_individual_target <- data.frame(target_w_zero_NPQ = target_list[tar])
    skipped_targets <- rbind(skipped_targets, skipped_individual_target)
 
    }

    if(as.vector(quantile(proteome_preATI_postR3_qt_no2523$NPQ, 0.50)) !=0)
    { 
  
    min_model_timepoint_no2523 <- lmer(data=proteome_preATI_postR3_qt_no2523, NPQ ~ (1|PID) + 
                                                                        timepoint)
    min_model_timepoint2_no2523 <- lmer(data=proteome_preATI_postR3_qt_no2523, NPQ ~ (1|PID) + timepoint + I(timepoint^2))

    full_model_timepoint_status_no2523 <- lmer(data=proteome_preATI_postR3_qt_no2523, NPQ ~ (1|PID) + 
                                                                                timepoint*status )
    full_model_timepoint2_status_no2523 <- lmer(data=proteome_preATI_postR3_qt_no2523, NPQ ~ (1|PID) +  timepoint*status + I(timepoint^2) + 
                                                                                 I(timepoint^2):status)
    anova_timepoint_no2523  <- anova(min_model_timepoint_no2523,  full_model_timepoint_status_no2523)
    anova_timepoint2_no2523<- anova(min_model_timepoint2_no2523,  full_model_timepoint2_status_no2523)
    
    anova_summary_no2523_df <- data.frame(range="preATI-postR3", 
                                   included="no_2523" ,
                                   variable="NPQ" ,
                                   target=target_list[tar], 
                                   best_fit=ifelse(as.numeric(anova(full_model_timepoint_status_no2523, full_model_timepoint2_status_no2523)$Pr[2]) <= 0.05 , "quadratic", "linear"),                                   
                                   estimate_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[4,5], digits=5),
                                   anova_linear_p=round(anova_timepoint_no2523[2,8], digits=5),
                                   estimate_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[5,5], digits=5),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[6,1], digits=3),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[6,2], digits=3),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[6,5], digits=5),
                                   anova_quadratic_p=round(anova_timepoint2_no2523[2,8], digits=5)
                                   )
                                   

  anova_summary_no2523 <- rbind(anova_summary_no2523, anova_summary_no2523_df)
  
  
  #####put the plot and change the tables
  #Plot quadratic
  if(anova_summary_no2523_df$anova_quadratic_p <=0.05){
    fixed_effects <- fixef(full_model_timepoint2_status_no2523)
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        choose_timepoints_2_predict, ##status:timepoint term
                                        choose_timepoints_2_predict^2 ##status:timepoint^2 term
    )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "timepoint2",
                                              "status_timepoint",
                                              "status_timepoint2")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                    choose_timepoints_2_predict, ##timepoint term
                                    rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                    choose_timepoints_2_predict^2, ##timepoint^2 term
                                    rep(0, length(choose_timepoints_2_predict)), ##status:timepoint term
                                    rep(0, length(choose_timepoints_2_predict)) ##status:timepoint^2 term
    )
    colnames(controller_pred_matrix) <- c("intercept",
                                          "timepoint",
                                          "statusnoncontroller",
                                          "timepoint2",
                                          "status_timepoint",
                                          "status_timepoint2")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint2_status_no2523)
    
    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                               NPQ = predicted_mean_controller,
                                               NPQ_se = predicted_se_controller,
                                               status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
    
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/quadratic/no_2523_preATI_postR3_associated_target_", str_split(gsub("/", "_", target_list[tar]), ";")[[1]][1],"_quadratic_lmer.pdf"))
    print(ggplot(proteome_preATI_postR3_qt_no2523, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
          
    )
    
    print(ggplot(proteome_preATI_postR3_qt_no2523, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    dev.off()
    
  }
  #Plot linear
  if(anova_summary_no2523_df$anova_linear_p <=0.05){
    fixed_effects <- fixef(full_model_timepoint_status_no2523)
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict ##status:timepoint term
    )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "status_timepoint")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                    choose_timepoints_2_predict, ##timepoint term
                                    rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                    rep(0, length(choose_timepoints_2_predict)) ##status:timepoint term
    )
    colnames(controller_pred_matrix) <- c("intercept",
                                          "timepoint",
                                          "statusnoncontroller",
                                          "status_timepoint")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint_status_no2523)
    
    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                               NPQ = predicted_mean_controller,
                                               NPQ_se = predicted_se_controller,
                                               status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
    
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/linear/no2523_preATI_postR3_associated_target_",str_split(gsub("/", "_", target_list[tar]), ";")[[1]][1],"_linear_lmer.pdf"))
    print(ggplot(proteome_preATI_postR3_qt_no2523, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
          
    )
    
    print(ggplot(proteome_preATI_postR3_qt_no2523, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    dev.off()
  }
  }
    #no2523 no 2887
    proteome_preATI_postR3_qt_no2523_no2887 <- 
    proteome_preATI_postR3_qt_no2523 %>% dplyr::filter(PID != "2887")

    if(as.vector(quantile(proteome_preATI_postR3_qt_no2523_no2887$NPQ, 0.50)) ==0)
    { 
    anova_summary_no2523_no2887 <- anova_summary_no2523_no2887
    print(paste0("skipped target ", target_list[tar]))
    skipped_individual_target <- data.frame(target_w_zero_NPQ = target_list[tar])
    skipped_targets <- rbind(skipped_targets, skipped_individual_target)
 
    }

    if(as.vector(quantile(proteome_preATI_postR3_qt_no2523_no2887$NPQ, 0.50)) !=0)
    { 
    min_model_timepoint_no2523_no2887 <- lmer(data=proteome_preATI_postR3_qt_no2523_no2887, NPQ ~ (1|PID) + 
                                                                        timepoint)
    min_model_timepoint2_no2523_no2887 <- lmer(data=proteome_preATI_postR3_qt_no2523_no2887, NPQ ~ (1|PID) + timepoint + I(timepoint^2))

    full_model_timepoint_status_no2523_no2887 <- lmer(data=proteome_preATI_postR3_qt_no2523_no2887, NPQ ~ (1|PID) + 
                                                                                timepoint*status )
    full_model_timepoint2_status_no2523_no2887 <- lmer(data=proteome_preATI_postR3_qt_no2523_no2887, NPQ ~ (1|PID) +  timepoint*status + I(timepoint^2) + 
                                                                                 I(timepoint^2):status)
    anova_timepoint_no2523_no2887  <- anova(min_model_timepoint_no2523_no2887,  full_model_timepoint_status_no2523_no2887)
    anova_timepoint2_no2523_no2887<- anova(min_model_timepoint2_no2523_no2887,  full_model_timepoint2_status_no2523_no2887)
    
    
    anova_summary_no2523_no2887_df  <- data.frame(range="preATI-postR3", 
                                   included="no_2523_no2887" ,
                                   variable="NPQ" ,
                                   target=target_list[tar], 
                                   best_fit=ifelse(as.numeric(anova(full_model_timepoint_status_no2523_no2887, full_model_timepoint2_status_no2523_no2887)$Pr[2]) <= 0.05 , "quadratic", "linear"),
                                   estimate_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[4,5], digits=5),
                                   anova_linear_p=round(anova_timepoint_no2523_no2887[2,8], digits=5),
                                   estimate_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[5,5], digits=5),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[6,1], digits=3),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[6,2], digits=3),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[6,5], digits=5),
                                   anova_quadratic_p=round(anova_timepoint2_no2523_no2887[2,8], digits=5)
                                   )
                         
      anova_summary_no2523_no2887 <- rbind(anova_summary_no2523_no2887, anova_summary_no2523_no2887_df)
      
          
######################################################################################################
#plot quadratic
if(anova_summary_no2523_no2887_df$anova_quadratic_p <=0.05){
    fixed_effects <- fixef(full_model_timepoint2_status_no2523_no2887)
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        choose_timepoints_2_predict, ##status:timepoint term
                                        choose_timepoints_2_predict^2 ##status:timepoint^2 term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "timepoint2",
                                              "status_timepoint",
                                              "status_timepoint2")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        rep(0, length(choose_timepoints_2_predict)), ##status:timepoint term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint^2 term
    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "timepoint2",
                                              "status_timepoint",
                                              "status_timepoint2")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint2_status_no2523_no2887)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/quadratic/no2523_no2887_preATI_postR3_associated_target_", str_split(gsub("/", "_", target_list[tar]), ";")[[1]][1],"_quadratic_lmer_no2523_no2887.pdf"))
    print(ggplot(proteome_preATI_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
      
    )
    
    print(ggplot(proteome_preATI_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    
    dev.off() 
    } 
######################################################################################################
######################################################################################################
#plot linear
    if(anova_summary_no2523_no2887_df$anova_linear_p <=0.05){
    fixed_effects <- fixef(full_model_timepoint_status_no2523_no2887)
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict ##status:timepoint term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "status_timepoint")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint term
                                    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "timepoint",
                                              "statusnoncontroller",
                                              "status_timepoint")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint_status_no2523_no2887)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/linear/no2523_no2887_preATI_postR3_associated_target_", str_split(gsub("/", "_", target_list[tar]), ";")[[1]][1],"_linear_lmer_no2523_no2887.pdf"))
    print(ggplot(proteome_preATI_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
            
    )
    
    print(ggplot(proteome_preATI_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    
    dev.off()

      
      
      }                               

  print(paste0(tar, " done"))
#   }
  
}#target
}
     anova_summary_adj <- arrange(anova_summary, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                         mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr) 
                         
    anova_summary_sorted <- anova_summary_adj %>% arrange(anova_linear_lfdr)
  
    write.csv(anova_summary_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/preATI_postR3_association_between_target_timepoint_timepoint2_by_status_p_adjusted_all.csv"), row.names=F, quote=F)

    anova_summary_no2523_adj <-  arrange(anova_summary_no2523, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                         mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr)     
    anova_summary_no2523_sorted <- anova_summary_no2523_adj %>% arrange(anova_linear_lfdr)
  
    write.csv(anova_summary_no2523_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/preATI_postR3_association_between_target_timepoint_timepoint2_by_status_p_adjusted_no2523.csv"), row.names=F, quote=F)

    anova_summary_no2523_no2887_adj <-  arrange(anova_summary_no2523_no2887, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                         mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr)     

anova_summary_no2523_no2887_sorted <- anova_summary_no2523_no2887_adj %>% arrange(anova_linear_lfdr)
  
    write.csv(anova_summary_no2523_no2887_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/preATI_postR3_association_between_target_timepoint_timepoint2_by_status_p_adjusted_no2523_no2887.csv"), row.names=F, quote=F)




    write.csv(unique(skipped_targets), file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/preATI_postR3_skipped_targets w_50perc_NPQ_equal_to_0_list.csv"), row.names=F, quote=F)

