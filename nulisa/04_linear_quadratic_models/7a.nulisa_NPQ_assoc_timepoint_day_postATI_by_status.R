rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(magrittr)
library(lme4)
library(lmerTest)
library(lmtest)
library(tidyr)
library(factoextra)
library(qvalue)
library(stringr)


setwd("~/Dropbox (Gladstone)/GB-JW-1447/data/")
plot_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/plots"
table_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/tables"
data_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/data"
if (!dir.exists(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/linear"))) {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/linear"), recursive=T)}  
if (!dir.exists(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/quadratic"))) {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/quadratic"), recursive=T)}  
if (!dir.exists(paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI"))) {dir.create(paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI"), recursive=T)}  

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
                             "ATI61",
                             "ATI62",
                             "ATI63"
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
              "postR3", "postR4")
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

nulisa_target_all_proteins<- data.frame()
anova_summary <- data.frame()
anova_summary_no2523 <- data.frame()
anova_summary_no2523_no2887 <- data.frame()

#for each target extract the timepoint
for (tar in 1:length(target_list)){
nulisa_target <- nulisa_all_t %>%  dplyr::filter(Target %in% target_list[tar])
nulisa_baseline0_target <-  nulisa_baseline0 %>% dplyr::filter(Target %in% target_list[tar])

 

  #adjust NPQ values by BL0
  nulisa_target_all_pid <- data.frame()

  for (pid in 1:length(pid_list)){
  covariatesForUnivariate_PID <- covariatesForUnivariate[covariatesForUnivariate$PID %in% pid_list[pid],]
  covariatesForUnivariate_PID$PID <- as.character(covariatesForUnivariate_PID$PID)
  covariatesForUnivariate_PID_subset <- covariatesForUnivariate_PID %>% dplyr::select(any_of(covariateToInclude) )
  nulisa_target_NPQ_pid <- nulisa_baseline0_target %>% 
                             dplyr::filter(PID %in% pid_list[pid]) %>% 
                             dplyr::select(NPQ) 
   
  
  nulisa_target_NPQ <-  
                nulisa_target %>% 
                dplyr::filter(PID %in% pid_list[pid]) %>% 
                 as.data.frame()
   
  

  
  nulisa_target_NPQ_df <- data.frame(target=target_list[tar], 
                                                      PID=as.character(pid_list[pid]), 
                                                      timepoint=nulisa_target_NPQ$timepoint,
                                                      status=nulisa_target_NPQ$status,
                                                      NPQ=nulisa_target_NPQ$NPQ
                                                       )
 #adding the covariates
  #adding the covariates
  nulisa_target_NPQ_dff_w_covar_pid <- full_join(nulisa_target_NPQ_df,covariatesForUnivariate_PID_subset, by="PID" )      
  #all timepoint for one given protein

  nulisa_target_all_pid  <- rbind(nulisa_target_all_pid,nulisa_target_NPQ_dff_w_covar_pid ) 
}#pid

  nulisa_target_all_proteins <- rbind(nulisa_target_all_proteins, nulisa_target_all_pid)




########################################################################################################################
     #including up to 70 days post ATI
    #uploading the day from ATI
    days_to_ATI <- read.csv(paste0(data_dir,"/time_to_ATI_day0.csv"))



    
   
        
    days_to_ATI_timepoint_range <- days_to_ATI  %>% arrange(PID)
    days_to_ATI_timepoint_range_long <- days_to_ATI_timepoint_range %>% pivot_longer(!c(PID), names_to = "timepoint", values_to = "days_from_ATI")%>% as.data.frame()
    days_to_ATI_timepoint_range_long <- na.omit(days_to_ATI_timepoint_range_long)
    days_to_ATI_timepoint_range_long_filtered <- days_to_ATI_timepoint_range_long %>% mutate(PID_timepoint = paste0(PID, "_", timepoint))  %>% 
                                                 dplyr::filter(days_from_ATI <= 70) %>%
                                                 dplyr::filter(timepoint != "BL0")
                                                  
    #include the PID_tiempoint selected above
    nulisa_timepoint_ranges_cov <- nulisa_target_all_pid %>% mutate(PID_timepoint = paste0(PID, "_", timepoint))%>%
                                       dplyr::filter(PID_timepoint %in%days_to_ATI_timepoint_range_long_filtered$PID_timepoint)  %>% arrange(PID)
 
   #merge days to corresponding nulisa timepoint 
   

    nulisa_days_post_ATI_cov <- left_join(nulisa_timepoint_ranges_cov, days_to_ATI_timepoint_range_long_filtered[-c(1,2)], by="PID_timepoint")
    
    
    nulisa_days_post_ATI_cov_qt <- 
    nulisa_days_post_ATI_cov%>% mutate(days_from_ATI = as.numeric(days_from_ATI ))
 
    #all samples
    min_model_days_from_ATI <- lmer(data=nulisa_days_post_ATI_cov_qt, NPQ ~ (1|PID) + 
                                                                        days_from_ATI)
    min_model_days_from_ATI2 <- lmer(data=nulisa_days_post_ATI_cov_qt, NPQ ~ (1|PID) + days_from_ATI + I(days_from_ATI^2))

    full_model_days_from_ATI_status <- lmer(data=nulisa_days_post_ATI_cov_qt, NPQ ~ (1|PID) + 
                                                                                days_from_ATI*status )
    full_model_days_from_ATI2_status <- lmer(data=nulisa_days_post_ATI_cov_qt, NPQ ~ (1|PID) +  days_from_ATI*status + I(days_from_ATI^2) + 
                                                                                 I(days_from_ATI^2):status)
    anova_days_from_ATI  <- anova(min_model_days_from_ATI,  full_model_days_from_ATI_status)
    anova_days_from_ATI2 <- anova(min_model_days_from_ATI2,  full_model_days_from_ATI2_status)
    
    anova_summary_df <- data.frame(range="70_days_postATI", 
                                   included="all_ind" ,
                                   variable="NPQ" ,
                                   target=target_list[tar], 
                                   best_fit=ifelse(as.numeric(anova(full_model_days_from_ATI_status, full_model_days_from_ATI2_status)$Pr[2]) <= 0.05 , "quadratic", "linear"),                                   
                                   estimate_status_in_linear=round(summary(full_model_days_from_ATI_status)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_days_from_ATI_status)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_days_from_ATI_status)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status)$coeff[4,5], digits=5),
                                   anova_linear_p=round(anova_days_from_ATI[2,8], digits=5),
                                   estimate_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[5,5], digits=5),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[6,1], digits=7),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[6,2], digits=7),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status)$coeff[6,5], digits=5),
                                   anova_quadratic_p=round(anova_days_from_ATI2[2,8], digits=5)

                                  )
                                   
    anova_summary <- rbind(anova_summary, anova_summary_df)
    
        
######################################################################################################
#plot quadratic
if(anova_summary_df$anova_quadratic_p <=0.05){
    fixed_effects <- fixef(full_model_days_from_ATI2_status)
    choose_timepoints_2_predict <- seq(min(nulisa_days_post_ATI_cov_qt$days_from_ATI), max(nulisa_days_post_ATI_cov_qt$days_from_ATI), by = 2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        choose_timepoints_2_predict, ##status:timepoint term
                                        choose_timepoints_2_predict^2 ##status:timepoint^2 term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "days_from_ATI2",
                                              "status_days_from_ATI",
                                              "status_days_from_ATI2")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        rep(0, length(choose_timepoints_2_predict)), ##status:timepoint term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint^2 term
    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "days_from_ATI2",
                                              "status_days_from_ATI",
                                              "status_days_from_ATI2")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_days_from_ATI2_status)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/quadratic/all_preATI_up_to70_days_postATI_associated_target_",target_list[tar],"_quadratic_lmer.pdf"))
    print(ggplot(nulisa_days_post_ATI_cov_qt, aes(days_from_ATI, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = days_from_ATI, y=NPQ), linewidth = 1)
      
    )
    dev.off() 
    } 
######################################################################################################
 ######################################################################################################
#plot linear
    if(anova_summary_df$anova_linear_p <=0.05){
    fixed_effects <- fixef(full_model_days_from_ATI_status)
    choose_timepoints_2_predict <- seq(min(nulisa_days_post_ATI_cov_qt$days_from_ATI), max(nulisa_days_post_ATI_cov_qt$days_from_ATI), by = 2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict ##status:timepoint term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "status_days_from_ATI")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint term
                                    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "status_days_from_ATI")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_days_from_ATI_status)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/linear/all_preATI_up_to70_days_postATI_associated_target_",target_list[tar],"_linear_lmer.pdf"))
    print(ggplot(nulisa_days_post_ATI_cov_qt, aes(days_from_ATI, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = days_from_ATI, y=NPQ), linewidth = 1)
            
    )
    dev.off()
    }
######################################################################################################
 
  
 ########################################################################################################################   
    
    #no2523
    nulisa_days_post_ATI_cov_no2523_qt <- 
    nulisa_days_post_ATI_cov%>% mutate(days_from_ATI = as.numeric(days_from_ATI )) %>% dplyr::filter(PID != "2523")

 
    #all samples
    min_model_days_from_ATI_no2523 <- lmer(data=nulisa_days_post_ATI_cov_no2523_qt, NPQ ~ (1|PID) + 
                                                                        days_from_ATI)
    min_model_days_from_ATI2_no2523  <- lmer(data=nulisa_days_post_ATI_cov_no2523_qt, NPQ ~ (1|PID) + days_from_ATI + I(days_from_ATI^2))

    full_model_days_from_ATI_status_no2523  <- lmer(data=nulisa_days_post_ATI_cov_no2523_qt, NPQ ~ (1|PID) + 
                                                                                days_from_ATI*status )
    full_model_days_from_ATI2_status_no2523  <- lmer(data=nulisa_days_post_ATI_cov_no2523_qt, NPQ ~ (1|PID) +  days_from_ATI*status + I(days_from_ATI^2) + 
                                                                                 I(days_from_ATI^2):status)
    anova_days_from_ATI_no2523   <- anova(min_model_days_from_ATI_no2523,  full_model_days_from_ATI_status_no2523)
    anova_days_from_ATI2_no2523  <- anova(min_model_days_from_ATI2_no2523,  full_model_days_from_ATI2_status_no2523)
    
    anova_summary_df_no2523 <- data.frame(range="70_days_postATI", 
                                   included="no2523" ,
                                   variable="NPQ" ,
                                   target=target_list[tar], 
                                   best_fit=ifelse(as.numeric(anova(full_model_days_from_ATI_status_no2523, full_model_days_from_ATI2_status_no2523)$Pr[2]) <= 0.05 , "quadratic", "linear"),                                   
                                   estimate_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523)$coeff[4,5], digits=5),
                                   anova_linear_p=round(anova_days_from_ATI_no2523[2,8], digits=5),
                                   estimate_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[5,5], digits=5),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[6,1], digits=7),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[6,2], digits=7),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523)$coeff[6,5], digits=5),
                                   anova_quadratic_p=round(anova_days_from_ATI2_no2523[2,8], digits=5)

                                  )
                                   
    anova_summary_no2523 <- rbind(anova_summary_no2523, anova_summary_df_no2523)
 ######################################################################################################
#plot quadratic
if(anova_summary_df_no2523$anova_quadratic_p <=0.05){
    fixed_effects <- fixef(full_model_days_from_ATI2_status_no2523)
    choose_timepoints_2_predict <- seq(min(nulisa_days_post_ATI_cov_no2523_qt$days_from_ATI), max(nulisa_days_post_ATI_cov_no2523_qt$days_from_ATI), by = 2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        choose_timepoints_2_predict, ##status:timepoint term
                                        choose_timepoints_2_predict^2 ##status:timepoint^2 term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "days_from_ATI2",
                                              "status_days_from_ATI",
                                              "status_days_from_ATI2")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        rep(0, length(choose_timepoints_2_predict)), ##status:timepoint term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint^2 term
    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "days_from_ATI2",
                                              "status_days_from_ATI",
                                              "status_days_from_ATI2")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_days_from_ATI2_status_no2523)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/quadratic/no2523_preATI_up_to70_days_postATI_associated_target_",target_list[tar],"_quadratic_lmer.pdf"))
    print(ggplot(nulisa_days_post_ATI_cov_no2523_qt, aes(days_from_ATI, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = days_from_ATI, y=NPQ), linewidth = 1)
      
    )
    dev.off() 
    } 
######################################################################################################
 ######################################################################################################
#plot linear
    if(anova_summary_df_no2523$anova_linear_p <=0.05){
    fixed_effects <- fixef(full_model_days_from_ATI_status_no2523)
    choose_timepoints_2_predict <- seq(min(nulisa_days_post_ATI_cov_no2523_qt$days_from_ATI), max(nulisa_days_post_ATI_cov_no2523_qt$days_from_ATI), by = 2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict ##status:timepoint term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "status_days_from_ATI")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint term
                                    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "status_days_from_ATI")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_days_from_ATI_status_no2523)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/linear/no_2523_preATI_up_to70_days_postATI_associated_target_",target_list[tar],"_linear_lmer.pdf"))
    print(ggplot(nulisa_days_post_ATI_cov_no2523_qt, aes(days_from_ATI, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = days_from_ATI, y=NPQ), linewidth = 1)
            
    )
    dev.off()
    }
######################################################################################################
 
 ########################################################################################################################   
 #no2523 no 2887

 nulisa_days_post_ATI_cov_no2523_no2887_qt <- 
    nulisa_days_post_ATI_cov%>% mutate(days_from_ATI = as.numeric(days_from_ATI )) %>% dplyr::filter(PID != "2523")  %>% dplyr::filter(PID != "2887")
 
    #all samples
    min_model_days_from_ATI_no2523_no2887 <- lmer(data=nulisa_days_post_ATI_cov_no2523_no2887_qt, NPQ ~ (1|PID) + 
                                                                        days_from_ATI)
    min_model_days_from_ATI2_no2523_no2887  <- lmer(data=nulisa_days_post_ATI_cov_no2523_no2887_qt, NPQ ~ (1|PID) + days_from_ATI + I(days_from_ATI^2))

    full_model_days_from_ATI_status_no2523_no2887  <- lmer(data=nulisa_days_post_ATI_cov_no2523_no2887_qt, NPQ ~ (1|PID) + 
                                                                                days_from_ATI*status )
    full_model_days_from_ATI2_status_no2523_no2887  <- lmer(data=nulisa_days_post_ATI_cov_no2523_no2887_qt, NPQ ~ (1|PID) +  days_from_ATI*status + I(days_from_ATI^2) + 
                                                                                 I(days_from_ATI^2):status)
    anova_days_from_ATI_no2523_no2887   <- anova(min_model_days_from_ATI_no2523_no2887,  full_model_days_from_ATI_status_no2523_no2887)
    anova_days_from_ATI2_no2523_no2887  <- anova(min_model_days_from_ATI2_no2523_no2887,  full_model_days_from_ATI2_status_no2523_no2887)
    
    anova_summary_df_no2523_no2887 <- data.frame(range="70_days_postATI", 
                                   included="no2523_no2887" ,
                                   variable="NPQ" ,
                                   target=target_list[tar], 
                                   best_fit=ifelse(as.numeric(anova(full_model_days_from_ATI_status_no2523_no2887, full_model_days_from_ATI2_status_no2523_no2887)$Pr[2]) <= 0.05 , "quadratic", "linear"),                                   
                                   estimate_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523_no2887)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523_no2887)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523_no2887)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523_no2887)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523_no2887)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_days_from_ATI_status_no2523_no2887)$coeff[4,5], digits=5),
                                   anova_linear_p=round(anova_days_from_ATI_no2523_no2887[2,8], digits=5),
                                   estimate_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[3,5], digits=5),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[5,5], digits=5),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[6,1], digits=7),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[6,2], digits=7),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_days_from_ATI2_status_no2523_no2887)$coeff[6,5], digits=5),
                                   anova_quadratic_p=round(anova_days_from_ATI2_no2523_no2887[2,8], digits=5)

                                  )
                                   
    anova_summary_no2523_no2887 <- rbind(anova_summary_no2523_no2887, anova_summary_df_no2523_no2887)
  ######################################################################################################
#plot quadratic
if(anova_summary_df_no2523_no2887$anova_quadratic_p <=0.05){
    fixed_effects <- fixef(full_model_days_from_ATI2_status_no2523_no2887)
    choose_timepoints_2_predict <- seq(min(nulisa_days_post_ATI_cov_no2523_no2887_qt$days_from_ATI), max(nulisa_days_post_ATI_cov_no2523_no2887_qt$days_from_ATI), by = 2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        choose_timepoints_2_predict, ##status:timepoint term
                                        choose_timepoints_2_predict^2 ##status:timepoint^2 term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "days_from_ATI2",
                                              "status_days_from_ATI",
                                              "status_days_from_ATI2")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict^2, ##timepoint^2 term
                                        rep(0, length(choose_timepoints_2_predict)), ##status:timepoint term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint^2 term
    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "days_from_ATI2",
                                              "status_days_from_ATI",
                                              "status_days_from_ATI2")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_days_from_ATI2_status_no2523_no2887)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/quadratic/no2523_no2887_preATI_up_to70_days_postATI_associated_target_",target_list[tar],"_quadratic_lmer.pdf"))
    print(ggplot(nulisa_days_post_ATI_cov_no2523_no2887_qt, aes(days_from_ATI, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = days_from_ATI, y=NPQ), linewidth = 1)
      
    )
    dev.off() 
    } 
######################################################################################################
 ######################################################################################################
#plot linear
    if(anova_summary_df_no2523_no2887$anova_linear_p <=0.05){
    fixed_effects <- fixef(full_model_days_from_ATI_status_no2523_no2887)
    choose_timepoints_2_predict <- seq(min(nulisa_days_post_ATI_cov_no2523_no2887_qt$days_from_ATI), max(nulisa_days_post_ATI_cov_no2523_no2887_qt$days_from_ATI), by = 2)
    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        choose_timepoints_2_predict ##status:timepoint term
                                        )
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "status_days_from_ATI")
    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        choose_timepoints_2_predict, ##timepoint term
                                        rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                        rep(0, length(choose_timepoints_2_predict)) ##status:timepoint term
                                    )
    colnames(controller_pred_matrix) <- c("intercept",
                                              "days_from_ATI",
                                              "statusnoncontroller",
                                              "status_days_from_ATI")
    
    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_days_from_ATI_status_no2523_no2887)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_noncontroller,
                                                  NPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(days_from_ATI = choose_timepoints_2_predict,
                                                  NPQ = predicted_mean_controller,
                                                  NPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/linear/no2523_no2887preATI_up_to70_days_postATI_associated_target_",target_list[tar],"_linear_lmer.pdf"))
    print(ggplot(nulisa_days_post_ATI_cov_no2523_no2887_qt, aes(days_from_ATI, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = days_from_ATI, y=NPQ), linewidth = 1)
            
    )
    dev.off()
    }
###
########################################################################################################################   



}#target
########################################################################################################################   

 
    nulisa_timepoint_range_all_target <-   nulisa_target_all_proteins %>% mutate(PID_timepoint = paste0(PID, "_", timepoint))  %>% 
                                         dplyr::filter(PID_timepoint %in%days_to_ATI_timepoint_range_long_filtered$PID_timepoint)  %>% arrange(PID) #transform timepoint in quantitative
    
   
  
    nulisa_all_proteins_cov_days_to_ATI <- left_join(nulisa_timepoint_range_all_target, days_to_ATI_timepoint_range_long_filtered[-c(1,2)], by="PID_timepoint")
    #for plots
    nulisa_all_proteins_cov_days_to_ATI_qt <-nulisa_all_proteins_cov_days_to_ATI%>% mutate(days_from_ATI = as.numeric(days_from_ATI ))


    head(nulisa_all_proteins_cov_days_to_ATI_qt)
    tail(nulisa_all_proteins_cov_days_to_ATI_qt)
    
    
    


    anova_summary_adj <-  arrange(anova_summary, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                         mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr)     
    anova_summary_sorted <- anova_summary_adj %>% arrange(anova_linear_lfdr)
    #anova_summary_sorted <- anova_summary_adj %>% arrange(anova_linear_p)
  
    write.csv(anova_summary_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/preATI_up_to70_days_postATI_association_between_target_timepoint_timepoint2_by_status_p_adjusted_no2523.csv"), row.names=F, quote=F)


    anova_summary_no2523_adj <-  arrange(anova_summary_no2523, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                         mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr)     
    anova_summary_no2523_sorted <- anova_summary_no2523_adj %>% arrange(anova_linear_lfdr)
  
    write.csv(anova_summary_no2523_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/preATI_up_to70_days_postATI_association_between_target_timepoint_timepoint2_by_status_p_adjusted_no2523.csv"), row.names=F, quote=F)

    anova_summary_no2523_no2887_adj <-  arrange(anova_summary_no2523_no2887, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                         mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr)     

    anova_summary_no2523_no2887_sorted <- anova_summary_no2523_no2887_adj %>% arrange(anova_linear_lfdr)
  
    write.csv(anova_summary_no2523_no2887_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/preATI_up_to70_days_postATI_association_between_target_timepoint_timepoint2_by_status_p_adjusted_no2523_no2887.csv"), row.names=F, quote=F)





