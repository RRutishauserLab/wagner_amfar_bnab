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
if (!dir.exists(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/quadratic"))) {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/quadratic"), recursive=T)}  
if (!dir.exists(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/linear"))) {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/linear"), recursive=T)}  
if (!dir.exists(paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal"))) {dir.create(paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal"), recursive=T)}  


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
nulisa_subset <- nulisa_subset  %>% arrange(PID)
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
#covariatesForUnivariate$ART_start <- as.factor(covariatesForUnivariate$ART_start)
#covariatesForUnivariate$ART_start <- relevel(covariatesForUnivariate$ART_start, ref="EARLY")
#covariatesForUnivariate$ART_start <- as.numeric(covariatesForUnivariate$ART_start)
#covariatesForUnivariate$ART_start <- covariatesForUnivariate$ART_start - 1

covariatesForUnivariate$Rebound_group <- as.factor(covariatesForUnivariate$Rebound_group)
covariatesForUnivariate$Rebound_group <- relevel(covariatesForUnivariate$Rebound_group, ref="fast")
covariatesForUnivariate$Rebound_group <- as.numeric(covariatesForUnivariate$Rebound_group)
covariatesForUnivariate$Rebound_group <- covariatesForUnivariate$Rebound_group - 1

#baseline value to subtract from each NQP




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

target_list <- unique(nulisa_subset$Target)
pid_list <- unique(nulisa_subset$PID)
  nulisa_target_all_proteins <- data.frame()
anova_summary <- data.frame()
anova_summary_no2523 <- data.frame()
anova_summary_no2523_no2887 <- data.frame()

#for each target extract the timepoint
for (tar in 1:length(target_list)){
nulisa_target <- nulisa_all_t %>%  dplyr::filter(Target %in% target_list[tar])

 
  nulisa_target_all_pid <- data.frame()

  for (pid in 1:length(pid_list)){
  covariatesForUnivariate_PID <- covariatesForUnivariate[covariatesForUnivariate$PID %in% pid_list[pid],]
  covariatesForUnivariate_PID$PID <- as.character(covariatesForUnivariate_PID$PID)
  covariatesForUnivariate_PID_subset <- covariatesForUnivariate_PID %>% dplyr::select(any_of(covariateToInclude) )
  nulisa_target_NPQ_pid <- nulisa_target %>% 
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
  nulisa_target_NPQ_dff_w_covar_pid <- full_join(nulisa_target_NPQ_df,covariatesForUnivariate_PID_subset, by="PID" )      
  #all timepoint for one given protein

  nulisa_target_all_pid  <- rbind(nulisa_target_all_pid,nulisa_target_NPQ_dff_w_covar_pid ) 
}#pid

  nulisa_target_all_proteins <- rbind(nulisa_target_all_proteins, nulisa_target_all_pid)




   #preR-postR3 - include all the individuals



    t_points <- c( "preR", "postR1", "postR2",  "postR3")    
 
 
    
    nulisa_preR_postR3 <- nulisa_target_all_pid %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)
    nulisa_preR_postR3_all_target <-   nulisa_target_all_proteins %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID) #transform timepoint in quantitative
    
    nulisa_preR_postR3_qt <- 
    nulisa_preR_postR3 %>% mutate(timepoint = as.factor(timepoint),
                                  timepoint = relevel(timepoint, ref = "preR") ,
                                  timepoint = as.numeric(timepoint),
                                  timepoint = timepoint - 1
                                  )

 
       
######################################################################################################
    
  
     #no2523
    nulisa_preR_postR3_qt_no2523 <- 
    nulisa_preR_postR3 %>% mutate(timepoint = as.factor(timepoint),
                                  timepoint = relevel(timepoint, ref = "preR") ,
                                  timepoint = as.numeric(timepoint),
                                  timepoint = timepoint - 1
                                  ) %>% dplyr::filter(PID != "2523")

    min_model_timepoint_no2523 <- lmer(data=nulisa_preR_postR3_qt_no2523, NPQ ~ (1|PID) + 
                                                                        timepoint)
    min_model_timepoint2_no2523 <- lmer(data=nulisa_preR_postR3_qt_no2523, NPQ ~ (1|PID) + timepoint + I(timepoint^2))

    full_model_timepoint_status_no2523 <- lmer(data=nulisa_preR_postR3_qt_no2523, NPQ ~ (1|PID) + 
                                                                                timepoint*status )
    full_model_timepoint2_status_no2523 <- lmer(data=nulisa_preR_postR3_qt_no2523, NPQ ~ (1|PID) +  timepoint*status + I(timepoint^2) + 
                                                                                 I(timepoint^2):status)
    anova_timepoint_no2523  <- anova(min_model_timepoint_no2523,  full_model_timepoint_status_no2523)
    anova_timepoint2_no2523<- anova(min_model_timepoint2_no2523,  full_model_timepoint2_status_no2523)
    
  anova_summary_no2523_df    <- data.frame(range="preR-postR3", 
                                   included="no_2523" ,
                                   variable="NPQ" ,
                                   target=target_list[tar], 
                                   best_fit=ifelse(as.numeric(anova(full_model_timepoint_status_no2523, full_model_timepoint2_status_no2523)$Pr[2]) <= 0.05 , "quadratic", "linear"),                                   
                                   estimate_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[3,5], digits=10),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523)$coeff[4,5], digits=10),
                                   anova_linear_p=round(anova_timepoint_no2523[2,8], digits=10),
                                   estimate_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[3,5], digits=10),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[5,5], digits=10),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[6,1], digits=3),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[6,2], digits=3),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523)$coeff[6,5], digits=10),
                                   anova_quadratic_p=round(anova_timepoint2_no2523[2,8], digits=10)
                                   )
                                   

  anova_summary_no2523 <- rbind(anova_summary_no2523, anova_summary_no2523_df)
  ######################################################################################################
#plot quadratic
if(anova_summary_no2523_df$anova_quadratic_p <=0.05){
    fixed_effects <- fixef(full_model_timepoint2_status_no2523)
    choose_timepoints_2_predict <- seq(0, 3, by = 0.2)
    
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
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/quadratic/no_2523_preR_postR3_associated_target_",target_list[tar],"_quadratic_lmer.pdf"))
    print(ggplot(nulisa_preR_postR3_qt_no2523, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
      
    )
    
       print(ggplot(nulisa_preR_postR3_qt_no2523, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
      
    )
    dev.off()  
    }
######################################################################################################
 ######################################################################################################
#plot linear
if(anova_summary_no2523_df$anova_linear_p <=0.05){
    fixed_effects <- fixef(full_model_timepoint_status_no2523)
    choose_timepoints_2_predict <- seq(0, 3, by = 0.2)
    
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
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/linear/no2523_preR_postR3_associated_target_",target_list[tar],"_linear_lmer.pdf"))
    print(ggplot(nulisa_preR_postR3_qt_no2523, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
            
    )
    
       print(ggplot(nulisa_preR_postR3_qt_no2523, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    dev.off()
    }
######################################################################################################
 
     
######################################################################################################
#no2523 no 2887

 nulisa_preR_postR3_qt_no2523_no2887 <- 
    nulisa_preR_postR3_qt_no2523 %>% dplyr::filter(PID != "2887")

    min_model_timepoint_no2523_no2887 <- lmer(data=nulisa_preR_postR3_qt_no2523_no2887, NPQ ~ (1|PID) + 
                                                                        timepoint)
    min_model_timepoint2_no2523_no2887 <- lmer(data=nulisa_preR_postR3_qt_no2523_no2887, NPQ ~ (1|PID) + timepoint + I(timepoint^2))

    full_model_timepoint_status_no2523_no2887 <- lmer(data=nulisa_preR_postR3_qt_no2523_no2887, NPQ ~ (1|PID) + 
                                                                                timepoint*status )
    full_model_timepoint2_status_no2523_no2887 <- lmer(data=nulisa_preR_postR3_qt_no2523_no2887, NPQ ~ (1|PID) +  timepoint*status + I(timepoint^2) + 
                                                                                 I(timepoint^2):status)
    anova_timepoint_no2523_no2887  <- anova(min_model_timepoint_no2523_no2887,  full_model_timepoint_status_no2523_no2887)
    anova_timepoint2_no2523_no2887<- anova(min_model_timepoint2_no2523_no2887,  full_model_timepoint2_status_no2523_no2887)
    
    
    anova_summary_no2523_no2887_df  <- data.frame(range="preR-postR3", 
                                   included="no_2523_no2887" ,
                                   variable="NPQ" ,
                                   target=target_list[tar], 
                                   best_fit=ifelse(as.numeric(anova(full_model_timepoint_status_no2523_no2887, full_model_timepoint2_status_no2523_no2887)$Pr[2]) <= 0.05 , "quadratic", "linear"),                                   
                                   estimate_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[3,1], digits=3),
                                   stderr_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[3,2], digits=3),
                                   p_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[3,5], digits=10),
                                   estimate_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[4,1], digits=3),
                                   stderr_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[4,2], digits=3),
                                   p_timepoint_by_status_in_linear=round(summary(full_model_timepoint_status_no2523_no2887)$coeff[4,5], digits=10),
                                   anova_linear_p=round(anova_timepoint_no2523_no2887[2,8], digits=10),
                                   estimate_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[3,1], digits=3),
                                   stderr_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[3,2], digits=3),
                                   p_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[3,5], digits=10),
                                   estimate_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[5,1], digits=3),
                                   stderr_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[5,2], digits=3),
                                   p_timepoint_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[5,5], digits=10),
                                   estimate_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[6,1], digits=3),
                                   stderr_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[6,2], digits=3),
                                   p_timepoint2_by_status_in_quadratic=round(summary(full_model_timepoint2_status_no2523_no2887)$coeff[6,5], digits=10),
                                   anova_quadratic_p=round(anova_timepoint2_no2523_no2887[2,8], digits=10)
                                   )
                         
    anova_summary_no2523_no2887 <- rbind(anova_summary_no2523_no2887, anova_summary_no2523_no2887_df)
  ######################################################################################################
#plot quadratic

if(anova_summary_no2523_no2887_df$anova_quadratic_p <=0.05){
    fixed_effects <- fixef(full_model_timepoint2_status_no2523_no2887)
    choose_timepoints_2_predict <- seq(0, 3, by = 0.2)
    
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
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/quadratic/no2523_no2887_preR_postR3_associated_target_",target_list[tar],"_quadratic_lmer.pdf"))
    print(ggplot(nulisa_preR_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
      
    )
    
        print(ggplot(nulisa_preR_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
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
    choose_timepoints_2_predict <- seq(0, 3, by = 0.2)
    
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
   
    pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/linear/no2523_no2887_preR_postR3_associated_target_",target_list[tar],"_linear_lmer.pdf"))
    print(ggplot(nulisa_preR_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)
            
    )
    
        print(ggplot(nulisa_preR_postR3_qt_no2523_no2887, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
  
        dev.off()
    }
######################################################################################################
 

}#target


######################################################################################################

    anova_summary_no2523_adj <-  arrange(anova_summary_no2523, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                        mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr)     
    anova_summary_no2523_sorted <- anova_summary_no2523_adj %>% arrange(anova_linear_lfdr)
    #anova_summary_no2523_sorted <- anova_summary_no2523_adj %>% arrange(anova_linear_p)
  
    write.csv(anova_summary_no2523_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/preR_postR3_association_between_target_timepoint_timepoint2_by_status_p_adjusted_no2523.csv"), row.names=F, quote=F)

    anova_summary_no2523_no2887_adj <-  arrange(anova_summary_no2523_no2887, anova_linear_p) %>% 
                         mutate(anova_linear_BH_p = p.adjust(anova_linear_p, method="BH")) %>%
                         mutate(anova_linear_q = qvalue(anova_linear_p)$qvalue)%>%
                         mutate(anova_linear_lfdr = qvalue(anova_linear_p)$lfdr)%>%
                         mutate(anova_quadratic_BH_p = p.adjust(anova_quadratic_p, method="BH")) %>%
                         mutate(anova_quadratic_q = qvalue(anova_quadratic_p)$qvalue) %>%
                         mutate(anova_quadratic_lfdr = qvalue(anova_quadratic_p)$lfdr)     

 anova_summary_no2523_no2887_sorted <- anova_summary_no2523_no2887_adj %>% arrange(anova_linear_lfdr)
  anova_summary_no2523_no2887_sorted <- anova_summary_no2523_no2887_adj %>% arrange(anova_linear_p)
  
    write.csv(anova_summary_no2523_no2887_sorted, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/preR_postR3_association_between_target_timepoint_timepoint2_by_status_p_adjusted_no2523_no2887.csv"), row.names=F, quote=F)



  
