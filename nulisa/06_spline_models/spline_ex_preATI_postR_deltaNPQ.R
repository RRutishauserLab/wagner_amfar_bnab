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
library(gam)


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
#covariatesForUnivariate$ART_start <- as.factor(covariatesForUnivariate$ART_start)
#covariatesForUnivariate$ART_start <- relevel(covariatesForUnivariate$ART_start, ref="EARLY")
#covariatesForUnivariate$ART_start <- as.numeric(covariatesForUnivariate$ART_start)
#covariatesForUnivariate$ART_start <- covariatesForUnivariate$ART_start - 1

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
anova_summary <- data.frame()
anova_summary_no2523 <- data.frame()
anova_summary_no2523_no2887 <- data.frame()
#for each target extract the timepoint

for (tar in 1:length(target_list)){

#tar <- which(target_list %in% c("IL12B"))
#tar <- which(target_list %in% c("TNFRSF1A"))
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



   #preATI-postR - include all the individuals


    t_points <- c("preATI", "ATI1", "ATI2",  "ATI3","ATI4", "preR", "postR1", "postR2",  "postR3")
    
    nulisa_preATI_postR3 <- nulisa_BL0_adj_NPQ_target_all_pid %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)
    nulisa_preATI_postR3_all_target <-   nulisa_BL0_adj_NPQ_all_proteins %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID) #transform timepoint in quantitative
    
    nulisa_preATI_postR3_qt <- 
    nulisa_preATI_postR3 %>% mutate(timepoint =  factor(timepoint, levels=c("preATI", "ATI1"  , "ATI2"  , "ATI3"   ,"ATI4"  , "preR","postR1" ,"postR2" ,"postR3"  )) ,
                                    timepoint = as.numeric(timepoint),
                                    timepoint = timepoint - 1
                                    )


#fit spline regression model

library(lme4)
library(splines)
library(mgcv)
nulisa_preATI_postR3_qt$timepoint <- as.factor(nulisa_preATI_postR3_qt$timepoint)
# Generate spline basis for predictor x
#nulisa_preATI_postR3_qt$spline_basis <- bs(nulisa_preATI_postR3_qt$timepoint, knots = c(0,5))
###cubic
nulisa_preATI_postR3_qt$spline_basis_cubic  <- bs(as.numeric(nulisa_preATI_postR3_qt$timepoint), degree = 3)
#head(nulisa_preATI_postR3_qt)
#  target  PID timepoint          status    deltaNPQ Age_enrollment Highest_viral_load_ever ART_start Rebound_group setpoint postR1 postR2 postR3 spline_basis.1 spline_basis.2 spline_basis.3

#nulisa_preATI_postR3_qt$spline_basis <- bs(as.numeric(nulisa_preATI_postR3_qt$timepoint), degree = 3, Boundary.knots=c(1,6))
#fit spline regression model
min_model_timepoint_cubic <- lmer(data=nulisa_preATI_postR3_qt, deltaNPQ ~ (1|PID) + 
                                                                                spline_basis_cubic )


full_model_timepoint_status_cubic <- lmer(data=nulisa_preATI_postR3_qt, deltaNPQ ~ (1|PID) + 
                                                                                spline_basis_cubic*status )

anova_timepoint_status_cubic  <- anova(min_model_timepoint_cubic,  full_model_timepoint_status_cubic)


#quadratic
nulisa_preATI_postR3_qt$spline_basis_quadratic <- bs(as.numeric(nulisa_preATI_postR3_qt$timepoint), degree = 2)
#nulisa_preATI_postR3_qt$spline_basis <- bs(as.numeric(nulisa_preATI_postR3_qt$timepoint), degree = 3, Boundary.knots=c(1,6))
#fit spline regression model
min_model_timepoint_quadratic <- lmer(data=nulisa_preATI_postR3_qt, deltaNPQ ~ (1|PID) + 
                                                                                spline_basis_quadratic )


full_model_timepoint_status_quadratic <- lmer(data=nulisa_preATI_postR3_qt, deltaNPQ ~ (1|PID) + 
                                                                                spline_basis_quadratic*status )

anova_timepoint_status_quadratic  <- anova(min_model_timepoint_quadratic,  full_model_timepoint_status_quadratic)

if((anova_timepoint_status_cubic$`Pr(>Chisq)`[2] <= 0.05) | (anova_timepoint_status_quadratic$`Pr(>Chisq)`[2] <= 0.05) )
{


AIC_models <- aictab(cand.set = list(full_model_timepoint_status_cubic, full_model_timepoint_status_quadratic), modnames = c("cubic", "quadratic"))
AIC_summary <- data.frame(Target = target_list[tar] , Models = AIC_models$Modnames, AIC = AIC_models$AICc, p_anova= rbind(anova_timepoint_status_quadratic$`Pr(>Chisq)`[2],anova_timepoint_status_cubic$`Pr(>Chisq)`[2] ))
AIC_quadratic <- AIC_summary %>% filter(Models == "quadratic") %>% select(AIC) %>% pull()
AIC_cubic <- AIC_summary %>% filter(Models == "cubic") %>% select(AIC) %>% pull()


write.csv(AIC_summary, file=paste0(table_dir,"/spline_models/deltaNPQ_preATI_postR3/all_preATI_postR3_associated_target_",target_list[tar],"_spline_lmer.csv"), quote=F)

if(AIC_quadratic <= AIC_cubic){
fixed_effects <- fixef(full_model_timepoint_status_quadratic)
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
     bs_res <- bs(choose_timepoints_2_predict, deg = 2)
non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        bs_res[,1], ##spline_basis1 term
                                        bs_res[,2], ##spline_basis2 term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        bs_res[,1], ##spline_basis1:timepoint term
                                        bs_res[,2] ##spline_basis2:timepoint term
                                        )
                                        
                                        
                                        
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "spline_basis_quadratic1",
                                              "spline_basis_quadratic2",
                                              "statusnoncontroller",
                                              "spline_basis_quadratic1:statusnon_controllers",
                                              "spline_basis_quadratic2:statusnon_controllers")
                                              

    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                    bs_res[,1], ##spline_basis1 term
                                    bs_res[,2], ##spline_basis2 term
                                    rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                    rep(0, length(choose_timepoints_2_predict)),##spline_basis1:timepoint term
                                    rep(0, length(choose_timepoints_2_predict))##spline_basis2:timepoint term
                                    )
    colnames(controller_pred_matrix) <-  c("intercept",
                                              "spline_basis_quadratic1",
                                              "spline_basis_quadratic2",
                                              "statusnoncontroller",
                                              "spline_basis_quadratic1:statusnon_controllers",
                                              "spline_basis_quadratic2:statusnon_controllers")
                                              

    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint_status_quadratic)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  deltaNPQ = predicted_mean_noncontroller,
                                                  deltaNPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  deltaNPQ = predicted_mean_controller,
                                                  deltaNPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/spline_models/deltaNPQ_preATI_postR3/quadratic/all_preATI_postR3_associated_target_",target_list[tar],"_spline_lmer.pdf"))
    g1<- ggplot(nulisa_preATI_postR3_qt, aes(timepoint, deltaNPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = deltaNPQ - deltaNPQ_se, ymax = deltaNPQ + deltaNPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=deltaNPQ), linewidth = 1)
            
    print(g1)
    
      
    g2 <- ggplot(nulisa_preATI_postR3_qt, aes(timepoint, deltaNPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000"))
    
 print(g2)
  
 dev.off()
 
 }
 
 if(AIC_cubic < AIC_quadratic)
 {
 
fixed_effects <- fixef(full_model_timepoint_status_cubic)
    choose_timepoints_2_predict <- seq(0, 9, by = 0.2)
     bs_res <- bs(choose_timepoints_2_predict, deg = 3)
non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        bs_res[,1], ##spline_basis1 term
                                        bs_res[,2], ##spline_basis2 termspline_test:

                                        bs_res[,3], ##spline_basis3 term
                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        bs_res[,1], ##spline_basis1:timepoint term
                                        bs_res[,2], ##spline_basis2:timepoint term
                                        bs_res[,3] ##spline_basis3:timepoint term
                                        )
                                        
                                        
                                        
    colnames(non_controller_pred_matrix) <- c("intercept",
                                              "spline_basis_cubic1",
                                              "spline_basis_cubic2",
                                              "spline_basis_cubic3",
                                              "statusnoncontroller",
                                              "spline_basis_cubic1:statusnon_controllers",
                                              "spline_basis_cubic2:statusnon_controllers",
                                              "spline_basis_cubic3:statusnon_controllers")
                                              

    
    controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                    bs_res[,1], ##spline_basis1 term
                                    bs_res[,2], ##spline_basis2 term
                                    bs_res[,3], ##spline_basis3 term
                                    rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                    rep(0, length(choose_timepoints_2_predict)),##spline_basis1:timepoint term
                                    rep(0, length(choose_timepoints_2_predict)),##spline_basis2:timepoint term
                                    rep(0, length(choose_timepoints_2_predict)) ##spline_basis3:timepoint term
                                    )
    colnames(controller_pred_matrix) <-  c("intercept",
                                              "spline_basis_cubic1",
                                              "spline_basis_cubic2",
                                              "spline_basis_cubic3",
                                              "statusnoncontroller",
                                              "spline_basis_cubic1:statusnon_controllers",
                                              "spline_basis_cubic2:statusnon_controllers",
                                              "spline_basis_cubic3:statusnon_controllers")
                                              

    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint_status_cubic)

    predicted_se_noncontroller <- vector(mode = "numeric")
    predicted_se_controller <- vector(mode = "numeric")
    for(t in 1:length(choose_timepoints_2_predict)) {
      predicted_se_noncontroller[t] <- sqrt(t(non_controller_pred_matrix[t,]) %*% varcor %*% (non_controller_pred_matrix[t,])) %>% as.numeric(.)
      predicted_se_controller[t] <- sqrt(t(controller_pred_matrix[t,]) %*% varcor %*% (controller_pred_matrix[t,])) %>% as.numeric(.)
    }
    
    predicted_values_noncontrollers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  deltaNPQ = predicted_mean_noncontroller,
                                                  deltaNPQ_se = predicted_se_noncontroller,
                                                  status = rep("non_controllers", length(choose_timepoints_2_predict)))
    
    predicted_values_controllers <- data.frame(timepoint = choose_timepoints_2_predict,
                                                  deltaNPQ = predicted_mean_controller,
                                                  deltaNPQ_se = predicted_se_controller,
                                                  status = rep("controllers", length(choose_timepoints_2_predict)))
    
    predicted_values <- rbind(predicted_values_noncontrollers,
                              predicted_values_controllers)
   
    pdf(paste0(plot_dir,"/spline_models/deltaNPQ_preATI_postR3/cubic/all_preATI_postR3_associated_target_",target_list[tar],"_spline_lmer.pdf"))
    g1 <- (ggplot(nulisa_preATI_postR3_qt, aes(timepoint, deltaNPQ, color = status)) +
            geom_point() +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = deltaNPQ - deltaNPQ_se, ymax = deltaNPQ + deltaNPQ_se), fill = "grey70", alpha = 0.5) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=deltaNPQ), linewidth = 1)
            
    )
    
    print(g1)
    g2 <- (ggplot(nulisa_preATI_postR3_qt, aes(timepoint, deltaNPQ, color = status, shape=as.factor(ART_start))) +
            geom_point() +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#56B4E9", "#Ff0000")) 
    )
    
    print(g2)
  
   dev.off()


}
 
} 
}

    