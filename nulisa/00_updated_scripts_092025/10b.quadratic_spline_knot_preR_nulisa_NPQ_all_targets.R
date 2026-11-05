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
library(stringr)



setwd("~/Dropbox (Gladstone)/GB-JW-1447/data/")
plot_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/plots_08292025"
table_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/tables_08292025"
data_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/data"


#table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/")
#plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/")
#if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/Manhattan_plots")))
# {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/Manhattan_plots"), recursive=T)}  

#table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/")
#plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/")
#if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/Manhattan_plots")))
#{dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/Manhattan_plots"), recursive=T)}  

#table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/")
#plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/")
#if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/Manhattan_plots")))
#{dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/Manhattan_plots"), recursive=T)}  

table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/")

setwd(paste0(table_workdir))
file_list <-  list.files(pattern = ".csv")
file_list <- file_list[!grepl("list.csv", file_list)]
file_list <- file_list[!grepl("knot", file_list)]
file_list <- file_list[grepl("preATI_postR3", file_list)]

f=2
#for (f in 1:length(file_list)){
results_lmer <- read.csv(file_list[f], header=T)
output_prefix <- gsub((str_split_i(file_list[f], "\\.", 2)), "",file_list[f] )

target_best_fit_quadratic <- results_lmer %>% 
                             # dplyr::filter(best_fit == "quadratic") %>%
                                 dplyr::select(target) %>%
                                 pull()
                                 



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
non_controllers<- c("35933", "42894", "31092")
controllers <- nulisa_subset$PID[!nulisa_subset$PID %in% non_controllers]
no_rebound <- nulisa_subset$PID[nulisa_subset$PID %in% "60610"]
controllers_rebound <- controllers[ !controllers %in% no_rebound]

nulisa_subset <- nulisa_subset %>%
   mutate(status = case_when(PID %in% non_controllers ~ 'non_controllers',
                           PID %in% controllers ~ 'controllers')) %>%
   mutate(rebound = case_when( PID %in% no_rebound ~ 'no' , PID %in% c(controllers_rebound,non_controllers) ~ 'yes'
                              ))

nulisa_subset$status <- relevel(factor(nulisa_subset$status), ref="non_controllers")
head(nulisa_subset)  

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

#for each target extract the timepoint
tar = 1
for (tar in 1:length(target_best_fit_quadratic)){
nulisa_target <- nulisa_all_t %>%  dplyr::filter(Target %in% target_best_fit_quadratic[tar])

 
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
   
  

  
  nulisa_target_NPQ_df <- data.frame(target=target_best_fit_quadratic[tar], 
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




   #preATI-postR3 - include all the individuals


   
t_points <- c("preATI", "ATI1", "ATI2",  "ATI3","ATI4", "preR", "postR1", "postR2",  "postR3")

 
    
    nulisa_preATI_postR3 <- nulisa_target_all_pid %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID)
    nulisa_preATI_postR3_all_target <-   nulisa_target_all_proteins %>% dplyr::filter(timepoint %in%t_points)  %>% arrange(PID) #transform timepoint in quantitative
    #diferent datasets
    if(f==1){
    nulisa_preATI_postR3_qt <- 
    nulisa_preATI_postR3 %>% mutate(timepoint =  factor(timepoint, levels=c("preATI", "ATI1"  , "ATI2"  , "ATI3"   ,"ATI4"  , "preR","postR1" ,"postR2" ,"postR3"  )) ,
                                    timepoint = as.numeric(timepoint)#,
                                    #timepoint = timepoint - 1
                                    )
                                    }
    
    if(f==2){
    nulisa_preATI_postR3_qt <- 
    nulisa_preATI_postR3 %>% mutate(timepoint =  factor(timepoint, levels=c("preATI", "ATI1"  , "ATI2"  , "ATI3"   ,"ATI4"  , "preR","postR1" ,"postR2" ,"postR3"  )) ,
                                    timepoint = as.numeric(timepoint)#,
                                    #timepoint = timepoint - 1
                                    ) %>% 
                                    dplyr::filter(PID != "48758") %>% 
                                    dplyr::filter(PID != "60610")
                                 }
   
    if(f==3){
    nulisa_preATI_postR3_qt <- 
    nulisa_preATI_postR3 %>% mutate(timepoint =  factor(timepoint, levels=c("preATI", "ATI1"  , "ATI2"  , "ATI3"   ,"ATI4"  , "preR","postR1" ,"postR2" ,"postR3"  )) ,
                                    timepoint = as.numeric(timepoint)#,
                                    #timepoint = timepoint - 1
                               ) %>% dplyr::filter(PID != "60610")
                                }     
    
      if(f==4){
    nulisa_preATI_postR3_qt <- 
    nulisa_preATI_postR3 %>% mutate(timepoint =  factor(timepoint, levels=c("preATI", "ATI1"  , "ATI2"  , "ATI3"   ,"ATI4"  , "preR","postR1" ,"postR2" ,"postR3"  )) ,
                                    timepoint = as.numeric(timepoint)#,
                                    #timepoint = timepoint - 1
                                )%>% 
                                    dplyr::filter(PID != "48758") 
                                    } 

#fit spline regression model

library(lme4)
library(splines)
library(mgcv)

nulisa_preATI_postR3_qt$timepoint <- as.factor(nulisa_preATI_postR3_qt$timepoint)
# Generate spline basis for predictor x
#nulisa_preATI_postR3_qt$spline_basis <- bs(nulisa_preATI_postR3_qt$timepoint, knots = c(0,5))

###quadratic with preR (timepoint=5) is the 6th timepoint and is set as internal knot
nulisa_preATI_postR3_qt$spline_basis <- bs(as.numeric(nulisa_preATI_postR3_qt$timepoint), degree = 2, knots = 6)
#head(nulisa_preATI_postR3_qt)
#  target  PID timepoint          status    deltaNPQ Age_enrollment Highest_viral_load_ever ART_start Rebound_group setpoint postR1 postR2 postR3 spline_basis.1 spline_basis.2 spline_basis.3
#preATI", "ATI1"  , "ATI2"  , "ATI3"   ,"ATI4"  , "preR","postR1" ,"postR2" ,"postR3"
#matplot(x, B2, type = "l", lty = 1, xlab="x", ylab="Basis functions")
#legend("topright", legend = paste0("B",1:ncol(B2)), col = 1:ncol(B2), lty=1)
#nulisa_preATI_postR3_qt$spline_basis <- bs(as.numeric(nulisa_preATI_postR3_qt$timepoint), degree = 3, Boundary.knots=c(1,6))
#fit spline regression model
min_model_timepoint <- lmer(data=nulisa_preATI_postR3_qt, NPQ ~ (1|PID) + spline_basis )


full_model_timepoint_status <- lmer(data=nulisa_preATI_postR3_qt, NPQ ~ (1|PID) + 
                                                                                spline_basis*status )

anova_timepoint_status  <- anova(min_model_timepoint,  full_model_timepoint_status)


#view summary of spline regression model
summary(full_model_timepoint_status)

 anova_summary_df <- data.frame(range="preATI-postR3", 
                                   included=ifelse(f==1 , "all", ifelse(f==2, "no60610_no48758", ifelse(f==3, "no60610", "no48758"))) ,
                                   variable="NPQ" ,
                                   reference=  "non_controllers",
                                   target=target_best_fit_quadratic[tar],
                                   AIC=anova_timepoint_status$AIC[2],
                                   estimate_status_in_spline=round(summary(full_model_timepoint_status)$coeff[5,1], digits=3),
                                   stderr_status_in_spline=round(summary(full_model_timepoint_status)$coeff[5,2], digits=3),
                                   p_status_in_spline=round(summary(full_model_timepoint_status)$coeff[5,5], digits=10),
                                   estimate_spline_basis1_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[6,1], digits=7),
                                   stderr_spline_basis1_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[6,2], digits=3),
                                   p_spline_basis1_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[6,5], digits=10),
                                   estimate_spline_basis2_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[7,1], digits=7),
                                   stderr_spline_basis2_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[7,2], digits=3),
                                   p_spline_basis2_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[7,5], digits=10),
                                   estimate_spline_basis3_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[8,1], digits=7),
                                   stderr_spline_basis3_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[8,2], digits=3),
                                   p_spline_basis3_statuscontrollers=round(summary(full_model_timepoint_status)$coeff[8,5], digits=10),
                                   anova_quadratic_spline_p=round(anova_timepoint_status[2,8], digits=10)
                                   )
                                   
anova_summary <- rbind(anova_summary, anova_summary_df)
    
 
 
 #calculate predictions using spline regression model
#x_lim <- range(as.numeric(nulisa_preATI_postR3_qt$timepoint))
#x_grid <- seq(x_lim[1], x_lim[2])
#preds <- predict(full_model_timepoint_status, newdata=list(x=seq(0, 9, by = 0.2)))

fixed_effects <- fixef(full_model_timepoint_status)
    choose_timepoints_2_predict <- seq(1, 9, by = 0.2)
   # choose_timepoints_2_predict <- seq(0, 8, by = 0.2)
    #bs_res <- bs(choose_timepoints_2_predict, deg = 3, Boundary.knots=c(1,6))
     bs_res <- bs(choose_timepoints_2_predict, deg = 2,  knots = 6)
controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                        bs_res[,1], ##spline_basis1 term
                                        bs_res[,2], ##spline_basis2 term
                                        bs_res[,3], ##spline_basis3 term
                                       # bs_res[,4], ##spline_basis4 term

                                        rep(1, length(choose_timepoints_2_predict)), ##non-controller status term
                                        bs_res[,1], ##spline_basis1:timepoint term
                                        bs_res[,2], ##spline_basis2:timepoint term
                                        bs_res[,3]#, ##spline_basis3:timepoint term
                                       # bs_res[,4] ##spline_basis4:timepoint term
                                        )
                                        
                                        
                                        
    colnames(controller_pred_matrix) <- c("intercept",
                                              "spline_basis1",
                                              "spline_basis2",
                                              "spline_basis3",
                                              #"spline_basis4",
                                              "statuscontroller",
                                              "spline_basis1:statusnon_controllers",
                                              "spline_basis2:statusnon_controllers",
                                              "spline_basis3:statusnon_controllers"#,
                                              #"spline_basis4:statusnon_controllers"
                                              )
                                              

    
    non_controller_pred_matrix <- cbind(rep(1, length(choose_timepoints_2_predict)), ##intercept term
                                    bs_res[,1], ##spline_basis1 term
                                    bs_res[,2], ##spline_basis2 term
                                    bs_res[,3], ##spline_basis3 term
                                   # bs_res[,4], ##spline_basis4 term
                                    rep(0, length(choose_timepoints_2_predict)), ##non-controller status term
                                    rep(0, length(choose_timepoints_2_predict)),##spline_basis1:timepoint term
                                    rep(0, length(choose_timepoints_2_predict)),##spline_basis2:timepoint term
                                    rep(0, length(choose_timepoints_2_predict))#, ##spline_basis3:timepoint term,
                                    #rep(0, length(choose_timepoints_2_predict)) ##spline_basis4:timepoint term
                                    )
    colnames(non_controller_pred_matrix) <-  c("intercept",
                                              "spline_basis1",
                                              "spline_basis2",
                                              "spline_basis3",
                                              #"spline_basis4",
                                              "statuscontroller",
                                              "spline_basis1:statusnon_controllers",
                                              "spline_basis2:statusnon_controllers",
                                              "spline_basis3:statusnon_controllers"#,
                                              #"spline_basis4:statusnon_controllers"
                                              )
                                              

    
    predicted_mean_noncontroller <- non_controller_pred_matrix %*% fixed_effects
    predicted_mean_controller <- controller_pred_matrix %*% fixed_effects
    
    varcor <- vcov(full_model_timepoint_status)

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
levels(nulisa_preATI_postR3_qt$timepoint) <- c(levels(nulisa_preATI_postR3_qt$timepoint),  "preATI", "ATI1"  , "ATI2"  , "ATI3"   ,"ATI4"  , "preR","postR1" ,"postR2" ,"postR3"  )

nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==1] <- "preATI"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==2] <- "ATI1"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==3] <- "ATI2"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==4] <- "ATI3"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==5] <- "ATI4"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==6] <- "preR"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==7] <- "postR1"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==8] <- "postR2"
nulisa_preATI_postR3_qt$timepoint[nulisa_preATI_postR3_qt$timepoint ==9] <- "postR3"
droplevels(nulisa_preATI_postR3_qt$timepoint)
  
  if (!dir.exists(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/spline_all_targets"))) {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/spline_all_targets"), recursive=T)}  

 
   pdf(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/spline_all_targets/",ifelse(f==1 , "all", ifelse(f==2, "no60610_no48758", ifelse(f==3, "no60610", "no48758"))),"_preATI_postR3_ALL_TARGETS_",target_best_fit_quadratic[tar],"_spline_lmer.pdf"), width=8, height=7)
    print(ggplot(nulisa_preATI_postR3_qt, aes(timepoint, NPQ, color = status)) +
            geom_point(size = 3) +
            scale_color_manual(values=c( "#B30000","#247EB3")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se,  fill = status), alpha = 0.4) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)+
            scale_fill_manual(values=c( "#247EB3","#B30000")) +
            theme_classic()  +
            theme(legend.position = "none")+
            theme(  axis.title = element_text(size = 18),   # axis title size
            axis.text  = element_text(size = 16)    # axis tick label size
  )               
    )
    
        print(ggplot(nulisa_preATI_postR3_qt, aes(timepoint, NPQ, color = status)) +
            geom_point(size = 3) +
            scale_color_manual(values=c( "#B30000","#247EB3")) + 
            geom_ribbon(data = predicted_values, mapping = aes(ymin = NPQ - NPQ_se, ymax = NPQ + NPQ_se,  fill = status), alpha = 0.4) +
            geom_line(data = predicted_values,  mapping = aes(x = timepoint, y=NPQ), linewidth = 1)+
            scale_fill_manual(values=c( "#247EB3","#B30000")) +
            theme_bw()    +
            theme(legend.position = "none")+
            theme(
            axis.title = element_text(size = 18),   # axis title size
            axis.text  = element_text(size = 16)    # axis tick label size
  )          
    )
    
    print(ggplot(nulisa_preATI_postR3_qt, aes(timepoint, NPQ, color = status, shape=as.factor(ART_start))) +
            geom_point(size = 3, stroke = 1.2) +
            scale_shape_manual(values = c(1, 2, 3)) +
            scale_color_manual(values=c("#B30000","#247EB3")) +
             theme_bw()+ 
            theme(legend.position = "none")+
            theme(
            axis.title = element_text(size = 18),   # axis title size
            axis.text  = element_text(size = 16)    # axis tick label size
  )          
    )
  
    dev.off()



} #tar
    print(unique(nulisa_preATI_postR3_qt$PID))   


 anova_summary_adj <-  anova_summary %>% mutate(anova_quadratic_spline_BH_p = p.adjust(anova_quadratic_spline_p, method="BH")) 
               
 if (!dir.exists(paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/spline_all_targets"))) {dir.create(paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/spline_all_targets"), recursive=T)}  
 
    write.csv(anova_summary_adj, file=paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/spline_all_targets/preATI_postR3_association_between_ALL_TARGETS_spline_preR_as_knot_p_adjusted_",ifelse(f==1 , "all", ifelse(f==2, "no60610_no48758", ifelse(f==3, "no60610", "no48758"))),".csv"), row.names=F, quote=F)



} #f





              
