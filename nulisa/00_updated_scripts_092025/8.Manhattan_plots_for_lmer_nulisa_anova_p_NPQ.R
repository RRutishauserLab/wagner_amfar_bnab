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
if (!dir.exists(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots"))) {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots"), recursive=T)}  


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
plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/")
if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots")))
{dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots"), recursive=T)}  

setwd(paste0(table_workdir))
file_list <-  list.files(pattern = ".csv")
file_list <- file_list[!grepl("list.csv", file_list)]
file_list <- file_list[!grepl("knot", file_list)]
for (f in 1:length(file_list)){
results_lmer <- read.csv(file_list[f], header=T)
output_prefix <- gsub((str_split_i(file_list[f], "\\.", 2)), "",file_list[f] )
strsplit(output_prefix, "_")

##visualization all points

##visualization all points#visualization all points
lmer_result_ord <- arrange(results_lmer, target) %>% mutate(rank = seq(1,nrow(.)))

pdf(paste0(plots_workdir,"Manhattan_plots/",output_prefix,"Manhattan.pdf"), height=8, width=20)
# Make the Manhattan plot on the  results dataset
lmer_result_ord$pvalue  <- "NA"   
lmer_result_ord$pvalue <- lmer_result_ord$anova_linear_p + 0.00000001

lmer_result_ord$dir <- "NA"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint_by_status_in_linear < 0 ] <- "POS_noncontroller_effect"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint_by_status_in_linear >= 0 ] <- "POS_controller_effect"
lmer_result_ord<- unique(lmer_result_ord)

#axis_set <- lmer_result_ord$rank

axis_set <- lmer_result_ord$rank

lmer_result_ord$group = "NA"

lmer_result_ord$group[(-log10(lmer_result_ord$pvalue)>=1.3 & (-log10(lmer_result_ord$p_timepoint_by_status_in_linear)>=1.3) & lmer_result_ord$dir %in% "POS_controller_effect"   & lmer_result_ord$anova_linear_BH_p >0.05)] <-  "POS_controllers_effect_p_0.05" 
lmer_result_ord$group[(-log10(lmer_result_ord$pvalue)>=1.3 & (-log10(lmer_result_ord$p_timepoint_by_status_in_linear)>=1.3) & lmer_result_ord$dir %in% "POS_noncontroller_effect"   & lmer_result_ord$anova_linear_BH_p >0.05)] <-  "POS_noncontrollers_effect_p_0.05" 
lmer_result_ord$group[(-log10(lmer_result_ord$pvalue)>=1.3 & (-log10(lmer_result_ord$p_timepoint_by_status_in_linear)>=1.3) & lmer_result_ord$dir %in% "POS_controller_effect"   & lmer_result_ord$anova_linear_BH_p <=0.05)] <-  "POS_controllers_effect_BH_0.05" 
lmer_result_ord$group[(-log10(lmer_result_ord$pvalue)>=1.3 & (-log10(lmer_result_ord$p_timepoint_by_status_in_linear)>=1.3) & lmer_result_ord$dir %in% "POS_noncontroller_effect"   & lmer_result_ord$anova_linear_BH_p <=0.05)] <-  "POS_noncontrollers_effect_BH_0.05" 
lmer_result_ord$group[(-log10(lmer_result_ord$pvalue)>=1.3 & (-log10(lmer_result_ord$p_timepoint_by_status_in_linear)<1.3) & lmer_result_ord$dir %in% "POS_controller_effect"   & lmer_result_ord$anova_linear_BH_p >0.05)] <-  "NS_POS_controllers_effect" 
lmer_result_ord$group[(-log10(lmer_result_ord$pvalue)>=1.3 & (-log10(lmer_result_ord$p_timepoint_by_status_in_linear)<1.3) & lmer_result_ord$dir %in% "POS_noncontroller_effect"   & lmer_result_ord$anova_linear_BH_p >0.05)] <-  "NS_POS_noncontrollers_effect" 
lmer_result_ord$group[(-log10(lmer_result_ord$pvalue)<1.3 )] <-  "NS_model" 



manhplot <- ggplot(lmer_result_ord,  aes(x = rank, y = -log10(pvalue), color = group)) +
  geom_point(aes(size=abs(estimate_timepoint_by_status_in_linear), shape=dir)) + #, alpha = 0.35)) +
  #geom_text(aes(label=ifelse((-log10(pvalue)>=1.3 & -log10(p_timepoint_by_status_in_linear)>=1.3) ,as.character(target),'')),hjust=0.5,vjust=-0.5)+
  geom_text(aes(label=ifelse((-log10(pvalue)>=1.3) ,as.character(target),'')),hjust=0.5,vjust=-0.9, show.legend=F)+
  scale_x_continuous(label = lmer_result_ord$target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(limits=c(0, max(-log10(lmer_result_ord$pvalue)) + 1), breaks = seq(0, max(-log10(lmer_result_ord$pvalue))+1 , by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("POS_noncontroller_effect" = 17,"POS_controller_effect"= 19))+
    scale_color_manual(values = c("POS_controllers_effect_p_0.05" = "#56B4E9", "POS_noncontrollers_effect_p_0.05" = "#F57965", "POS_controllers_effect_BH_0.05" = "#2C5985", "POS_noncontrollers_effect_BH_0.05" =  "#D63A3A", "NS_model" = "grey81", "NS_POS_controllers_effect" =  "grey61", "NS_POS_noncontrollers_effect" =  "grey61")) + 


    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("linear models - timepoint:status effects")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
         theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), axis.title=element_text(size=20), axis.text.y = element_text(size = 15))+
            theme(legend.position="bottom") 
 
   print(manhplot)
   
dev.off()

 }
