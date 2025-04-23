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
plot_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/plots"
table_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/tables"
data_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/data"

#table_workdir <- paste0(table_dir,"/lmer_quadratic_models/timepoint_in_days/days_from_R/")
#plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models/timepoint_in_days/days_from_R/")
table_workdir <- paste0(table_dir,"/lmer_quadratic_models/timepoint_ordinal/")
plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models/timepoint_ordinal/")

setwd(paste0(table_workdir))
file_list <-  list.files(pattern = ".csv")
for (f in 1:length(file_list)){
results_lmer <- read.csv(file_list[f], header=T)
output_prefix <- gsub((str_split_i(file_list[f], "\\.", 2)), "",file_list[f] )


##visualization all points
lmer_result_ord <- arrange(results_lmer, target) %>% mutate(rank = seq(1,nrow(.)))



pdf(paste0(plots_workdir,"Manhattan_plots/",output_prefix,"Manhattan.pdf"), height=10, width=15)
# Make the Manhattan plot on the  results dataset
lmer_result_ord$pvalue  <- "NA"   
lmer_result_ord$pvalue <- lmer_result_ord$p_status_in_linear + 0.00000001

lmer_result_ord$dir <- "NA"
lmer_result_ord$dir[lmer_result_ord$estimate_status_in_linear < 0 ] <- "NEG_noncontroller_effect"
lmer_result_ord$dir[lmer_result_ord$estimate_status_in_linear >= 0 ] <- "POS_noncontroller_effect"

axis_set <- lmer_result_ord$rank


manhplot <- ggplot(lmer_result_ord,  aes(x = rank, y = -log10(pvalue), color = target)) +
  geom_point(aes(size=abs(estimate_status_in_linear), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse((-log10(pvalue)>=1.3 & -log10(anova_linear_p)>=1.3) ,as.character(target),'')),hjust=0.5,vjust=0)+
  scale_x_continuous(label = lmer_result_ord$target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(breaks = seq(0, 2, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("NEG_noncontroller_effect" = 17,"POS_noncontroller_effect"= 19))+

    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("significant status linear effect by target")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), axis.title=element_text(size=20), axis.text.y = element_text(size = 20))+
            theme(legend.position="none") 
 
   print(manhplot)
   
   
lmer_result_ord$pvalue  <- "NA"   
lmer_result_ord$pvalue <- lmer_result_ord$p_timepoint_by_status_in_linear + 0.00000001

lmer_result_ord$dir <- "NA"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint_by_status_in_linear < 0 ] <- "NEG_noncontroller_effect"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint_by_status_in_linear >= 0 ] <- "POS_noncontroller_effect"

axis_set <- lmer_result_ord$rank


manhplot <- ggplot(lmer_result_ord,  aes(x = rank, y = -log10(pvalue), color = target)) +
  geom_point(aes(size=abs(estimate_timepoint_by_status_in_linear), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse((-log10(pvalue)>=1.3 & -log10(anova_linear_p)>=1.3) ,as.character(target),'')),hjust=0.5,vjust=0)+
  scale_x_continuous(label = lmer_result_ord$target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(breaks = seq(0, 2, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("NEG_noncontroller_effect" = 17,"POS_noncontroller_effect"= 19))+

    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("significant timepoint:status linear effect by target")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), axis.title=element_text(size=20), axis.text.y = element_text(size = 20))+
            theme(legend.position="none") 
 
   print(manhplot)
 # Make the Manhattan plot on the  results dataset
lmer_result_ord$pvalue  <- "NA"   
lmer_result_ord$pvalue <- lmer_result_ord$p_status_in_quadratic + 0.00000001

lmer_result_ord$dir <- "NA"
lmer_result_ord$dir[lmer_result_ord$estimate_status_in_quadratic < 0 ] <- "NEG_noncontroller_effect"
lmer_result_ord$dir[lmer_result_ord$estimate_status_in_quadratic >= 0 ] <- "POS_noncontroller_effect"

axis_set <- lmer_result_ord$rank


manhplot <- ggplot(lmer_result_ord,  aes(x = rank, y = -log10(pvalue), color = target)) +
  geom_point(aes(size=abs(estimate_status_in_quadratic), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse((-log10(pvalue)>=1.3 & -log10(anova_quadratic_p)>=1.3),as.character(target),'')),hjust=0.5,vjust=0)+
  scale_x_continuous(label = lmer_result_ord$target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(limits=c(0, max(-log10(lmer_result_ord$pvalue))),breaks = seq(0, 2, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("NEG_noncontroller_effect" = 17,"POS_noncontroller_effect"= 19))+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("significant status quadratic effect by target")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), axis.title=element_text(size=20), axis.text.y = element_text(size = 20))+
            theme(legend.position="none") 
 
   print(manhplot)

 # Make the Manhattan plot on the  results dataset
lmer_result_ord$pvalue  <- "NA"   
lmer_result_ord$pvalue <- lmer_result_ord$p_timepoint_by_status_in_quadratic + 0.00000001

lmer_result_ord$dir <- "NA"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint_by_status_in_quadratic < 0 ] <- "NEG_noncontroller_effect"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint_by_status_in_quadratic >= 0 ] <- "POS_noncontroller_effect"

axis_set <- lmer_result_ord$rank


manhplot <- ggplot(lmer_result_ord,  aes(x = rank, y = -log10(pvalue), color = target)) +
  geom_point(aes(size=abs(estimate_timepoint_by_status_in_quadratic), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse((-log10(pvalue)>=1.3 & -log10(anova_quadratic_p)>=1.3),as.character(target),'')),hjust=0.5,vjust=0)+
  scale_x_continuous(label = lmer_result_ord$target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(limits=c(0, max(-log10(lmer_result_ord$pvalue))),breaks = seq(0, 2, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("NEG_noncontroller_effect" = 17,"POS_noncontroller_effect"= 19))+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("significant timepoint:status quadratic effect by target")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), axis.title=element_text(size=20), axis.text.y = element_text(size = 20))+
            theme(legend.position="none") 
 
   print(manhplot)

 # Make the Manhattan plot on the  results dataset
lmer_result_ord$pvalue  <- "NA"   
lmer_result_ord$pvalue <- lmer_result_ord$p_timepoint2_by_status_in_quadratic + 0.00000001

lmer_result_ord$dir <- "NA"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint2_by_status_in_quadratic < 0 ] <- "NEG_noncontroller_effect"
lmer_result_ord$dir[lmer_result_ord$estimate_timepoint2_by_status_in_quadratic >= 0 ] <- "POS_noncontroller_effect"

axis_set <- lmer_result_ord$rank


manhplot <- ggplot(lmer_result_ord,  aes(x = rank, y = -log10(pvalue), color = target)) +
  geom_point(aes(size=abs(estimate_timepoint2_by_status_in_quadratic), shape=dir, alpha = 0.35)) +
  geom_text(aes(label=ifelse((-log10(pvalue)>=1.3 & -log10(anova_quadratic_p)>=1.3),as.character(target),'')),hjust=0.5,vjust=0)+
  scale_x_continuous(label = lmer_result_ord$target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(limits=c(0, max(-log10(lmer_result_ord$pvalue))),breaks = seq(0, 2, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("NEG_noncontroller_effect" = 17,"POS_noncontroller_effect"= 19))+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("significant timepoint^2:status quadratic effect by target")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), axis.title=element_text(size=20), axis.text.y = element_text(size = 20))+
            theme(legend.position="none") 
 
   print(manhplot)
 dev.off()
 }
  
