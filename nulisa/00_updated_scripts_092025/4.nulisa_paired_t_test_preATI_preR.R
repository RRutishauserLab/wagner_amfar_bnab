#install.packages("readxl")

#For the comparisons between 2 timepoints, Rachel was pointing out that the VL for 1292 is quite a bit higher than anyone else at the postReb1 timepoint. So we wanted to make sure it wasn't just them driving preR vs postR1 differences?
#Is there a way to include in the graphs which participants/timepoints are being used for each comparison (just so we know for sure it is the correct ones)? Are you basing this off the table in the power point I had sent in January? I'
#m sorry if this messes anything up but we had originally planned to exclude some (8106 preATI, 8092 preATI) cause we thought there was a sample mix-up but it turns out there wasn't. So for the comparisons between 2 timepoints, the only participants/timepoints that should be excluded are 8106 ATI1 and 2523 preR + postReb1



rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(magrittr)

setwd("~/Dropbox (Gladstone)/GB-JW-1447/data/")
plot_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/plots_08292025"
table_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/tables_08292025"
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


nulisa_subset$status <- relevel(factor(nulisa_subset$status), ref="non_controllers")
head(nulisa_subset)  

#comparing timepoints between controllers and non controllers gene by gene
#Comparisons: Baseline0 vs 1w, Baseline0 vs preATI, preATI vs ATI1, preATI vs preR, preR vs postR1
target <-unique(nulisa_subset$Target)

t1="preATI"
t2="preR"

paired.w.test <- data.frame()
paired.t.test <- data.frame()
linear_models <- data.frame()
paired.t.test_noOut <- data.frame()
linear_models_noOut <- data.frame()
for (i in 1:length(target)){
nulisa_subset_paired <- nulisa_subset %>% dplyr::filter(timepoint %in% c(t1,t2)) %>% 
dplyr::filter(Target %in% target[i]) 

#if missing data - remove the pair
#nulisa_subset_paired %>% dplyr::filter(duplicated(PID))


nulisa_subset_paired <- nulisa_subset_paired %>% dplyr::filter(!(PID %in% (((nulisa_subset_paired %>%  group_by(PID) %>%
    dplyr::filter(n() ==1))$PID)) ))

t1_table <- nulisa_subset_paired %>% dplyr::filter(timepoint %in% t1) %>% arrange(PID)
t2_table <- nulisa_subset_paired %>% dplyr::filter(timepoint %in% t2) %>% arrange(PID)

differences <- (t1_table$NPQ - t2_table$NPQ)

#perform shapiro-wilk test for normality on this vector of values
#if(shapiro.test(differences)$p.value <=0.05){
pvalue <- t.test(t1_table$NPQ, t2_table$NPQ, paired = TRUE)$p.value
estimate <- t.test(t1_table$NPQ, t2_table$NPQ, paired = TRUE)$estimate
CI_L <- t.test(t1_table$NPQ, t2_table$NPQ, paired = TRUE)$conf.int[1]
CI_U <- t.test(t1_table$NPQ, t2_table$NPQ, paired = TRUE)$conf.int[2]
pvalue_summary <- data.frame(Target= target[i], t1_mean=mean(t1_table$NPQ ), t2_mean=mean(t2_table$NPQ ) , t1_t2_mean_diff = estimate , t1_t2_mean_diff_95CI_L=CI_L , t1_t2_mean_diff_95CI_U=CI_U ,  pvalue=pvalue, t1_t2_pair= paste0(t1,"_",t2), test="paired_t_test")
paired.t.test <- rbind(paired.t.test, pvalue_summary)
#}
#if(shapiro.test(differences)$p.value >0.05){
#pvalue <- wilcox.test(nulisa_subset_paired$NPQ[nulisa_subset_paired$timepoint == t1], nulisa_subset_paired$NPQ[nulisa_subset_paired$timepoint == t2], paired = TRUE)$p.value

#pvalue_summary <- data.frame(Target= target[i], pvalue=pvalue, pair= paste0(t1,"_",t2), test="paired_wilcoxon_test")
#paired.w.test <- rbind(paired.w.test, pvalue_summary)
#}

###check for controller and non controller status effects using a linear model modeling the differences between the timepoints 
#lm(diff ~ status)
diff_table = data.frame(PID=t1_table$PID, status=t1_table$status, diff=differences)
linear_model_target <- lm(data=diff_table, diff ~ status )

linear_model_target_df <-  data.frame(Target= target[i], 
t1_t2_diff= paste0(t1,"_", t2, "_diff"),
intercept_beta= summary(linear_model_target)$coeff[1,1], 
intercept_p= summary(linear_model_target)$coeff[1,4],
contr_vs_noncontrollers_beta= summary(linear_model_target)$coeff[2,1], 
contr_vs_noncontrollers_p= summary(linear_model_target)$coeff[2,4], 
test="linear_model",
ref="non_controllers")
linear_models <- rbind(linear_models, linear_model_target_df)




###removing sample 2887

nulisa_subset_paired_noOut <- nulisa_subset %>% dplyr::filter(timepoint %in% c(t1,t2)) %>% 
dplyr::filter(Target %in% target[i]) %>% dplyr::filter(!PID == "2887")

#if missing data - remove the pair
#nulisa_subset_paired %>% dplyr::filter(duplicated(PID))


nulisa_subset_paired_noOut <- nulisa_subset_paired_noOut %>% dplyr::filter(!(PID %in% (((nulisa_subset_paired_noOut %>%  group_by(PID) %>%
    dplyr::filter(n() ==1))$PID)) ))

t1_table_noOut <- nulisa_subset_paired_noOut %>% dplyr::filter(timepoint %in% t1) %>% arrange(PID)
t2_table_noOut <- nulisa_subset_paired_noOut %>% dplyr::filter(timepoint %in% t2) %>% arrange(PID)

differences_noOut <- (t1_table_noOut$NPQ - t2_table_noOut$NPQ)

#perform shapiro-wilk test for normality on this vector of values
#if(shapiro.test(differences)$p.value <=0.05){
pvalue_noOut <- t.test(t1_table_noOut$NPQ, t2_table_noOut$NPQ, paired = TRUE)$p.value
estimate_noOut <- t.test(t1_table_noOut$NPQ, t2_table_noOut$NPQ, paired = TRUE)$estimate
CI_L_noOut <- t.test(t1_table_noOut$NPQ, t2_table_noOut$NPQ, paired = TRUE)$conf.int[1]
CI_U_noOut <- t.test(t1_table_noOut$NPQ, t2_table_noOut$NPQ, paired = TRUE)$conf.int[2]
pvalue_summary_noOut <- data.frame(Target= target[i], t1_mean=mean(t1_table_noOut$NPQ ), t2_mean=mean(t2_table_noOut$NPQ ),t1_t2_mean_diff = estimate_noOut , t1_t2_mean_diff_95CI_L=CI_L_noOut , t1_t2_mean_diff_95CI_U=CI_U_noOut ,  pvalue=pvalue_noOut, t1_t2_pair= paste0(t1,"_",t2), test="paired_t_test")
paired.t.test_noOut <- rbind(paired.t.test_noOut, pvalue_summary_noOut)
###check for controller and non controller status effects using a linear model modeling the differences between the timepoints 
#lm(diff ~ status)
diff_table_noOut = data.frame(PID=t1_table_noOut$PID, status=t1_table_noOut$status, diff=differences_noOut)
linear_model_target_noOut <- lm(data=diff_table_noOut, diff ~ status )

linear_model_target_df_noOut <-  data.frame(Target= target[i], 
t1_t2_diff= paste0(t1,"_", t2, "_diff"),
intercept_beta= summary(linear_model_target_noOut)$coeff[1,1], 
intercept_p= summary(linear_model_target_noOut)$coeff[1,4],
contr_vs_noncontrollers_beta= summary(linear_model_target_noOut)$coeff[2,1], 
contr_vs_noncontrollers_p= summary(linear_model_target_noOut)$coeff[2,4],
test="linear_model",
ref="non_controllers")
linear_models_noOut <- rbind(linear_models_noOut, linear_model_target_df_noOut)

}

#with outlier
#using only t test for the small sample size
combined_par_nonpar_test <- rbind(paired.t.test, paired.w.test  )
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
 mutate(p_adj_BH = p.adjust(pvalue, method="BH"))
write.csv(combined_par_nonpar_test_ord_noOut, file=paste0(table_dir, "/paired_t_test/overall/paired_t_test_targets_",t1,"_",t2,"_noOut_PID_2887.csv"), quote=F, row.names=F)


#linear models 
linear_models_ord_noOut <- arrange(linear_models_noOut, intercept_p) %>% mutate(p_adj_BH = p.adjust(pvalue, method="BH"))
write.csv(linear_models_ord_noOut, file=paste0(table_dir, "/paired_t_test/by_status/paired_t_test_targets_",t1,"_",t2,"_status_effect_noOut_PID_2887.csv"), quote=F, row.names=F)

#####

##visualization all points
combined_par_nonpar_test_ord<-read.csv(paste0(table_dir, "/paired_t_test/overall/paired_t_test_targets_",t1,"_",t2,".csv"), header=T)
mean_diff_ord <- arrange(combined_par_nonpar_test_ord, Target) %>% mutate(rank = seq(1,nrow(.)))

# Make the Manhattan plot on the  results dataset
mean_diff_ord$pvalue <- mean_diff_ord$pvalue + 0.00000001

mean_diff_ord$dir <- "NA"
mean_diff_ord$dir[mean_diff_ord$t1_t2_mean_diff < 0 ] <- "t2>t1"
mean_diff_ord$dir[mean_diff_ord$t1_t2_mean_diff >= 0 ] <- "t1>t2"

axis_set <- mean_diff_ord$rank

mean_diff_ord$group = "NA"

mean_diff_ord$group[(-log10(mean_diff_ord$pvalue)>=1.3 & mean_diff_ord$dir %in% "t1>t2"   & mean_diff_ord$p_adj_BH >0.05)] <-  "t1>t2_p_0.05" 
mean_diff_ord$group[(-log10(mean_diff_ord$pvalue)>=1.3 & mean_diff_ord$dir %in% "t2>t1"   & mean_diff_ord$p_adj_BH >0.05)] <-  "t2>t1_p_0.05" 
mean_diff_ord$group[(-log10(mean_diff_ord$pvalue)>=1.3 & mean_diff_ord$dir %in% "t1>t2"   & mean_diff_ord$p_adj_BH <=0.05)] <-  "t1>t2_BH_0.05" 
mean_diff_ord$group[(-log10(mean_diff_ord$pvalue)>=1.3 & mean_diff_ord$dir %in% "t2>t1"   & mean_diff_ord$p_adj_BH <=0.05)] <-  "t2>t1_BH_0.05" 
mean_diff_ord$group[(-log10(mean_diff_ord$pvalue)<1.3 )] <-  "not_sign" 


# Add a conditional group column


pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,".pdf"), height=5, width=22)
manhplot <- ggplot(mean_diff_ord,  aes(x = rank, y = -log10(pvalue), color = group)) +
  geom_point(aes(size=abs(t1_t2_mean_diff), shape=dir)) + #, alpha = 0.35
  geom_text(aes(label=ifelse(-log10(pvalue)>=1.3,as.character(Target),'')),hjust=0.5,vjust=-0.9, show.legend = FALSE)+
  scale_x_continuous(label = mean_diff_ord$Target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(limits=c(0, max(-log10(mean_diff_ord$pvalue)) + 1), breaks = seq(0, max(-log10(mean_diff_ord$pvalue)) + 1, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("t2>t1" = 17,"t1>t2"= 19))+
  #scale_color_manual(values = c("POS_mean_diff_p_0.05" = "#Ff0000", "NEG_mean_diff_p_0.05" = "#56B4E9", "POS_mean_diff_BH_0.05" = "#D63A3A", "NEG_mean_diff_BH_0.05" =  "#2C5985", "not_sign" = "grey81")) + 
  scale_color_manual(values = c("t1>t2_p_0.05" = "grey51", "t2>t1_p_0.05" = "grey51", "t1>t2_BH_0.05" = "black", "t2>t1_BH_0.05" =  "black", "not_sign" = "grey91")) + 
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0(t1,"-",t2, " mean differences")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
         theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), axis.title=element_text(size=20), axis.text.y = element_text(size = 15))+
        theme(legend.position =  "bottom")
   print(manhplot)
  dev.off()
  
#
##visualization noOut 
combined_par_nonpar_test_ord_noOut<-read.csv(paste0(table_dir, "/paired_t_test/overall/paired_t_test_targets_",t1,"_",t2,"_noOut_PID_2887.csv"), header=T)
mean_diff_ord_noOut <- arrange(combined_par_nonpar_test_ord_noOut, Target) %>% mutate(rank = seq(1,nrow(.)))

# Make the Manhattan plot on the  results dataset
mean_diff_ord_noOut$pvalue <- mean_diff_ord_noOut$pvalue + 0.00000001

mean_diff_ord_noOut$dir <- "NA"
mean_diff_ord_noOut$dir[mean_diff_ord_noOut$t1_t2_mean_diff < 0 ] <- "t2>t1"
mean_diff_ord_noOut$dir[mean_diff_ord_noOut$t1_t2_mean_diff >= 0 ] <- "t1>t2"

axis_set <- mean_diff_ord_noOut$rank
mean_diff_ord_noOut$group = "NA"

mean_diff_ord_noOut$group[(-log10(mean_diff_ord_noOut$pvalue)>=1.3 & mean_diff_ord_noOut$dir %in% "t2>t1"   & mean_diff_ord_noOut$p_adj_BH >0.05)] <-  "t1>t2_p_0.05" 
mean_diff_ord_noOut$group[(-log10(mean_diff_ord_noOut$pvalue)>=1.3 & mean_diff_ord_noOut$dir %in% "t1>t2"   & mean_diff_ord_noOut$p_adj_BH >0.05)] <-  "t2>t1_p_0.05" 
mean_diff_ord_noOut$group[(-log10(mean_diff_ord_noOut$pvalue)>=1.3 & mean_diff_ord_noOut$dir %in% "t2>t1"   & mean_diff_ord_noOut$p_adj_BH <=0.05)] <-  "t1>t2_BH_0.05" 
mean_diff_ord_noOut$group[(-log10(mean_diff_ord_noOut$pvalue)>=1.3 & mean_diff_ord_noOut$dir %in% "t1>t2"   & mean_diff_ord_noOut$p_adj_BH <=0.05)] <-  "t2>t1_BH_0.05" 
mean_diff_ord_noOut$group[(-log10(mean_diff_ord_noOut$pvalue)<1.3 )] <-  "not_sign" 



pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,"_noOut_PID_2887.pdf"), height=5, width=22)
manhplot <- ggplot(mean_diff_ord_noOut,  aes(x = rank, y = -log10(pvalue), color = group)) +
  geom_point(aes(size=abs(t1_t2_mean_diff), shape=dir)) + #, alpha = 0.35)) +
  geom_text(aes(label=ifelse(-log10(pvalue)>=1.3,as.character(Target),'')),hjust=0.5,vjust=-0.9, show.legend = FALSE)+
  scale_x_continuous(label = mean_diff_ord_noOut$Target, breaks = seq(from = 1, by = 1, length.out = length(axis_set))
  ) +
  scale_y_continuous(limits=c(0, max(-log10(mean_diff_ord_noOut$pvalue)) + 1), breaks = seq(0, max(-log10(mean_diff_ord_noOut$pvalue)) + 1, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("t2>t1" = 17,"t1>t2"= 19))+
  #scale_color_manual(values = c("POS_mean_diff_p_0.05" = "#Ff0000", "NEG_mean_diff_p_0.05" = "#56B4E9", "POS_mean_diff_BH_0.05" = "#D63A3A", "NEG_mean_diff_BH_0.05" =  "#2C5985", "not_sign" = "grey81")) + 
  scale_color_manual(values = c("t1>t2_p_0.05" = "grey51", "t2>t1_p_0.05" = "grey51", "t1>t2_BH_0.05" = "black", "t2>t1_BH_0.05" =  "black", "not_sign" = "grey91")) + 
      theme(legend.position="none") + 
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0(t1,"-",t2, " mean differences")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), axis.title=element_text(size=20), axis.text.y = element_text(size = 15))+
        theme(legend.position="bottom") 
   print(manhplot)
  dev.off()
  
  
###linear models effects STATUS
linear_models_ord <- read.csv(paste0(table_dir, "/paired_t_test/by_status/paired_t_test_targets_",t1,"_",t2,"_status_effect.csv"))
non_controller_effect_ord <- arrange(linear_models_ord, Target) %>% mutate(rank = seq(1,nrow(.)))

# Make the Manhattan plot on the  results dataset
non_controller_effect_ord$contr_vs_noncontrollers_p <- non_controller_effect_ord$contr_vs_noncontrollers_p + 0.00000001

non_controller_effect_ord$dir <- "NA"
non_controller_effect_ord$dir[non_controller_effect_ord$contr_vs_noncontrollers_beta < 0 ] <- "POS_noncontrollers_effect"
non_controller_effect_ord$dir[non_controller_effect_ord$contr_vs_noncontrollers_beta >= 0 ] <- "POS_controllers_effect"

axis_set <- non_controller_effect_ord$rank

non_controller_effect_ord$group = "NA"

non_controller_effect_ord$group[(-log10(non_controller_effect_ord$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord$dir %in% "POS_controllers_effect"   & non_controller_effect_ord$p_adj_BH >0.05)] <-  "POS_controllers_effect_p_0.05" 
non_controller_effect_ord$group[(-log10(non_controller_effect_ord$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord$dir %in% "POS_noncontrollers_effect"   & non_controller_effect_ord$p_adj_BH >0.05)] <-  "POS_noncontrollers_effect_p_0.05" 
non_controller_effect_ord$group[(-log10(non_controller_effect_ord$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord$dir %in% "POS_controllers_effect"   & non_controller_effect_ord$p_adj_BH <=0.05)] <-  "POS_controllers_effect_BH_0.05" 
non_controller_effect_ord$group[(-log10(non_controller_effect_ord$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord$dir %in% "POS_noncontrollers_effect"   & non_controller_effect_ord$p_adj_BH <=0.05)] <-  "POS_noncontrollers_effect_BH_0.05" 
non_controller_effect_ord$group[(-log10(non_controller_effect_ord$contr_vs_noncontrollers_p)<1.3 )] <-  "not_sign" 



pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,"_by_status.pdf"), height=5, width=22)
manhplot <- ggplot(non_controller_effect_ord,  aes(x = rank, y = -log10(contr_vs_noncontrollers_p), color = group)) +
  geom_point(aes(size=abs(contr_vs_noncontrollers_beta), shape=dir))+ #, alpha = 0.35)) +
  geom_text(aes(label=ifelse(-log10(contr_vs_noncontrollers_p)>=1.3,as.character(Target),'')),hjust=0.5,vjust=-0.9, show.legend = FALSE)+
  scale_x_continuous(label = non_controller_effect_ord$Target, breaks = seq(from = 1, by = 1, length.out = length(non_controller_effect_ord$Target))
  ) +
  scale_y_continuous(limits=c(0, max(-log10(non_controller_effect_ord$contr_vs_noncontrollers_p)) + 1), breaks = seq(0, max(-log10(non_controller_effect_ord$contr_vs_noncontrollers_p)) + 1, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("POS_noncontrollers_effect" = 17,"POS_controllers_effect"= 19))+
    scale_color_manual(values = c("POS_controllers_effect_p_0.05" = "#56B4E9", "POS_noncontrollers_effect_p_0.05" = "#E67E7E", "POS_controllers_effect_BH_0.05" = "#2C5985", "POS_noncontrollers_effect_BH_0.05" =  "darkred2", "not_sign" = "grey81")) + 

    #"#Ff0000"
    
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0(t1,"-",t2, " mean differences - controllers vs non_controllers effects")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
         theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), axis.title=element_text(size=20), axis.text.y = element_text(size = 15)) +
            theme(legend.position="bottom") 
 
   print(manhplot)
  dev.off()


###linear models effects no Out
linear_models_ord_noOut <- read.csv(paste0(table_dir, "/paired_t_test/by_status/paired_t_test_targets_",t1,"_",t2,"_status_effect_noOut_PID_2887.csv"))
non_controller_effect_ord_noOut <- arrange(linear_models_ord_noOut, Target) %>% mutate(rank = seq(1,nrow(.)))

# Make the Manhattan plot on the  results dataset
non_controller_effect_ord_noOut$contr_vs_noncontrollers_p <- non_controller_effect_ord_noOut$contr_vs_noncontrollers_p + 0.00000001

non_controller_effect_ord_noOut$dir <- "NA"
non_controller_effect_ord_noOut$dir[non_controller_effect_ord_noOut$contr_vs_noncontrollers_beta < 0 ] <- "POS_noncontrollers_effect"
non_controller_effect_ord_noOut$dir[non_controller_effect_ord_noOut$contr_vs_noncontrollers_beta >= 0 ] <- "POS_controllers_effect"

axis_set <- non_controller_effect_ord_noOut$rank
non_controller_effect_ord$group = "NA"

non_controller_effect_ord_noOut$group[(-log10(non_controller_effect_ord_noOut$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord_noOut$dir %in% "POS_controllers_effect"   & non_controller_effect_ord_noOut$p_adj_BH >0.05)] <-  "POS_controllers_effect_p_0.05" 
non_controller_effect_ord_noOut$group[(-log10(non_controller_effect_ord_noOut$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord_noOut$dir %in% "POS_noncontrollers_effect"   & non_controller_effect_ord_noOut$p_adj_BH >0.05)] <-  "POS_noncontrollers_effect_p_0.05" 
non_controller_effect_ord_noOut$group[(-log10(non_controller_effect_ord_noOut$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord_noOut$dir %in% "POS_controllers_effect"   & non_controller_effect_ord_noOut$p_adj_BH <=0.05)] <-  "POS_controllers_effect_BH_0.05" 
non_controller_effect_ord_noOut$group[(-log10(non_controller_effect_ord_noOut$contr_vs_noncontrollers_p)>=1.3 & non_controller_effect_ord_noOut$dir %in% "POS_noncontrollers_effect"   & non_controller_effect_ord_noOut$p_adj_BH <=0.05)] <-  "POS_noncontrollers_effect_BH_0.05" 
non_controller_effect_ord_noOut$group[(-log10(non_controller_effect_ord_noOut$contr_vs_noncontrollers_p)<1.3 )] <-  "not_sign" 


pdf(paste0(plot_dir,"/paired_t_test/Manhattan_mean_diff_",t1,"_",t2,"_noOut_PID_2887_by_status.pdf"), height=5, width=22)
manhplot <- ggplot(non_controller_effect_ord_noOut,  aes(x = rank, y = -log10(contr_vs_noncontrollers_p), color = group)) +
  geom_point(aes(size=abs(contr_vs_noncontrollers_beta), shape=dir)) + #, alpha = 0.35)) +
  geom_text(aes(label=ifelse(-log10(contr_vs_noncontrollers_p)>=1.3,as.character(Target),'')),hjust=0.5,vjust=-0.9, show.legend = FALSE)+
  scale_x_continuous(label = non_controller_effect_ord_noOut$Target, breaks = seq(from = 1, by = 1, length.out = length(non_controller_effect_ord_noOut$Target))
  ) +
    scale_color_manual(values = c("POS_controllers_effect_p_0.05" = "#56B4E9", "POS_noncontrollers_effect_p_0.05" = "#E67E7E", "POS_controllers_effect_BH_0.05" = "#2C5985", "POS_noncontrollers_effect_BH_0.05" =  "darkred2", "not_sign" = "grey81")) + 
  scale_y_continuous(limits=c(0, max(-log10(non_controller_effect_ord_noOut$contr_vs_noncontrollers_p)) + 1), breaks = seq(0, max(-log10(non_controller_effect_ord_noOut$contr_vs_noncontrollers_p)) + 1, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("POS_noncontrollers_effect" = 17,"POS_controllers_effect"= 19))+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0(t1,"-",t2, "  mean differences - controllers vs non controllers effects")) +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), axis.title=element_text(size=20), axis.text.y = element_text(size = 15)) +
        theme(legend.position="bottom")  
   print(manhplot)
  dev.off()
