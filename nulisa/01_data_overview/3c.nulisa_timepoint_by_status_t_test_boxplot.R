#install.packages("readxl")
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(tidyr)

setwd("~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/data/")
plot_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/nulisa/plots"
table_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/nulisa/tables"
data_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/data"

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
timepoint_of_interest="BL0"
nulisa_subset <- nulisa %>% dplyr::filter(timepoint %in% timepoint_of_interest)
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

#comparing one timepoint between controllers and non controllers overall and target by target
t_test_table <- t.test(data=nulisa_subset, NPQ ~ status )
p <- ggplot(nulisa_subset, aes(x=status, y=NPQ)) + 
  geom_boxplot(aes(fill=as.character(status)))
p <- p + guides(fill="none") + 
  annotate("text", label=paste0("t.test.p = ", round(t.test(data=nulisa_subset, NPQ ~ status )$p.value, digits=4)),x=1.5, y=15 )
print(p) %>%   ggsave(., filename = paste0(plot_dir,"/boxplot/nulisa_seq_",timepoint_of_interest,"_by_status_across_target.jpg"), height=5, width=5)


#all targets one by one
#for each target t-test 
target_list <- unique(nulisa_subset$Target)

t_test_all_nulisa <- data.frame()
#for each target extract the timepoint
for (tar in 1:length(target_list)){
nulisa_subset_tar <- nulisa_subset %>% dplyr::filter(Target %in% target_list[tar])
t_test_one_target<- data.frame(Target=target_list[tar], mean_controllers=round(t.test(data=nulisa_subset_tar, NPQ ~ status )$estimate[1], digits=4), mean_noncontrollers=round(t.test(data=nulisa_subset_tar, NPQ ~ status )$estimate[2], digits=4) ,t_test_p=round(t.test(data=nulisa_subset_tar, NPQ ~ status )$p.value, digits=6) )
rownames(t_test_one_target)<- tar
t_test_all_nulisa <- rbind(t_test_all_nulisa,t_test_one_target)
}
#include the different adjustment
t_test_all_nulisa_ord<- t_test_all_nulisa %>% mutate(adj_t_test_p= p.adjust(t_test_p)) %>% arrange(t_test_p)
write.csv(t_test_all_nulisa_ord, file=paste0(table_dir,"/timepoint_by_status/nulisa_seq_t_test_",timepoint_of_interest,"_by_status_by_target.csv"), row.names=F, quote=F)
#A plot for each target by status in multipanel and multipages
p <- ggplot(nulisa_subset, aes(x=as.character(status), y=NPQ, fill=Target, color=as.character(status))) + 
  geom_boxplot(aes(fill=as.character(status)))
p <- p + facet_wrap_paginate(~Target, ncol = 4, nrow = 4)+
         theme(text = element_text(size = 30),
         axis.text.x = element_text(angle = 90, hjust = 1)) 
p <- p + guides(fill="none")

for(i in 1:n_pages(p)){
  p_save <-  p + 
    facet_wrap_paginate(~ Target, ncol = 4, nrow = 4, page = i) 
  ggsave(plot = p_save, filename = paste0(plot_dir,"/boxplot/nulisa_seq_",timepoint_of_interest,"_by_status_by_target_page_", i, ".jpg"), height=30, width=30)
}

t_test_all_nulisa_sign <- t_test_all_nulisa_ord %>% 
                           dplyr::filter(t_test_p <=0.05)%>% 
                            dplyr::select(Target) %>%
                            pull()
                            
nulisa_subset_sign <-nulisa_subset  %>% dplyr::filter(Target %in% t_test_all_nulisa_sign)                         
p2 <- ggplot(nulisa_subset_sign, aes(x=as.character(status), y=NPQ, fill=Target, color=as.character(status))) + 
  geom_boxplot(aes(fill=as.character(status)))
p2 <- p2 + facet_wrap_paginate(~Target, ncol = 4, nrow = 4)+
         theme(text = element_text(size = 30),
         axis.text.x = element_text(angle = 90, hjust = 1)) 
p2 <- p2 + guides(fill="none")

for(i in 1:n_pages(p2)){
  p2_save <-  p2 + 
    facet_wrap_paginate(~ Target, ncol = 4, nrow = 4, page = i) 
  ggsave(plot = p2_save, filename = paste0(plot_dir,"/boxplot/nulisa_seq_",timepoint_of_interest,"_by_status_by_target_significant_t_test_", i, ".jpg"), height=30, width=30)
}


print(p)
