#install.packages("readxl")
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(tidyr)




setwd("~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/data/")
plot_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/nulisa/plots_08292025"
table_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/nulisa/tables_08292025"
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

nulisa_subset$status <- relevel(factor(nulisa_subset$status), ref="non_controllers")
head(nulisa_subset)  


#comparing one timepoint between controllers and non controllers overall and target by target
t_test_table <- t.test(data=nulisa_subset, NPQ ~ status )
p <- ggplot(nulisa_subset, aes(x=status, y=NPQ)) + 
  geom_boxplot(aes(fill=as.factor(status)))+
  scale_fill_manual(values = c("#Ff0000", "#56B4E9"))
p <- p + guides(fill="none") + 
  annotate("text", label=paste0("t.test.p = ", round(t.test(data=nulisa_subset, NPQ ~ status )$p.value, digits=4)),x=1.5, y=15 )
p %>%   ggsave(., filename = paste0(plot_dir,"/boxplot/nulisa_seq_",timepoint_of_interest,"_by_status_across_target.jpg"), height=5, width=5)


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
t_test_all_nulisa_ord<- t_test_all_nulisa %>% mutate(adj_t_test_p= p.adjust(t_test_p, method="BH")) %>% arrange(t_test_p)
write.csv(t_test_all_nulisa_ord, file=paste0(table_dir,"/timepoint_by_status/nulisa_seq_t_test_",timepoint_of_interest,"_by_status_by_target.csv"), row.names=F, quote=F)
#A plot for each target by status in multipanel and multipages
p <- ggplot(nulisa_subset, aes(x=as.factor(status), y=NPQ, fill=Target)) + 
  geom_boxplot(aes(fill=as.factor(status)))+
  scale_fill_manual(values = c("#Ff0000", "#56B4E9"))+
  scale_color_manual(values = c("#Ff0000", "#56B4E9"))
p <- p + facet_wrap_paginate(~Target, ncol = 4, nrow = 4)+
         theme(text = element_text(size = 30),
         axis.text.x = element_text(angle = 90, hjust = 1)) 
p <- p + guides(fill="none")

for(i in 1:n_pages(p)){
  p_save <-  p + 
    facet_wrap_paginate(~ Target, ncol = 4, nrow = 4, page = i) 
  ggsave(plot = p_save, filename = paste0(plot_dir,"/boxplot/nulisa_seq_",timepoint_of_interest,"_by_status_by_target_page_", i, ".jpg"), height=30, width=30)
}

label_nulisa_sign <- t_test_all_nulisa_ord %>% 
                           dplyr::filter(t_test_p <=0.05)
t_test_all_nulisa_sign <- t_test_all_nulisa_ord %>% 
                           dplyr::filter(t_test_p <=0.05)%>% 
                            dplyr::select(Target) %>%
                            pull()



nulisa_subset_sign <-nulisa_subset  %>% dplyr::filter(Target %in% t_test_all_nulisa_sign)                         
rownames(nulisa_subset_sign) <- seq(1,nrow(nulisa_subset_sign))



                     
 p2 <- ggplot(nulisa_subset_sign, aes(x=as.factor(status), y=NPQ, fill=Target)) + 
  geom_boxplot(aes(fill=as.factor(status)))+scale_fill_manual(values = c("#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9"))+
  scale_color_manual(values = c("#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9","#Ff0000", "#56B4E9"))
p2 <- p2 + facet_wrap_paginate(~Target, ncol = 3, nrow = 2)+
         theme(text = element_text(size = 30),
         axis.text.x = element_text(angle = 90, hjust = 1)) + 
         labs(x = " ") +
         geom_text(size = 7,
         data    = label_nulisa_sign,
         mapping = aes(x = -Inf, y = -Inf, label = paste0("t-test p =", round(t_test_p, digits=3))),
         hjust   = -0.1,
         vjust   = -1) +
         guides(fill="none")

p2 %>%
  ggsave(., filename = paste0(plot_dir,"/boxplot/nulisa_seq_",timepoint_of_interest,"_by_status_by_target_significant_t_test.jpg"), height=12, width=12)



#print(p)
