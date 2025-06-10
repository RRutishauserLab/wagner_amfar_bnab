#install.packages("readxl")
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)

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
nulisa_subset <- nulisa %>% dplyr::filter(timepoint %in% "BL0")
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

#comparing BL0 between controllers and non controllers overall and gene by gene
t_test_table <- t.test(data=nulisa_subset, NPQ ~ status )
p <- ggplot(nulisa_subset, aes(x=status, y=NPQ)) + 
  geom_boxplot(aes(fill=as.character(status)))
p <- p + guides(fill="none") + 
  annotate("text", label=paste0("t.test.p = ", round(t.test(data=nulisa_subset, NPQ ~ status )$p.value, digits=4)),x=1.5, y=15 )
print(p) %>%   ggsave(., filename = paste0(plot_dir,"/boxplot/nulisa_seq_BL0_controllers_noncontrollers_across_target_boxplot.jpg"), height=5, width=5)

##BL0 is sign different in controllers vs non controllers, non controllers showed higher mean overall across targets, BL0 should be included in downstream models.

#all targets



#for each target extract the timepoint


p <- ggplot(nulisa_subset, aes(x=as.character(status), y=NPQ, fill=Target, color=as.character(status))) + 
  geom_boxplot(aes(fill=as.character(status)))
p <- p + facet_wrap(~Target)
p <- p + guides(fill="none")

for(i in 1:n_pages(p)){
  p_save <-  p + 
    facet_wrap_paginate(~ Target, ncol = 4, nrow = 4, page = i) 
  ggsave(plot = p_save, filename = paste0(plot_dir,"/boxplot/nulisa_seq_BL0_controllers_noncontrollers_by_target_boxplot_page_", i, ".jpg"), height=20, width=20)
}


print(p)
