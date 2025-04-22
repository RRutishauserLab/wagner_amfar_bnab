# Install and load readxl package
#install.packages("readxl")
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)

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
t_points <- c("BL0", "BL1", "1w", "preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","ATI10",
              "ATI7",  "ATI11",    "ATI8" ,"ATI9", "preR", "postR1", "postR2",
              "postR3", "postR4")
t_points             

#dataframe
nulisa <- data.frame(nulisa)
head(nulisa)
#basic R
#nulisa_subset <- nulisa[nulisa$timepoint %in% t_points,]
#tidyverse
nulisa_subset <- nulisa %>% dplyr::filter(timepoint %in% t_points)
head(nulisa_subset)

dim(nulisa)
dim(nulisa_subset)

#write subset
write.table(nulisa_subset, file=paste0(table_dir,"/nulisa_seq_timepoints_of_interest.txt"), sep="\t", quote=F, row.names=F)
#there are commas in the Protein names - don't use csv
#write.csv(nulisa_subset, file=paste0(table_dir,"/nulisa_seq_timepoints_of_interest.csv"), quote=F)

#re upload table
#nulisa_subset <- read.table(file=paste0(table_dir,"/nulisa_seq_timepoints_of_interest.txt"),  sep="\t", header=T)

head(nulisa_subset)
nulisa_subset$PID_timepoint <- paste(nulisa_subset$PID, nulisa_subset$timepoint, sep="_")


#example 
nulisa_subset_AGER <- nulisa_subset %>% dplyr::filter(Target == "AGER")
head(nulisa_subset_AGER)
#boxplots
#pdf(paste0(plot_dir,"/nulisa_seq_by_target_across_samples.pdf"))
p <- ggplot(nulisa_subset_AGER, aes(x=as.character(PID), y=NPQ)) + 
  geom_boxplot(aes(fill=Target))
p <- p + facet_wrap( ~ Target, ncol = 3)
p <- p + guides(fill="none")
print(p)

#dev.off()


length(unique(nulisa_subset$Target))
#boxplots

p <- ggplot(nulisa_subset, aes(x=as.character(PID), y=NPQ)) + 
  geom_boxplot(aes(fill=Target))
p <- p + facet_wrap_paginate(~ Target, ncol = 4, nrow = 4)
p <- p + guides(fill="none")

for(i in 1:n_pages(p)){
  p_save <-  p + 
    facet_wrap_paginate(~ Target, ncol = 4, nrow = 4, page = i) 
  ggsave(plot = p_save, filename = paste0(plot_dir,"/boxplot/nulisa_seq_by_target_across_samples_boxplot_page_", i, ".jpg"), height=20, width=20)
}





#boxplots across targets and PID_timepoints
pdf(paste0(plot_dir,"/boxplot/nulisa_seq_all_target_across_PID_timepoints.pdf"), height=5, width=15)
p <- ggplot(nulisa_subset, aes(x=PID_timepoint, y=NPQ)) + 
  geom_boxplot(aes(fill=as.character(PID)))
p <- p + guides(fill="none") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
print(p)
dev.off()

pdf(paste0(plot_dir,"/boxplot/nulisa_seq_all_target_across_PID.pdf"), height=5, width=15)
p2 <- ggplot(nulisa_subset, aes(x=as.character(PID), y=NPQ)) + 
  geom_boxplot(aes(fill=as.character(PID)))
p2 <- p2 + guides(fill="none") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
print(p2)
dev.off()

#2887

potential_outlier <- nulisa_subset %>% dplyr::filter(PID == unique(controllers)[1])
table(potential_outlier$SampleQC)

#removing WARN
potential_outlier %>% dplyr::filter(!SampleQC == "WARN") %>% dim()

nulisa_subset_noOut <- nulisa_subset %>% dplyr::filter(!SampleQC == "WARN")

#boxplots across targets and PID_timepoints
pdf(paste0(plot_dir,"/boxplot/nulisa_seq_all_target_across_PID_timepoints_noOut.pdf"), height=5, width=15)
p3 <- ggplot(nulisa_subset_noOut, aes(x=PID_timepoint, y=NPQ)) + 
  geom_boxplot(aes(fill=as.character(PID)))
p3 <- p3 + guides(fill="none") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
print(p3)
dev.off()


pdf(paste0(plot_dir,"/boxplot/nulisa_seq_all_target_across_PID_timepoints_2887_OutOnly.pdf"), height=5, width=15)
nulisa_subset_2887_Out <- nulisa_subset %>% dplyr::filter(!(PID == "2887" & SampleQC == "PASS"))

p4 <- ggplot(nulisa_subset_2887_Out, aes(x=PID_timepoint, y=NPQ)) + 
  geom_boxplot(aes(fill=as.character(PID)))
p4 <- p4 + guides(fill="none") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
print(p4)
dev.off()

