# Install and load readxl package
#install.packages("readxl")
#remove variables
rm( list=ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(factoextra)
library(stringr)

setwd("~/Dropbox (Gladstone)/GB-JW-1447/data/")
plot_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/plots"
table_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/nulisa/tables"
data_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/data"

nulisa_subset <- read.table(file=paste0(table_dir,"/nulisa_seq_timepoints_of_interest.txt"),  sep="\t", header=T)
head(nulisa_subset)
unique(nulisa_subset$timepoint)

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
unique(nulisa_subset$timepoint)
unique(nulisa_subset$PID_timepoint)
length(unique(nulisa_subset$PID_timepoint))
#
nulisa_subset_t <- nulisa_subset %>% 
                   dplyr::select(c(PID_timepoint, NPQ, Target, PID, timepoint, status, rebound)) %>% 
                   pivot_wider(names_from = Target, values_from = NPQ) %>%
                   as.data.frame()

?pivot_wider()

dim(nulisa_subset_t)
head(nulisa_subset_t)

#variance across ssample and timepoints across 250 proteins
target_var <- data.frame(Target=unique(nulisa_subset$Target), var=diag(var(na.omit(nulisa_subset_t[,-c(1:5)]))))
target_var_rank <- target_var[order(target_var$var, decreasing=T),]
head(target_var)
head(target_var_rank)
median_target_var_rank <- median(target_var_rank$var)
median_target_var_rank

perc3_target_var_rank <- quantile(target_var_rank$var, 0.75)
perc3_target_var_rank



#filtering the target above the threshold of var

Target_75perc<- target_var_rank %>% 
  dplyr::filter(var >= perc3_target_var_rank) %>% 
  dplyr::select(Target)


most_var_target_75th <- dplyr::select(nulisa_subset_t ,one_of(Target_75perc$Target))
head(most_var_target_75th)

#PCA
res.pca.norm_75th<- prcomp(most_var_target_75th,  scale=T, center = T, rank.=10)
?prcomp()
res.pca.norm_75th


PCA_75th <- data.frame(PID_timepoint=nulisa_subset_t$PID_timepoint, 
                       res.pca.norm_75th$x,
                       PID=nulisa_subset_t$PID, 
                       timepoint=nulisa_subset_t$timepoint,
                       status=nulisa_subset_t$status,
                       rebound=nulisa_subset_t$rebound)
                       

#plot
pdf(paste0(plot_dir, "/pca/nulisa_PCA_var75th_PID_timepoint.pdf"), height=5, width=15)

g1<-fviz_pca_ind(res.pca.norm_75th, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_75th), aes(color = as.character(PID))) #+
  #geom_label(aes(label = PCA_75th$status))
#hjust = .5, vjust = 0, label.size = 0, fill = NA,
#label.padding = unit(0, "pt")
#) 
print(g1)

g2<-fviz_pca_ind(res.pca.norm_75th, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_75th), aes(color = as.character(timepoint))) #+
  #geom_label(aes(label = PCA_75th$status))
#hjust = .5, vjust = 0, label.size = 0, fill = NA,
#label.padding = unit(0, "pt")
#) 
print(g2)

g3<-fviz_pca_ind(res.pca.norm_75th, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_75th), aes(color = as.character(status))) #+
  #geom_label(aes(label = PCA_75th$status))
#hjust = .5, vjust = 0, label.size = 0, fill = NA,
#label.padding = unit(0, "pt")
#) 
print(g3)

dev.off()



perc95_target_var_rank <- quantile(target_var_rank$var, 0.95)
perc95_target_var_rank



#filtering the target above the threshold of var

Target_95perc<- target_var_rank %>% 
  dplyr::filter(var >= perc95_target_var_rank) %>% 
  dplyr::select(Target)


most_var_target_95th <- dplyr::select(nulisa_subset_t ,one_of(Target_95perc$Target))
head(most_var_target_95th)

#PCA
res.pca.norm_95th<- prcomp(most_var_target_95th,  scale=T, center = T, rank.=10)
?prcomp()
res.pca.norm_95th


PCA_95th <- data.frame(PID_timepoint=nulisa_subset_t$PID_timepoint, 
                       res.pca.norm_95th$x,
                       PID=nulisa_subset_t$PID, 
                       timepoint=nulisa_subset_t$timepoint,
                       status=nulisa_subset_t$status,
                       rebound=nulisa_subset_t$rebound)
                       

#plot
pdf(paste0(plot_dir, "/pca/nulisa_PCA_var95th_PID_timepoint.pdf"), height=5, width=8)

g1<-fviz_pca_ind(res.pca.norm_95th, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_95th), aes(color = as.character(PID))) #+
  #geom_label(aes(label = PCA_75th$status))
#hjust = .5, vjust = 0, label.size = 0, fill = NA,
#label.padding = unit(0, "pt")
#) 
print(g1)

g2<-fviz_pca_ind(res.pca.norm_95th, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_95th), aes(color = as.character(timepoint))) #+
  #geom_label(aes(label = PCA_75th$status))
#hjust = .5, vjust = 0, label.size = 0, fill = NA,
#label.padding = unit(0, "pt")
#) 
print(g2)

g3<-fviz_pca_ind(res.pca.norm_95th, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_95th), aes(color = as.character(status))) #+
  #geom_label(aes(label = PCA_75th$status))
#hjust = .5, vjust = 0, label.size = 0, fill = NA,
#label.padding = unit(0, "pt")
#) 
print(g3)


g4<-fviz_pca_ind(res.pca.norm_95th, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_95th), aes(color = as.character(rebound))) #+
#geom_label(aes(label = PCA_75th$status))
#hjust = .5, vjust = 0, label.size = 0, fill = NA,
#label.padding = unit(0, "pt")
#) 
print(g4)

dev.off()



