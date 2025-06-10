#install.packages("readxl")

#Removes all objects from workspace
rm(list=ls())

#Load relevant packages into environment
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(tidyr)
library(factoextra)

#setwd("~/Dropbox (Gladstone)/GB-JW-1447/data/")
#plot_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/proteomics/plots"
#table_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/results/proteomics/tables"
#data_dir <- "~/Dropbox (Gladstone)/GB-JW-1447/data"

#setwd sets working directory where R reads and writes from
#plot_dir and table_dir stores where plot and table outputs are saved
#data_dir specifies where data is found

setwd("/Users/alanxie/Library/CloudStorage/OneDrive-Personal/Documents/2025/Gladstone/GB-JW-1447")
plot_dir <- "/Users/alanxie/Library/CloudStorage/OneDrive-Personal/Documents/2025/Gladstone/GB-JW-1447/Results/Plots"
table_dir <- "/Users/alanxie/Library/CloudStorage/OneDrive-Personal/Documents/2025/Gladstone/GB-JW-1447/Results/Tables"
data_dir <- "/Users/alanxie/Library/CloudStorage/OneDrive-Personal/Documents/2025/Gladstone/GB-JW-1447/Data"

#uploading the proteomics seq dataset using the function that read a csv file, the columns with a number in the name are not automatically read by R - check.names=F overcome the issue
proteome <- read.csv(paste0(data_dir,"/imputed_log2_transformed_protein_intensities.csv"), check.names=F)

#head() displays the first six rows/elements of object
#summary() displays information about each row (length, class, mode -->how an object is stored) 
#and statistical info (min., max., median, mean, 1st and 3rd quartile)

head(proteome)
summary(proteome)
colnames(proteome)

proteome_long <- proteome %>% pivot_longer(
    cols = `3555696563_8092_ATI1`:`3555696339_8012_ATI5`, 
    names_to = "id_PID_timepoint",
    values_to = "log2_quantification") %>%
    as.data.frame()

#splitting the id_PID_timepoint in the column and create a new variable PID_timepoint
proteome_long <- proteome_long %>% 
                      separate(id_PID_timepoint, c('id', 'PID', 'timepoint')) %>% 
                      mutate(PID_timepoint= paste0(PID,"_", timepoint ))

#Unneeded:We already have log transformed
#Create a new column with a function to have log2 transformed values
#proteome_long  <- proteome_long %>% mutate(log2_quantification=log2(transformed_quantification))
#unique(proteome_long$timepoint)



#create a vector of timepoints of interest
t_points <- c("BL0", "BL1", "1w", "preATI", "ATI1", "ATI2",  "ATI3","ATI4", "ATI5","ATI6","ATI10",
              "ATI7",  "ATI11",    "ATI8" ,"ATI9", "preR", "postR1", "postR2",
              "postR3", "postR4")
t_points             

#the table is already a dataframe - check
#proteome <- data.frame(proteome)
head(str(proteome_long))
head(proteome_long)
#basic R
#proteome_subset <- proteome[proteome$timepoint %in% t_points,]

#tidyverse - filter rows
#function we need is filter()
proteome_subset <- proteome_long %>% dplyr::filter(timepoint %in% t_points)


head(proteome_subset)

dim(proteome)
dim(proteome_subset)

#create a column of controllers / non controllers
non_controllers<- c("1292", "1408", "8106")
controllers <- proteome_subset$PID[!proteome_subset$PID %in% non_controllers]
no_rebound <- proteome_subset$PID[proteome_subset$PID %in% "2523"]
controllers_rebound <- controllers[ !controllers %in% no_rebound]

proteome_subset <- proteome_subset %>%
   mutate(status = case_when(PID %in% non_controllers ~ 'non_controllers',
                           PID %in% controllers ~ 'controllers')) %>%
   mutate(rebound = case_when( PID %in% no_rebound ~ 'no' , PID %in% c(controllers_rebound,non_controllers) ~ 'yes'
                              ))

head(proteome_subset)  

#boxplots across targets and PID_timepoints

pdf(paste0(plot_dir,"/boxplot/mass_spec_seq_all_target_across_PID_timepoints.pdf"), height=5, width=15)
p <- ggplot(proteome_subset, aes(x=as.factor(PID_timepoint), y=log2_quantification)) + 
  geom_boxplot(aes(fill=as.factor(PID)))
p <- p + guides(fill="none") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
print(p)
dev.off()

pdf(paste0(plot_dir,"/boxplot/mass_spec__seq_all_target_across_PID.pdf"), height=5, width=15)
p2 <- ggplot(proteome_subset, aes(x=as.character(PID), y=log2_quantification)) + 
  geom_boxplot(aes(fill=as.character(PID)))
p2 <- p2 + guides(fill="none") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
print(p2)
dev.off()

###PCA
##variance across sample and timepoints across 4000 proteins
proteome_subset_t <- proteome_subset %>% 
  dplyr::select(c(PID_timepoint, log2_quantification, Genes)) %>% 
  pivot_wider(names_from = Genes, values_from = log2_quantification) %>%
  as.data.frame()

#Extra step to identify samples with multiple entries
proteome_subset  |>
  dplyr::summarise(n = dplyr::n(), .by = c(PID_timepoint, Genes)) |>
  dplyr::filter(n > 1L) 

dim(proteome_subset_t)
head(proteome_subset_t)

proteome_subset %>% dplyr::filter(Genes == "KIF2A")%>% dplyr::filter(PID_timepoint == "8092_ATI1")

proteome_subset_t <-  proteome_subset %>% 
  mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions )) %>%
  dplyr::select(c(PID_timepoint, log2_quantification, Genes_protein_accessions)) %>% 
  pivot_wider(names_from = Genes_protein_accessions, values_from = log2_quantification) %>%
  as.data.frame()

proteome_subset_t[1:10,1:10]

gene_protein_accession <- proteome_subset %>% 
  mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions )) %>%
  dplyr::select(c( Genes_protein_accessions)) 

target_variance <- data.frame()

for (col in 2:ncol(proteome_subset_t)){
  target_var <- data.frame(Target=colnames(proteome_subset_t[col]), var=var(na.omit(proteome_subset_t[col]))[1,1])
  target_variance <- rbind(target_variance,target_var )
}


target_var_rank <- target_variance[order(target_variance$var, decreasing=T),]
head(target_variance)
head(target_var_rank)
perc3_target_var_rank <- quantile(target_var_rank$var, 0.90)
#perc3_target_var_rank


#filtering the target above the threshold of var
Target_Nperc<- target_var_rank %>% 
  dplyr::filter(var >= perc3_target_var_rank) %>% 
  dplyr::select(Target)


most_var_target_nth <- dplyr::select(proteome_subset_t ,one_of(Target_Nperc$Target))
head(most_var_target_nth)
length(most_var_target_nth)

#Write most_var_target_nth into R 
#write.csv(most_var_target_nth, file = "most_var_target_nth.csv", na = "NA", row.names = FALSE)
#Isolate columns with missing values and output into csv
#missing_var_vals <- most_var_target_nth[,colSums(is.na(most_var_target_nth)) > 0]
#write.csv(missing_var_vals, file = "var_target_nth_missing.csv", na = "NA", row.names = FALSE)

#PCA
res.pca.norm_nth<- prcomp(most_var_target_nth,  scale=T, center = T, rank.=10)
#res.pca.norm_nth

PCA_nth <- data.frame(PID_timepoint=proteome_subset_t$PID_timepoint, res.pca.norm_nth$x)
PCA_nth$PID <- substr(PCA_nth$PID_timepoint, 1, 4)

# three plots colored by status, PID and timepoint
pdf(paste0(plot_dir, "/pca/proteome_PCA_varNth_PID_timepoint.pdf"), height=5, width=15)

#Use fviz_pca_ind, which is a function from factoextra library to create a PCA plot
g1<-fviz_pca_ind(res.pca.norm_nth, geom = "blank") +
  geom_point(data = ~ bind_cols(.x, PCA_nth), aes(color = PID)) 

print(g1)
dev.off()




