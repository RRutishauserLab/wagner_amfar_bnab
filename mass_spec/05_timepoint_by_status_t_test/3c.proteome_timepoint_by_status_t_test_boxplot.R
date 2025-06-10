##Removes all objects from workspace
rm(list=ls())
#install.packages("readxl")
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(tidyr)

setwd("~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/data/")
plot_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/proteomics/plots"
table_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/results/proteomics/tables"
data_dir <- "~/Gladstone Dropbox/Michela Traglia/GB-JW-1447/data"

#uploading the proteome seq dataset
proteome <- read.csv(paste0(data_dir,"/imputed_log2_transformed_protein_intensities_v2.csv"), check.names=F)

proteome_long <- proteome %>% pivot_longer(
  cols = `3555696563_8092_ATI1`:`3555696339_8012_ATI5`, 
  names_to = "id_PID_timepoint",
  values_to = "log2_quantification") %>%
  as.data.frame()

#splitting the id_PID_timepoint in the column and create a new variable PID_timepoint
proteome_long <- proteome_long %>% 
  separate(id_PID_timepoint, c('id', 'PID', 'timepoint')) %>% 
  mutate(PID_timepoint= paste0(PID,"_", timepoint ))

unique(proteome_long$timepoint)

# Make proteome_long into data frame
# proteome_long <- data.frame(proteome_long)
# head(proteome_long)
timepoint_of_interest = "BL0"
proteome_subset <- proteome_long %>% dplyr::filter(timepoint %in% timepoint_of_interest)
head(proteome_subset)

#Create a column of controllers / non controllers
non_controllers<- c("1292", "1408", "8106")
controllers <- proteome_subset$PID[!proteome_subset$PID %in% non_controllers]
no_rebound <- proteome_subset$PID[proteome_subset$PID %in% "2523"]
controllers_rebound <- controllers[ !controllers %in% no_rebound]

proteome_subset <- proteome_subset %>%
  mutate(status = case_when(PID %in% non_controllers ~ 'non_controllers',
                            PID %in% controllers ~ 'controllers')) %>%
  mutate(rebound = case_when( PID %in% no_rebound ~ 'no' , PID %in% c(controllers_rebound,non_controllers) ~ 'yes'
  ))
proteome_subset <-  proteome_subset %>% 
  mutate(Genes_protein_accessions = paste0(Genes, "_",protein_accessions ))



#Double-check that new columns were added successfully
head(proteome_subset)  

#comparing BL0 between controllers and non controllers overall and gene by gene
t_test_table <- t.test(data=proteome_subset, log2_quantification ~ status )
p <- ggplot(proteome_subset, aes(x=status, y=log2_quantification)) + 
  geom_boxplot(aes(fill=as.character(status)))
p <- p + guides(fill="none") + 
  annotate("text", label=paste0("p-value= ", round(t.test(data=proteome_subset, log2_quantification ~ status )$p.value, digits=15)),x=1.5, y=15 )
print(p) %>%   ggsave(., filename = paste0(plot_dir,"/boxplot/mass_spec_",timepoint_of_interest,"_by_status_across_target.jpg"), height=5, width=5)


#all targets one by one
#for each target t-test 
t_test_all_proteome <- data.frame()

target_list <- unique(proteome_subset$Genes_protein_accessions)
skipped_targets<- data.frame()
for (tar in 1:length(target_list)){
proteome_subset_tar <- proteome_subset %>% dplyr::filter(Genes_protein_accessions %in% target_list[tar])

      if(var(proteome_subset_tar$log2_quantification) ==0)
    { 
    print(paste0("skipped protein ", target_list[tar]))
    skipped_individual_target <- data.frame(protein_w_zero_variance = target_list[tar])
    skipped_targets <- rbind(skipped_targets, skipped_individual_target)
 
    }
  if(var(proteome_subset_tar$log2_quantification) > 0){

t_test_one_target<- data.frame(Genes_protein_accessions=target_list[tar], mean_controllers=round(t.test(data=proteome_subset_tar, log2_quantification ~ status )$estimate[1], digits=4), mean_noncontrollers=round(t.test(data=proteome_subset_tar, log2_quantification ~ status )$estimate[2], digits=4) ,t_test_p=round(t.test(data=proteome_subset_tar, log2_quantification ~ status )$p.value, digits=6) )
rownames(t_test_one_target)<- tar
t_test_all_proteome <- rbind(t_test_all_proteome,t_test_one_target)
}
}
#include the different adjustment
t_test_all_proteome_ord<- t_test_all_proteome %>% mutate(adj_t_test_p= p.adjust(t_test_p)) %>% arrange(t_test_p)
write.csv(t_test_all_proteome_ord, file=paste0(table_dir,"/timepoint_by_status/proteome_seq_t_test_",timepoint_of_interest,"_by_status_by_target.csv"), row.names=F, quote=F)

t_test_all_proteome_sign <- t_test_all_proteome_ord %>% 
                           dplyr::filter(t_test_p <=0.05)%>% 
                            dplyr::select(Genes_protein_accessions) %>%
                            pull()
                            
proteome_subset_sign <-proteome_subset  %>% dplyr::filter(Genes_protein_accessions %in% t_test_all_proteome_sign)                         
p2 <- ggplot(proteome_subset_sign, aes(x=as.character(status), y=log2_quantification, fill=Genes_protein_accessions, color=as.character(status))) + 
  geom_boxplot(aes(fill=as.character(status)))
p2 <- p2 + facet_wrap_paginate(~Genes_protein_accessions, ncol = 5, nrow = 5)+
         theme(text = element_text(size = 30),legend.position = "none",
         axis.text.x = element_text(angle = 90, hjust = 1)) 
p2 <- p2 + guides(fill="none") 

for(i in 1:n_pages(p2)){
  p2_save <-  p2 + 
    facet_wrap_paginate(~ Genes_protein_accessions, ncol = 5, nrow = 5, page = i) 
  ggsave(plot = p2_save, filename = paste0(plot_dir,"/boxplot/proteome_seq_",timepoint_of_interest,"_by_status_by_target_significant_t_test_", i, ".jpg"), height=30, width=40)
}


print(p2)



