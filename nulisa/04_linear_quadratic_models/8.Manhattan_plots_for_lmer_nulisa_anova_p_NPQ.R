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
if (!dir.exists(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots"))) {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots"), recursive=T)}  


table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/")
plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/")
if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/Manhattan_plots")))
 {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_postATI/Manhattan_plots"), recursive=T)}  

#table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/")
#plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/")
#if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/Manhattan_plots")))
#{dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/days_from_R/Manhattan_plots"), recursive=T)}  

#table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/")
#plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/")
#if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/Manhattan_plots")))
# {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_in_days/56days_beforeR_to_70days_postR/Manhattan_plots"), recursive=T)}  

#table_workdir <- paste0(table_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/")
#plots_workdir <- paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/")
#if (!dir.exists(paste0(plots_workdir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots")))
# {dir.create(paste0(plot_dir,"/lmer_quadratic_models_NPQ/timepoint_ordinal/Manhattan_plots"), recursive=T)}  

setwd(paste0(table_workdir))
file_list <-  list.files(pattern = ".csv")
file_list <- file_list[!grepl("list.csv", file_list)]
for (f in 1:length(file_list)){
results_lmer <- read.csv(file_list[f], header=T)
output_prefix <- gsub((str_split_i(file_list[f], "\\.", 2)), "",file_list[f] )


##visualization all points

##visualization all points
lmer_result_ord <- arrange(results_lmer, target) %>% mutate(rank = seq(1,nrow(.)))
lmer_result_ord$tar <- NA
for( m in 1:nrow(lmer_result_ord)){
lmer_result_ord$tar[m] <- strsplit(lmer_result_ord$target, "_")[[m]][1]
}

#gene location
# Install and load the biomaRt package
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define your list of genes
gene_list <- unique(lmer_result_ord$tar)

# Query gene locations
gene_locations <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = gene_list,
  mart = ensembl
)

gene_locations <- gene_locations[!grepl("H", gene_locations$chromosome_name),]

gene_locations$chromosome_name[gene_locations$chromosome_name == "X"] <- 23
gene_locations$chromosome_name[gene_locations$chromosome_name == "Y"] <- 24
gene_locations$chromosome_name[gene_locations$chromosome_name == "MT"] <- 25





pdf(paste0(plots_workdir,"Manhattan_plots/",output_prefix,"Manhattan.pdf"), height=5, width=40)
# Make the Manhattan plot on the  results dataset
lmer_result_ord$pvalue  <- "NA"   
lmer_result_ord$pvalue <- lmer_result_ord$anova_linear_p + 0.00000001

lmer_result_ord$dir <- "NA"
lmer_result_ord$dir[lmer_result_ord$estimate_status_in_linear < 0 ] <- "NEG_noncontroller_effect"
lmer_result_ord$dir[lmer_result_ord$estimate_status_in_linear >= 0 ] <- "POS_noncontroller_effect"
lmer_result_ord<- unique(lmer_result_ord)


lmer_result_ord_chr_location <- merge(lmer_result_ord, gene_locations, by.x="tar", by.y="hgnc_symbol")
#axis_set <- lmer_result_ord$rank

axis_set <- lmer_result_ord_chr_location %>% 
  group_by(chromosome_name) %>% 
  summarize(center = mean(start_position)) %>% 
  arrange(as.numeric(chromosome_name)) %>% 
  as.data.frame()

lmer_result_ord_chr_location$chromosome_name <- as.numeric(lmer_result_ord_chr_location$chromosome_name)
lmer_result_ord_chr_location  <- lmer_result_ord_chr_location %>% arrange(as.numeric(chromosome_name))

manhplot <- ggplot(lmer_result_ord_chr_location,  aes(x = start_position, y = -log10(pvalue), color = as.factor(chromosome_name))) +
  geom_point(aes(size=abs(estimate_status_in_linear), shape=dir, alpha = 0.35)) +
  geom_text(aes( size = 0.5,label=ifelse((-log10(pvalue)>=1.3 | -log10(anova_linear_p)>=1.3) ,as.character(tar),'')),hjust=0.5,vjust=0, angle = 15)+

  scale_y_continuous(breaks = seq(0, 12, by = 1))+
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
    scale_shape_manual(values=c("NEG_noncontroller_effect" = 17,"POS_noncontroller_effect"= 19))+

    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("significant anova linear effect by target")) +
    facet_grid(cols = vars(chromosome_name),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
    labs(x = NULL, 
       y = "-log10(P)" , size=10)+
        theme_minimal() +
            theme(legend.position="none", panel.grid = element_blank(),
             panel.spacing = unit(0.2, "cm")) + 
             scale_x_continuous(label = axis_set$chromosome_name,  breaks = NULL)
            


 
   print(manhplot)
   
dev.off()

 }
