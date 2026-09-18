source("/krummellab/data1/DSCoLab/AMFAR/analysis/src/amfar/src/10x/load_object_v2.R")

# Removing 2523_UD1 from RR2 pool, because this sample is in two pools (RR2 & RR10), and there are 2000+ more cells in RR10
amfar = subset( amfar, cells = Cells(amfar)[! (amfar$pid_tp == "2523_UD1" & amfar$orig.ident=="AMFAR-POOL-RR2-SCG1")])

# Keep only those sample x celltype combination where at least 10 cells were captured.
cell_count = as.data.frame(table(amfar$pid_tp, amfar$fine_annotation))
cell_count = cell_count[ cell_count$Freq >= 10, ] %>% mutate(key = paste0(Var1, "#", Var2))
amfar = subset(amfar, cells = Cells(amfar)[ paste0(amfar$pid_tp, "#", amfar$fine_annotation) %in% cell_count$key ] )

celltype_levels = sort(unique(amfar$fine_annotation))
ct_df = data.frame(ori=celltype_levels, new=gsub("_","-",celltype_levels))
amfar$key = paste0(match_UDs_to_otherTPs(amfar$pid_tp), "#", amfar$fine_annotation)
pb = AggregateExpression(amfar, assays='RNA', slot="counts", group.by="key")


# AggregateExpression stupidly converts all underscores to dashes. And adds 'g' in the beginning since pid is a number. Converting them back.
cn = colnames(pb$RNA)
colnames(pb$RNA) = paste0( 
    gsub("g","",gsub("-","_",str_extract(cn, "^[^#]+"))), 
    "#", 
    ct_df$ori[ match(str_extract(cn, "[^#]+$"), ct_df$new) ])


pb_list = list()
for(ctype in celltype_levels) {
        pb_list[[ctype]] = pb$RNA %>% as.data.frame() %>% dplyr::select(ends_with(paste0("#",ctype))) %>% setNames(gsub(paste0("#", gsub("\\+","\\\\\\+",ctype) ),"",colnames(.)))
}
print("Saving the pseudobulked expression")
saveRDS(pb_list, "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/pb_exp_all_genes_fine_annot.rds")



# Pseudobulk for broad populations
library(tidyverse)
library(Seurat)
amfar = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/5.merged_annotation_final_object_v1/sobj_final_v1.rds")

# Removing 2523_UD1 from RR2 pool, because this sample is in two pools (RR2 & RR10), and there are 2000+ more cells in RR10
amfar = subset( amfar, cells = Cells(amfar)[! (amfar$pid_tp == "2523_UD1" & amfar$orig.ident=="AMFAR-POOL-RR2-SCG1")])

# Keep only those sample x celltype combination where at least 10 cells were captured.
cell_count = as.data.frame(table(amfar$pid_tp, amfar$broad_annotation_new))
cell_count = cell_count[ cell_count$Freq >= 10, ] %>% mutate(key = paste0(Var1, "#", Var2))
amfar = subset(amfar, cells = Cells(amfar)[ paste0(amfar$pid_tp, "#", amfar$broad_annotation_new) %in% cell_count$key ] )

celltype_levels = sort(unique(amfar$broad_annotation_new))
ct_df = data.frame(ori=celltype_levels, new=gsub("_","-",celltype_levels))
amfar$key = paste0(match_UDs_to_otherTPs(amfar$pid_tp), "#", amfar$broad_annotation_new)
pb = AggregateExpression(amfar, assays='RNA', slot="counts", group.by="key")


# AggregateExpression stupidly converts all underscores to dashes. And adds 'g' in the beginning since pid is a number. Converting them back.
cn = colnames(pb$RNA)
colnames(pb$RNA) = paste0( 
    gsub("g","",gsub("-","_",str_extract(cn, "^[^#]+"))), 
    "#", 
    ct_df$ori[ match(str_extract(cn, "[^#]+$"), ct_df$new) ])


pb_list = list()
for(ctype in celltype_levels) {
        pb_list[[ctype]] = pb$RNA %>% as.data.frame() %>% dplyr::select(ends_with(paste0("#",ctype))) %>% setNames(gsub(paste0("#", gsub("\\+","\\\\\\+",ctype) ),"",colnames(.)))
}
print("Saving the pseudobulked expression")
saveRDS(pb_list, "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/pb_exp_all_genes_broad_annot.rds")

