library(Seurat)
library(tidyverse)


# V3 clustering of the full object - removing Isotype control Abs from the analysis
# Res 0.5 did not match well with v2 res 0.5. So performed clustering using res 0.6 which matched better
sobjs_all = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/3.rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3/sobjs_all_final.rds")
sobjs_all = FindClusters(sobjs_all, graph.name = "wsnn_harmony", algorithm = 3, resolution = 0.6, verbose = TRUE, n.start = 10)

sub9 = subset(sobjs_all, subset = wsnn_harmony_res.0.6 == 9)
sub12 = subset(sobjs_all, subset = wsnn_harmony_res.0.6 == 12)

sub9 = FindClusters(sub9, graph.name = "wsnn_harmony", algorithm = 3, resolution = 0.2, verbose = TRUE, n.start = 10, cluster.name="wsnn_harmony_c9_subclust_r0.2")
sub12 = FindClusters(sub12, graph.name = "wsnn_harmony", algorithm = 3, resolution = 0.3, verbose = TRUE, n.start = 10, cluster.name="wsnn_harmony_c12_subclust_r0.3")


wsnn_harmony_res.0.6_sub = as.character(sobjs_all$wsnn_harmony_res.0.6) %>% setNames(Cells(sobjs_all))
wsnn_harmony_res.0.6_sub[ Cells(sub9) ] = paste0("9_", as.character(sub9$wsnn_harmony_c9_subclust_r0.2))
wsnn_harmony_res.0.6_sub[ Cells(sub12) ] = paste0("12_", as.character(sub12$wsnn_harmony_c12_subclust_r0.3))
sobjs_all$wsnn_harmony_res.0.6_sub = factor(wsnn_harmony_res.0.6_sub)

#9_0 is CD4+ and 9_1 is CD8+
#12_0 is CD8+ and 12_1 CD4+, 12_2 and 12_3 are NKs
saveRDS(sobjs_all, "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/3.rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3/sobjs_all_final_wRes0.6_sub9_sub12.rds")

