library(Seurat)
library(tidyverse)
sobjs_all = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/3.rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3/sobjs_all_final_wRes0.6_sub9_sub12.rds")
annot=read.csv("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/3.rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3/annotation_res0.6.csv")
sobjs_all$broad_annotation = annot$broad_annotation[ match(sobjs_all$wsnn_harmony_res.0.6_sub, annot$cluster) ]


DefaultAssay(sobjs_all) = 'RNA'
all_cd3r = FetchData(sobjs_all, c("CD3D","CD3E","CD3G"))
all_cd4r = FetchData(sobjs_all, c("CD4"))
all_cd8r = FetchData(sobjs_all, c("CD8A","CD8B"))
all_br = FetchData(sobjs_all, c("CD19","MS4A1"))
all_mr = FetchData(sobjs_all, c("CD14","ITGAX"))

DefaultAssay(sobjs_all) = 'integrated.ADT'
all_cd3a = FetchData(sobjs_all, c("Hu.CD3-UCHT1"))
all_cd4a = FetchData(sobjs_all, c("Hu.CD4-RPA.T4"))
all_cd8a = FetchData(sobjs_all, c("Hu.CD8"))
all_ba = FetchData(sobjs_all, c("Hu.CD19","Hu.CD20-2H7"))
all_ma = FetchData(sobjs_all, c("Hu.CD14-M5E2","Hu.CD11c"))


pivot_longer(all_cd3r, names_to="feat", values_to="exp", cols=colnames(all_cd3r)) %>% ggplot(aes(x=exp, color=feat)) + geom_density() + theme_classic() + geom_vline(xintercept=0.5)


plot(density(all_cd4a[,1]))
abline(v=23)

plot(density(all_cd8a[,1]))
abline(v=37)


sub = subset(sobjs_all, subset = broad_annotation %in% c("CD4","CD8","otherT","Proliferating T"))
DefaultAssay(sub) = 'RNA'
cd3r = FetchData(sub, c("CD3D","CD3E","CD3G"))
cd4r = FetchData(sub, c("CD4"))
cd8r = FetchData(sub, c("CD8A","CD8B"))
br = FetchData(sub, c("CD19","MS4A1"))
mr = FetchData(sub, c("CD14","ITGAX"))

DefaultAssay(sub) = 'integrated.ADT'
cd3a = FetchData(sub, c("Hu.CD3-UCHT1"))
cd4a = FetchData(sub, c("Hu.CD4-RPA.T4"))
cd8a = FetchData(sub, c("Hu.CD8"))
ba = FetchData(sub, c("Hu.CD19","Hu.CD20-2H7"))
ma = FetchData(sub, c("Hu.CD14-M5E2","Hu.CD11c"))



pivot_longer(cd3r, names_to="feat", values_to="exp", cols=colnames(cd3r)) %>% ggplot(aes(x=exp, color=feat)) + geom_density() + theme_classic() + geom_vline(xintercept=0.5)

plot(density(cd4a[,1]))
abline(v=30)
abline(v=23)

plot(density(cd8a[,1]))
abline(v=32)
abline(v=37)


cd4 = (cd3a > 25 | rowMeans(cd3r > 0.5) > 0) & (rowSums(cd4a > 30) > 0 | rowSums(cd4r > 0.5) > 0) & (ba[,1] < 45 & ba[,2] < 38) & (ma[,1] < 32 & ma[,2] < 45)
cd8 = (cd3a > 25 | rowMeans(cd3r > 0.5) > 0) & (rowSums(cd8a > 32) > 0 | rowSums(cd8r > 0.5) > 0) & (ba[,1] < 45 & ba[,2] < 38) & (ma[,1] < 32 & ma[,2] < 45)

s_cd8 = subset(sub, cells=Cells(sub)[cd8 & !cd4])
s_cd4 = subset(sub, cells=Cells(sub)[cd4 & !cd8])

saveRDS(s_cd4, "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/manualSubset_CD4/sobjs_all_final.rds")
saveRDS(s_cd8, "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/manualSubset_CD8/sobjs_all_final.rds")

