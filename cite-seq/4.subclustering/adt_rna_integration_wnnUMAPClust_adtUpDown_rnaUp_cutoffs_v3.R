library(tidyverse)
library(Seurat)
library(future)
library(harmony)
library(SeuratWrappers)
library(Azimuth)

print_message <- function(message) {
  cat("[", format(Sys.time()), "]", message, "\n")
}

args = commandArgs(trailingOnly=TRUE)
ctype = args[1]
print_message(paste0("Processing for ", ctype))

print_message("Reading in the Seurat object.")
# Read in the seurat object and add broad annotation
sobjs_all = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/3.rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3/sobjs_all_final_wRes0.6_sub9_sub12.rds")
sobjs_all[['raw.ADT']] = NULL
#sobjs_all$prelim_cluster = ifelse( sobjs_all$wsnn_harmony_res.0.5 != 8, as.character(sobjs_all$wsnn_harmony_res.0.5), ifelse( sobjs_all$wsnn_harmony_res.4 %in% c(21,31,43), paste0("8_", as.character(sobjs_all$wsnn_harmony_res.4)), "8_other" ) )
sobjs_all$prelim_cluster = sobjs_all$wsnn_harmony_res.0.6_sub

annot = read.csv("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/3.rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3/annotation_res0.6.csv")
sobjs_all$broad_annotation = annot$broad_annotation[ match(sobjs_all$prelim_cluster, annot$cluster) ]
ctype_levels = sort(unique(sobjs_all$broad_annotation))


if(! ctype %in% ctype_levels) {
    stop(paste0("Can't find ", ctype, ".\n"))
}

print_message(paste0("Processing for ", ctype, "..."))
out_dir = paste0("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/",ctype)
dir.create(out_dir)
print_message(paste("output directory: ", out_dir))

print_message(paste("Subsetting data for ", ctype))
sobjs_all = subset(sobjs_all, subset = broad_annotation == ctype)

DefaultAssay(sobjs_all) <- "ADT"
sobjs.list <- SplitObject(sobjs_all, split.by = "orig.ident")

features <- rownames(sobjs_all[["ADT"]])

# Remove isotype control Abs
features = grep("Isotype", features, value=T, invert=T)

# GC
rm(sobjs_all)
gc()


## For "rpca", Data needs to be scaled and PCA needs to be calculated for each individual object.
print_message("Scaling and PCAing each library")
# For Proliferating T and HSPCs npcs=50 in PCA calculation for ADT data caused an error since some pools/batches had less than 50 cells, e.g. 17 cells for one batch for Proliferating T. So, I am setting the npcs=min(ncells)-1 for these two cell types.
npcs = 50
ncells = unlist(lapply(sobjs.list, ncol))
if(any(ncells <= 50) ) {
  npcs = round((min(ncells) - 1)/2)
}
sobjs.list <- lapply(X = sobjs.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE, npcs=npcs)
})

print_message("Calculating anchors")
options(future.globals.maxSize = 650 * 1024^3) # Setting maxSize to 100Gb. https://github.com/satijalab/seurat/issues/1845
plan("multisession", workers = 10)
anchors <- FindIntegrationAnchors(object.list = sobjs.list, anchor.features = features, reduction = "rpca", dims=1:npcs )
print_message("Saving anchors")
saveRDS(anchors, file = file.path(out_dir,paste0("anchors_v1.rds")))


# GC
rm(sobjs.list)
gc()

# this command creates an 'integrated' data assay
print_message("Performing integration")
sobj_all_adt_corr <- IntegrateData(anchorset = anchors, new.assay.name = "integratedADT", dims=1:npcs)
print_message("Saving ADT-integrated data")
saveRDS(sobj_all_adt_corr, file = file.path(out_dir,paste0("sobj_adt_corr_v1.rds")))

# GC
rm(anchors)
gc()

# Calculate S and G2M scores.
# Preparing cell-cycle genes. cc.genes is loaded with Seurat.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(sobj_all_adt_corr) = 'RNA'
sobj_all_adt_corr = CellCycleScoring(sobj_all_adt_corr, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)


# Calculate PCA for the integratedADT
DefaultAssay(sobj_all_adt_corr) <- "integratedADT"
sobj_all_adt_corr <- ScaleData(sobj_all_adt_corr, vars.to.regress = c("nCount_ADT","nFeature_ADT","S.Score", "G2M.Score") ) %>%
                RunPCA(npcs = 30, verbose = FALSE, reduction.name = "int.apca")
saveRDS(sobj_all_adt_corr, file = file.path(out_dir,paste0("sobj_adt_corr_v1.1.rds")))



sobj_all_adt_corr = readRDS(file.path(out_dir,paste0("sobj_adt_corr_v1.1.rds")))

sobjs_all = sobj_all_adt_corr
rm(sobj_all_adt_corr)
gc()


## RNA integration
remove_XY_genes = FALSE
remove_TCR_BCR_genes = TRUE
xygenes_file = "/krummellab/data1/rpatel5/data/10X_ref_v5.0.0/chrXY_genes.tsv"
TCR_gene_regex = "^TR[ABDG][VJC]"
BCR_gene_regex = "^IG[HKL][VJC]|IGH[ADEGM]"



DefaultAssay(sobjs_all) = 'RNA'
sobjs_all <- RunAzimuth(sobjs_all, reference = "pbmcref")

print_message("RNA normalization, finding var. feat., scaling, PCA and scVI-based batch-correction")
sobjs_all[["RNA5"]] <- as(object = sobjs_all[["RNA"]], Class = "Assay5")
DefaultAssay(sobjs_all) = 'RNA5'
sobjs_all[["RNA5"]] <- split(sobjs_all[["RNA5"]], f = sobjs_all$orig.ident)

sobjs_all <- NormalizeData(sobjs_all) %>%
                FindVariableFeatures(selection.method = "vst", nfeatures = 2000)


var_feat = VariableFeatures(sobjs_all)

xygenes = read.csv(xygenes_file, header=F, sep="\t")
print("Variable features: chrX genes")
print(var_feat[ var_feat %in% xygenes$V2[ xygenes$V1 == "chrX"] ])
print("Variable features: chrY genes")
print(var_feat[ var_feat %in% xygenes$V2[ xygenes$V1 == "chrY"] ])
print("Variable features: mitochondrial genes")
print(grep("^MT-", var_feat, value=T))
print("Variable features: ribosomal genes")
print(grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", var_feat, value=T))
print("Variable features: starting with HB")
print(grep("^HB", var_feat, value=T))
print("Variable features: starting with IG")
print(grep("^IG", var_feat, value=T))
print("Variable features: TCR genes")
print(grep(TCR_gene_regex, var_feat, value=T))
print("Variable features: BCR genes")
print(grep(BCR_gene_regex, var_feat, value=T))


if(remove_XY_genes) {
  print_message("Removing chrX/chrY genes from the list of variable features")
  print(paste0("Number of variable features before: ", length(var_feat)))
  var_feat = var_feat[ ! var_feat %in% xygenes$V2  ]
  VariableFeatures(sobjs_all) = var_feat
  print(paste0("Number of variable features after: ", length(var_feat)))
}

if(remove_TCR_BCR_genes) {
  print_message("Removing T/BCR genes from the list of variable features")
  print(paste0("Number of variable features before: ", length(var_feat)))
  var_feat = var_feat[ ! ( grepl(TCR_gene_regex, var_feat) | grepl(BCR_gene_regex, var_feat) )  ]
  VariableFeatures(sobjs_all) = var_feat
  print(paste0("Number of variable features after: ", length(var_feat)))
}


sobjs_all <- sobjs_all %>%
                ScaleData(vars.to.regress = c("percent.mt","percent.ribo","nCount_RNA","nFeature_RNA","S.Score", "G2M.Score") ) %>%
                RunPCA(features = VariableFeatures(object = .), npcs=50)

Sys.setenv("RETICULATE_PYTHON"="/krummellab/data1/rpatel5/software/miniconda3/envs/scvi/bin/python")
sobjs_all <- IntegrateLayers(
 object = sobjs_all, method = scVIIntegration,
 new.reduction = "integrated.scvi",
 assay = "RNA5",
 ndims = 50,
 conda_env = "/krummellab/data1/rpatel5/software/miniconda3/envs/scvi/", verbose = FALSE, group="orig.ident"
)

print_message("Saving object after batch-correction")
saveRDS(sobjs_all, file = file.path(out_dir,"sobjs_all_vSCVI.rds"))

sobjs_all <- IntegrateLayers(
  object = sobjs_all, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  assay = "RNA5",
  verbose = FALSE
)
print_message("Saving object after harmony batch-correction")
saveRDS(sobjs_all, file = file.path(out_dir,"sobjs_all_vHarmony.rds"))


########### WNN
print("Multi-modal data integration: scvi + rpca")
sobjs_all <- FindMultiModalNeighbors(
  sobjs_all, reduction.list = list("integrated.scvi", "int.apca"),
  knn.graph.name = "wknn_scvi",
  snn.graph.name = "wsnn_scvi",
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight_scvi", weighted.nn.name = "weighted.nn_scvi", prune.SNN = 1/20
)
sobjs_all$integratedADT.weight_scvi = sobjs_all$integratedADT.weight
sobjs_all$RNA.weight_scvi = sobjs_all$RNA.weight

print("Multi-modal data integration: harmony + rpca")
sobjs_all <- FindMultiModalNeighbors(
  sobjs_all, reduction.list = list("harmony", "int.apca"),
  knn.graph.name = "wknn_harmony",
  snn.graph.name = "wsnn_harmony",
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight_harmony", weighted.nn.name = "weighted.nn_harmony", prune.SNN = 1/20
)
sobjs_all$integratedADT.weight_harmony = sobjs_all$integratedADT.weight
sobjs_all$RNA.weight_harmony = sobjs_all$RNA5.weight

print_message("Save object")
saveRDS(sobjs_all, file = file.path(out_dir,"sobjs_all_wnn.rds"))




####### Clustering
print_message("Clustering")
ress = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3,1.5,2)
for(res in ress) {
  print_message(paste0("resolution: ", res))
  sobjs_all = FindClusters(sobjs_all, graph.name = "wsnn_scvi", algorithm = 3, resolution = res, verbose = TRUE, n.start = 10, cluster.name = paste0("wsnn_scvi_", ctype ,"_res.", res) )
  sobjs_all = FindClusters(sobjs_all, graph.name = "wsnn_harmony", algorithm = 3, resolution = res, verbose = TRUE, n.start = 10, cluster.name = paste0("wsnn_harmony_", ctype ,"_res.", res) )
  gc()
}
print_message("Saving object after clustering")
saveRDS(sobjs_all, file = file.path(out_dir,"sobjs_all_clust.rds"))

###### UMAP
print_message("UMAP calculations")
#sobjs_all <- RunUMAP(sobjs_all, reduction = "harmony", reduction.name = "harmony.umap", reduction.key = "harmony.UMAP_", dims = 1:30)
sobjs_all <- RunUMAP(sobjs_all, reduction = "int.apca", reduction.name = "int.aumap", reduction.key = "int.aUMAP_", dims = 1:18)
sobjs_all <- RunUMAP(sobjs_all, nn.name = "weighted.nn_scvi", reduction.name = "wnn.umap_scvi", reduction.key = "wnnUMAPSCVI_")
sobjs_all <- RunUMAP(sobjs_all, nn.name = "weighted.nn_harmony", reduction.name = "wnn.umap_harmony", reduction.key = "wnnUMAPHarmony_")
saveRDS(sobjs_all, file = file.path(out_dir,"sobjs_all_clustUmap.rds"))


DefaultAssay(sobjs_all) = 'RNA5'
sobjs_all <- JoinLayers(sobjs_all)


# Add clrADT
DefaultAssay(sobjs_all) <- "clrADT"
sobjs_all <- NormalizeData(sobjs_all, normalization.method = 'CLR', margin = 2)
saveRDS(sobjs_all, file = file.path(out_dir,"sobjs_all_final.rds"))


stop()

# Find Markers
for(res in c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) {
  print(res)
  Idents(sobjs_all) <- paste0("wsnn_harmony_",ctype,"_res.", res)
  if(file.exists(file.path(out_dir, paste0("marker_adt_res",res,".csv")))) {
    marker_adt = read.csv(file.path(out_dir, paste0("marker_adt_res",res,".csv")), row.names=1)
  } else {
    marker_adt = FindAllMarkers(sobjs_all, assay = "clrADT", test.use = "MAST", only.pos = T,  max.cells.per.ident = 1000, random.seed = 1234, latent.vars = "orig.ident")
    write.csv(marker_adt, file=file.path(out_dir, paste0("marker_adt_res",res,".csv")) )
  }

  if(file.exists(file.path(out_dir, paste0("marker_rna_res",res,".csv")))) {
    marker_rna = read.csv(file.path(out_dir, paste0("marker_rna_res",res,".csv")), row.names=1)
  } else {
    marker_rna = FindAllMarkers(sobjs_all, assay = "RNA", test.use = "MAST", only.pos = T,  max.cells.per.ident = 1000, random.seed = 1234, latent.vars = "orig.ident")
    write.csv(marker_rna, file=file.path(out_dir, paste0("marker_rna_res",res,".csv")) )
  }



  markers = c("MS4A1", "CD79A", "CD27", "IGHD","CD19", "IGHM","GNLY", "NKG7", "KLRF1",  "KLRB1", "CD3E","GZMA", "GZMB", "GZMH", "FGFBP2", "CD8A", "CD8B", "CCR7", "PTPRC","CD4", "IL7R", "S100A4","LEF1", "SELL","TNFRSF4", "TIMP1", "USP10", "LGALS3","FOXP3", "IL2RA", "CCR10","SLC4A10", "TRAV1-2","LYZ", "CD14", "FCGR3A","FCER1A",  "CST3","CLEC4C", "MZB1","CX3CR1", "MS4A7",  "IGHG1", "TSPAN13","PPBP", "PF4","HBB") # From Kim's list.

  top_5_genes = marker_rna %>% arrange(p_val_adj, -avg_log2FC) %>% group_by(cluster) %>% dplyr::filter(row_number() %in% 1:10) %>% arrange(cluster) %>% pull(gene) %>% unique()

  pdf(file.path(out_dir, paste0("marker_rna_dotplots_res",res,".pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = unique(c(top_5_genes,markers)) , cols = c("blue","red"), assay = "RNA") + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()

  pdf(file.path(out_dir, paste0("marker_rna_dotplots_res",res,"_clust.pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = unique(c(top_5_genes,markers)) , cols = c("blue","red"), assay = "RNA", cluster.idents = T) + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()


  top_5_genes = marker_adt %>% arrange(p_val_adj, -avg_log2FC) %>% group_by(cluster) %>% dplyr::filter(row_number() %in% 1:10) %>% arrange(cluster) %>% pull(gene) %>% unique()

  pdf(file.path(out_dir, paste0("marker_adt_dotplots_res",res,".pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = top_5_genes, cols = c("blue","red"), assay = "clrADT") + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()

  pdf(file.path(out_dir, paste0("marker_adt_dotplots_res",res,"_clust.pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = top_5_genes, cols = c("blue","red"), assay = "clrADT", cluster.idents = T) + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()




  # SCVI
  print(res)
  Idents(sobjs_all) <- paste0("wsnn_scvi_",ctype,"_res.", res)
  if(file.exists(file.path(out_dir, paste0("marker_adt_scvi_res",res,".csv")))) {
    marker_adt = read.csv(file.path(out_dir, paste0("marker_adt_scvi_res",res,".csv")), row.names=1)
  } else {
    marker_adt = FindAllMarkers(sobjs_all, assay = "clrADT", test.use = "MAST", only.pos = T,  max.cells.per.ident = 1000, random.seed = 1234, latent.vars = "orig.ident")
    write.csv(marker_adt, file=file.path(out_dir, paste0("marker_adt_scvi_res",res,".csv")) )
  }

  if(file.exists(file.path(out_dir, paste0("marker_rna_scvi_res",res,".csv")))) {
    marker_rna = read.csv(file.path(out_dir, paste0("marker_rna_scvi_res",res,".csv")), row.names=1)
  } else {
    marker_rna = FindAllMarkers(sobjs_all, assay = "RNA", test.use = "MAST", only.pos = T,  max.cells.per.ident = 1000, random.seed = 1234, latent.vars = "orig.ident")
    write.csv(marker_rna, file=file.path(out_dir, paste0("marker_rna_scvi_res",res,".csv")) )
  }



  markers = c("MS4A1", "CD79A", "CD27", "IGHD","CD19", "IGHM","GNLY", "NKG7", "KLRF1",  "KLRB1", "CD3E","GZMA", "GZMB", "GZMH", "FGFBP2", "CD8A", "CD8B", "CCR7", "PTPRC","CD4", "IL7R", "S100A4","LEF1", "SELL","TNFRSF4", "TIMP1", "USP10", "LGALS3","FOXP3", "IL2RA", "CCR10","SLC4A10", "TRAV1-2","LYZ", "CD14", "FCGR3A","FCER1A",  "CST3","CLEC4C", "MZB1","CX3CR1", "MS4A7",  "IGHG1", "TSPAN13","PPBP", "PF4","HBB") # From Kim's list.

  top_5_genes = marker_rna %>% arrange(p_val_adj, -avg_log2FC) %>% group_by(cluster) %>% dplyr::filter(row_number() %in% 1:10) %>% arrange(cluster) %>% pull(gene) %>% unique()

  pdf(file.path(out_dir, paste0("marker_rna_dotplots_scvi_res",res,".pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = unique(c(top_5_genes,markers)) , cols = c("blue","red"), assay = "RNA") + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()

  pdf(file.path(out_dir, paste0("marker_rna_dotplots_scvi_res",res,"_clust.pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = unique(c(top_5_genes,markers)) , cols = c("blue","red"), assay = "RNA", cluster.idents = T) + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()


  top_5_genes = marker_adt %>% arrange(p_val_adj, -avg_log2FC) %>% group_by(cluster) %>% dplyr::filter(row_number() %in% 1:10) %>% arrange(cluster) %>% pull(gene) %>% unique()

  pdf(file.path(out_dir, paste0("marker_adt_dotplots_scvi_res",res,".pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = top_5_genes, cols = c("blue","red"), assay = "clrADT") + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()

  pdf(file.path(out_dir, paste0("marker_adt_dotplots_scvi_res",res,"_clust.pdf")), width=25, height=8)
  print(DotPlot(sobjs_all, features = top_5_genes, cols = c("blue","red"), assay = "clrADT", cluster.idents = T) + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + scale_color_distiller(palette = "RdBu"))
  dev.off()
}

