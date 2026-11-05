# The data used for v2 was where the upper bounds were removed for nCount_RNA and nFeature_RNA in qc_and_addModalities step.

library(tidyverse)
library(Seurat)
library(dsb)
args = commandArgs(trailingOnly=TRUE)
s_idx = as.numeric(args[1])
proc_dir = "/krummellab/data1/immunox/AMFAR/data/single_cell_GEX/processed/"
out_dir = "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/1.dsb_normalization_adtUpDown_rnaUp_cutoffs_v2/"
pools = paste0("AMFAR-POOL-RR", 1:15, "-SCG1")
s = pools[s_idx]
sobjs = list()
stat = list()

# Combine sample-to-fml cluster mapping, colab ids, and sample metadata 
sample_map = read.csv("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/sample_mapping/sample_mapping.csv")
sample_map$key = paste0(sample_map$orig.ident, "#", sample_map$BEST.GUESS)
# sample_map$pid_tp = paste0(sample_map$participant_id, "_", sample_map$timepoint)
# Add participant metadata
participant_meta = read.csv("/krummellab/data1/DSCoLab/AMFAR/analysis/sample_metadata/participant_metadata_v2.csv")
sample_meta = cbind(sample_map, participant_meta[ match(sample_map$participant_id, participant_meta$pid) ,])
sample_meta$pid = NULL
sample_meta$new_pid = NULL
# Add colab ids
colab_ids = read.csv("/krummellab/data1/DSCoLab/AMFAR/analysis/sample_metadata/sample_collection_shifted_dates_and_ids.csv") %>% dplyr::select(c("tp","pid_tp","pid","new_pid","days_from_bl","colabs_participant","colabs_tp","colabs_10x_sample","colabs_gex","colabs_pool"))
sample_meta = inner_join(
	sample_meta %>% mutate(pid_tp_pool = paste0(participant_id, "_", timepoint, orig.ident) ), 
	colab_ids %>% mutate(pid_tp_pool = paste0(pid_tp, colabs_pool) ), 
	by = "pid_tp_pool"
	) %>% mutate( participant_id = NULL, timepoint = NULL, pid_tp_pool = NULL )

#final_junk_markers = c("CD202b--Tie2--Tek","CD338--ABCG2","mast-cell-tryptase","CD9.1","CD140b","CD140a","CD326--EpCAM","CD138--Syndecan-1","XCR1.1","TCR-G--TCR-D","CD79a--IGA","C5L2","CD336--NKp44","CD1a","CD61","TCR-V-A-24-J-A-18--iNKT cell","CD135","CD102","CD49f","CD18","CD29","CD45","CD112--Nectin-2","CD52.1","CD48.1","CD47.1","CD90--Thy1","CD269--BCMA")

pdf(paste0(out_dir,"/plots_", s, ".pdf"), width=10, height=10)
print(paste0("Processing data for ", s))
stat_local = list()

sobj_local = readRDS( file.path(proc_dir, s, paste0("/qc_and_addModalities/",s,"_seurat_object_qcAndAddModalities.rds")) )

# Removing platelets and RBCs.
pl4genes_sum = colSums(sobj_local[['RNA']]@counts[c("PPBP","PF4","CAVIN2","GNG11"), ])
hbb3genes_sum = colSums(sobj_local[['RNA']]@counts[c("HBB","HBA1","HBA2"), ])
sobj_local = subset(sobj_local, cells=Cells(sobj_local)[ ! (hbb3genes_sum > 1 | pl4genes_sum > 1) ])


DefaultAssay(sobj_local) = 'ADT'


meta_local = sobj_local@meta.data %>% mutate(key = paste0(orig.ident, "#", BEST.GUESS))
sample_meta_idx = match(meta_local$key, sample_meta$key)
sobj_local = AddMetaData( sobj_local,
  sample_meta[sample_meta_idx,] %>% `rownames<-`(Cells(sobj_local))
  )

# Remove cells that did not map to any sample or mapped to a wrong sample (in case of RR14).
sobj_local = subset(sobj_local, cells = Cells(sobj_local)[ ! is.na(sobj_local$pid) ] )

# Perform dsb normalization for ADT data.
# Select 'background' droplets
raw_data <- Read10X_h5( file.path(proc_dir, s, "/cellranger/raw_feature_bc_matrix.h5") )
nCount_RNA = colSums(raw_data$`Gene Expression`)
nCount_ADT = colSums(raw_data$`Antibody Capture`)

rna_upper_cutoffs = c(315,130,315,100,215,215,215,245,235,315,185,100,215,150,170)
adt_lower_cutoffs = c(22,55,10,30,15,15,15,17,17,14,30,18,19,10,11)
adt_upper_cutoffs = c(200,550,300,400,400,200,200,400,150,200,300,250,250,100,300)

rna_upper = rna_upper_cutoffs[s_idx]
adt_lower = adt_lower_cutoffs[s_idx]
adt_upper = adt_upper_cutoffs[s_idx]
# Make plot of nCount_ADT vs nCount_RNA.
# 2d density in log space
p = data.frame(nCount_RNA = nCount_RNA, nCount_ADT = nCount_ADT) %>% ggplot(aes(log10(nCount_ADT+1), log10(nCount_RNA+1))) + geom_point(size=0.01) + geom_hex(aes(fill = stat(log(count))), bins=100) + scale_fill_distiller(palette="RdYlBu") + theme_classic() + geom_hline(yintercept=log10(rna_upper + 1)) + geom_vline(xintercept=c(log10(adt_lower + 1), log10(adt_upper + 1)))
ggExtra::ggMarginal( ggrastr::rasterize(p), type = "histogram", bins=50)

# 2d density in linear space
# p = data.frame(nCount_RNA = nCount_RNA, nCount_ADT = nCount_ADT) %>% ggplot(aes(log10(nCount_ADT+1), log10(nCount_RNA+1))) + geom_point(size=0.01) + geom_hex(bins=100) + scale_fill_distiller(palette="RdYlBu") + theme_classic() + geom_hline(yintercept=log10(101)) + geom_vline(xintercept=log10(31))
# ggExtra::ggMarginal( ggrastr::rasterize(p), type = "histogram", bins=50)


empty_drops <- colnames(raw_data$`Gene Expression`)[nCount_RNA < rna_upper & nCount_ADT > adt_lower & nCount_ADT < adt_upper ]
empty_drops <- empty_drops[! empty_drops %in% Cells(sobj_local)]
print(paste0("Number of empty droplets: ", length(empty_drops)))
ADT_background <- raw_data$`Antibody Capture`[, empty_drops]
rownames(ADT_background) = gsub("_", "-", rownames(ADT_background))
ADT_counts <- sobj_local[['ADT']]@counts

#Preserving the raw ADT data
sobj_local[['raw.ADT']] = sobj_local[['ADT']]

stat_local[[ 'empty_drops_count' ]] <- length(empty_drops)


adt_norm = DSBNormalizeProtein(
  cell_protein_matrix = ADT_counts,
  empty_drop_matrix = ADT_background,
  denoise.counts = TRUE,
  use.isotype.control = TRUE,
  isotype.control.name.vec = grep("Isotype", rownames(ADT_counts), value = T, ignore.case = T)
)

adt_norm = apply(adt_norm,
                 2,
                 function(x){ ifelse(test = x < -10, yes = -10, no = x)})

sobj_local[['ADT']] = CreateAssayObject(data = adt_norm)


p = VlnPlot(sobj_local, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0, group.by="colabs_10x_sample")
print(p)

print(dim(sobj_local))
stat_local[[ "dim_raw" ]] <- dim(sobj_local)
stat_local[[ "droplet.types" ]] = table(sobj_local$DROPLET.TYPE.FINAL)

stat[[ s ]] = stat_local
sobjs[[ s ]] = sobj_local
dev.off()
save(sobjs,stat, file= file.path(out_dir,paste0("objs_",s,".rdata")))

