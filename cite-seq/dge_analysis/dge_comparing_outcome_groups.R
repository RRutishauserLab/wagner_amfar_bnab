library(tidyverse)
library(edgeR)
library(variancePartition)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)

match_UDs_to_otherTPs <- function(tp) {
  tp = as.character(tp)
  tp = case_when(
    grepl("UD1", tp) ~ gsub("UD1", "preR", tp),
    grepl("UD2", tp) ~ gsub("UD2", "postR1", tp),
    TRUE ~ tp
    )
  return(tp)
}


meta_tmp = read.csv("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/sample_mapping/sample_mapping.csv")
p_meta = read.csv("/krummellab/data1/DSCoLab/AMFAR/analysis/sample_metadata/participant_metadata_v2.csv")
# Convert NR to C
p_meta[ p_meta$outcome == "NR", "outcome"] = "C"

meta = data.frame(pid_tp=paste0(meta_tmp$participant_id, "_", match_UDs_to_otherTPs(meta_tmp$timepoint)), pid=meta_tmp$participant_id, tp=match_UDs_to_otherTPs(meta_tmp$timepoint), orig.ident=meta_tmp$orig.ident, outcome = p_meta$outcome[ match(meta_tmp$participant_id, p_meta$pid) ]) %>% na.omit() %>% filter(!(orig.ident == "AMFAR-POOL-RR2-SCG1" & pid == "60610") & !(pid == "35933" & tp == "postR1") & !(pid == "60610" & tp == "postR1") & !(pid == "60610" & tp == "preR") ) 

pb_list = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/pb_exp_all_genes_fine_annot.rds")

# preATI vs vaccinal
tp = args[1] # 

outdir = paste0("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/dge_analysis/CvsNC_", tp, "/")
if(! dir.exists(outdir)) {
  dir.create(outdir)
}

print(paste0("Comparing C and NC at ", tp, ": log2FC(C/NC)"))

ctypes = names(pb_list)
# param <- SnowParam(60, "SOCK", progressbar = TRUE)
for(ctype in ctypes) {
  print(ctype)
  cd = select(pb_list[[ctype]], contains( intersect( paste0(p_meta$pid[ p_meta$outcome %in% c("C","NC")], "_", tp), meta$pid_tp ) ))
  meta_tmp = meta[ match(colnames(cd), meta$pid_tp), ] %>% remove_rownames() %>% column_to_rownames("pid_tp")
  if( length(unique(meta_tmp$outcome)) < 2) {
    print( paste0("Skipping ", ctype, ", because it does not have both C and NC samples."))
    print(colnames(cd))
    next
  }
  if( ! any(table(meta_tmp$outcome) > 1) ) {
    print( paste0("Skipping ", ctype, ", because neither groups have more than one sample."))
    print(colnames(cd))
    next
  }

  #feats=names(which(rowMeans(cpm(cd) > 1) > 0.4))
  # Keep features that are expressed with CPM > 1 in more than 50% of samples in one or both outcome groups.
  feats = names(which(apply(cpm(cd) > 1, 1, function(x) { any(tapply(x, meta_tmp$outcome, mean) > 0.5) } )))
  print(length(feats))
  cd = cd[feats,]

  dds <- DESeqDataSetFromMatrix(cd, colData=meta_tmp, design = ~ outcome)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("outcome","C","NC")) %>% as.data.frame() %>% arrange(padj)
  res = rename(res, "log2FoldChange"="logFC", "baseMean"="AveExpr", "pvalue"="P.Value", "padj"="adj.P.Val")
  write.csv(res, paste0(outdir, "/v3_DESeq2woBatch_", ctype, ".csv"))
}

