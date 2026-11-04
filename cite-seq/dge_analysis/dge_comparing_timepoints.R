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
meta = data.frame(pid_tp=paste0(meta_tmp$participant_id, "_", match_UDs_to_otherTPs(meta_tmp$timepoint)), pid=meta_tmp$participant_id, tp=match_UDs_to_otherTPs(meta_tmp$timepoint), orig.ident=meta_tmp$orig.ident) %>% na.omit() %>% filter(!(orig.ident == "AMFAR-POOL-RR2-SCG1" & pid == "60610") & !(pid == "35933" & tp == "postR1") & !(pid == "60610" & tp == "postR1") & !(pid == "60610" & tp == "preR") ) 

pb_list = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/pb_exp_all_genes_fine_annot.rds")

# preATI vs vaccinal
tp1 = args[1] # This will be denominator in log2FC calculations
tp2 = args[2] # This will be numerator in log2FC calculations

outdir = paste0("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/dge_analysis/", tp2, "_vs_", tp1, "/")
if(!dir.exists(outdir)) {
  dir.create(outdir)
}

print(paste0("Comparing ", tp2, " and ", tp1, ": log2FC(", tp2, "/", tp1, ")"))

ctypes = names(pb_list)
param <- SnowParam(60, "SOCK", progressbar = TRUE)
for(ctype in ctypes) {
  print(ctype)
  cd = select(pb_list[[ctype]], contains(c(tp1, tp2))) %>% select(contains(meta$pid_tp))
  feats=names(which(rowMeans(cpm(cd) > 1) > 0.4))
  print(length(feats))
  cd = cd[feats,]
  meta_tmp = meta[ match(colnames(cd), meta$pid_tp), ] %>% remove_rownames() %>% column_to_rownames("pid_tp")
  meta_tmp$tp = factor(meta_tmp$tp, levels=c(tp1, tp2))

  if( length(unique(meta_tmp$tp)) < 2 ) {
    print( paste0("Skipping ", ctype, ", because it does not have both ", tp1, " and ", tp2, " samples."))
    print(colnames(cd))
    next
  }
  if( ! any(table(meta_tmp$tp) > 1) ) {
    print( paste0("Skipping ", ctype, ", because neither groups have more than one sample."))
    print(colnames(cd))
    next
  }
  
  y <- DGEList(counts = cd)
  y <- calcNormFactors(y)

  # model: time fixed, patient random
  form <- ~ tp + (1|pid) + (1|tp:orig.ident)

  result <- try({
    v <- voomWithDreamWeights(y, form, data = meta_tmp, BPPARAM=param)
    fit <- dream(v, form, data = meta_tmp, BPPARAM=param)
  }, silent = TRUE)

  if (inherits(result, "try-error")) {
    print(paste0("There was a following error for ", ctype))
    print(result)
    print("Using DESeq2 without batch correction")
    dds <- DESeqDataSetFromMatrix(cd, colData=meta_tmp, design = ~tp)
    dds <- DESeq(dds)
    contrast = c("tp",tp2,tp1)
    res <- results(dds, contrast = contrast,alpha = 0.05)
    res = rename(res, "log2FoldChange"="logFC", "baseMean"="AveExpr", "pvalue"="P.Value", "padj"="adj.P.Val")
    write.csv(res, paste0(outdir, "/v2_nestedBatchRandom_", ctype, "_deseq.csv"))
  } else {
    #fit <- dream(v, form, data = meta_tmp, BPPARAM=param)

    # extract results
    tt <- topTable(fit, coef=paste0("tp",tp2), "Inf")

    write.csv(tt, paste0(outdir,"/v2_nestedBatchRandom_", ctype, ".csv"))
  }
}

