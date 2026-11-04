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

new_days_from_LEF = c("35933"=7,"42894"=104,"10227"=7,"60610"=28,"47019"=162,"48758"=7,"71190"=163,"72210"=11,"30988"=140,"31092"=126)
p_meta$days_lastLEF_ATI_old = p_meta$days_lastLEF_ATI_old
p_meta$days_lastLEF_ATI = new_days_from_LEF[as.character(p_meta$new_pid)]

meta = data.frame(
  pid_tp=paste0(meta_tmp$participant_id, "_", match_UDs_to_otherTPs(meta_tmp$timepoint)), 
  pid=meta_tmp$participant_id, 
  tp=match_UDs_to_otherTPs(meta_tmp$timepoint), 
  orig.ident=meta_tmp$orig.ident, 
  outcome = p_meta$outcome[ match(meta_tmp$participant_id, p_meta$pid) ],
  daysLEF = p_meta$days_lastLEF_ATI[ match(meta_tmp$participant_id, p_meta$pid) ],
  dosesLEF = p_meta$LEF_doses[ match(meta_tmp$participant_id, p_meta$pid) ]
) %>% 
  na.omit() %>% 
  filter(!(orig.ident == "AMFAR-POOL-RR2-SCG1" & pid == "60610") & !(pid == "35933" & tp == "postR1") & !(pid == "60610" & tp == "postR1") & !(pid == "60610" & tp == "preR") )

pb_list = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/pb_exp_all_genes_fine_annot.rds")

# preATI vs vaccinal
tp = "preATI"

#outdir = paste0("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/dge_analysis/CvsNC_", tp, "/")
#if(! dir.exists(outdir)) {
#  dir.create(outdir)
#}

print(paste0("Comparing C and NC at ", tp, ": log2FC(C/NC)"))

ctypes = c("cMo", "cDC2")
res_list = list()
for(ctype in ctypes) {
  print(ctype)
  cd = select(pb_list[[ctype]], contains( intersect( paste0(p_meta$pid[ p_meta$outcome %in% c("C","NC")], "_", tp), meta$pid_tp ) ))
  meta_tmp = meta[ match(colnames(cd), meta$pid_tp), ] %>% remove_rownames() %>% column_to_rownames("pid_tp")
  meta_tmp$outcome = factor(meta_tmp$outcome)
  meta_tmp$daysLEF = scale(log(meta_tmp$daysLEF))
  meta_tmp$dosesLEF = scale(meta_tmp$dosesLEF)
  
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
  
  dds <- DESeqDataSetFromMatrix(cd, colData=meta_tmp, design = ~ outcome + daysLEF + dosesLEF )
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("outcome","C","NC")) %>% as.data.frame() %>% arrange(padj)
  res = rename(res, "log2FoldChange"="logFC", "baseMean"="AveExpr", "pvalue"="P.Value", "padj"="adj.P.Val")
  res_list[[ctype]] = res
  #write.csv(res, paste0(outdir, "/v3_DESeq2woBatch_", ctype, ".csv"))
}






GSEA = function(gene_list, GO_file, padj, useGAGE=TRUE, plot_title, nSigPath=20) {   # https://bioinformaticsbreakdown.com/how-to-gsea/
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes_full <- fgsea::fgsea(pathways = myGO,
                             stats = gene_list,
                             minSize=15,
                             maxSize=600) %>%
    #nperm=10000) %>%
    as.data.frame()
  fgRes <- fgRes_full %>%
    dplyr::filter(padj < !!padj)
  print(dim(fgRes))
  #fgRes_full = fgRes
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
  
  ups = as.data.frame(gaRes$greater) %>%
    tibble::rownames_to_column("Pathway") %>%
    dplyr::filter(!is.na(p.geomean) & q.val < padj ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>%
    tibble::rownames_to_column("Pathway") %>%
    dplyr::filter(!is.na(p.geomean) & q.val < padj ) %>%
    dplyr::select("Pathway")
  
  ## Define up / down pathways which are significant in both tests
  if(useGAGE) {
    keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
    keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  } else {
    keepups = fgRes[fgRes$NES > 0, ]
    keepdowns = fgRes[fgRes$NES < 0, ]
  }
  
  print(dim(rbind(keepups,keepdowns)))
  if( (nrow(keepups)+nrow(keepdowns)) == 0 ) {
    return(list("results_full"=fgRes_full))
  }
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway,
                              c( keepups$pathway, keepdowns$pathway))), ] %>%
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = head(fgRes, n = nSigPath)
  if(nrow(fgRes) > nSigPath & nrow(fgRes) < 2*nSigPath) {
    n = nrow(fgRes) - nSigPath
    filtRes = rbind(filtRes, tail(fgRes, n = n ) )
  } else if(nrow(fgRes) >= 2*nSigPath) {
    filtRes = rbind(filtRes, tail(fgRes, n = nSigPath ) )
  }
  
  
  upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
  downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
  colos = c(upcols, downcols)
  names(colos) = 1:length(colos)
  filtRes$Index = as.factor(1:nrow(filtRes))
  
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=plot_title) +
    theme_minimal()
  
  
  output = list("results" = fgRes, "results_full" = fgRes_full, "plot" = g)
  return(output)
}

gsea_dir = "/krummellab/data1/rpatel5/data/msigdb/20240320/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/"

gsea_res = list()
for(ctype in ctypes) {
  ct_lca1_vs_lca2 = res_list[[ctype]]
  ct_lca1_vs_lca2 = ct_lca1_vs_lca2 %>% arrange(-logFC)
  
  gene_rank = ct_lca1_vs_lca2$logFC %>% setNames(rownames(ct_lca1_vs_lca2))
  hallm_file = file.path(gsea_dir, "h.all.v2023.2.Hs.symbols.gmt")
  padj_cutoff = 0.1
  hallm_res = GSEA(gene_rank, hallm_file, padj_cutoff, useGAGE = FALSE, plot_title="Hallmark")
  gsea_res[[ctype]] = hallm_res
}





####### GSEA running ES plot
library(genekitr)

#' Reads GMT file and converts it to a data.frame containing two columns, first being the pathway names, and second the name of features that are part of the corresponding pathway.
#'
#' @param file: GMT file path.
gmtFile_to_Term2Gene <- function(file) {
  gmt = fgsea::gmtPathways(file)
  TERM2GENE = data.frame()
  for(p in names(gmt)) {
    TERM2GENE = rbind(TERM2GENE, data.frame(description=p, gene_name=gmt[[p]]) )
  }
  return(TERM2GENE)
}


my_geneset = gmtFile_to_Term2Gene(hallm_file)

fgres = gsea_res2$cDC2$results %>% filter(pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
fgres_new = fgres[c(2,2,8,6,7,3,4,4)-1] %>% setNames(c("ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue"))
fgres_new$rank = 1
fgres_new$leading_edge = 1
fgres_new$geneID = unlist(lapply(fgres$leadingEdge, function(x) { paste0(x, collapse = "/") } ))
fgres_new$Count = 1

res = res_list2$cDC2 %>% rownames_to_column("X")
my_genelist = data.frame(ID=res$X , logfc=res$logFC) %>% arrange(-logfc)

my_gse = list(gsea_df = fgres_new, geneset = my_geneset, genelist = my_genelist, exponent = 1, org = "hsapiens")

paths = list(
  "a" = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
)
plotlist_cd4 = list()
for(p in names(paths)) {
  plotlist_cd4[[p]] = plotGSEA(my_gse, plot_type = "classic", show_pathway = paths[[p]], wrap_length = 20 )
}
