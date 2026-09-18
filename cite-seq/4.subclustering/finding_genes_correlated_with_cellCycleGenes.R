# Finding genes correlated with with cell-cycle genes as described in this paper. https://pmc.ncbi.nlm.nih.gov/articles/PMC11215410/
# To remove the effects of proliferation on cell-type identification, genes correlated (P > 0.2) with either Mki67, Top2a, Hist1h2ap, or Hist1h1e were excluded from the set of highly variable genes before running all dimensionality reduction and clustering algorithms.
library(Seurat)
library(tidyverse)
s_cd4 = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/manualSubset_CD4/sobjs_all_final.rds")
sobj = s_cd4
expr_matrix <- GetAssayData(sobj, slot = "data", assay='RNA5')
ccgene = c("MKI67","TOP2A","H2AC11","H1-4")
ccgene_expr <- expr_matrix[ccgene, ]

cor_results <- apply(expr_matrix, 1, function(gene_expr) {
  apply(ccgene_expr, 1, function(ccg) {
    cor.test(gene_expr,ccg, method = "spearman")
  })
})

# Extract correlation coefficients and p-values
cor_dfs = list()
for(ccg in ccgene) {
  tmp_cor_df <- data.frame(
    gene = names(cor_results),
    ccgene = ccg,
    cor = sapply(cor_results, function(x) x[[ccg]]$estimate),
    pval = sapply(cor_results, function(x) x[[ccg]]$p.value)
  )
  # Adjust p-values for multiple testing (FDR)
  tmp_cor_df$padj <- p.adjust(tmp_cor_df$pval, method = "fdr")
  cor_dfs[[ccg]] = tmp_cor_df
}
saveRDS(cor_dfs, "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/manualSubset_CD4/ccgene_cor_genes.rds")



s_cd8 = readRDS("/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/manualSubset_CD8/sobjs_all_final.rds")
sobj = s_cd8
expr_matrix <- GetAssayData(sobj, slot = "data", assay='RNA5')
ccgene = c("MKI67","TOP2A","H2AC11","H1-4")
ccgene_expr <- expr_matrix[ccgene, ]

cor_results <- apply(expr_matrix, 1, function(gene_expr) {
  apply(ccgene_expr, 1, function(ccg) {
    cor.test(gene_expr,ccg, method = "spearman")
  })
})

# Extract correlation coefficients and p-values
cor_dfs = list()
for(ccg in ccgene) {
  tmp_cor_df <- data.frame(
    gene = names(cor_results),
    ccgene = ccg,
    cor = sapply(cor_results, function(x) x[[ccg]]$estimate),
    pval = sapply(cor_results, function(x) x[[ccg]]$p.value)
  )
  # Adjust p-values for multiple testing (FDR)
  tmp_cor_df$padj <- p.adjust(tmp_cor_df$pval, method = "fdr")
  cor_dfs[[ccg]] = tmp_cor_df
}
saveRDS(cor_dfs, "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/4.subclustering_adtUpDown_rnaUp_cutoffs_v3/manualSubset_CD8/ccgene_cor_genes.rds")




