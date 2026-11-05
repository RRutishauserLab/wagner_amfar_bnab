library(tidyverse)
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


args <- commandArgs(trailingOnly = TRUE)

#gsea_dir = "/krummellab/data1/rpatel5/projects/202112_DexCOMET/gsea_pathways"
gsea_dir = "/krummellab/data1/rpatel5/data/msigdb/20240320/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/"
data_dir = args[1] # "/krummellab/data1/DSCoLab/AMFAR/analysis/10x/dge_analysis/postR1_vs_preR/"
#dge_files = list.files(data_dir, pattern="csv$", full.name=T)
dge_files = list.files(data_dir, pattern=".csv$", full.name=T, recursive=F) %>% grep("gsea", ., invert=T, value=T)

#dge_files = list.files(data_dir, pattern="_dgea.csv$", full.name=T) %>% grep("v4.*CD[48]|v4.*cDC",., value=T)
#dge_files = list.files(data_dir, pattern="_dgea.csv$", full.name=T) %>% grep("v5",., value=T)
for(f in dge_files) {
  print(f)
  if(grepl("gsea.csv",f)) {
    next
  }
  out_f = gsub(".csv","_gsea.pdf",f)
  out_f_csv = gsub(".csv","_gsea.csv",f)
  if(file.exists(out_f)) {
    #print(paste0("Skipping ", f))
    #next
  }
  print(out_f)
  print(out_f_csv)
  ct_lca1_vs_lca2 = read.csv(f)
  ct_lca1_vs_lca2 = ct_lca1_vs_lca2 %>% column_to_rownames("X") %>% arrange(-logFC)

  gene_rank = ct_lca1_vs_lca2$logFC %>% setNames(rownames(ct_lca1_vs_lca2))
  customDB_file = file.path("/krummellab/data1/DSCoLab/AMFAR/analysis/src/amfar/src/10x/dge_analysis/20210705_GeneLists_Tetpaper.gmt")
  customDB2_file = file.path("/krummellab/data1/DSCoLab/AMFAR/analysis/src/amfar/src/10x/dge_analysis/rutishauser_lab_cd8_gene_lists_for_cohn_lab.gmt")
  biocarta_file = file.path(gsea_dir, "c2.cp.biocarta.v2023.2.Hs.symbols.gmt")
  kegg_file = file.path(gsea_dir, "c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")
  reactome_file = file.path(gsea_dir, "c2.cp.reactome.v2023.2.Hs.symbols.gmt")
  tft_file = file.path(gsea_dir, "c3.tft.v2023.2.Hs.symbols.gmt")
  c7_file = file.path(gsea_dir, "c7.all.v2023.2.Hs.symbols.gmt")
  c5_file = file.path(gsea_dir, "c5.go.bp.v2023.2.Hs.symbols.gmt")
  hallm_file = file.path(gsea_dir, "h.all.v2023.2.Hs.symbols.gmt")
  padj_cutoff = 0.1
  customDB_res = GSEA(gene_rank, customDB_file, padj_cutoff, useGAGE = FALSE, plot_title="customDB")
  customDB2_res = GSEA(gene_rank, customDB2_file, padj_cutoff, useGAGE = FALSE, plot_title="customDB2")
  bioc_res = GSEA(gene_rank, biocarta_file, padj_cutoff, useGAGE = FALSE, plot_title="Biocarta")
  kegg_res = GSEA(gene_rank, kegg_file, padj_cutoff, useGAGE = FALSE, plot_title="KEGG")
  reac_res = GSEA(gene_rank, reactome_file, padj_cutoff, useGAGE = FALSE, plot_title="Reactome")
  tft_res = GSEA(gene_rank, tft_file, padj_cutoff, useGAGE = FALSE, plot_title="TF-target")
  c7_res = GSEA(gene_rank, c7_file, padj_cutoff, useGAGE = FALSE, plot_title="Immune-C7")
  c5_res = GSEA(gene_rank, c5_file, padj_cutoff, useGAGE = FALSE, plot_title="C5-BP")
  hallm_res = GSEA(gene_rank, hallm_file, padj_cutoff, useGAGE = FALSE, plot_title="Hallmark")
  gsea_res = list("customDB1" = customDB_res, "customDB2" = customDB2_res, "biocarta" = bioc_res, "kegg" = kegg_res, "reactome" = reac_res, "tft" = tft_res, "c7" = c7_res, "c5bp" = c5_res, "hallmark" = hallm_res)
  all_res = data.frame()
  pdf(out_f, width=15, height=10)
  for(d in names(gsea_res)) {
    if(! is.null(gsea_res[[d]]$plot) ) {
      all_res = rbind(all_res, gsea_res[[d]]$results_full %>% mutate(db = d))
      print(gsea_res[[d]]$plot)
    }
  }
  dev.off()
  all_res$leadingEdge = unlist(lapply(all_res$leadingEdge, paste0, collapse=", "))
  write.csv(all_res, out_f_csv)
}


