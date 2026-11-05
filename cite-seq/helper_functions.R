modify_annot <- function(fine_annotation_ori) {
  fine_annotation_ori = as.character(fine_annotation_ori)
  fine_annot = case_when(
    fine_annotation_ori == "CD4_memory-like_CTL" ~ "CD4_CTL_1",
    fine_annotation_ori == "CD4_CTL" ~ "CD4_CTL_2",
    fine_annotation_ori == "CD4_IFN_responsive" ~ "CD4_IFN",
    fine_annotation_ori == "CD8_IFN_responsive" ~ "CD8_IFN",
    fine_annotation_ori == "NK_IFN_stim" ~ "NK_IFN",
    fine_annotation_ori == "CD4_memory_CCR6_Th17-like" ~ "CD4_Th17-like",
    fine_annotation_ori == "CD8_Tcm_vs_Ttm" ~ "CD8_Tcm",
    fine_annotation_ori == "CD8_effector_memory" ~ "CD8_Tem_1",
    fine_annotation_ori == "CD8_effector_memory_2" ~ "CD8_Tem_2",
    fine_annotation_ori == "CD8_A6B1_Memory" ~ "CD8_mem_A6B1",
    fine_annotation_ori == "CD8_activated_1" ~ "CD8_Tem_activated",
    fine_annotation_ori == "CD8_activated_2" ~ "CD8_Tcm_activated",
    fine_annotation_ori == "NKT" ~ "gdT",
    fine_annotation_ori == "CD8_NKT" ~ "gdT_CD8",
    fine_annotation_ori == "NK1" ~ "NK_mature_cytotoxic",
    fine_annotation_ori == "NK2_1" ~ "NK_immature_1",
    fine_annotation_ori == "NK2_2" ~ "NK_immature_2",
    fine_annotation_ori == "NK2" ~ "NK_adaptive_3",
    #    fine_annotation_ori == "NKT" ~ "gdT_1",
    fine_annotation_ori == "Treg" ~ "CD4_Treg",
    #    fine_annotation_ori == "CD8_NKT" ~ "gdT_2",
    TRUE ~ fine_annotation_ori
  )
  return(fine_annot)
}

rename_tp <- function(tp) {
  tp = gsub("vaccinal","bNAb",tp)
  tp = gsub("postR1","postR",tp)
  return(tp)
}


tp_order = c("preATI", "bNAb", "preR", "postR", "postR2")


all_ct=c(
  "cMo","ncMo_cd16mid","ncMo_cd16high","cDC1","cDC2","pDC",
  "B_naive","B_IgM+IgD+","B_mem_CD11c","B_mem_not_switched","B_mem_switched","Plasmablasts",
  "CD4_naive","CD4_naive_SOX4","CD4_memory_NOS","CD4_Treg","CD4_CTL_1","CD4_CTL_2","CD4_Th17-like","CD4_cycling","CD4_IFN",
  "CD8_naive","CD8_Tcm","CD8_Tem_1","CD8_Tem_2","CD8_Tcm_activated","CD8_Tem_activated","CD8_mem_A6B1","CD8_Trm-like","CD8_IFN",
  "MAIT",
  "gdT","gdT_CD8",
  "NKT_like",
  "NK_mature_cytotoxic", "NK_immature_1","NK_immature_2","NK_adaptive_1","NK_adaptive_2","NK_adaptive_3","NK_IFN",
  "HSPC"
)


get_pid_from_pidtp = function(pid_tp) {
  pid = str_extract(pid_tp, "^\\d+")
  return(pid)
}

