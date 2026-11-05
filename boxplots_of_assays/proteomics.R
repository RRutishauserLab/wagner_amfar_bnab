#' Make publication-quality plots of the proteomics data.
#' 
#' This script requires the "NULISA_Data.xlsx" file, containing NPQ values
#' for each protein. This script will output separate plots for each protein
#' in the output folder.
#' 
#' @author George Gruenhagen and Julia Wagner
#' @dateCreated 02/06/2026

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper Functions ============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Setup ======================================================================

#' Check for required libraries
#' 
#' @return TRUE if all required libraries are installed
checkForRequiredLibraries = function(required_libs) {
  is_installed = find.package(required_libs)
  if (length(is_installed) != length(required_libs)) {
      num_missing = length(required_libs) - length(is_installed)
      stop(paste0("Please install the ", num_missing, " missing libraries listed and try again."))
  }
  return(TRUE)
}

## Stats ======================================================================

#' Run a generalized linear mixed effects model on the data for a given population
#' 
#' @param pop Population to model
#' @param raw_df Raw data frame
#' 
#' @return Data frame with model results
modelData = function(pop, raw_df, full_formula_str, update_formula_str) {
  # Select data for population
  raw_df = raw_df[which(raw_df$Target == pop),]
  
  #
  full_formula = as.formula(full_formula_str)
  update_formula = as.formula(update_formula_str)

  # Fit models
  full_model    = suppressMessages(suppressWarnings( lme4::glmer(full_formula, data=raw_df) ))
  reduced_model = suppressMessages(suppressWarnings( update(full_model, update_formula) ))
  kr_compare = pbkrtest::KRmodcomp(full_model, reduced_model)
  
  # Extract results
  p_value = kr_compare$stats$p.value
  f_stat = kr_compare$stats$Fstat
  ddf = kr_compare$stats$ddf
  full_singular = ifelse(paste0(full_model@optinfo$conv$lme4$messages, collapse="; ") == "boundary (singular) fit: see help('isSingular')", T, F)
  reduced_singular = ifelse(paste0(reduced_model@optinfo$conv$lme4$messages, collapse="; ") == "boundary (singular) fit: see help('isSingular')", T, F)
  full_messages = paste0(full_model@optinfo$message, collapse="; ")
  full_warnings = paste0(full_model@optinfo$warnings, collapse="; ")
  full_messages = ifelse(full_messages == "NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached.", "", full_messages)

  res_df = data.frame(pop=pop, p=p_value, f=f_stat, ddf=ddf, full_messages=full_messages, full_warnings=full_warnings, full_singular=full_singular, reduced_singular=reduced_singular, stringsAsFactors=F)
  
  return(res_df)
}

modelDataAtTimepoints = function(timepoints, long_df, n_cores, full_formula_str="NPQ ~ timepoint + (1|PID_New)", update_formula_str=" .~. - timepoint") {
  # Make a nice message that lists the timepoints being modeled
  is_all_timpoints = all(levels(long_df$timepoint) %in% timepoints)
  timepoints_str = paste0("between ", paste0(timepoints, collapse=" and "))
  timepoints_str = ifelse(is_all_timpoints, "at all timepoints", timepoints_str)
  message(paste0("Modeling data ", timepoints_str, "..."))

  # Subset data by specified timepoints
  df = long_df[which(long_df$timepoint %in% timepoints),]
  df$timepoint = droplevels(df$timepoint)

  # Model the data
  res_df_list = bettermc::mclapply(unique(df$Target), modelData, raw_df=df, full_formula_str=full_formula_str, update_formula_str=update_formula_str, mc.stdout='ignore', mc.messages='ignore', mc.allow.fatal=TRUE, mc.allow.error=TRUE, mc.fail.early=TRUE, mc.dump.frames="no", mc.progress=TRUE)
  res_df = as.data.frame(data.table::rbindlist(res_df_list))

  # Rename columns to indicate timepoints
  timepoints_str = paste0(tolower(timepoints), collapse="_")
  timepoints_str = ifelse(is_all_timpoints, "all", timepoints_str)
  colnames(res_df) = paste0(colnames(res_df), "_", timepoints_str)
  colnames(res_df)[1] = "pop"

  return(res_df)
}

#' Convert p-values to significance stars
#' 
#' @param p Vector of p-values
#' 
#' @return Vector of significance stars
getSigStars = function(p) {
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,  1),
                      symbols = c("****", "***", "**", "*",  "ns"))
  symnum.args$x <- p
  pvalue.signif <- do.call(stats::symnum, symnum.args) %>% as.character()
  return(pvalue.signif)
}

#' Adjust p-values for multiple testing using multiple methods (q, Benjamini-Hochber, Bonferroni)
#' 
#' @param res_df Data frame with p-values
#' 
#' @return Data frame with adjusted p-values
pvalueAdjustment = function(res_df) {
  qvalue_res = tryCatch(                  { qvalue::qvalue(res_df$p)        }, 
                        error=function(e) { qvalue::qvalue(res_df$p, pi0=1) }) 
  res_df$q_pi0      = qvalue_res$pi0
  res_df$q          = qvalue_res$qvalues
  res_df$q_format   = prettyNum(res_df$q, digits=3, eps=0.001)
  res_df$q_stars    = getSigStars(res_df$q)
  res_df$bh         = p.adjust(res_df$p, method="BH")            # Benjamini-Hochberg correction (FDR)
  res_df$bh_format  = prettyNum(res_df$bh,  digits=3, eps=0.001)
  res_df$bh_stars   = getSigStars(res_df$bh)
  res_df$bon        = p.adjust(res_df$p, method="bonferroni")    # Bonferroni correction
  res_df$bon_format = prettyNum(res_df$bon, digits=3, eps=0.001)
  res_df$bon_stars  = getSigStars(res_df$bon)
  return(res_df)
}

categorizeAndAdjust = function(res_df, narrow_adpt_pops=narrow_adpt_pops, narrow_innt_pops=narrow_innt_pops, base_col_name="all") {
  # Keep a copy of the input dataframe
  res_df_input = res_df

  # Get the p-value column
  p_value_col = paste0("p_", base_col_name)
  if (!(p_value_col %in% colnames(res_df))) { stop(paste0("Column ", p_value_col, " not found in res_df")) }
  res_df$p = res_df[, p_value_col]
  res_df$p_format = prettyNum(res_df$p, digits=3, eps=0.001)
  res_df$p_stars = getSigStars(res_df$p)

  # Adjust the p-values population category separately
  res_df = pvalueAdjustment(res_df)

  # Rename added columns to indicate base column name
  cols_added = setdiff(colnames(res_df), colnames(res_df_input))
  colnames(res_df)[which(colnames(res_df) %in% cols_added)] = paste0(cols_added, "_", base_col_name)

  return(res_df)
}

## Plotting ===================================================================

# #' Custom ggplot2 theme for publication-quality figures
# theme_classic_pub = function(title_font_size=6, title_font_face="plain", axis_font_size=6, axis_font_face="plain", axis_line_width=0.5, legend_shape_text_space=0.25) {
#   theme_classic(base_size=6, base_line_size=0.25, base_family='Helvetica') %+replace% theme(panel.grid.major=element_line(linetype="dashed", color='grey90', linewidth=axis_line_width*lwd_converter), panel.grid.minor=element_blank(), axis.title=element_text(size=axis_font_size, face=axis_font_face), axis.text=element_text(size=axis_font_size), plot.title=element_text(size=title_font_size, face=title_font_face, margin=margin(0,0,3,0)), legend.title=element_text(size=6), legend.text=element_text(size=6, margin=margin(l=legend_shape_text_space, unit="pt")), legend.background=element_rect(fill='transparent', color='transparent'), axis.line=element_line(linewidth=axis_line_width*lwd_converter, lineend='round'), plot.background=element_blank(), panel.background=element_blank(), strip.text=element_text(size=6))
# }

# #' Function to create a single plot for one protein
# singlePlot = function(protein, df, pid_col_df, base_name, dpi=900, file_ext='png') {
#   # Subset the data to the input protein
#   df_sub = df[which(df$Target == protein),]
  
#   # Set the beeswarmed coordinates of the points
#   df_sub$timepoint = droplevels(df_sub$timepoint)
#   df_sub$timepoint_num = as.numeric(df_sub$timepoint) - 1
#   df_sub$x = NA
#   for (time_pt in unique(df_sub$timepoint_num)) {
#     df_sub_tmp = data.frame(x=as.numeric(df_sub$timepoint_num[which(df_sub$timepoint_num == time_pt)]), y=df_sub[which(df_sub$timepoint_num == time_pt), "NPQ"])
#     df_sub$x[which(df_sub$timepoint_num == time_pt)] = ggbeeswarm:::offset_beeswarm(df_sub_tmp, yLim.expand=1, xRange=2, yRange=max(df_sub_tmp$y)*0.6, cex=5)$x + 1
#   }
#   df_paired = df_sub %>% group_by(PID_New) %>% filter(n()==2) %>% tidyr::pivot_wider(id_cols=PID_New, names_from=timepoint_num, values_from=c(x, NPQ)) %>% ungroup() %>% as.data.frame()

#   # Subset the data to the PIDs that should be filled vs not filled
#   df_sub_pid_fill    = df_sub[which(df_sub$PID_New %in% pid_col_df$PID[ pid_col_df$do_fill]),]
#   df_sub_pid_no_fill = df_sub[which(df_sub$PID_New %in% pid_col_df$PID[!pid_col_df$do_fill]),]

#   # Set plotting area size
#   width_per_time = 0.4
#   num_time = length(levels(df_sub$timepoint))
#   plotting_area_width = width_per_time*num_time
#   plotting_area_height = 1.25

#   # Plot title
#   plot_title = protein
#   plot_title = gsub("TNFRSF4", "OX40", plot_title)
#   plot_title = gsub("TNFRSF13C", "BAFF-R", plot_title)

#   # Signficance brackets
  
#   # get the max value per timepoint in wide format
#   # df_sub_sig = df_sub %>% group_by(timepoint) %>% summarise(NPQ = max(NPQ)) %>% ungroup() %>% tidyr::pivot_wider(df_sub_sig, names_from=timepoint, values_from=max_NPQ)
#   df_sub$timepoint2 = droplevels(df_sub$timepoint)
#   df_sub_sig = data.frame(max=max(df_sub$NPQ*1.015), Pre=levels(df_sub$timepoint2)[1], Post=levels(df_sub$timepoint2)[2])
#   df_sub_sig$stat_label = "*"
#   if (plot_title %in% c("LAMP3", "BAFF-R")) { df_sub_sig$stat_label = "**"  }
#   if (plot_title %in% c("ICOSLG"))          { df_sub_sig$stat_label = "***" }
#   # df_sub_sig = data.frame(timepoint=levels(df_sub$timepoint2), NPQ=tapply(df_sub$NPQ, df_sub$timepoint2, max))

#   # Make the plot
#   p = ggplot(df_sub, aes(x=timepoint, y=NPQ)) + 
#         geom_boxplot(color='gray30', fill="#ffffff99", alpha=0.75, outliers=F, width=0.5, linewidth=0.5*lwd_converter) +
#         geom_segment(data=df_paired, aes(x=x_0, xend=x_1, y=NPQ_0, yend=NPQ_1, color=as.character(PID_New)), linetype='solid', alpha=0.4, linewidth=0.5*lwd_converter) + 
#         geom_point(data=df_sub_pid_fill,    aes(x=x, fill=PID_New,  shape=PID_New), alpha=0.9, size=1.25, stroke=0.5*0.753, color="black") + 
#         geom_point(data=df_sub_pid_no_fill, aes(x=x, color=PID_New, shape=PID_New), alpha=0.9, size=1.25, stroke=0.5*0.753) + 
#         ggpubr::geom_bracket(data=df_sub_sig, aes(xmin=Pre, xmax=Post, y.position=max, label=stat_label), color='black', size=0.5*lwd_converter, label.size=6/.pt, vjust=0.5) + 
#         # scale_y_continuous(limits = c(10.5, 15.9)) + 
#         scale_fill_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
#         scale_color_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
#         scale_shape_manual(values=setNames(pid_col_df$shape, pid_col_df$PID)) +
#         ggh4x::force_panelsizes(cols=unit(plotting_area_width, "in"), rows=unit(plotting_area_height, "in")) +
#         labs(x=NULL, y="NPQ", title=plot_title) + 
#         theme_classic_pub() + theme(panel.grid.major.x=element_blank(), legend.position='none', panel.grid.minor.y=element_line(linetype="dashed", color='grey90', linewidth=0.5*lwd_converter))
#         # guides(y=guide_axis(minor.ticks=T))

#   # Save the plot
#   ggsave(paste0(base_name, protein, ".", file_ext), plot=p, width=plotting_area_width+0.4, height=plotting_area_height+0.4, dpi=my_dpi)
#   return(p)
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main ========================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Setup ======================================================================

# Set file paths
in_dir = "/c4/home/ggruenh/re/julia/amfar/proteomics/data/"
out_dir = "/c4/home/ggruenh/re/julia/amfar/proteomics/results/"
file_ext = "pdf"
dirs_to_check = c(out_dir, paste0(out_dir, "/plots_bh"))
for (dir_i in dirs_to_check) {
  if (!dir.exists(dir_i)) { dir.create(dir_i, recursive=TRUE) }
}

# Set number of cores for parallelization
n_cores = 5

# Check for required libraries
is_installed = checkForRequiredLibraries(required_libs=c("bettermc", "dplyr", "ggh4x", "ggplot2", "plyr", "tidyr"))

# Load libraries
suppressWarnings(suppressMessages( library('dplyr',    quietly=T, warn.conflicts=F, verbose=F) ))
suppressWarnings(suppressMessages( library('ggplot2',  quietly=T, warn.conflicts=F, verbose=F) ))

## Data Processing ============================================================

# Read data
nulisa = as.data.frame(readxl::read_excel(paste0(in_dir, "NULISA_Data.xlsx"), sheet="NPQ Counts"))

# Add new columns
nulisa$PID = as.character(nulisa$PID)
nulisa$PID_timepoint = paste0(nulisa$PID, "_", nulisa$timepoint)
nulisa$PID_New = plyr::revalue(nulisa$PID, c("35933"="35933", "42894"="42894", "10227"="10227", "60610"="60610", "47019"="47019", "71190"="71190", "72210"="72210", "30988"="30988", "31092"="31092"))
nulisa$PID_New = factor(nulisa$PID_New, levels=c("31092", "35933", "42894", "10227", "30988", "47019", "71190", "72210", "60610"))
nulisa$C = ! nulisa$PID_New %in% c("31092", "35933", "42894")
nulisa$C = plyr::revalue(as.character(nulisa$C), c("TRUE"="C", "FALSE"="NC"))
nulisa$C = factor(nulisa$C, levels=c("C", "NC"))

# Subset data
timepoints_to_use = c("preATI", "ATI1", "ATI2", "preR")
# proteins_to_use = c("LAMP3","IL18R1", "CCL21", "ICOSLG", "TNFRSF13C", "TNFRSF4")
nulisa = nulisa[which(nulisa$timepoint %in% timepoints_to_use),]                 # subset to timepoints of interest
# nulisa = nulisa[which(nulisa$Target %in% proteins_to_use),]                    # subset to proteins of interest
nulisa = nulisa[which(nulisa$PID != "48758"),]                                   # remove PID 48758
nulisa = nulisa[which(nulisa$PID_timepoint != "72210_ATI2"),]                    # remove PID 72210 at ATI2
nulisa = nulisa[which(nulisa$PID_timepoint != "60610_preR"),]                    # remove PID 60610 at preR
nulisa = nulisa[which( !(nulisa$timepoint == "ATI1" & nulisa$PID != "72210") ),] # remove the ATI1 timepoint for all PIDs except 72210

# Recode timepoints
nulisa$timepoint[which(nulisa$timepoint %in% c("ATI1", "ATI2"))] = "bNAb"
nulisa$timepoint = factor(nulisa$timepoint, levels=c("preATI", "bNAb", "preR"))
long_df = nulisa[,c("PID", "PID_New", "C", "PID_timepoint", "timepoint", "Target", "NPQ")]
timepoints_to_use = c("preATI", "bNAb", "preR")

## Stats time =================================================================

# Model data across all timepoints
res_df = modelDataAtTimepoints(timepoints=timepoints_to_use, long_df=long_df, n_cores=n_cores)
res_df = categorizeAndAdjust(res_df, base_col_name="all")

# Model data b/w BL and bNAb
res_df_preati_bnab = modelDataAtTimepoints(timepoints=c("preATI", "bNAb"), long_df=long_df, n_cores=n_cores)
res_df_preati_bnab = categorizeAndAdjust(res_df_preati_bnab, base_col_name="preati_bnab")
res_df[,colnames(res_df_preati_bnab)[-1]] = res_df_preati_bnab[match(res_df$pop, res_df_preati_bnab$pop), colnames(res_df_preati_bnab)[-1]]

# Model data b/w bNAb and PreR
res_df_bnab_prer = modelDataAtTimepoints(timepoints=c("bNAb", "preR"), long_df=long_df, n_cores=n_cores)
res_df_bnab_prer = categorizeAndAdjust(res_df_bnab_prer, base_col_name="bnab_prer")
res_df[,colnames(res_df_bnab_prer)[-1]] = res_df_bnab_prer[match(res_df$pop, res_df_bnab_prer$pop), colnames(res_df_bnab_prer)[-1]]

# Model data b/w PreATI and PreR
res_df_preati_prer = modelDataAtTimepoints(timepoints=c("preATI", "preR"), long_df=long_df, n_cores=n_cores)
res_df_preati_prer = categorizeAndAdjust(res_df_preati_prer, base_col_name="preati_prer")
res_df[,colnames(res_df_preati_prer)[-1]] = res_df_preati_prer[match(res_df$pop, res_df_preati_prer$pop), colnames(res_df_preati_prer)[-1]]

# Calculate means
res_df$mean_NPQ_preati = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preATI")], na.rm=T) })
res_df$mean_NPQ_bnab   = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "bNAb")],   na.rm=T) })
res_df$mean_NPQ_prer   = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preR")],   na.rm=T) })
res_df$median_NPQ_preati = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preATI")], na.rm=T) })
res_df$median_NPQ_bnab   = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "bNAb")],   na.rm=T) })
res_df$median_NPQ_prer   = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preR")],   na.rm=T) })
res_df = res_df[, c("pop", "mean_NPQ_preati", "mean_NPQ_bnab", "mean_NPQ_prer", "median_NPQ_preati", "median_NPQ_bnab", "median_NPQ_prer", setdiff(colnames(res_df), c("pop", "mean_NPQ_preati", "mean_NPQ_bnab", "mean_NPQ_prer", "median_NPQ_preati", "median_NPQ_bnab", "median_NPQ_prer")))]

## Save time csv ==============================================================

# # Remove bh from the results because really the adjustment should've been done on all proteins
# res_df = res_df[,!grepl("bh_", colnames(res_df))]
# res_df = res_df[,!grepl("q_", colnames(res_df))]
# res_df = res_df[,!grepl("bon_", colnames(res_df))]

# Make separate dataframes for results with any p or bh significance
p_cols_to_keep = c("pop", "mean_NPQ_preati", "mean_NPQ_bnab", "mean_NPQ_prer", "median_NPQ_preati", "median_NPQ_bnab", "median_NPQ_prer", "p_all", "p_preati_bnab", "p_bnab_prer", "p_preati_prer")
bh_cols_to_keep = c("pop", "mean_NPQ_preati", "mean_NPQ_bnab", "mean_NPQ_prer", "median_NPQ_preati", "median_NPQ_bnab", "median_NPQ_prer", "bh_all", "bh_preati_bnab", "bh_bnab_prer", "bh_preati_prer")
res_any_p_sig  = res_df[which(res_df$p_all  < 0.05 | res_df$p_preati_bnab  < 0.05 | res_df$p_bnab_prer  < 0.05 | res_df$p_preati_prer  < 0.05), p_cols_to_keep]
res_any_bh_sig = res_df[which(res_df$bh_all < 0.05 | res_df$bh_preati_bnab < 0.05 | res_df$bh_bnab_prer < 0.05 | res_df$bh_preati_prer < 0.05), bh_cols_to_keep]

# Write results
write.csv(res_df, file=paste0(out_dir, "stats_detailed_time.csv"), row.names=FALSE)
write.csv(res_any_p_sig, file=paste0(out_dir, "stats_any_p_sig_time.csv"), row.names=FALSE)
write.csv(res_any_bh_sig, file=paste0(out_dir, "stats_any_bh_sig_time.csv"), row.names=FALSE)

## Stats interaction ==========================================================

# Model data across all timepoints
res_df = modelDataAtTimepoints(timepoints=timepoints_to_use, long_df=long_df, n_cores=n_cores, full_formula_str="NPQ ~ timepoint*C + (1|PID_New)", update_formula_str=" .~. - timepoint:C")
res_df = categorizeAndAdjust(res_df, base_col_name="all")

# Model data b/w BL and bNAb
res_df_preati_bnab = modelDataAtTimepoints(timepoints=c("preATI", "bNAb"), long_df=long_df, n_cores=n_cores, full_formula_str="NPQ ~ timepoint*C + (1|PID_New)", update_formula_str=" .~. - timepoint:C")
res_df_preati_bnab = categorizeAndAdjust(res_df_preati_bnab, base_col_name="preati_bnab")
res_df[,colnames(res_df_preati_bnab)[-1]] = res_df_preati_bnab[match(res_df$pop, res_df_preati_bnab$pop), colnames(res_df_preati_bnab)[-1]]

# Model data b/w bNAb and PreR
res_df_bnab_prer = modelDataAtTimepoints(timepoints=c("bNAb", "preR"), long_df=long_df, n_cores=n_cores, full_formula_str="NPQ ~ timepoint*C + (1|PID_New)", update_formula_str=" .~. - timepoint:C")
res_df_bnab_prer = categorizeAndAdjust(res_df_bnab_prer, base_col_name="bnab_prer")
res_df[,colnames(res_df_bnab_prer)[-1]] = res_df_bnab_prer[match(res_df$pop, res_df_bnab_prer$pop), colnames(res_df_bnab_prer)[-1]]

# Model data b/w PreATI and PreR
res_df_preati_prer = modelDataAtTimepoints(timepoints=c("preATI", "preR"), long_df=long_df, n_cores=n_cores, full_formula_str="NPQ ~ timepoint*C + (1|PID_New)", update_formula_str=" .~. - timepoint:C")
res_df_preati_prer = categorizeAndAdjust(res_df_preati_prer, base_col_name="preati_prer")
res_df[,colnames(res_df_preati_prer)[-1]] = res_df_preati_prer[match(res_df$pop, res_df_preati_prer$pop), colnames(res_df_preati_prer)[-1]]

# Calculate means
res_df$mean_C_NPQ_preati   = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop   & long_df$timepoint == "preATI" & long_df$C == "C")], na.rm=T) })
res_df$mean_C_NPQ_bnab     = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop   & long_df$timepoint == "bNAb"   & long_df$C == "C")],   na.rm=T) })
res_df$mean_C_NPQ_prer     = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop   & long_df$timepoint == "preR"   & long_df$C == "C")],   na.rm=T) })
res_df$median_C_NPQ_preati = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preATI" & long_df$C == "C")], na.rm=T) })
res_df$median_C_NPQ_bnab   = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "bNAb"   & long_df$C == "C")],   na.rm=T) })
res_df$median_C_NPQ_prer   = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preR"   & long_df$C == "C")],   na.rm=T) })
res_df$mean_NC_NPQ_preati   = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop   & long_df$timepoint == "preATI" & long_df$C == "NC")], na.rm=T) })
res_df$mean_NC_NPQ_bnab     = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop   & long_df$timepoint == "bNAb"   & long_df$C == "NC")],   na.rm=T) })
res_df$mean_NC_NPQ_prer     = sapply(res_df$pop, function(pop) { mean(long_df$NPQ[which(long_df$Target == pop   & long_df$timepoint == "preR"   & long_df$C == "NC")],   na.rm=T) })
res_df$median_NC_NPQ_preati = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preATI" & long_df$C == "NC")], na.rm=T) })
res_df$median_NC_NPQ_bnab   = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "bNAb"   & long_df$C == "NC")],   na.rm=T) })
res_df$median_NC_NPQ_prer   = sapply(res_df$pop, function(pop) { median(long_df$NPQ[which(long_df$Target == pop & long_df$timepoint == "preR"   & long_df$C == "NC")],   na.rm=T) })
mean_and_median_cols = c("pop", "mean_C_NPQ_preati", "mean_NC_NPQ_preati", "mean_C_NPQ_bnab", "mean_NC_NPQ_bnab", "mean_C_NPQ_prer", "mean_NC_NPQ_prer", "median_C_NPQ_preati", "median_NC_NPQ_preati", "median_C_NPQ_bnab", "median_NC_NPQ_bnab", "median_C_NPQ_prer", "median_NC_NPQ_prer")
res_df = res_df[, c(mean_and_median_cols, setdiff(colnames(res_df), mean_and_median_cols))]

## Save interaction csv =======================================================

# Make separate dataframes for results with any p or bh significance
p_cols_to_keep = c("pop", "p_all", "p_preati_bnab", "p_bnab_prer", "p_preati_prer")
bh_cols_to_keep = c("pop", "bh_all", "bh_preati_bnab", "bh_bnab_prer", "bh_preati_prer")
res_any_p_sig  = res_df[which(res_df$p_all  < 0.05 | res_df$p_preati_bnab  < 0.05 | res_df$p_bnab_prer  < 0.05 | res_df$p_preati_prer  < 0.05), p_cols_to_keep]
res_any_bh_sig = res_df[which(res_df$bh_all < 0.05 | res_df$bh_preati_bnab < 0.05 | res_df$bh_bnab_prer < 0.05 | res_df$bh_preati_prer < 0.05), bh_cols_to_keep]

# Write results
write.csv(res_df, file=paste0(out_dir, "stats_detailed_time-C-interact.csv"), row.names=FALSE)
write.csv(res_any_p_sig, file=paste0(out_dir, "stats_any_p_sig_time-C-interact.csv"), row.names=FALSE)
write.csv(res_any_bh_sig, file=paste0(out_dir, "stats_any_bh_sig_time-C-interact.csv"), row.names=FALSE)
