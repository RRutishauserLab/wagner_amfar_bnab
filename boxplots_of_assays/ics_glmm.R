#' 
#' @author George Gruenhagen
#' @dateCreated 02/03/2026

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
modelData = function(pop, raw_df) {
  # Select data for population
  raw_df = raw_df[which(raw_df$pop == pop),]

  #
  raw_df$prop = raw_df$pct / 100
  raw_df$WEIGHTS = raw_df$COUNT_P1
  if (grepl("CD4|CD8|B", pop)) { raw_df$WEIGHTS = raw_df$COUNT_P2 }
  
  # Fit models
  full_model    = suppressMessages(suppressWarnings( lme4::glmer(pct ~ timepoint + (1|PID_New), data=raw_df) ))
  reduced_model = suppressMessages(suppressWarnings( update(full_model, .~. - timepoint) ))
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

modelDataAtTimepoints = function(timepoints, long_df, n_cores) {
  # Make a nice message that lists the timepoints being modeled
  is_all_timpoints = all(levels(long_df$timepoint) %in% timepoints)
  timepoints_str = paste0("between ", paste0(timepoints, collapse=" and "))
  timepoints_str = ifelse(is_all_timpoints, "at all timepoints", timepoints_str)
  message(paste0("Modeling data ", timepoints_str, "..."))

  # Subset data by specified timepoints
  df = long_df[which(long_df$timepoint %in% timepoints),]
  df$timepoint = droplevels(df$timepoint)

  # Model the data
  res_df_list = bettermc::mclapply(unique(df$pop), modelData, raw_df=df, mc.stdout='ignore', mc.messages='ignore', mc.allow.fatal=TRUE, mc.allow.error=TRUE, mc.fail.early=TRUE, mc.dump.frames="no", mc.progress=TRUE)
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

#' Custom ggplot2 theme for publication-quality figures
theme_classic_pub = function(title_font_size=6, title_font_face="plain", axis_font_size=6, axis_font_face="plain", axis_line_width=0.5, legend_shape_text_space=0.25) {
  theme_classic(base_size=6, base_line_size=0.25, base_family='Helvetica') %+replace% theme(panel.grid.major=element_line(linetype="dashed", color='grey90', linewidth=axis_line_width*lwd_converter), panel.grid.minor=element_blank(), axis.title=element_text(size=axis_font_size, face=axis_font_face), axis.text=element_text(size=axis_font_size), plot.title=element_text(size=title_font_size, face=title_font_face, margin=margin(0,0,3,0)), legend.title=element_text(size=6), legend.text=element_text(size=6, margin=margin(l=legend_shape_text_space, unit="pt")), legend.background=element_rect(fill='transparent', color='transparent'), axis.line=element_line(linewidth=axis_line_width*lwd_converter, lineend='round'), plot.background=element_blank(), panel.background=element_blank(), strip.text=element_text(size=6))
}

#' Plot data for a given population
plotData = function(pop, res_df, raw_df, out_path, dpi=dpi, stat_to_display="bh_stars") {
  # Subset dataframes on the population
  res_df_i = res_df[which(res_df$pop == pop),]
  df = raw_df[which(raw_df$pop == pop),]
  
  # Set plotting area size
  width_per_time = 0.3
  num_time = length(levels(df$timepoint))
  plotting_area_width = width_per_time*num_time
  plotting_area_height = 1

  # Set the beeswarmed coordinates of the points
  df$x = NA
  for (time_pt in levels(df$timepoint)) {
    df_tmp = data.frame(x=as.numeric(df$timepoint[which(df$timepoint == time_pt)]), y=df[which(df$timepoint == time_pt), "pct"])
    df$x[which(df$timepoint == time_pt)] = ggbeeswarm:::offset_beeswarm(df_tmp, yLim.expand=1, xRange=2, yRange=max(df_tmp$y)*0.6, cex=7)$x
  }

  # Get the paired lines data
  df_paired_bnab_preati = df %>% group_by(PID_New) %>% filter(timepoint %in% c("preATI", "bNAb")) %>% filter(n() == 2) %>% reshape2::dcast(PID_New ~ timepoint, value.var="pct") %>% ungroup()
  df_paired_bnab_preati$Pre_x  = df$x[match(paste0(df_paired_bnab_preati$PID_New, "_preATI"), paste0(df$PID_New, "_", df$timepoint))]
  df_paired_bnab_preati$Post_x = df$x[match(paste0(df_paired_bnab_preati$PID_New, "_bNAb"),   paste0(df$PID_New, "_", df$timepoint))]
  
  df_paired_prer_bnab = df %>% group_by(PID_New) %>% filter(timepoint %in% c("bNAb", "preR")) %>% filter(n() == 2) %>% reshape2::dcast(PID_New ~ timepoint, value.var="pct") %>% ungroup()
  df_paired_prer_bnab$Pre_x  = df$x[match(paste0(df_paired_prer_bnab$PID_New, "_bNAb"), paste0(df$PID_New, "_", df$timepoint))]
  df_paired_prer_bnab$Post_x = df$x[match(paste0(df_paired_prer_bnab$PID_New, "_preR"), paste0(df$PID_New, "_", df$timepoint))]
  
  df_paired_prer_preati = df %>% group_by(PID_New) %>% filter(timepoint %in% c("preATI", "bNAb", "preR")) %>% filter(n() == 2) 
  if (nrow(df_paired_prer_preati) > 0) { 
    df_paired_prer_preati= df_paired_prer_preati %>% reshape2::dcast(PID_New ~ timepoint, value.var="pct") %>% filter(is.na(bNAb)) %>% ungroup()
    df_paired_prer_preati$Pre_x  = df$x[match(paste0(df_paired_prer_preati$PID_New, "_BL"),     paste0(df$PID_New, "_", df$timepoint))]
    df_paired_prer_preati$Post_x = df$x[match(paste0(df_paired_prer_preati$PID_New, "_preR"),   paste0(df$PID_New, "_", df$timepoint))]
  } else ( df_paired_prer_preati = data.frame(PID_New=character(0), Pre_x=character(0), Post_x=character(0), preATI=numeric(0), preR=numeric(0)) )

  # Define colors and shapes for PIDs
  pid_col_df = data.frame(
    PID=c('31092', '35933', '42894', '60610', '10227', '30988', '47019', '48758', '71190', '72210'),
    old_pid=c('8106', '1292', '1408', '2523', '2251', '8092', '2594', '2887', '2908', '8012'),
    group=c('NC', 'NC', 'NC', "C", 'C', 'C', 'C', 'C', 'C', 'C'),
    do_fill=c(T, T, T, T, T, T, T, T, T, F),
    color=c('#FF7e79', '#8d1000', '#FF2600','#929292', '#137ee3', '#4fafdd', '#005493', '#0433ff', '#5c5f8b', '#08b1ff'),
    shape=c(24,21,21,24,22,24,22,22,21,0)
  )

  # Get the stats for significance brackets
  base_multiplaction_factor = 0.065
  pct_max = max(df$pct, na.rm=T)
  sig_df_pair = data.frame(tp_1 = c("preATI", "bNAb", "preATI"),
                           tp_2 = c("bNAb",   "preR", "preR"),
                           stat_label=unname(t(res_df_i[,paste0(stat_to_display, c("_preati_bnab", "_bnab_prer", "_preati_prer"))])))
  # sig_df_all = data.frame(tp_1=c("BL"), tp_2=c("PreR"), stat_label=as.vector(res_df_i[,paste0(stat_to_display, c("_all"))]))
  sig_df_pair = sig_df_pair[which(sig_df_pair$stat_label != "ns"), ]
  # sig_df_all = sig_df_all[which(sig_df_all$stat_label != "ns"), ]
  sig_df_pair$row_num = seq_len(nrow(sig_df_pair))
  sig_df_pair$multiplaction_factor = 1 + sig_df_pair$row_num*base_multiplaction_factor
  sig_df_pair$y = pct_max * sig_df_pair$multiplaction_factor
  # if (nrow(sig_df_all) > 0) {
  #   pair_max = ifelse(nrow(sig_df_pair) > 0, max(sig_df_pair$y), pct_max)
  #   sig_df_all$y = pair_max + (base_multiplaction_factor * pct_max)
  # }

  # Clean up the population name for plotting
  pop_clean = pop
  pop_first  = strsplit(pop_clean, "_")[[1]][1]
  pop_middle = strsplit(pop_clean, "_")[[1]][2]
  pop_last   = strsplit(pop_clean, "_")[[1]][3]
  if (pop_first == "IFNG") {
    pop_clean = paste0(pop_last, "-specific\nof ", pop_middle, "+ T cells (%)")
    # pop_clean = paste0("%", pop_last, "-specific (IFN\u03B3+)\nof ", pop_middle, "+ T cells")
    # pop_clean = bquote(atop(NA, atop(textstyle(paste("%", .(pop_last), "-specific (IFN", gamma, "+)")), textstyle(paste("of ", .(pop_middle), "+ T cells")))))
    pop_clean = bquote(atop(NA, textstyle(paste("%", .(pop_last), "-specific (IFN", gamma, "+)"))))
    plot_title = paste0(pop_middle, "+ T cells")
  } else { plot_title = "Gag-specific CD8+ T cells"; pop_clean = "%PD-1+" }

  # Determine the plot title and y_axis title
  y_axis_title = pop_clean

  # Subset the data to the PIDs that should be filled vs not filled
  df_pid_fill    = df[which(df$PID_New %in% pid_col_df$PID[ pid_col_df$do_fill]),]
  df_pid_no_fill = df[which(df$PID_New %in% pid_col_df$PID[!pid_col_df$do_fill]),]

  # Make the plot
  p = ggplot(df, ggplot2::aes(x=timepoint, y=10^pct)) +
        geom_boxplot(color='gray30', fill="#ffffff99", alpha=0.75, outliers=F, width=0.5, linewidth=0.5*lwd_converter) +
        geom_segment(data=df_paired_bnab_preati, aes(x=Pre_x, xend=Post_x, y=10^preATI, yend=10^bNAb, color=as.character(PID_New)), linetype='solid', alpha=0.4, linewidth=0.5*lwd_converter) + 
        geom_segment(data=df_paired_prer_bnab,   aes(x=Pre_x, xend=Post_x, y=10^bNAb,   yend=10^preR, color=as.character(PID_New)), linetype='solid', alpha=0.4, linewidth=0.5*lwd_converter) +
        geom_segment(data=df_paired_prer_preati, aes(x=Pre_x, xend=Post_x, y=10^preATI, yend=10^preR, color=as.character(PID_New)), linetype='solid', alpha=0.4, linewidth=0.5*lwd_converter) +
        geom_point(data=df_pid_fill,    aes(x=x, fill=as.character(PID_New),  shape=as.character(PID_New)), alpha=0.8, size=1, stroke=0.5*0.753, color='black') +
        geom_point(data=df_pid_no_fill, aes(x=x, color=as.character(PID_New), shape=as.character(PID_New)), alpha=0.8, size=1, stroke=0.5*0.753) +
        scale_fill_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
        scale_color_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
        scale_shape_manual(values=setNames(pid_col_df$shape, pid_col_df$PID)) +
        ggh4x::force_panelsizes(cols=unit(plotting_area_width, "in"), rows=unit(plotting_area_height, "in")) +
        coord_cartesian(clip='off') +
        theme_classic_pub() + theme(panel.grid.major.x=element_blank(), legend.position="none") +
        labs(title=plot_title, x=NULL, y=y_axis_title)

  if (pop_first == "IFNG") {
    p = p + geom_hline(yintercept=1e-3, linetype='solid', color='gray60', linewidth=0.5*lwd_converter) +
    geom_text(data=data.frame(timepoint=3.4, pct=-3), label="LOD", size=6/.pt, vjust=-0.3, color="gray60") +
    scale_y_log10() +
    annotation_logticks(sides="l", outside=TRUE, short=unit(.25,"mm"), mid=unit(0.5,"mm"), long=unit(0.75,"mm"), linewidth=0.5*lwd_converter)
  }
  
  # Add significance brackets
  if (nrow(sig_df_pair) > 0) {
    p = p + ggpubr::geom_bracket(data=sig_df_pair, aes(xmin=tp_1, xmax=tp_2, y.position=y, label=stat_label), color='black', size=0.5*lwd_converter, label.size=6/.pt, vjust=0.5)
  }
  # if (nrow(sig_df_all) > 0) {
  #   p = p + ggpubr::geom_bracket(data=sig_df_all,  aes(xmin=tp_1, xmax=tp_2, y.position=y, label=stat_label), color='black', size=1.0*lwd_converter, label.size=6/.pt, vjust=0.5, tip.length=0)
  # }
  
  # Save the plot
  ggsave(filename=paste0(out_path, pop, "_boxplot.", "png"), plot=p, width=plotting_area_width+0.8, height=plotting_area_height+0.4, units="in", dpi=dpi)
  ggsave(filename=paste0(out_path, pop, "_boxplot.", "pdf"), plot=p, width=plotting_area_width+0.8, height=plotting_area_height+0.4, units="in", dpi=dpi)
  return(TRUE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main ========================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Setup ======================================================================

# Start timer
start_time = proc.time()[3]

# Check for required libraries
is_installed = checkForRequiredLibraries(required_libs=c("bettermc", "dplyr", "ggh4x", "ggplot2", "ggpubr", "lme4", "pbkrtest", "plyr", "reshape2", "survey"))

# Load libraries
message("Loading libraries...")
pacman::p_load(dplyr, ggplot2, survey)

# Constants
in_file = "/c4/home/ggruenh/re/julia/amfar/ics_lm/data/20260212_ICS_Data_for_GG.csv"
out_path = "/c4/home/ggruenh/re/julia/amfar/ics_lm/results/log10/"
# out_path = "/c4/home/ggruenh/re/julia/amfar/ics_lm/results/raw/" # not doing raw since they want logticks
n_cores = 5
use_log10 = TRUE
dirs_to_check = c(out_path, paste0(out_path, "/plots_p"), paste0(out_path, "/plots_bh"))
for (dir_i in dirs_to_check) {
  if (!dir.exists(dir_i)) { dir.create(dir_i, recursive=TRUE) }
}

# Read in data
raw_df = data.table::fread(in_file, data.table=F)
raw_df = raw_df[,which(!colnames(raw_df) %in% c("PostR_Info", "Date", "Notes"))]
long_df = reshape2::melt(raw_df, id.vars=c("PID", "PID_NEW", "Timepoint"), variable.name="pop", value.name="pct")
timepoints_to_use = c("preATI", "bNAb", "preR")
long_df = long_df[which(long_df$Timepoint %in% timepoints_to_use),]
long_df$Timepoint = factor(long_df$Timepoint, levels=timepoints_to_use)
long_df$timepoint = long_df$Timepoint
long_df$PID = as.character(long_df$PID)
long_df$PID_New = as.character(long_df$PID_NEW)
long_df = long_df[, c("PID", "PID_New", "timepoint", "pop", "pct")]
long_df$pct_raw = long_df$pct
long_df = long_df[which(!is.na(long_df$pct)),]
# table(long_df$PID_New, long_df$timepoint) # check that we have data for all PIDs at all timepoints

# Take log10 of the percentages if specified
if (use_log10) {
  message("Taking log10 of the percentages...")
  long_df$pct = log10(long_df$pct)
}

# # Exclude populations that are <0.5% in >10% of samples at any timepoint
# pops_to_exclude = c()
# for (pop_i in unique(long_df$pop)) {
#   pop_i_df = long_df[which(long_df$pop == pop_i),]
#   prop_low = sapply(levels(pop_i_df$timepoint), function(tp) { mean(pop_i_df$pct[which(pop_i_df$timepoint == tp)] < 0.5, na.rm=T) })
#   if (any(prop_low > 0.1)) {
#     pops_to_exclude = c(pops_to_exclude, pop_i)
#   }
# }

## Stats ======================================================================

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
res_df$mean_pct_preati = sapply(res_df$pop, function(pop) { mean(long_df$pct_raw[which(long_df$pop == pop & long_df$timepoint == "preATI")], na.rm=T) })
res_df$mean_pct_bnab   = sapply(res_df$pop, function(pop) { mean(long_df$pct_raw[which(long_df$pop == pop & long_df$timepoint == "bNAb")],   na.rm=T) })
res_df$mean_pct_prer   = sapply(res_df$pop, function(pop) { mean(long_df$pct_raw[which(long_df$pop == pop & long_df$timepoint == "preR")],   na.rm=T) })
res_df$median_pct_preati = sapply(res_df$pop, function(pop) { median(long_df$pct_raw[which(long_df$pop == pop & long_df$timepoint == "preATI")], na.rm=T) })
res_df$median_pct_bnab   = sapply(res_df$pop, function(pop) { median(long_df$pct_raw[which(long_df$pop == pop & long_df$timepoint == "bNAb")],   na.rm=T) })
res_df$median_pct_prer   = sapply(res_df$pop, function(pop) { median(long_df$pct_raw[which(long_df$pop == pop & long_df$timepoint == "preR")],   na.rm=T) })
res_df = res_df[, c("pop", "mean_pct_preati", "mean_pct_bnab", "mean_pct_prer", "median_pct_preati", "median_pct_bnab", "median_pct_prer", setdiff(colnames(res_df), c("pop", "mean_pct_preati", "mean_pct_bnab", "mean_pct_prer", "median_pct_preati", "median_pct_bnab", "median_pct_prer")))]

## Save csv ===================================================================

# Make separate dataframes for results with any p or bh significance
p_cols_to_keep = c("pop", "mean_pct_preati", "mean_pct_bnab", "mean_pct_prer", "median_pct_preati", "median_pct_bnab", "median_pct_prer", "p_all", "p_preati_bnab", "p_bnab_prer", "p_preati_prer")
bh_cols_to_keep = c("pop", "mean_pct_preati", "mean_pct_bnab", "mean_pct_prer", "median_pct_preati", "median_pct_bnab", "median_pct_prer", "bh_all", "bh_preati_bnab", "bh_bnab_prer", "bh_preati_prer")
res_any_p_sig  = res_df[which(res_df$p_all  < 0.05 | res_df$p_preati_bnab  < 0.05 | res_df$p_bnab_prer  < 0.05 | res_df$p_preati_prer  < 0.05), p_cols_to_keep]
res_any_bh_sig = res_df[which(res_df$bh_all < 0.05 | res_df$bh_preati_bnab < 0.05 | res_df$bh_bnab_prer < 0.05 | res_df$bh_preati_prer < 0.05), bh_cols_to_keep]

# Write results
write.csv(res_df, file=paste0(out_path, "stats_detailed.csv"), row.names=FALSE)
write.csv(res_any_p_sig, file=paste0(out_path, "stats_any_p_sig.csv"), row.names=FALSE)
write.csv(res_any_bh_sig, file=paste0(out_path, "stats_any_bh_sig.csv"), row.names=FALSE)

## Plot =======================================================================

# Turn showtext on for Helvetica font in png's
my_dpi = 900
# showtext::showtext_auto()
# showtext::showtext_opts(dpi=my_dpi)
lwd_converter = 15/32 # multiple the number you want in pts by this in order to get the correct size in ggplot2

# Generate plots
message("Generating plots with p-values...")
plot_success = bettermc::mclapply(as.vector(res_df$pop), plotData, res_df=res_df, raw_df=long_df, out_path=paste0(out_path, "/plots_p/"), dpi=my_dpi, stat_to_display="p_stars", mc.cores=n_cores, mc.stdout='ignore', mc.messages='ignore', mc.allow.fatal=TRUE, mc.allow.error=TRUE, mc.fail.early=TRUE, mc.dump.frames="no", mc.progress=TRUE)

message("Generating plots with BH-adjusted p-values...")
plot_success = bettermc::mclapply(as.vector(res_df$pop), plotData, res_df=res_df, raw_df=long_df, out_path=paste0(out_path, "/plots_bh/"), dpi=my_dpi, stat_to_display="bh_stars", mc.cores=n_cores, mc.stdout='ignore', mc.messages='ignore', mc.allow.fatal=TRUE, mc.allow.error=TRUE, mc.fail.early=TRUE, mc.dump.frames="no", mc.progress=TRUE)

uploadFile = function(source_path, dest_path, regex=NULL) {
  rclonev = "/c4/home/ggruenh/bin/miniforge3/envs/scpy/bin/rclone copy -v -P --no-check-certificate"
  if (!is.null(regex)) {
    rclonev = paste0(rclonev, " --include ", '"', regex, '"')
  }
  system(paste0("bash -c '", rclonev, " ", source_path, " ", dest_path, "'"))
}
uploadFile("/c4/home/ggruenh/re/julia/amfar/ics_lm/results/", "box:1.Rutishauser_Lab/_Rutishauser_Lab_Projects/HIV/HIV_amfAR_trial/Analysis/JW_bNAb\\ Paper\\ Analyses/Paper_Figures/Data\\ for\\ GG/ics/v3_04292026/")
