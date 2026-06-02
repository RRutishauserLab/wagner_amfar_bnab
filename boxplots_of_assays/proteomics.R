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

#' Custom ggplot2 theme for publication-quality figures
theme_classic_pub = function(title_font_size=6, title_font_face="plain", axis_font_size=6, axis_font_face="plain", axis_line_width=0.5, legend_shape_text_space=0.25) {
  theme_classic(base_size=6, base_line_size=0.25, base_family='Helvetica') %+replace% theme(panel.grid.major=element_line(linetype="dashed", color='grey90', linewidth=axis_line_width*lwd_converter), panel.grid.minor=element_blank(), axis.title=element_text(size=axis_font_size, face=axis_font_face), axis.text=element_text(size=axis_font_size), plot.title=element_text(size=title_font_size, face=title_font_face, margin=margin(0,0,3,0)), legend.title=element_text(size=6), legend.text=element_text(size=6, margin=margin(l=legend_shape_text_space, unit="pt")), legend.background=element_rect(fill='transparent', color='transparent'), axis.line=element_line(linewidth=axis_line_width*lwd_converter, lineend='round'), plot.background=element_blank(), panel.background=element_blank(), strip.text=element_text(size=6))
}

#' Function to create a single plot for one protein
singlePlot = function(protein, df, pid_col_df, base_name, dpi=900, file_ext='png') {
  # Subset the data to the input protein
  df_sub = df[which(df$Target == protein),]
  
  # Set the beeswarmed coordinates of the points
  df_sub$timepoint = droplevels(df_sub$timepoint)
  df_sub$timepoint_num = as.numeric(df_sub$timepoint) - 1
  df_sub$x = NA
  for (time_pt in unique(df_sub$timepoint_num)) {
    df_sub_tmp = data.frame(x=as.numeric(df_sub$timepoint_num[which(df_sub$timepoint_num == time_pt)]), y=df_sub[which(df_sub$timepoint_num == time_pt), "NPQ"])
    df_sub$x[which(df_sub$timepoint_num == time_pt)] = ggbeeswarm:::offset_beeswarm(df_sub_tmp, yLim.expand=1, xRange=2, yRange=max(df_sub_tmp$y)*0.6, cex=5)$x + 1
  }
  df_paired = df_sub %>% group_by(PID_New) %>% filter(n()==2) %>% tidyr::pivot_wider(id_cols=PID_New, names_from=timepoint_num, values_from=c(x, NPQ)) %>% ungroup() %>% as.data.frame()

  # Subset the data to the PIDs that should be filled vs not filled
  df_sub_pid_fill    = df_sub[which(df_sub$PID_New %in% pid_col_df$PID[ pid_col_df$do_fill]),]
  df_sub_pid_no_fill = df_sub[which(df_sub$PID_New %in% pid_col_df$PID[!pid_col_df$do_fill]),]

  # Set plotting area size
  width_per_time = 0.4
  num_time = length(levels(df_sub$timepoint))
  plotting_area_width = width_per_time*num_time
  plotting_area_height = 1.25

  # Plot title
  plot_title = protein
  plot_title = gsub("TNFRSF4", "OX40", plot_title)
  plot_title = gsub("TNFRSF13C", "BAFF-R", plot_title)

  # Signficance brackets
  
  # get the max value per timepoint in wide format
  # df_sub_sig = df_sub %>% group_by(timepoint) %>% summarise(NPQ = max(NPQ)) %>% ungroup() %>% tidyr::pivot_wider(df_sub_sig, names_from=timepoint, values_from=max_NPQ)
  df_sub$timepoint2 = droplevels(df_sub$timepoint)
  df_sub_sig = data.frame(max=max(df_sub$NPQ*1.015), Pre=levels(df_sub$timepoint2)[1], Post=levels(df_sub$timepoint2)[2])
  df_sub_sig$stat_label = "*"
  if (plot_title %in% c("LAMP3", "BAFF-R")) { df_sub_sig$stat_label = "**"  }
  if (plot_title %in% c("ICOSLG"))          { df_sub_sig$stat_label = "***" }
  # df_sub_sig = data.frame(timepoint=levels(df_sub$timepoint2), NPQ=tapply(df_sub$NPQ, df_sub$timepoint2, max))

  # Make the plot
  p = ggplot(df_sub, aes(x=timepoint, y=NPQ)) + 
        geom_boxplot(color='gray30', fill="#ffffff99", alpha=0.75, outliers=F, width=0.5, linewidth=0.5*lwd_converter) +
        geom_segment(data=df_paired, aes(x=x_0, xend=x_1, y=NPQ_0, yend=NPQ_1, color=as.character(PID_New)), linetype='solid', alpha=0.4, linewidth=0.5*lwd_converter) + 
        geom_point(data=df_sub_pid_fill,    aes(x=x, fill=PID_New,  shape=PID_New), alpha=0.9, size=1.25, stroke=0.5*0.753, color="black") + 
        geom_point(data=df_sub_pid_no_fill, aes(x=x, color=PID_New, shape=PID_New), alpha=0.9, size=1.25, stroke=0.5*0.753) + 
        ggpubr::geom_bracket(data=df_sub_sig, aes(xmin=Pre, xmax=Post, y.position=max, label=stat_label), color='black', size=0.5*lwd_converter, label.size=6/.pt, vjust=0.5) + 
        # scale_y_continuous(limits = c(10.5, 15.9)) + 
        scale_fill_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
        scale_color_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
        scale_shape_manual(values=setNames(pid_col_df$shape, pid_col_df$PID)) +
        ggh4x::force_panelsizes(cols=unit(plotting_area_width, "in"), rows=unit(plotting_area_height, "in")) +
        labs(x=NULL, y="NPQ", title=plot_title) + 
        theme_classic_pub() + theme(panel.grid.major.x=element_blank(), legend.position='none', panel.grid.minor.y=element_line(linetype="dashed", color='grey90', linewidth=0.5*lwd_converter))
        # guides(y=guide_axis(minor.ticks=T))

  # Save the plot
  ggsave(paste0(base_name, protein, ".", file_ext), plot=p, width=plotting_area_width+0.4, height=plotting_area_height+0.4, dpi=my_dpi)
  return(p)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main ========================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Setup ======================================================================

# Set file paths
# in_dir = "/Users/jwagner1/Desktop/"
# out_dir = "/Users/jwagner1/Desktop/"
in_dir = "/c4/home/ggruenh/re/julia/amfar/proteomics/data/"
out_dir = "/c4/home/ggruenh/re/julia/amfar/proteomics/results/"
file_ext = "pdf"

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
nulisa$PID_New = plyr::revalue(nulisa$PID, c("1292"="35933", "1408"="42894", "2251"="10227", "2523"="60610", "2594"="47019", "2908"="71190", "8012"="72210", "8092"="30988", "8106"="31092"))
nulisa$PID_New = factor(nulisa$PID_New, levels=c("31092", "35933", "42894", "10227", "30988", "47019", "71190", "72210", "60610"))

# Subset data
timepoints_to_use = c("preATI", "ATI1", "ATI2", "preR")
proteins_to_use = c("LAMP3","IL18R1", "CCL21", "ICOSLG", "TNFRSF13C", "TNFRSF4")
nulisa = nulisa[which(nulisa$timepoint %in% timepoints_to_use),]                # subset to timepoints of interest
nulisa = nulisa[which(nulisa$Target %in% proteins_to_use),]                     # subset to proteins of interest
nulisa = nulisa[which(nulisa$PID != "2887"),]                                   # remove PID 2887
nulisa = nulisa[which(nulisa$PID_timepoint != "8012_ATI2"),]                    # remove PID 8012 at ATI2
nulisa = nulisa[which( !(nulisa$timepoint == "ATI1" & nulisa$PID != "8012") ),] # remove the ATI1 timepoint for all PIDs except 8012
# table(nulisa$PID, nulisa$timepoint) # double-check that all the correct samples were removed

# Recode timepoints
nulisa$timepoint[which(nulisa$timepoint %in% c("ATI1", "ATI2"))] = "bNAb"
nulisa$timepoint = factor(nulisa$timepoint, levels=c("preATI", "bNAb", "preR"))

# Make separate dataframes for preATI vs bNAb and preATI vs preR
nulisa_preati_vs_bnab = nulisa[which(nulisa$timepoint %in% c("preATI", "bNAb")),]
nulisa_preati_vs_bnab = nulisa_preati_vs_bnab[which(nulisa_preati_vs_bnab$PID_New != "31092"),] # remove PID 31092 from preATI vs bNAb plot
nulisa_preati_vs_bnab = nulisa_preati_vs_bnab[which(nulisa_preati_vs_bnab$Target %in% c("LAMP3", "IL18R1", "CCL21")),] # only plot these three proteins in preATI vs bNAb plot
nulisa_preati_vs_prer = nulisa[which(nulisa$timepoint %in% c("preATI", "preR")),]
nulisa_preati_vs_prer = nulisa_preati_vs_prer[which(nulisa_preati_vs_prer$PID_New != "60610"),] # remove PID 60610 from preATI vs preR plot
nulisa_preati_vs_prer = nulisa_preati_vs_prer[which(nulisa_preati_vs_prer$Target %in% c("ICOSLG", "TNFRSF13C", "TNFRSF4")),] # only plot these three proteins in preATI vs preR plot

## Plotting ===================================================================

# Turn showtext on for Helvetica font in png's
my_dpi = 900
showtext::showtext_auto()
showtext::showtext_opts(dpi=my_dpi)
lwd_converter = 15/32 # multiple the number you want in pts by this in order to get the correct size in ggplot2

# Set the colors and shapes for the PIDs
pid_col_df = data.frame(
  PID=c('31092', '35933', '42894', '60610', '10227', '30988', '47019', '48758', '71190', '72210'),
  old_pid=c('8106', '1292', '1408', '2523', '2251', '8092', '2594', '2887', '2908', '8012'),
  group=c('NC', 'NC', 'NC', "C", 'C', 'C', 'C', 'C', 'C', 'C'),
  do_fill=c(T, T, T, T, T, T, T, T, T, F),
  color=c('#FF7e79', '#8d1000', '#FF2600','#929292', '#137ee3', '#4fafdd', '#005493', '#0433ff', '#5c5f8b', '#08b1ff'),
  shape=c(24,21,21,24,22,24,22,22,21,0)
)

# Make preATI vs bNAb plots
plot_list_preati_vs_bnab = bettermc::mclapply(unique(nulisa_preati_vs_bnab$Target), singlePlot, df=nulisa_preati_vs_bnab, pid_col_df=pid_col_df, base_name=paste0(out_dir, "preATI_vs_bnab_protein_"), dpi=my_dpi, file_ext=file_ext, mc.cores=n_cores, mc.stdout='ignore', mc.messages='ignore', mc.allow.fatal=TRUE, mc.allow.error=TRUE, mc.fail.early=TRUE, mc.dump.frames="no", mc.progress=TRUE)
p_preati_vs_bnab = cowplot::plot_grid(plotlist=plot_list_preati_vs_bnab, ncol=3)
ggsave(paste0(out_dir, "preATI_vs_bnab_protein_all.", file_ext), plot=p_preati_vs_bnab, width=4, height=1.5, dpi=my_dpi)

# Make preATI vs preR plots
plot_list_preati_vs_prer = bettermc::mclapply(unique(nulisa_preati_vs_prer$Target), singlePlot, df=nulisa_preati_vs_prer, pid_col_df=pid_col_df, base_name=paste0(out_dir, "preati_vs_prer_protein_"), dpi=my_dpi, file_ext=file_ext, mc.cores=n_cores, mc.stdout='ignore', mc.messages='ignore', mc.allow.fatal=TRUE, mc.allow.error=TRUE, mc.fail.early=TRUE, mc.dump.frames="no", mc.progress=TRUE)
p_preati_vs_prer = cowplot::plot_grid(plotlist=plot_list_preati_vs_prer, ncol=3)
ggsave(paste0(out_dir, "preati_vs_prer_protein_all.", file_ext), plot=p_preati_vs_prer, width=4, height=1.5, dpi=my_dpi)

# # Make just one plot as test
# singlePlot("LAMP3", df=nulisa_preati_vs_bnab, pid_col_df=pid_col_df, base_name=paste0(out_dir, "preATI_vs_bnab_protein_"), dpi=my_dpi)
# singlePlot("CLEC4A", df=nulisa_preati_vs_prer, pid_col_df=pid_col_df, base_name=paste0(out_dir, "preati_vs_prer_protein_"), dpi=my_dpi)

uploadFile = function(source_path, dest_path, regex=NULL) {
  rclonev = "/c4/home/ggruenh/bin/miniforge3/envs/scpy/bin/rclone copy -v -P --no-check-certificate"
  if (!is.null(regex)) {
    rclonev = paste0(rclonev, " --include ", '"', regex, '"')
  }
  system(paste0("bash -c '", rclonev, " ", source_path, " ", dest_path, "'"))
}
uploadFile(out_dir, "box:1.Rutishauser_Lab/_Rutishauser_Lab_Projects/HIV/HIV_amfAR_trial/Analysis/JW_bNAb\\ Paper\\ Analyses/Paper_Figures/Data\\ for\\ GG/Proteomics/v3_03162026/")
