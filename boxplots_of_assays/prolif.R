#' Paired Wilcox test of preR vs preATI CD4 & CD8 proliferation data
#' 
#' @author George Gruenhagen
#' @dateCreated 05/12/2026

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
singlePlot = function(ct, df, pid_col_df, base_name, dpi=900, file_ext='pdf') {
  # Subset the data to the input protein
  df_sub = df[which(df$celltype == ct),]
  
  # Set the beeswarmed coordinates of the points
  df_sub$timepoint_num = as.numeric(df_sub$timepoint) - 1
  df_sub$x = NA
  for (time_pt in unique(df_sub$timepoint_num)) {
    df_sub_tmp = data.frame(x=as.numeric(df_sub$timepoint_num[which(df_sub$timepoint_num == time_pt)]), y=df_sub[which(df_sub$timepoint_num == time_pt), "value"])
    df_sub$x[which(df_sub$timepoint_num == time_pt)] = ggbeeswarm:::offset_beeswarm(df_sub_tmp, yLim.expand=1, xRange=2, yRange=max(df_sub_tmp$y)*0.6, cex=5)$x + 1
  }
  df_paired = df_sub %>% group_by(PID_New) %>% filter(n()==2) %>% tidyr::pivot_wider(id_cols=PID_New, names_from=timepoint_num, values_from=c(x, value)) %>% ungroup() %>% as.data.frame()

  # Subset the data to the PIDs that should be filled vs not filled
  df_sub_pid_fill    = df_sub[which(df_sub$PID_New %in% pid_col_df$PID[ pid_col_df$do_fill]),]
  df_sub_pid_no_fill = df_sub[which(df_sub$PID_New %in% pid_col_df$PID[!pid_col_df$do_fill]),]

  # Set plotting area size
  width_per_time = 0.4
  num_time = length(levels(df_sub$timepoint))
  plotting_area_width = width_per_time*num_time
  plotting_area_height = 1

  # Plot title
  plot_title = paste0(ct, "+ T cells")
  # y_axis_title = bquote(atop(NA, atop(textstyle('%CTV'^lo*' of non-na\u00efve'), atop(textstyle(.(ct)*'+ T cells'), textstyle('(Gag+Pol+Nef+Env)')))))
  y_axis_title = bquote(atop(NA, atop(textstyle('%CTV'^lo*' of non-na\u00efve'), textstyle('(Gag+Pol+Nef+Env)'))))

  # Make the plot
  p = ggplot(df_sub, aes(x=timepoint, y=value)) + 
        geom_boxplot(color='gray30', fill="#ffffff99", alpha=0.75, outliers=F, width=0.5, linewidth=0.5*lwd_converter) +
        geom_segment(data=df_paired, aes(x=x_0, xend=x_1, y=value_0, yend=value_1, color=as.character(PID_New)), linetype='solid', alpha=0.4, linewidth=0.5*lwd_converter) + 
        geom_point(data=df_sub_pid_fill,    aes(x=x, fill=PID_New,  shape=PID_New), alpha=0.9, size=1.25, stroke=0.5*0.753, color="black") + 
        geom_point(data=df_sub_pid_no_fill, aes(x=x, color=PID_New, shape=PID_New), alpha=0.9, size=1.25, stroke=0.5*0.753) + 
        scale_fill_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
        scale_color_manual(values=setNames(pid_col_df$color, pid_col_df$PID)) +
        scale_shape_manual(values=setNames(pid_col_df$shape, pid_col_df$PID)) +
        ggh4x::force_panelsizes(cols=unit(plotting_area_width, "in"), rows=unit(plotting_area_height, "in")) +
        labs(x=NULL, y=y_axis_title, title=plot_title) + 
        theme_classic_pub() + theme(panel.grid.major.x=element_blank(), legend.position='none')

  # Save the plot
  ggsave(paste0(base_name, ct, ".", file_ext), plot=p, width=plotting_area_width+0.7, height=plotting_area_height+0.4, dpi=my_dpi)
  return(p)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main ========================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Setup ======================================================================

# Set file paths
in_dir = "/c4/home/ggruenh/re/julia/amfar/prolif/data/"
out_dir = "/c4/home/ggruenh/re/julia/amfar/prolif/results/"
file_ext = "pdf"

# Set number of cores for parallelization
n_cores = 0

# Check for required libraries
is_installed = checkForRequiredLibraries(required_libs=c("bettermc", "dplyr", "ggh4x", "ggplot2", "plyr", "tidyr"))

# Load libraries
suppressWarnings(suppressMessages( library('dplyr',    quietly=T, warn.conflicts=F, verbose=F) ))
suppressWarnings(suppressMessages( library('ggplot2',  quietly=T, warn.conflicts=F, verbose=F) ))

## Data Processing ============================================================

# Read data
prolif = read.csv(paste0(in_dir, "20260512_RR_Prolif_Data_for_GG.csv"))

# Add new columns
prolif$PID = as.character(prolif$PID)
prolif$timepoint = factor(prolif$timepoint, levels=c("preATI", "preR"))
prolif$PID_New = factor(prolif$PID_New, levels=c("31092", "35933", "42894", "10227", "30988", "47019", "48758", "71190", "72210", "60610"))

# Melt the data to long format for plotting
prolif_melt = prolif %>% tidyr::pivot_longer(cols=c("CD4", "CD8"), names_to="celltype", values_to="value") %>% as.data.frame()

## Check for Normality ========================================================

# Check for normality
norm_test = shapiro.test(prolif$CD4) # p = 0.0003348, so not normal
norm_test = shapiro.test(prolif$CD8) # p = 0.06327,   so trending to not normal
norm_test = nortest::ad.test(prolif$CD4)      # p = 0.0003459, so not normal
norm_test = nortest::ad.test(prolif$CD8)      # p = 0.0553,    so trending to not normal

## Analyze ====================================================================

# Perform paired Wilcoxon test
res = prolif_melt %>% group_by(celltype) %>% summarise(p=wilcox.test(value ~ timepoint, paired=T)$p.value, mean_preATI=mean(value[which(timepoint == "preATI")]), mean_preR=mean(value[which(timepoint == "preR")]), median_preATI=median(value[which(timepoint == "preATI")]), median_preR=median(value[which(timepoint == "preR")])) %>% as.data.frame()

# Save results
write.csv(res, paste0(out_dir, "preR_vs_preATI_stats.csv"), row.names=F)

## Plotting ===================================================================

# Turn showtext on for Helvetica font in png's
my_dpi = 900
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

# Make CD4 plot
p = singlePlot(ct="CD4", df=prolif_melt, pid_col_df=pid_col_df, base_name=paste0(out_dir, "preR_vs_preATI_"), dpi=my_dpi, file_ext=file_ext)

# Make CD8 plot
p = singlePlot(ct="CD8", df=prolif_melt, pid_col_df=pid_col_df, base_name=paste0(out_dir, "preR_vs_preATI_"), dpi=my_dpi, file_ext=file_ext)

uploadFile = function(source_path, dest_path, regex=NULL) {
  rclonev = "/c4/home/ggruenh/bin/miniforge3/envs/scpy/bin/rclone copy -v -P --no-check-certificate"
  if (!is.null(regex)) {
    rclonev = paste0(rclonev, " --include ", '"', regex, '"')
  }
  system(paste0("bash -c '", rclonev, " ", source_path, " ", dest_path, "'"))
}
uploadFile(out_dir, "box:1.Rutishauser_Lab/_Rutishauser_Lab_Projects/HIV/HIV_amfAR_trial/Analysis/JW_bNAb\\ Paper\\ Analyses/Paper_Figures/Data\\ for\\ GG/Proliferation/v2_05132026/")
