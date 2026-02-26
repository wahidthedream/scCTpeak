##############################################################################################################################################
### Figure 2 A, B, D, E
### Peak Count and Peak Width Analysis for Real HumanPBMC data With and With Input Case 
### for without input case we need to ~/result_frip/HumanPBMC_combined_peak_frip_summary_withot_input.csv and save to corressponding figures
##############################################################################################################################################


# Load required packages
library(readr)      # For reading CSV
library(dplyr)      # For data manipulation
library(knitr)      # For kable
library(kableExtra) # For styling table
library(ggplot2)
library(patchwork)
library(ggsci)

# Read the CSV file
data <- read_csv("~/result_frip/HumanPBMC_combined_peak_frip_summary_with_input.csv")

# Replace NA with 0
data[is.na(data)] <- 0

# Round numeric columns for better display
data <- data %>%
  mutate(
    Number_of_Peaks = round(Number_of_Peaks),
    Mean_Peak_Width = round(Mean_Peak_Width),
    Total_Reads = round(Total_Reads),
    Reads_in_Peaks = round(Reads_in_Peaks),
    FRiP = round(FRiP, 4)
  )


data$log10_Peaks <- log10(data$Number_of_Peaks + 1)  # +1 to avoid log(0)
data$log10_Peak_Width <- log10(data$Mean_Peak_Width + 1)  
# Calculate comprehensive statistics
median_data <- data %>%
  group_by(Celltype, Method, Histone) %>%
  summarize(
    median_log10_peaks = median(log10_Peaks, na.rm = TRUE),
    median_log10_Peak_Width = median(log10_Peak_Width, na.rm = TRUE),
    mean_log10_peaks = mean(log10_Peaks, na.rm = TRUE),
    mean_log10_Peak_Width = mean(log10_Peak_Width, na.rm = TRUE),
    sd_log10_peaks = sd(log10_Peaks, na.rm = TRUE),
    sd_log10_Peak_Width = sd(log10_Peak_Width, na.rm = TRUE),
    se_log10_peaks = sd_log10_peaks / sqrt(n()),
    se_log10_Peak_Width = sd_log10_Peak_Width / sqrt(n()),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Celltype = factor(Celltype),
    Method = factor(Method),
    Histone = factor(Histone)
  )

# Ensure all cell types are properly ordered
all_celltypes <- unique(median_data$Celltype)
median_data <- median_data %>%
  mutate(Celltype = factor(Celltype, levels = all_celltypes))

# Get all unique methods and create a consistent color mapping
all_methods <- unique(data$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)

# Create a named color vector for consistent method coloring
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

# Calculate consistent y-axis limits for both plots
y_min_peaks <- min(min(median_data$median_log10_peaks - median_data$se_log10_peaks, na.rm = TRUE),
                   min(data$log10_Peaks, na.rm = TRUE))
y_max_peaks <- max(max(median_data$median_log10_peaks + median_data$se_log10_peaks, na.rm = TRUE),
                   max(data$log10_Peaks, na.rm = TRUE))

y_min_width <- min(min(median_data$median_log10_Peak_Width - median_data$se_log10_Peak_Width, na.rm = TRUE),
                   min(data$log10_Peak_Width, na.rm = TRUE))
y_max_width <- max(max(median_data$median_log10_Peak_Width + median_data$se_log10_Peak_Width, na.rm = TRUE),
                   max(data$log10_Peak_Width, na.rm = TRUE))

# Add a small buffer to the y-axis limits
y_buffer_peaks <- (y_max_peaks - y_min_peaks) * 0.05
y_limits_peaks <- c(y_min_peaks - y_buffer_peaks, y_max_peaks + y_buffer_peaks)

y_buffer_width <- (y_max_width - y_min_width) * 0.05
y_limits_width <- c(y_min_width - y_buffer_width, y_max_width + y_buffer_width)

# Enhanced Nature Methods Theme with better readability
nature_methods_enhanced_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "Helvetica", color = "black", size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, 
                                margin = margin(b = 12), color = "#2C3E50", lineheight = 1.1),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "#566573",
                                  margin = margin(b = 18), lineheight = 1.2),
      axis.title = element_text(face = "bold", size = 14, color = "#2C3E50"),
      axis.title.x = element_text(margin = margin(t = 12), vjust = -0.5, size = 13),
      axis.title.y = element_text(margin = margin(r = 12), vjust = 2, size = 13),
      axis.text = element_text(size = 11, color = "#34495E"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                margin = margin(t = 5), size = 11, face = "bold"),
      axis.text.y = element_text(margin = margin(r = 5), size = 11, face = "bold"),
      axis.line = element_line(color = "#2C3E50", linewidth = 0.8),
      axis.ticks = element_line(color = "#2C3E50", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(face = "bold", size = 13, color = "#2C3E50",
                                 margin = margin(b = 8)),
      legend.text = element_text(size = 12, color = "#34495E",
                                margin = margin(r = 15)),
      legend.position = "bottom",
      legend.key.size = unit(0.6, "cm"),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width = unit(1.2, "cm"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA),
      legend.margin = margin(t = -5, b = 8),
      legend.spacing.x = unit(0.3, "cm"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(1.5, "lines"),
      panel.border = element_blank(),
      strip.background = element_rect(fill = "#F8F9F9", color = "#D5D8DC",
                                     linewidth = 0.8),
      strip.text = element_text(face = "bold", size = 12, color = "#2C3E50",
                               margin = margin(8, 8, 8, 8)),
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(fill = "white", color = NA),
      plot.caption = element_text(size = 11, color = "#7F8C8D", 
                                 margin = margin(t = 15), hjust = 1, lineheight = 1.1)
    )
}

# Create p1 - Peak Counts (Line plot) with FIXED Y-AXIS LABELS
p1 <- ggplot(median_data, aes(x = Celltype, y = median_log10_peaks)) +
  geom_line(aes(group = Method, color = Method), 
            linewidth = 1.4, alpha = 0.9, show.legend = TRUE,
            lineend = "round", linejoin = "round") +
  geom_point(aes(fill = Method), 
             shape = 21, size = 5, color = "white", stroke = 1.2,
             alpha = 0.95) +
  geom_errorbar(aes(ymin = median_log10_peaks - se_log10_peaks,
                    ymax = median_log10_peaks + se_log10_peaks,
                    color = Method),
                width = 0.3, linewidth = 1.0, alpha = 0.8) +
  facet_wrap(~ Histone, nrow = 1, scales = "free_x") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = method_colors, 
                     name = "Peak-Calling Method",
                     guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_fill_manual(values = method_colors,
                    name = "Peak-Calling Method",
                    guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_y_continuous(
    limits = c(0, 6.5), #y_limits_peaks
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title = "Comparative Performance of Peak Calling Methods Across Cell-Types",
    #subtitle = "Median log-transformed peak counts with standard error bars across histone modifications",
    x = "Cell-Type",
    y = expression(paste("log"[10], "(Peak Count + 1)"))#,  # FIXED: Using expression()
    #caption = paste("Data analysis: n =", format(nrow(data), big.mark = ","), 
    #               "total samples | Error bars represent standard error of the mean (SEM)")
  ) +
  nature_methods_enhanced_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                              margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold")
  )

# Create p2 - Peak Counts (Box plot) with FIXED Y-AXIS LABELS
p2 <- ggplot(data, aes(x = Method, y = log10_Peaks, fill = Method)) +
  geom_boxplot(outlier.size = 1.5, outlier.shape = 21, outlier.fill = "white", 
               outlier.color = "black", lwd = 0.4, alpha = 0.8) +
  facet_wrap(~ Histone, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    limits = c(0, 6.5), #y_limits_peaks,
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title = "Distribution of Peak Count by Peak Calling Method",
    #subtitle = "Boxplots showing variability in log-transformed peak counts across methods",
    x = "Peak Calling Method",
    y = expression(paste("log"[10], "(Peak Count + 1)"))  # FIXED: Using expression()
  ) +
  nature_methods_enhanced_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                              margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold")
  )

# Create p3 - Peak Width (Line plot) with FIXED Y-AXIS LABELS
p3 <- ggplot(median_data, aes(x = Celltype, y = median_log10_Peak_Width)) +
  geom_line(aes(group = Method, color = Method), 
            linewidth = 1.4, alpha = 0.9, show.legend = TRUE,
            lineend = "round", linejoin = "round") +
  geom_point(aes(fill = Method), 
             shape = 21, size = 5, color = "white", stroke = 1.2,
             alpha = 0.95) +
  geom_errorbar(aes(ymin = median_log10_Peak_Width - se_log10_Peak_Width,
                    ymax = median_log10_Peak_Width + se_log10_Peak_Width,
                    color = Method),
                width = 0.3, linewidth = 1.0, alpha = 0.8) +
  facet_wrap(~ Histone, nrow = 1, scales = "free_x") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = method_colors, 
                     name = "Peak-calling Method",
                     guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_fill_manual(values = method_colors,
                    name = "Peak-calling Method",
                    guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_y_continuous(
    limits = c(0, 4.0), #y_limits_width
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title = "Comparative Performace of Peak Calling Methods Across Cell-Types",
    #subtitle = "Median log-transformed peak width with standard error bars across histone modifications",
    x = "Cell Type",
    y = expression(paste("log"[10], "(Peak Width + 1)"))#,  # FIXED: Using expression()
    #caption = paste("Data analysis: n =", format(nrow(data), big.mark = ","), 
    #               "total samples | Error bars represent standard error of the mean (SEM)")
  ) +
  nature_methods_enhanced_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                              margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold")
  )

# Create p4 - Peak Width (Violin plot) with FIXED Y-AXIS LABELS
p4 <- ggplot(data, aes(x = Method, y = log10_Peak_Width, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.8, linewidth = 0.3, width=1) +
  geom_boxplot(width = 0.3, fill = "white", alpha = 0.7, outlier.size = 0.8, 
               outlier.color = "black", linewidth = 0.3) +
  facet_wrap(~ Histone, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    limits =c(0, 4.0),# y_limits_width,
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title = "Peak Width Distribution by Peak Calling Method",
    x = "Peak Calling Method",
    y = expression(paste("log"[10], "(Peak Width + 1)"))
  ) +
  nature_methods_enhanced_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                              margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold")
  )


# Combine plots for peak counts
# COMPACT VERSION FOR PEAK COUNTS (p1 + p2)
compact_peaks <- p1 / p2 +
  plot_annotation(
    title = "PEAK COUNT ANALYSIS",
    subtitle = "A: Performance across cell types | B: Distribution by method",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, 
                               margin = margin(b = 5), color = "#2C3E50"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10),
                                  color = "#566573"),
      plot.margin = margin(10, 10, 10, 10)
    )
  ) &
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

# COMPACT VERSION FOR PEAK WIDTH (p3 + p4)
compact_width <- p3 / p4 +
  plot_annotation(
    title = "PEAK WIDTH ANALYSIS",
    subtitle = "A: Performance across cell types | B: Distribution by method",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, 
                               margin = margin(b = 5), color = "#2C3E50"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10),
                                  color = "#566573"),
      plot.margin = margin(10, 10, 10, 10)
    )
  ) &
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Display compact plots
print(compact_peaks)
print(compact_width)

# Save compact versions
ggsave("~/result_frip/Figure/HumanPBMC_PeakCount_Analysis_With_Input_compact.pdf", plot = compact_peaks, width = 16, height = 12, device = cairo_pdf,bg = "white")

ggsave("~/result_frip/Figure/HumanPBMC_PeakWidth_Analysis_With_Input_compact.pdf", plot = compact_width, width = 16, height = 12,device = cairo_pdf, bg = "white")


###############################################################################
## Figure 2 C: PEAK Count SCATTER PLOTS human with input
###############################################################################

## Create scatter plot with median and standard deviation of LOG10-TRANSFORMED peak counts
# Calculate statistics from log10-transformed peak counts (using pre-calculated log10_Peaks)
scatter_data_log <- data %>%
  group_by(Method, Histone) %>%
  summarize(
    median_log_peak_count = median(log10_Peaks, na.rm = TRUE),  # Use pre-calculated log10_Peaks
    sd_log_peak_count = sd(log10_Peaks, na.rm = TRUE),
    # Also keep raw counts for reference
    median_raw_peak_count = median(Number_of_Peaks, na.rm = TRUE),
    sd_raw_peak_count = sd(Number_of_Peaks, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(Method),
    Histone = factor(Histone)
  )

# Create shape mapping for Histones
all_histones <- unique(scatter_data_log$Histone)
histone_shapes <- c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10)
histone_shape_map <- setNames(histone_shapes[1:length(all_histones)], all_histones)

# Create scatter plot using LOG-TRANSFORMED data


# Alternative version without text labels with custom axis limits
scatter_plot_clean <- ggplot(scatter_data_log, 
                             aes(x = sd_log_peak_count,
                                 y = median_log_peak_count,
                                 color = Method,
                                 shape = Histone)) +
  geom_point(size = 7, alpha = 0.9, stroke = 1.5) +
  scale_color_manual(values = method_colors, 
                     name = "Peak-Calling Method") +
  scale_shape_manual(values = histone_shape_map,
                     name = "Histone Modification") +
  scale_x_continuous(
    labels = scales::comma_format(),
    limits = c(0, 2.0),  # Set X-axis limits (0 to 1.5)
    breaks = seq(0, 2.0, by = 0.5),  # Set breaks every 0.25
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_y_continuous(
    labels = scales::comma_format(),
    limits = c(3, 6),    # Set Y-axis limits (3 to 6)
    breaks = seq(3, 6, by = 0.5),  # Set breaks every 0.5
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title = "Standard Deviation vs Median ",
    x = "SD of log10(Peak Count + 1)",
    y = "Median of log10(Peak Count + 1)"
  ) +
  nature_methods_enhanced_theme() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.3, "cm")
  )

print(scatter_plot_clean)

ggsave("~/result_frip/Figure/HumanPBMC_SD_vs_Median_Scatter_With_Input.pdf", plot = scatter_plot_clean, width = 12, height = 9, device = cairo_pdf, bg = "white")

###############################################################################
## Figure 2 F:PEAK WIDTHS SCATTER PLOTS human with input
## 
###############################################################################

# Calculate statistics from log10-transformed peak widths (using pre-calculated log10_Peak_Width)
scatter_data_width_log <- data %>%
  group_by(Method, Histone) %>%
  summarize(
    median_log_peak_width = median(log10_Peak_Width, na.rm = TRUE),  # Use pre-calculated log10_Peak_Width
    sd_log_peak_width = sd(log10_Peak_Width, na.rm = TRUE),
    # Also keep raw widths for reference
    median_raw_peak_width = median(Mean_Peak_Width, na.rm = TRUE),
    sd_raw_peak_width = sd(Mean_Peak_Width, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(Method),
    Histone = factor(Histone)
  )

# Clean version for peak widths with custom axis limits
scatter_plot_width_clean <- ggplot(scatter_data_width_log, 
                                   aes(x = sd_log_peak_width,
                                       y = median_log_peak_width,
                                       color = Method,
                                       shape = Histone)) +
  geom_point(size = 7, alpha = 0.9, stroke = 1.5) +
  scale_color_manual(values = method_colors, 
                     name = "Peak-Calling Method") +
  scale_shape_manual(values = histone_shape_map,
                     name = "Histone Modification") +
  scale_x_continuous(
    labels = scales::comma_format(),
    limits = c(0, 1),  # Set X-axis limits (0 to 0.6)
    breaks = seq(0, 1, by = 0.1),  # Set breaks every 0.1
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_y_continuous(
    labels = scales::comma_format(),
    limits = c(1.5, 4.5),  # Set Y-axis limits (2.5 to 3.5)
    breaks = seq(1.5, 4.5, by = 0.5),  # Set breaks every 0.2
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title = "Standard Deviation vs Median of log10(Peak Width + 1)",
    x = "SD of log10(Peak Width + 1)",
    y = "Median of log10(Peak Width + 1)"
  ) +
  nature_methods_enhanced_theme() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.3, "cm")
  )

print(scatter_plot_width_clean)

ggsave("~/result_frip/Figure/HumanPBMC_SD_vs_Median_Scatter_With_Input.pdf", plot = scatter_plot_width_clean, width = 12, height = 9, device = cairo_pdf, bg = "white")

#########################################################################################################
## Figure 2 G
## Functional Analysis of scCUT&Tag (HumanPBMC without) 
## ChIPseeker assignment of peaks to regulatory elements 
#########################################################################################################
# ------ Package Setup ------
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(org.Hs.eg.db)
library(RColorBrewer)
library(ggpubr)

# ------ Data Preparation ------
# Define base directory and methods
base_dir <- "~/HumanPBMC_peakbed"
methods <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2",  "SEACR",  "SICER2")
histone_marks <- c("H3K27ac-b", "H3K27ac-s",  "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3") # Add other marks if needed

# Find all BED files for each method and histone mark
peak_files <- list()
for (method in methods) {
  method_files <- list.files(
    path = paste0(base_dir, method, "_peakbed"),
    pattern = paste0(method, "_.*\\.bed$"),
    full.names = TRUE
  )
  
  # Only keep files for specified histone marks
  method_files <- method_files[grepl(paste(histone_marks, collapse = "|"), method_files)]
  
  if (length(method_files) > 0) {
    peak_files[[method]] <- method_files
  }
}

# Load TxDb (Human hg38)
#
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ------ Peak Annotation ------
# Process all files with progress tracking
peak_anno <- list()
for (method in names(peak_files)) {
  message("Processing ", method, " (", length(peak_files[[method]]), " files)")
  
  peak_anno[[method]] <- lapply(peak_files[[method]], function(peak_file) {
    tryCatch({
      peak <- readPeakFile(peak_file)
      annotatePeak(peak, 
                  tssRegion = c(-3000, 3000),
                  TxDb = txdb,
                  annoDb = "org.Hs.eg.db")
    }, error = function(e) {
      message("Error processing ", basename(peak_file), ": ", e$message)
      return(NULL)
    })
  })
  
  # Name each annotation by cell type
  names(peak_anno[[method]]) <- gsub(paste0("^", method, "_|.bed$"), "", basename(peak_files[[method]]))
}

# ------ Visualization ------
# Load required packages (ensure they're loaded)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

# Create consolidated data frame for plotting
plot_data <- map_dfr(names(peak_anno), function(method) {
  map_dfr(names(peak_anno[[method]]), function(sample) {
    if (!is.null(peak_anno[[method]][[sample]])) {
      anno <- peak_anno[[method]][[sample]]@annoStat
      anno$Method <- method
      anno$Sample <- sample
      # Extract histone mark from sample name
      anno$HistoneMark <- str_extract(sample, paste(histone_marks, collapse = "|"))
      return(anno)
    }
  })
}) %>%
  mutate(Feature = factor(Feature, 
                         levels = c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)",
                                   "5' UTR", "3' UTR", "1st Exon", "Other Exon",
                                   "1st Intron", "Other Intron", 
                                   "Downstream (<1kb)", "Distal Intergenic")))

# Create a robust color palette with 11 distinct colors
feature_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62"
)
names(feature_colors) <- levels(plot_data$Feature)


#####################################################################################################

## Merge All figure together on a single page
#####################################################################################################

# ------ All Methods Parallel Visualization (A4 Page) ------
# Create a plot with all methods arranged horizontally
parallel_methods_plot <- ggplot(plot_data, aes(x = Frequency, y = Sample, fill = Feature)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.7) +
  scale_fill_manual(values = feature_colors) +
  labs(x = "Percentage (%)", 
       y = "",
       fill = "Genomic Feature") +
  facet_grid(. ~ Method, scales = "free_x", space = "free_x") +  # Horizontal faceting by method
  theme_pubr(base_size = 10) +
  theme(
    plot.title = element_blank(),  # Removed main title
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "top",  # Moved legend to bottom
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),  # Reduced margins
    strip.text.x = element_text(angle = 0, size = 9, face = "bold"),  # Method names as headers
    strip.background = element_rect(fill = "grey95"),
    legend.box.margin = margin(t = -10)  # Move legend closer to plot
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

# Calculate dimensions for A4 page (21.0 x 29.7 cm)
# Width will be full A4 width (21cm), height adjusted based on content
n_methods <- length(unique(plot_data$Method))
plot_width <- 21  # A4 width in cm
plot_height <- min(29.7, 5 + n_methods * 2)  # Cap at A4 height, adjust as needed
  

# ------ Save Output ------
output_dir <- "~/result_regulatory-elements_ChIPseeker/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save the parallel plot on A4 page
if (FALSE) {  # Change to TRUE if you want a title
  parallel_methods_plot <- parallel_methods_plot +
    labs(title = "Peak Feature Distribution Across All Methods") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
}

#  Histone mark version if multiple marks exist
if (length(unique(plot_data$HistoneMark)) > 1) {
  histone_facet_plot <- parallel_methods_plot +
    facet_grid(HistoneMark ~ Method, scales = "free", space = "free") +
    theme(strip.text.y = element_text(size = 8, face = "bold", angle = 0))
  
  ggsave(file.path(output_dir, "Figure_2G.pdf"), 
         histone_facet_plot, 
         width = plot_width, 
         height = min(29.7, plot_height * 1.5),  # Increased height for histone marks
         units = "cm", 
         dpi = 300, 
         limitsize = FALSE)
}

write_csv(plot_data, file.path(output_dir, "HumanPBMC_all_peak_annotation_stats_without_input.csv"))



###############################################################################
#Figure 2 H 
# COMBINED PROMOTER ANALYSIS FROM SAVED DATA
# Using MEDIAN with Standard Deviation (Recommended)
###############################################################################

# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)

# Read the saved CSV file
plot_data <- read_csv("~/result_regulatory-elements_ChIPseeker/HumanPBMC_all_peak_annotation_stats_without_input.csv")

# Define Nature Methods color palette for methods
all_methods <- unique(plot_data$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

###############################################################################
# DATA PROCESSING: Extract and combine promoter categories
###############################################################################

# Extract and combine all promoter categories
promoter_data <- plot_data %>%
  filter(Feature %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")) %>%
  mutate(Promoter_Type = Feature,
         Combined_Promoters = "All Promoters (≤3kb)")

# Calculate total promoter percentage for each sample-method combination
combined_promoter_summary <- promoter_data %>%
  group_by(Method, Sample, HistoneMark) %>%
  summarize(
    Total_Promoter_Percent = sum(Frequency, na.rm = TRUE),
    n_promoter_types = n_distinct(Promoter_Type),
    .groups = 'drop'
  )

###############################################################################
#  MEDIAN with STANDARD DEVIATION 
###############################################################################

# Calculate statistics using Median with SD
promoter_stats_median <- combined_promoter_summary %>%
  group_by(Method, HistoneMark) %>%
  summarize(
    median_promoter = median(Total_Promoter_Percent, na.rm = TRUE),
    mean_promoter = mean(Total_Promoter_Percent, na.rm = TRUE),
    sd_promoter = sd(Total_Promoter_Percent, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  )

# Create plot with Median and SD
median_promoter_plot <- ggplot(promoter_stats_median,
                               aes(x = HistoneMark,
                                   y = median_promoter,
                                   color = Method,
                                   group = Method)) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = pmax(median_promoter - sd_promoter, 0),
                    ymax = pmin(median_promoter + sd_promoter, 100)),
                width = 0.4,
                position = position_dodge(width = 0.7),
                alpha = 0.7,
                linewidth = 0.8) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title = "Median Promoter-Associated Peaks with Standard Deviation",
    subtitle = "Combined promoter regions within 3kb of TSS",
    x = "Histone Modification",
    y = "Median Percentage of Peaks in Promoters (%)",
    caption = "Error bars represent standard deviation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(median_promoter_plot)

# Save the plot
ggsave("~/result_regulatory-elements_ChIPseeker/HumanPBMC_Median_Promoters_with_SD_Without_Input.pdf",
       plot = median_promoter_plot,
       width = 16,
       height = 10,
       device = "pdf",
       bg = "white")



###############################################################################
# For DISTAL INTERGENIC ANALYSIS FROM SAVED DATA
# 
###############################################################################

# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)

# Read the saved CSV file
plot_data <- read_csv("~/result_regulatory-elements_ChIPseeker/HumanPBMC_all_peak_annotation_stats_without_input.csv")

# Define Nature Methods color palette for methods
all_methods <- unique(plot_data$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

###############################################################################
# DATA PROCESSING: Extract Distal Intergenic data
###############################################################################

# Extract Distal Intergenic data
intergenic_data <- plot_data %>%
  filter(Feature == "Distal Intergenic") %>%
  mutate(Feature_Type = "Distal Intergenic")

# Rename for consistency with analysis
distal_intergenic_summary <- intergenic_data %>%
  rename(Distal_Intergenic_Percent = Frequency) %>%
  select(Method, Sample, HistoneMark, Distal_Intergenic_Percent)

###############################################################################
#  MEDIAN with STANDARD DEVIATION 
###############################################################################

# Calculate statistics using Median with SD
intergenic_stats_median <- distal_intergenic_summary %>%
  group_by(Method, HistoneMark) %>%
  summarize(
    median_intergenic = median(Distal_Intergenic_Percent, na.rm = TRUE),
    mean_intergenic = mean(Distal_Intergenic_Percent, na.rm = TRUE),
    sd_intergenic = sd(Distal_Intergenic_Percent, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  )

# Create plot with Median and SD
median_intergenic_plot <- ggplot(intergenic_stats_median,
                                 aes(x = HistoneMark,
                                     y = median_intergenic,
                                     color = Method,
                                     group = Method)) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = pmax(median_intergenic - sd_intergenic, 0),
                    ymax = pmin(median_intergenic + sd_intergenic, 100)),
                width = 0.4,
                position = position_dodge(width = 0.7),
                alpha = 0.7,
                linewidth = 0.8) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title = "Median Distal Intergenic Peaks with Standard Deviation",
    subtitle = "Analysis of peaks in distal intergenic regions",
    x = "Histone Modification",
    y = "Median Percentage of Peaks in Distal Intergenic Regions (%)",
    caption = "Error bars represent standard deviation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(median_intergenic_plot)

# Save the plot
ggsave("~/result_regulatory-elements_ChIPseeker/HumanPBMC_Median_Distal_Intergenic_with_SD_With_Input.pdf",
       plot = median_intergenic_plot,
       width = 16,
       height = 10,
       device = "pdf",
       bg = "white")



