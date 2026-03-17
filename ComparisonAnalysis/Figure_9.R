###################################################################################
### Figure 9 A B
###################################################################################

# =============================================================================
# Add normalized log10 peak count (norm_Log10_peak) as Log10_peak_100 - Log10_peak
# Data: /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/all_combined_peak.tsv
# Output: all_combined_peak_with_norm.tsv
# =============================================================================

# Load required libraries
library(dplyr)
library(tidyr)

# ----------------------------------------------------------------------
# 1. Read the data (tab-separated)
# ----------------------------------------------------------------------
df <- read.table("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/all_combined_peak.tsv",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ----------------------------------------------------------------------
# 2. (Optional) Rename columns if needed – they are already fine.
# ----------------------------------------------------------------------
# The file has columns: Percentage, Caller, Peak, Histone, Type, Log10_peak

# ----------------------------------------------------------------------
# 3. Extract the 100% Log10_peak for each histone and caller
# ----------------------------------------------------------------------
log10_at_100 <- df %>%
  filter(Percentage == 100) %>%
  select(Histone, Caller, Log10_peak_100 = Log10_peak)

# ----------------------------------------------------------------------
# 4. Join with original data and compute normalized value (deficit)
# ----------------------------------------------------------------------
df_norm <- df %>%
  left_join(log10_at_100, by = c("Histone", "Caller")) %>%
  mutate(
    norm_Log10_peak = Log10_peak-Log10_peak_100  # positive deficit
  ) %>%
  select(-Log10_peak_100)   # remove helper column (optional)

# ----------------------------------------------------------------------
# 5. View the first few rows to check
# ----------------------------------------------------------------------
head(df_norm)

# ----------------------------------------------------------------------
# 6. Save the result to a TSV file
# ----------------------------------------------------------------------
output_file <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/all_combined_peak_with_norm_subtract.tsv"
write.table(df_norm, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("✅ Normalized data saved to:", output_file, "\n")



# =============================================================================
# Two‑panel figure for normalized direct peak counts (norm_Log10_peak)
# Data: /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/all_combined_peak_with_norm.tsv
# Output: HumanPBMC_peakcount_norm_two_panel.pdf
# =============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(patchwork)   # for combining plots

# ----------------------------------------------------------------------
# 1. Read the normalized data (tab‑separated)
# ----------------------------------------------------------------------
df <- read.table("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/all_combined_peak_with_norm_subtract.tsv",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ----------------------------------------------------------------------
# 2. Keep only Direct peaks and prepare factors
# ----------------------------------------------------------------------
df <- df %>%
  filter(Type == "Direct") %>%
  mutate(
    # Create cell percentage factor with correct order
    celltype_pct = factor(paste0(Percentage, "%"),
                          levels = paste0(sort(unique(Percentage)), "%")),
    # Ensure Histone is a factor preserving original order
    Histone = factor(Histone, levels = unique(Histone)),
    # Ensure Caller is a factor
    Caller = factor(Caller)
  )

# ----------------------------------------------------------------------
# 3. Compute summary statistics for the line plot (median ± SE)
# ----------------------------------------------------------------------
summary_df <- df %>%
  group_by(Caller, Histone, celltype_pct) %>%
  summarise(
    median_norm = median(norm_Log10_peak, na.rm = TRUE),
    sd_norm     = sd(norm_Log10_peak, na.rm = TRUE),
    n           = n(),
    se_norm     = sd_norm / sqrt(n),
    .groups     = "drop"
  )

# ----------------------------------------------------------------------
# 4. Define tool colors (Nature‑Methods palette, extended as needed)
# ----------------------------------------------------------------------
all_callers <- levels(df$Caller)
nature_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
caller_colors <- setNames(nature_colors[1:length(all_callers)], all_callers)

# ----------------------------------------------------------------------
# 5. Enhanced Nature Methods theme (from your earlier code)
# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
# 6. Create the line plot (median ± SE of norm_Log10_peak)
# ----------------------------------------------------------------------
p1 <- ggplot(summary_df, aes(x = celltype_pct, y = median_norm,
                             group = Caller, color = Caller)) +
  geom_line(linewidth = 1.4, alpha = 0.9, lineend = "round", linejoin = "round") +
  geom_point(aes(fill = Caller), shape = 21, size = 4, color = "white", stroke = 1.2,
             alpha = 0.95) +
  geom_errorbar(aes(ymin = median_norm - se_norm, ymax = median_norm + se_norm,
                    color = Caller),
                width = 0.3, linewidth = 1.0, alpha = 0.8) +
  facet_wrap(~ Histone, nrow = 1, scales = "free_x") +
  scale_color_manual(values = caller_colors, name = "Tool") +
  scale_fill_manual(values = caller_colors, name = "Tool") +
  labs(
    x = "Percentage of cells",
    y = "Normalized log10(peak count)"
  ) +
  nature_methods_enhanced_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                               margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold"),
    legend.position = "bottom"
  )

# ----------------------------------------------------------------------
# 7. Create the boxplot (distribution of norm_Log10_peak)
# ----------------------------------------------------------------------
p2 <- ggplot(df, aes(x = Caller, y = norm_Log10_peak, fill = Caller)) +
  geom_boxplot(outlier.size = 1.5, outlier.shape = 21, outlier.fill = "white",
               outlier.color = "black", lwd = 0.4, alpha = 0.8) +
  facet_wrap(~ Histone, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = caller_colors) +
  labs(
    x = "Tool",
    y = "Normalized log10(peak count)"
  ) +
  nature_methods_enhanced_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                               margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold")
  )

# ----------------------------------------------------------------------
# 8. Combine plots with patchwork and add overall title & tags
# ----------------------------------------------------------------------
combined_plot <- (p1 / p2) +
  plot_annotation(
    title = "hPBMC: Normalized peak counts of saturation across the histone marks",

    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, 
                                margin = margin(b = 5), color = "#2C3E50"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10),
                                   color = "#566573")
    )
  ) &
  # Apply consistent smaller font sizes for the combined figure
  theme(
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ----------------------------------------------------------------------
# 9. Display and save the figure
# ----------------------------------------------------------------------
print(combined_plot)

ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/HumanPBMC_peakcount_norm_two_panel_subtract_real.pdf",
       plot = combined_plot,
width = 16,
       height = 12,    # inches (adjust height as needed)
       device = cairo_pdf,
       bg = "white")



# =============================================================================
# Two‑panel figure for normalized direct peak counts (norm_Log10_peak)
# Data: /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/all_combined_peak_with_norm.tsv
# Output: HumanPBMC_peakcount_norm_two_panel.pdf
# =============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(patchwork)   # for combining plots

# ----------------------------------------------------------------------
# 1. Read the normalized data (tab‑separated)
# ----------------------------------------------------------------------
df <- read.table("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/all_combined_peak_with_norm_subtract.tsv",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ----------------------------------------------------------------------
# 2. Keep only Direct peaks and prepare factors
# ----------------------------------------------------------------------
df <- df %>%
  filter(Type == "Direct") %>%
  mutate(
    # Create cell percentage factor with correct order
    celltype_pct = factor(paste0(Percentage, "%"),
                          levels = paste0(sort(unique(Percentage)), "%")),
    # Ensure Histone is a factor preserving original order
    Histone = factor(Histone, levels = unique(Histone)),
    # Ensure Caller is a factor
    Caller = factor(Caller)
  )


# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# 3. Summarise by tool and cell percentage (across histones)
# ----------------------------------------------------------------------
df_cell_summary <- df %>%
  group_by(Caller, celltype_pct) %>%
  summarise(
    median_log = median(norm_Log10_peak, na.rm = TRUE),
    sd_log     = sd(norm_Log10_peak, na.rm = TRUE),   # variability across histones
    n          = n(),
    .groups    = "drop"
  ) %>%
  mutate(
    tool = factor(Caller)
  )

# ----------------------------------------------------------------------
# 4. Nature‑Methods theme (same as before)
# ----------------------------------------------------------------------
nature_methods_enhanced_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "Helvetica", color = "black", size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, 
                                margin = margin(b = 12), color = "#2C3E50"),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "#566573",
                                   margin = margin(b = 18)),
      axis.title = element_text(face = "bold", size = 14, color = "#2C3E50"),
      axis.title.x = element_text(margin = margin(t = 12), vjust = -0.5),
      axis.title.y = element_text(margin = margin(r = 12), vjust = 2),
      axis.text = element_text(size = 11, color = "#34495E"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                 margin = margin(t = 5), face = "bold"),
      axis.text.y = element_text(margin = margin(r = 5), face = "bold"),
      axis.line = element_line(color = "#2C3E50", linewidth = 0.8),
      axis.ticks = element_line(color = "#2C3E50", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(face = "bold", size = 13, color = "#2C3E50",
                                  margin = margin(b = 8)),
      legend.text = element_text(size = 12, color = "#34495E",
                                 margin = margin(r = 15)),
      legend.position = "bottom",
      legend.key.size = unit(0.6, "cm"),
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
                                  margin = margin(t = 15), hjust = 1)
    )
}

# ----------------------------------------------------------------------
# 5. Define tool colors (Nature‑Methods palette)
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# 5. Define tool colors (Nature‑Methods palette)
# ----------------------------------------------------------------------
all_tools <- unique(df$Caller)   # corrected: use Caller, not tool
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
tool_colors <- setNames(nature_methods_colors[1:length(all_tools)], all_tools)
# ----------------------------------------------------------------------
# 6. Set y‑axis limits automatically (include all error bars)
# ----------------------------------------------------------------------
y_min <- min(df_cell_summary$median_log - df_cell_summary$sd_log, na.rm = TRUE)
y_max <- max(df_cell_summary$median_log + df_cell_summary$sd_log, na.rm = TRUE)
y_pad <- 0.05 * (y_max - y_min)
y_limits <- c(y_min - y_pad, y_max + y_pad)

# ----------------------------------------------------------------------
# 7. Create the plot: median ± SD across histones, per cell percentage
# ----------------------------------------------------------------------
p_cell <- ggplot(df_cell_summary, aes(x = celltype_pct, y = median_log, 
                                       group = tool, color = tool)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = median_log - sd_log, 
                    ymax = median_log + sd_log),
                width = 0.2, size = 0.8, alpha = 0.5) +
  scale_color_manual(values = tool_colors, name = "Peak caller") +
  scale_y_continuous(limits = y_limits) +
  labs(
    title = "hPBMC: Saturation Summary Percentage by tools",
    x = "Cell percentage",
    y = "Normalized log10(peak count)"
  ) +
  nature_methods_enhanced_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# ----------------------------------------------------------------------
# 8. Save the plot
# ----------------------------------------------------------------------
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional_final/HumanPBMC_saturation_summary_subtract_real.pdf",
       plot = p_cell,
       width = 7, height = 7,      # adjust as needed
       device = cairo_pdf,
       bg = "white")

# ----------------------------------------------------------------------
# 9. Display the plot (optional)
# ----------------------------------------------------------------------
print(p_cell)


###################################################################################
### Figure 9 C
###################################################################################

##############################################################################################################
###
### COMPACT PUBLICATION VERSION for Nature Methods
### Percentage-wise Peak Saturation Analysis - All Histone Modifications Combined
### Author: Md Wahiduzzaman
### Date: 2025-10-17
### VERSION: Compact Publication with Full Page Figure 4
### MODIFICATION: Y-axis for direct peak counts uses log(direct peak count + 1)
### MODIFICATION 2: All panels in one row for Fig1
##############################################################################################################

#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(scales)
})

cat("📊 PEAK CALLING SATURATION ANALYSIS - COMPACT PUBLICATION VERSION\n")
cat("=================================================================\n")
cat("📈 MODIFICATION: Y-axis uses log(direct peak count + 1)\n")
cat("📈 MODIFICATION: All panels in one row for Fig1\n\n")

# Configuration
output_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_saturation_analysis_additional2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", output_dir, "\n")
}

# Define constants
histones <- c("H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3","H3K9me3")
peak_callers <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
percentages <- seq(10, 100, 10)

# Nature Methods color palette (Tableau 10)
nature_colors <- c(
  "DROMPAplus" = "#0072B2",    # Blue
  "Genrich" = "#D55E00",       # Vermillion
  "GoPeaks" = "#009E73",       # Bluish green
  "HOMER" = "#CC79A7",         # Reddish purple
  "MACS2" = "#F0E442",         # Yellow
  "SEACR" = "#56B4E9",         # Sky blue
  "SICER2" = "#E69F00"         # Orange
)

# Compact Nature Methods theme for regular figures - UPDATED FOR SINGLE ROW
theme_nature_compact <- function(base_size = 8, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black", family = base_family),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5,
                               margin = margin(b = 5)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 2)),
      axis.text.y = element_text(margin = margin(r = 2)),
      axis.line = element_line(color = "black", linewidth = 0.25),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = unit(1.5, "pt"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.margin = margin(t = -5, b = 0),
      legend.spacing.x = unit(0.2, "cm"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.15),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(5, 5, 5, 5),
      strip.text = element_text(size = base_size, face = "bold", 
                               margin = margin(t = 3, b = 3)),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}

# Enhanced theme for full-page Figure 4 - UPDATED FOR SINGLE ROW
theme_nature_fullpage <- function(base_size = 9, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black", family = base_family),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5,
                               margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40",
                                  margin = margin(b = 15)),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 3)),
      axis.text.y = element_text(margin = margin(r = 3)),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.position = "bottom",
      legend.key.size = unit(0.5, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.margin = margin(t = 5, b = 5),
      legend.spacing.x = unit(0.3, "cm"),
      legend.box.margin = margin(t = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(15, 15, 15, 15),
      strip.text = element_text(size = base_size + 1, face = "bold", 
                               margin = margin(t = 5, b = 5)),
      strip.background = element_rect(fill = "grey93", color = "grey70", linewidth = 0.3)
    )
}

# Function to extract peaks from BED files
extract_peaks <- function(bed_file) {
  if (!file.exists(bed_file)) {
    return(data.frame(chr = character(), start = numeric(), end = numeric()))
  }
  
  data <- tryCatch({
    read.table(bed_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", 
               comment.char = "#", fill = TRUE, quote = "")
  }, error = function(e) {
    return(data.frame(V1 = character()))
  })
  
  if (nrow(data) == 0 || ncol(data) < 3) {
    return(data.frame(chr = character(), start = numeric(), end = numeric()))
  }
  
  coords <- data.frame(
    chr = as.character(data[, 1]),
    start = suppressWarnings(as.numeric(data[, 2])),
    end = suppressWarnings(as.numeric(data[, 3]))
  )
  
  # Filter valid entries
  coords <- coords[grepl("^chr", coords$chr), ]
  coords <- coords[!is.na(coords$start) & !is.na(coords$end), ]
  coords <- coords[coords$start < coords$end, ]
  
  if (nrow(coords) > 0) {
    return(unique(coords))
  } else {
    return(data.frame(chr = character(), start = numeric(), end = numeric()))
  }
}

# Count peaks in BED file
count_peaks <- function(bed_file) {
  peaks <- extract_peaks(bed_file)
  return(nrow(peaks))
}

# Merge peak coordinates
merge_peak_coords <- function(coord_list) {
  if (length(coord_list) == 0) return(0)
  
  all_coords <- bind_rows(coord_list)
  if (nrow(all_coords) == 0) return(0)
  
  unique_peaks <- distinct(all_coords, chr, start, end)
  return(nrow(unique_peaks))
}

# Main analysis function - MODIFIED TO RETURN ACTUAL COUNTS
analyze_histone_saturation <- function(hist) {
  cat("Processing:", hist, "\n")
  
  # Define peak directory
  peakdir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed_percentagewise_corrected/Allpeakbed/", hist, "_peakbed")
  
  if (!dir.exists(peakdir)) {
    cat("  Directory not found:", peakdir, "\n")
    return(NULL)
  }
  
  # Initialize results data frames
  direct_results <- data.frame(Percentage = percentages)
  cumulative_results <- data.frame(Percentage = percentages)
  
  # Process each peak caller
  for (caller in peak_callers) {
    cat("  ", caller, "... ", sep = "")
    
    direct_peaks <- numeric(length(percentages))
    cumulative_peaks <- numeric(length(percentages))
    all_coords_list <- list()
    
    for (i in 1:length(percentages)) {
      pct <- percentages[i]
      bed_file <- paste0(peakdir, "/", caller, "_", hist, "_", pct, "cell.bed")
      
      # Count direct peaks
      direct_peaks[i] <- count_peaks(bed_file)
      
      # Add to cumulative list
      current_coords <- extract_peaks(bed_file)
      if (nrow(current_coords) > 0) {
        all_coords_list[[length(all_coords_list) + 1]] <- current_coords
      }
      
      # Count cumulative unique peaks
      cumulative_peaks[i] <- merge_peak_coords(all_coords_list)
    }
    
    direct_results[[caller]] <- direct_peaks
    cumulative_results[[caller]] <- cumulative_peaks
    
    cat("Done\n")
  }
  
  # Calculate percentages for cumulative (but keep direct as actual counts)
  cumulative_percentage <- cumulative_results
  
  for (caller in peak_callers) {
    max_cumulative <- max(cumulative_results[[caller]], na.rm = TRUE)
    
    if (!is.na(max_cumulative) && max_cumulative > 0) {
      cumulative_percentage[[caller]] <- (cumulative_results[[caller]] / max_cumulative) * 100
    } else {
      cumulative_percentage[[caller]] <- 0
    }
  }
  
  # Save results
  tryCatch({
    write.table(direct_results, paste0(peakdir, "/saturation_direct_results.tsv"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(cumulative_results, paste0(peakdir, "/saturation_cumulative_results.tsv"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
  }, error = function(e) {
    cat("  Error saving results:", e$message, "\n")
  })
  
  # Prepare data for plotting - DIRECT as log(count+1), CUMULATIVE as percentage
  direct_long <- direct_results %>%
    pivot_longer(cols = all_of(peak_callers), 
                names_to = "Caller", values_to = "Value") %>%
    mutate(
      Histone = hist,
      Type = "Direct",
      # Apply log transformation for plotting
      LogValue = log10(Value + 1)
    )
  
  cumulative_long <- cumulative_percentage %>%
    pivot_longer(cols = all_of(peak_callers), 
                names_to = "Caller", values_to = "Value") %>%
    mutate(Histone = hist, Type = "Cumulative")
  
  return(list(
    direct = direct_long,
    cumulative = cumulative_long,
    direct_counts = direct_results,
    cumulative_counts = cumulative_results
  ))
}

# Create compact plot for regular figures - MODIFIED FOR SINGLE ROW
create_compact_plot <- function(data_combined, plot_type = "Direct") {
  
  if (plot_type == "Direct") {
    y_lab <- expression(bold("Direct Peak Count [log"[10]*"(count + 1)]"))
    title_text <- "Direct Peak Counts (Log Scale)"
  } else {
    y_lab <- "Unique Peaks (% of total)"
    title_text <- "Cumulative Peak Discovery"
  }
  
  # Filter data
  if (plot_type == "Direct") {
    # Use LogValue for direct plots
    plot_data <- data_combined %>% 
      filter(Type == plot_type) %>%
      mutate(Value = LogValue)  # Use the pre-calculated log values
  } else {
    plot_data <- data_combined %>% filter(Type == plot_type)
  }
  
  # Determine number of columns for facets
  # For Direct plot: 7 histones in one row (ncol = 7)
  # For Cumulative plot: Keep original layout (ncol = 4)
  if (plot_type == "Direct") {
    ncol_facets <- 7  # All 7 histones in one row
  } else {
    ncol_facets <- 4  # Original layout for cumulative
  }
  
  ggplot(plot_data, aes(x = Percentage, y = Value, color = Caller, group = Caller)) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_point(size = 1, shape = 16, alpha = 0.8) +
    facet_wrap(~ Histone, ncol = ncol_facets) +  # MODIFIED: ncol based on plot type
    scale_color_manual(values = nature_colors) +
    scale_x_continuous(
      breaks = c(10, 30, 50, 70, 90, 100),
      labels = c("10%", "30%", "50%", "70%", "90%", "100%"),
      limits = c(5, 105),
      expand = expansion(mult = 0.01)
    ) +
    labs(
      x = "Percentage of Cells",
      y = y_lab,
      title = title_text
    ) +
    theme_nature_compact() +
    theme(
      strip.text = element_text(size = 7, face = "bold", margin = margin(t = 3, b = 3))
    ) +
    guides(color = guide_legend(nrow = 1, 
                               override.aes = list(linewidth = 1, alpha = 1, size = 1.5)))
}

# Create full-page Figure 4 with all histones - MODIFIED FOR SINGLE ROW
create_fullpage_figure4 <- function(data_combined) {
  
  # Prepare data with appropriate transformations
  plot_data <- data_combined %>%
    mutate(
      # Use log values for direct, percentage for cumulative
      PlotValue = ifelse(Type == "Direct", LogValue, Value),
      Type = factor(Type, levels = c("Direct", "Cumulative"),
                   labels = c("A. Direct Peak Counts", "B. Cumulative Unique Peaks")),
      Histone = factor(Histone, levels = histones)
    )
  
  ggplot(plot_data, aes(x = Percentage, y = PlotValue, color = Caller, group = Caller)) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    geom_point(size = 1.2, shape = 16, alpha = 0.85) +
    facet_grid(Type ~ Histone, scales = "free_y") +  # All 7 histones in one row
    scale_color_manual(
      name = "Peak Caller",
      values = nature_colors,
      breaks = peak_callers,
      labels = peak_callers
    ) +
    scale_x_continuous(
      breaks = c(10, 30, 50, 70, 90, 100),
      labels = c("10%", "30%", "50%", "70%", "90%", "100%"),
      limits = c(5, 105),
      expand = expansion(mult = 0.01)
    ) +
    labs(
      x = "Percentage of Cells",
      y = NULL,
      title = "Comprehensive Peak Calling Saturation Analysis Across Histone Modifications",
      subtitle = "Comparison of seven peak callers across increasing cell percentages"
    ) +
    theme_nature_fullpage() +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.spacing.x = unit(1.0, "lines"),  # Reduced spacing for single row
      panel.spacing.y = unit(1.5, "lines"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(margin = margin(r = 10)),
      plot.title = element_text(margin = margin(b = 10)),
      plot.subtitle = element_text(margin = margin(b = 20))
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(linewidth = 1.2, size = 2, alpha = 1)
      )
    )
}

# Create an alternative full-page figure with separate y-axis labels
create_fullpage_figure4_alt <- function(data_combined) {
  
  # Create separate y-axis labels for each row
  plot_data <- data_combined %>%
    mutate(
      PlotValue = ifelse(Type == "Direct", LogValue, Value),
      Type = factor(Type, levels = c("Direct", "Cumulative"),
                   labels = c("A. Direct Peak Counts [log10(count + 1)]", 
                             "B. Cumulative Unique Peaks (% of total)")),
      Histone = factor(Histone, levels = histones)
    )
  
  ggplot(plot_data, aes(x = Percentage, y = PlotValue, color = Caller, group = Caller)) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    geom_point(size = 1.2, shape = 16, alpha = 0.85) +
    facet_grid(Type ~ Histone, scales = "free_y") +  # All 7 histones in one row
    scale_color_manual(
      name = "Peak Caller",
      values = nature_colors,
      breaks = peak_callers,
      labels = peak_callers
    ) +
    scale_x_continuous(
      breaks = c(10, 30, 50, 70, 90, 100),
      labels = c("10%", "30%", "50%", "70%", "90%", "100%"),
      limits = c(5, 105),
      expand = expansion(mult = 0.01)
    ) +
    labs(
      x = "Percentage of Cells",
      title = "Comprehensive Peak Calling Saturation Analysis Across Histone Modifications",
      subtitle = "Comparison of seven peak callers across increasing cell percentages"
    ) +
    theme_nature_fullpage() +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.spacing.x = unit(1.0, "lines"),  # Reduced spacing for single row
      panel.spacing.y = unit(1.5, "lines"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(margin = margin(r = 10)),
      plot.title = element_text(margin = margin(b = 10)),
      plot.subtitle = element_text(margin = margin(b = 20))
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(linewidth = 1.2, size = 2, alpha = 1)
      )
    )
}

# Main analysis
cat("\n🚀 STARTING ANALYSIS...\n")
cat("=======================\n")

# Process all histones
all_results <- list()
for (hist in histones) {
  hist_results <- analyze_histone_saturation(hist)
  if (!is.null(hist_results)) {
    all_results[[hist]] = hist_results
  }
}

# Combine all data
cat("\n📊 COMBINING DATA...\n")

direct_combined <- bind_rows(lapply(all_results, function(x) x$direct))
cumulative_combined <- bind_rows(lapply(all_results, function(x) x$cumulative))

all_combined <- bind_rows(direct_combined, cumulative_combined)

# Create plots
cat("\n🎨 CREATING PLOTS...\n")

# Individual plots (Fig 1-3)
# For Direct plot: All 7 panels in one row
p_direct <- create_compact_plot(all_combined, "Direct") +
  scale_y_continuous(
    name = expression(bold("Direct Peak Count [log"[10]*"(count + 1)]"))
  )

# For Cumulative plot: Keep original layout (4 columns)
p_cumulative <- create_compact_plot(all_combined, "Cumulative") +
  scale_y_continuous(
    name = "Unique Peaks (% of total)",
    labels = function(x) paste0(round(x), "%"),
    limits = c(0, 105),
    expand = expansion(mult = c(0, 0.02))
  )

# Save individual plots
cat("  Saving individual plots (Fig 1-3)...\n")

# For Direct plot with one row: Wider figure
ggsave(
  file.path(output_dir, "Fig1_direct_peaks_log_single_row.tiff"),
  p_direct,
  width = 11, height = 5, units = "in",  # Wider for 7 panels in one row
  dpi = 600, compression = "lzw"
)

ggsave(
  file.path(output_dir, "Fig1_direct_peaks_log_single_row.pdf"),
  p_direct,
  width = 11, height = 5, units = "in",  # Wider for 7 panels in one row
  device = cairo_pdf
)

# For comparison, also save the original layout
p_direct_original <- p_direct + facet_wrap(~ Histone, ncol = 4)
ggsave(
  file.path(output_dir, "Fig1_direct_peaks_log_original.tiff"),
  p_direct_original,
  width = 7.2, height = 4.5, units = "in",
  dpi = 600, compression = "lzw"
)

ggsave(
  file.path(output_dir, "Fig1_direct_peaks_log_original.pdf"),
  p_direct_original,
  width = 7.2, height = 4.5, units = "in",
  device = cairo_pdf
)

# Cumulative plot (original layout)
ggsave(
  file.path(output_dir, "Fig2_cumulative_peaks.tiff"),
  p_cumulative,
  width = 7.2, height = 4.5, units = "in",
  dpi = 600, compression = "lzw"
)

ggsave(
  file.path(output_dir, "Fig2_cumulative_peaks.pdf"),
  p_cumulative,
  width = 7.2, height = 4.5, units = "in",
  device = cairo_pdf
)

# Create combined figure (Fig 3) with single row for direct plots
cat("  Creating combined figure (Fig 3)...\n")

legend <- get_legend(p_direct)

# Create direct plot with single row (no legend)
p_direct_noleg <- p_direct + 
  theme(legend.position = "none") +
  theme(plot.margin = margin(5, 5, 5, 5))

# Create cumulative plot (no legend, original layout)
p_cumulative_noleg <- p_cumulative + 
  theme(legend.position = "none") +
  theme(plot.margin = margin(5, 5, 5, 5))

# Create combined plot with single row for direct
combined_plot_single_row <- grid.arrange(
  arrangeGrob(
    p_direct_noleg,
    p_cumulative_noleg,
    ncol = 1,  # Stack vertically
    heights = c(1, 1.2)  # Slightly more space for cumulative (4 rows vs 1 row)
  ),
  legend,
  nrow = 2,
  heights = c(10, 1),
  top = textGrob("Peak Calling Saturation Analysis", 
                gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica"))
)

ggsave(
  file.path(output_dir, "Fig3_combined_saturation_log_single_row.tiff"),
  combined_plot_single_row,
  width = 10, height = 8, units = "in",  # Adjusted for new layout
  dpi = 600, compression = "lzw"
)

ggsave(
  file.path(output_dir, "Fig3_combined_saturation_log_single_row.pdf"),
  combined_plot_single_row,
  width = 10, height = 8, units = "in",  # Adjusted for new layout
  device = cairo_pdf
)

# Also create original combined plot
p_direct_noleg_original <- p_direct_original + theme(legend.position = "none")
combined_plot_original <- grid.arrange(
  arrangeGrob(
    p_direct_noleg_original,
    p_cumulative_noleg,
    ncol = 2,
    widths = c(1, 1)
  ),
  legend,
  nrow = 2,
  heights = c(10, 1),
  top = textGrob("Peak Calling Saturation Analysis", 
                gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica"))
)

ggsave(
  file.path(output_dir, "Fig3_combined_saturation_log_original.tiff"),
  combined_plot_original,
  width = 9, height = 5, units = "in",
  dpi = 600, compression = "lzw"
)

ggsave(
  file.path(output_dir, "Fig3_combined_saturation_log_original.pdf"),
  combined_plot_original,
  width = 9, height = 5, units = "in",
  device = cairo_pdf
)

# Create FULL-PAGE Figure 4 with log scale for direct counts
cat("  Creating FULL-PAGE Figure 4 with log scale (single row)...\n")

# Create a custom full-page figure with proper y-axis labeling
fig4_fullpage <- create_fullpage_figure4_alt(all_combined)

# Save as full-page figure (Nature Methods full page is typically ~7.2" wide × 9.8" height)
ggsave(
  file.path(output_dir, "Fig4_fullpage_comprehensive_log_single_row.tiff"),
  fig4_fullpage,
  width = 11, height = 6, units = "in",  # Wider for single row
  dpi = 600, compression = "lzw"
)

ggsave(
  file.path(output_dir, "Fig4_fullpage_comprehensive_log_single_row.pdf"),
  fig4_fullpage,
  width = 11, height = 6, units = "in",  # Wider for single row
  device = cairo_pdf
)

# Also create a landscape version for supplementary
fig4_landscape <- fig4_fullpage + 
  theme(
    strip.text.y = element_text(angle = 90, hjust = 0.5),
    panel.spacing.x = unit(0.8, "lines"),  # Even tighter for landscape
    panel.spacing.y = unit(1.2, "lines")
  )

ggsave(
  file.path(output_dir, "FigS1_landscape_comprehensive_log_single_row.tiff"),
  fig4_landscape,
  width = 14, height = 8.5, units = "in",  # Wider landscape
  dpi = 600, compression = "lzw"
)

ggsave(
  file.path(output_dir, "FigS1_landscape_comprehensive_log_single_row.pdf"),
  fig4_landscape,
  width = 14, height = 8.5, units = "in",  # Wider landscape
  device = cairo_pdf
)

# Generate summary statistics with log values
cat("\n📈 GENERATING SUMMARY STATISTICS...\n")

# First, compile actual counts for summary
direct_counts_list <- lapply(all_results, function(x) x$direct_counts)
cumulative_counts_list <- lapply(all_results, function(x) x$cumulative_counts)

# Combine all direct counts into one dataframe
all_direct_counts <- data.frame(Histone = character(),
                               Caller = character(),
                               Percentage = numeric(),
                               Direct_Count = numeric(),
                               Log_Direct_Count = numeric())

for (hist in histones) {
  if (!is.null(direct_counts_list[[hist]])) {
    hist_counts <- direct_counts_list[[hist]] %>%
      pivot_longer(cols = all_of(peak_callers),
                  names_to = "Caller",
                  values_to = "Direct_Count") %>%
      mutate(Histone = hist,
             Log_Direct_Count = log10(Direct_Count + 1))
    all_direct_counts <- bind_rows(all_direct_counts, hist_counts)
  }
}

summary_stats <- all_direct_counts %>%
  group_by(Histone, Caller) %>%
  summarize(
    Max_Direct_Count = max(Direct_Count, na.rm = TRUE),
    Max_Log_Direct_Count = max(Log_Direct_Count, na.rm = TRUE),
    Direct_Count_at_50pct = Direct_Count[Percentage == 50],
    Log_Direct_Count_at_50pct = Log_Direct_Count[Percentage == 50],
    .groups = "drop"
  ) %>%
  arrange(Histone, desc(Max_Direct_Count))

# Save summary statistics
write.csv(summary_stats, file.path(output_dir, "saturation_summary_statistics_log.csv"), 
          row.names = FALSE)

# Create a concise summary table
top_performers <- summary_stats %>%
  group_by(Histone) %>%
  slice_max(order_by = Max_Direct_Count, n = 2) %>%
  arrange(Histone, desc(Max_Direct_Count))

write.csv(top_performers, file.path(output_dir, "top_performers_summary_log.csv"), 
          row.names = FALSE)

# Create performance heatmap for supplementary
performance_matrix <- summary_stats %>%
  select(Histone, Caller, Max_Log_Direct_Count) %>%
  pivot_wider(names_from = Caller, values_from = Max_Log_Direct_Count) %>%
  column_to_rownames("Histone")

write.csv(performance_matrix, file.path(output_dir, "performance_matrix_log.csv"))

# Print final summary
cat("\n✅ ANALYSIS COMPLETE!\n")
cat("====================\n\n")

cat("Output files saved to:", output_dir, "\n\n")

cat("📁 MAIN FIGURES (SINGLE ROW LAYOUT):\n")
cat("1. Fig1_direct_peaks_log_single_row.tiff/pdf - Direct peaks in ONE ROW (10 x 4 in)\n")
cat("2. Fig1_direct_peaks_log_original.tiff/pdf   - Direct peaks original (7.2 x 4.5 in)\n")
cat("3. Fig2_cumulative_peaks.tiff/pdf           - Cumulative (7.2 x 4.5 in)\n")
cat("4. Fig3_combined_saturation_log_single_row.tiff/pdf - Combined with single row (10 x 8 in)\n")
cat("5. Fig3_combined_saturation_log_original.tiff/pdf   - Combined original (9 x 5 in)\n")
cat("6. Fig4_fullpage_comprehensive_log_single_row.tiff/pdf - FULL PAGE single row (11 x 6 in) ★\n\n")

cat("📁 SUPPLEMENTARY FIGURES:\n")
cat("1. FigS1_landscape_comprehensive_log_single_row.tiff/pdf - Landscape single row (14 x 8.5 in)\n\n")

cat("📊 DATA FILES:\n")
cat("1. saturation_summary_statistics_log.csv - Full statistics with log values\n")
cat("2. top_performers_summary_log.csv        - Top 2 performers per histone\n")
cat("3. performance_matrix_log.csv            - Log performance matrix\n\n")

cat("📋 SUMMARY OF TOP PERFORMERS (BY PEAK COUNT):\n")
print(top_performers %>% select(Histone, Caller, Max_Direct_Count, Max_Log_Direct_Count))

cat("\n🎯 KEY METRICS:\n")
cat("- Fig1: All 7 histone panels in ONE ROW (ncol = 7 in facet_wrap)\n")
cat("- All direct peak counts transformed using log10(count + 1)\n")
cat("- Adjusted figure dimensions for better aspect ratio\n")
cat("- Enhanced typography and spacing for readability\n")
cat("- High-resolution TIFF (600 DPI) and vector PDF output\n")

# Visual representation of figure layout
cat("\n📐 FIGURE 1 LAYOUT (SINGLE ROW):\n")
cat("┌─────────────────────────────────────────────────────────────────────────┐\n")
cat("│                     Direct Peak Counts (Log Scale)                     │\n")
cat("├───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┤\n")
cat("│H3K│H3K│H3K│H3K│H3K│H3K│H3K│   │H3K│H3K│H3K│H3K│H3K│H3K│H3K│   │H3K│H3K│\n")
cat("│27a│27a│27m│4me│4me│4me│9me│   │27a│27a│27m│4me│4me│4me│9me│   │27a│27a│\n")
cat("│c-b│c-s│e3 │1  │2  │3  │3  │   │c-b│c-s│e3 │1  │2  │3  │3  │   │c-b│c-s│\n")
cat("├───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┴───┤\n")
cat("│                                                                        │\n")
cat("│                     [7 peak callers - colored lines]                   │\n")
cat("└─────────────────────────────────────────────────────────────────────────┘\n")

cat("\n📐 FIGURE 4 LAYOUT (SINGLE ROW):\n")
cat("┌─────────────────────────────────────────────────────────────────────────┐\n")
cat("│                    Comprehensive Peak Calling Saturation Analysis      │\n")
cat("│                    Comparison of seven peak callers...                │\n")
cat("├─────────────────────────────────────────────────────────────────────────┤\n")
cat("│ A. Direct Peak Counts [log10(count + 1)]                              │\n")
cat("│ ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┐                          │\n")
cat("│ │H3K27│H3K27│H3K27│H3K4m│H3K4m│H3K4m│H3K9m│                          │\n")
cat("│ │ac-b │ac-s │me3  │e1   │e2   │e3   │e3   │                          │\n")
cat("│ └─────┴─────┴─────┴─────┴─────┴─────┴─────┘                          │\n")
cat("│                                                                        │\n")
cat("│ B. Cumulative Unique Peaks (% of total)                               │\n")
cat("│ ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┐                          │\n")
cat("│ │H3K27│H3K27│H3K27│H3K4m│H3K4m│H3K4m│H3K9m│                          │\n")
cat("│ │ac-b │ac-s │me3  │e1   │e2   │e3   │e3   │                          │\n")
cat("│ └─────┴─────┴─────┴─────┴─────┴─────┴─────┘                          │\n")
cat("├─────────────────────────────────────────────────────────────────────────┤\n")
cat("│                     LEGEND (7 peak callers)                           │\n")
cat("└─────────────────────────────────────────────────────────────────────────┘\n")

# Save session info
sink(file.path(output_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n📝 Session info saved to: session_info.txt\n")

###################################################################################
### Figure 9 D
###################################################################################

# =============================================================================
# COMPLETE SCRIPT: Running time of peak callers (log10 scale)
# Combined figure: line plot (across percentages) + box plot (per tool)
# =============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(patchwork)      # for combining plots

# ----------------------------------------------------------------------
# 1. Read and clean data
# ----------------------------------------------------------------------
df <- read.csv("/home/wahid/project_scHMTF/GSE195725_processed_data/saturation_with_scCTPEAK/HumanPBMC_running_time.csv",
               header = TRUE,
               fileEncoding = "UTF-8-BOM",
               strip.white = TRUE)
df <- df[, -5]                # drop empty last column
names(df)[1] <- "tool"         # ensure first column is named "tool"

# ----------------------------------------------------------------------
# 2. Transform celltype to percentage labels and proper order
# ----------------------------------------------------------------------
df <- df %>%
  mutate(
    percent = as.integer(gsub("cell", "", celltype)),
    celltype_pct = factor(paste0(percent, "%"),
                          levels = paste0(sort(unique(percent)), "%")),
    # Add log10-transformed time (time > 0, so no +1 needed)
    log10_time = log10(time)
  )

# ----------------------------------------------------------------------
# 3. Compute summary statistics on log10_time (for line plot)
# ----------------------------------------------------------------------
summary_df <- df %>%
  group_by(tool, histone, celltype_pct) %>%
  summarise(
    median_log10_time = median(log10_time, na.rm = TRUE),
    mean_log10_time   = mean(log10_time, na.rm = TRUE),
    sd_log10_time     = sd(log10_time, na.rm = TRUE),
    n                 = n(),
    se_log10_time     = sd_log10_time / sqrt(n),
    .groups           = "drop"
  ) %>%
  mutate(
    tool    = factor(tool),
    histone = factor(histone)
  )

# ----------------------------------------------------------------------
# 4. Define a consistent colour palette for the tools
# ----------------------------------------------------------------------
all_tools <- levels(factor(df$tool))
nature_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                   "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
                   "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39")
tool_colors <- setNames(nature_colors[1:length(all_tools)], all_tools)

# ----------------------------------------------------------------------
# 5. Enhanced Nature Methods theme (used for both plots)
# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
# 6. Create the line plot (median log10(time) + SE) – without title
# ----------------------------------------------------------------------
p1_time <- ggplot(summary_df, aes(x = celltype_pct, y = median_log10_time)) +
  geom_line(aes(group = tool, color = tool), 
            linewidth = 1.4, alpha = 0.9, show.legend = TRUE,
            lineend = "round", linejoin = "round") +
  geom_point(aes(fill = tool), 
             shape = 21, size = 5, color = "white", stroke = 1.2,
             alpha = 0.95) +
  geom_errorbar(aes(ymin = median_log10_time - se_log10_time,
                    ymax = median_log10_time + se_log10_time,
                    color = tool),
                width = 0.3, linewidth = 1.0, alpha = 0.8) +
  facet_wrap(~ histone, nrow = 1, scales = "free_x") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = tool_colors, name = "Tool") +
  scale_fill_manual(values = tool_colors, name = "Tool") +
  labs(
    title = NULL,                                # remove title for combined figure
    x = "Percentage of cells",
    y = expression(log[10](Time) ~ "(seconds)")
  ) +
  nature_methods_enhanced_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                               margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold")
  )

# ----------------------------------------------------------------------
# 7. Create the box plot (distribution of log10(time)) – without title
# ----------------------------------------------------------------------
p2_time <- ggplot(df, aes(x = tool, y = log10_time, fill = tool)) +
  geom_boxplot(outlier.size = 1.5, outlier.shape = 21, outlier.fill = "white", 
               outlier.color = "black", lwd = 0.4, alpha = 0.8) +
  facet_wrap(~ histone, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = tool_colors) +
  labs(
    title = NULL,                                # remove title
    x = "Tool",
    y = expression(log[10](Time) ~ "(seconds)")
  ) +
  nature_methods_enhanced_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13,
                               margin = margin(t = 8), face = "bold"),
    axis.text.y = element_text(size = 13, margin = margin(r = 8), face = "bold")
  )

# ----------------------------------------------------------------------
# 8. Combine plots with patchwork and add overall title & tags
# ----------------------------------------------------------------------
compact_peaks <- (p1_time / p2_time) +
  plot_annotation(
    title = "hPBMC Running time of peak callers",
    subtitle = "A: Performance across cell types | B: Distribution by method",
    tag_levels = "A",                            # automatically adds "A" and "B"
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, 
                                margin = margin(b = 5), color = "#2C3E50"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10),
                                   color = "#566573")
    )
  ) &
  # Apply consistent smaller font sizes for the combined figure
  theme(
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ----------------------------------------------------------------------
# 9. Display the combined plot
# ----------------------------------------------------------------------
print(compact_peaks)

# ----------------------------------------------------------------------
# 10. Save the figure as PDF with correct dimensions
#     (double‑column width ≈ 7.2 inches, height proportionally scaled)
# ----------------------------------------------------------------------
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/saturation_with_scCTPEAK/HumanPBMC_saturation_data_running_time_with_percentage.pdf",
       plot = compact_peaks,
width = 16,
       height = 12,    # inches (adjust height as needed)
       device = cairo_pdf,
       bg = "white")

# Optionally, save as high‑resolution TIFF for journals that require it
# ggsave("/path/to/figure.tiff", plot = compact_peaks,
#        width = 7.2, height = 5, dpi = 300, compression = "lzw")

###################################################################################
### Figure 9 E
###################################################################################
# =============================================================================
# Running time by cell percentage for each tool (aggregated across histones)
# Data: /home/wahid/project_scHMTF/GSE195725_processed_data/saturation_with_scCTPEAK/HumanPBMC_running_time.csv
# Outputs: 
#   - HumanPBMC_running_time_by_celltype_log.pdf   (log10 scale)
#   - HumanPBMC_running_time_by_celltype_raw.pdf   (raw seconds)
# =============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------------
# 1. Read and clean data
# ----------------------------------------------------------------------
df <- read.csv("/home/wahid/project_scHMTF/GSE195725_processed_data/saturation_with_scCTPEAK/HumanPBMC_running_time.csv",
               header = TRUE,
               fileEncoding = "UTF-8-BOM",
               strip.white = TRUE)
df <- df[, -5]                # drop empty last column
names(df)[1] <- "tool"         # ensure first column is named "tool"

# ----------------------------------------------------------------------
# 2. Create cell percentage factor and log10 time
# ----------------------------------------------------------------------
df <- df %>%
  mutate(
    percent = as.integer(gsub("cell", "", celltype)),
    celltype_pct = factor(paste0(percent, "%"),
                          levels = paste0(sort(unique(percent)), "%")),
    log10_time = log10(time)
  )

# ----------------------------------------------------------------------
# 3. Summarise by tool and cell percentage (across histones) – for log and raw
# ----------------------------------------------------------------------
df_cell_summary_log <- df %>%
  group_by(tool, celltype_pct) %>%
  summarise(
    median = median(log10_time, na.rm = TRUE),
    sd     = sd(log10_time, na.rm = TRUE),   # variability across histones
    n      = n(),
    .groups = "drop"
  ) %>%
  mutate(tool = factor(tool))

df_cell_summary_raw <- df %>%
  group_by(tool, celltype_pct) %>%
  summarise(
    median = median(time, na.rm = TRUE),
    sd     = sd(time, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  ) %>%
  mutate(tool = factor(tool))

# ----------------------------------------------------------------------
# 4. Nature‑Methods theme (same as before)
# ----------------------------------------------------------------------
nature_methods_enhanced_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "Helvetica", color = "black", size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, 
                                margin = margin(b = 12), color = "#2C3E50"),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "#566573",
                                   margin = margin(b = 18)),
      axis.title = element_text(face = "bold", size = 14, color = "#2C3E50"),
      axis.title.x = element_text(margin = margin(t = 12), vjust = -0.5),
      axis.title.y = element_text(margin = margin(r = 12), vjust = 2),
      axis.text = element_text(size = 11, color = "#34495E"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                 margin = margin(t = 5), face = "bold"),
      axis.text.y = element_text(margin = margin(r = 5), face = "bold"),
      axis.line = element_line(color = "#2C3E50", linewidth = 0.8),
      axis.ticks = element_line(color = "#2C3E50", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(face = "bold", size = 13, color = "#2C3E50",
                                  margin = margin(b = 8)),
      legend.text = element_text(size = 12, color = "#34495E",
                                 margin = margin(r = 15)),
      legend.position = "bottom",
      legend.key.size = unit(0.6, "cm"),
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
                                  margin = margin(t = 15), hjust = 1)
    )
}

# ----------------------------------------------------------------------
# 5. Define tool colors (Nature‑Methods palette)
# ----------------------------------------------------------------------
all_tools <- unique(df$tool)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
tool_colors <- setNames(nature_methods_colors[1:length(all_tools)], all_tools)

# ========================= LOG10 VERSION ================================
# ----------------------------------------------------------------------
# 6a. Set y‑axis limits automatically for log10 version
# ----------------------------------------------------------------------
y_min_log <- min(df_cell_summary_log$median - df_cell_summary_log$sd, na.rm = TRUE)
y_max_log <- max(df_cell_summary_log$median + df_cell_summary_log$sd, na.rm = TRUE)
y_pad_log <- 0.05 * (y_max_log - y_min_log)
y_limits_log <- c(y_min_log - y_pad_log, y_max_log + y_pad_log)

# ----------------------------------------------------------------------
# 7a. Create the log10 plot
# ----------------------------------------------------------------------
p_cell_log <- ggplot(df_cell_summary_log, aes(x = celltype_pct, y = median, 
                                               group = tool, color = tool)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = median - sd, ymax = median + sd),
                width = 0.2, size = 0.8, alpha = 0.5) +
  scale_color_manual(values = tool_colors, name = "Peak caller") +
  scale_y_continuous(limits = y_limits_log) +
  labs(
    title = "hPBMC: Running time by cell percentage for each tool",
    x = "Cell percentage",
    y = expression(log[10]("Time (seconds)"))
  ) +
  nature_methods_enhanced_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# ----------------------------------------------------------------------
# 8a. Save the log10 plot
# ----------------------------------------------------------------------
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/saturation_with_scCTPEAK/HumanPBMC_running_time_by_celltype_log.pdf",
       plot = p_cell_log,
       width = 7, height = 7,
       device = cairo_pdf,
       bg = "white")

print(p_cell_log)

# =============================================================================
# Running time by cell percentage for each tool (aggregated across histones)
# Data: /home/wahid/project_scHMTF/GSE195725_processed_data/saturation_with_scCTPEAK/HumanPBMC_running_time.csv
# Outputs: 
#   - HumanPBMC_running_time_by_celltype_log.pdf   (log10 scale)
#   - HumanPBMC_running_time_by_celltype_raw.pdf   (raw seconds)
# =============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------------
# 1. Read and clean data
# ----------------------------------------------------------------------
df <- read.csv("/home/wahid/project_scHMTF/GSE157637_processed_data/saturation_with_scCTPEAK/MouseBrain_running_time.csv",
               header = TRUE,
               fileEncoding = "UTF-8-BOM",
               strip.white = TRUE)
df <- df[, -5]                # drop empty last column
names(df)[1] <- "tool"         # ensure first column is named "tool"

# ----------------------------------------------------------------------
# 2. Create cell percentage factor and log10 time
# ----------------------------------------------------------------------
df <- df %>%
  mutate(
    percent = as.integer(gsub("cell", "", celltype)),
    celltype_pct = factor(paste0(percent, "%"),
                          levels = paste0(sort(unique(percent)), "%")),
    log10_time = log10(time)
  )

# ----------------------------------------------------------------------
# 3. Summarise by tool and cell percentage (across histones) – for log and raw
# ----------------------------------------------------------------------
df_cell_summary_log <- df %>%
  group_by(tool, celltype_pct) %>%
  summarise(
    median = median(log10_time, na.rm = TRUE),
    sd     = sd(log10_time, na.rm = TRUE),   # variability across histones
    n      = n(),
    .groups = "drop"
  ) %>%
  mutate(tool = factor(tool))

df_cell_summary_raw <- df %>%
  group_by(tool, celltype_pct) %>%
  summarise(
    median = median(time, na.rm = TRUE),
    sd     = sd(time, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  ) %>%
  mutate(tool = factor(tool))

# ----------------------------------------------------------------------
# 4. Nature‑Methods theme (same as before)
# ----------------------------------------------------------------------
nature_methods_enhanced_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "Helvetica", color = "black", size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, 
                                margin = margin(b = 12), color = "#2C3E50"),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "#566573",
                                   margin = margin(b = 18)),
      axis.title = element_text(face = "bold", size = 14, color = "#2C3E50"),
      axis.title.x = element_text(margin = margin(t = 12), vjust = -0.5),
      axis.title.y = element_text(margin = margin(r = 12), vjust = 2),
      axis.text = element_text(size = 11, color = "#34495E"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                 margin = margin(t = 5), face = "bold"),
      axis.text.y = element_text(margin = margin(r = 5), face = "bold"),
      axis.line = element_line(color = "#2C3E50", linewidth = 0.8),
      axis.ticks = element_line(color = "#2C3E50", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(face = "bold", size = 13, color = "#2C3E50",
                                  margin = margin(b = 8)),
      legend.text = element_text(size = 12, color = "#34495E",
                                 margin = margin(r = 15)),
      legend.position = "bottom",
      legend.key.size = unit(0.6, "cm"),
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
                                  margin = margin(t = 15), hjust = 1)
    )
}

# ----------------------------------------------------------------------
# 5. Define tool colors (Nature‑Methods palette)
# ----------------------------------------------------------------------
all_tools <- unique(df$tool)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
tool_colors <- setNames(nature_methods_colors[1:length(all_tools)], all_tools)

# ========================= LOG10 VERSION ================================
# ----------------------------------------------------------------------
# 6a. Set y‑axis limits automatically for log10 version
# ----------------------------------------------------------------------
y_min_log <- min(df_cell_summary_log$median - df_cell_summary_log$sd, na.rm = TRUE)
y_max_log <- max(df_cell_summary_log$median + df_cell_summary_log$sd, na.rm = TRUE)
y_pad_log <- 0.05 * (y_max_log - y_min_log)
y_limits_log <- c(y_min_log - y_pad_log, y_max_log + y_pad_log)

# ----------------------------------------------------------------------
# 7a. Create the log10 plot
# ----------------------------------------------------------------------
p_cell_log <- ggplot(df_cell_summary_log, aes(x = celltype_pct, y = median, 
                                               group = tool, color = tool)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = median - sd, ymax = median + sd),
                width = 0.2, size = 0.8, alpha = 0.5) +
  scale_color_manual(values = tool_colors, name = "Peak caller") +
  scale_y_continuous(limits = y_limits_log) +
  labs(
    title = "mBrain: Running time by cell percentage for each tool",
    x = "Cell percentage",
    y = expression(log[10]("Time (seconds)"))
  ) +
  nature_methods_enhanced_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# ----------------------------------------------------------------------
# 8a. Save the log10 plot
# ----------------------------------------------------------------------
ggsave("/home/wahid/project_scHMTF/GSE157637_processed_data/saturation_with_scCTPEAK/MouseBrain_running_time_by_celltype_log.pdf",
       plot = p_cell_log,
       width = 7, height = 7,
       device = cairo_pdf,
       bg = "white")


###################################################################################
### Figure 9 F
###################################################################################
# =============================================================================
# Combined script: two heatmaps with equal total height
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)

# ----------------------------------------------------------------------
# Helper: Nature‑Methods theme (identical for both plots)
# ----------------------------------------------------------------------
nature_heatmap_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                 size = 12, face = "bold", color = "#2C3E50"),
      axis.text.y = element_text(size = 12, face = "bold", color = "#2C3E50"),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(face = "bold", size = 11, color = "#2C3E50"),
      legend.text = element_text(size = 10, color = "#34495E"),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.5, "cm"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                                margin = margin(b = 10), color = "#2C3E50"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ----------------------------------------------------------------------
# 1. FIRST HEATMAP: peak‑caller ranks (rank_tools.csv)
# ----------------------------------------------------------------------
df_tools <- read_csv("/home/wahid/project_scHMTF/GSE195725_processed_data/rank_figure/rank_tools.csv",
                     col_names = TRUE)

tools <- colnames(df_tools)[-1]
n_tools <- length(tools)

df_tools_long <- df_tools %>%
  pivot_longer(cols = all_of(tools), names_to = "Tool", values_to = "Rank") %>%
  rename(Metric = 1) %>%
  mutate(
    Rank = ifelse(Rank == 0, NA, Rank),
    Metric = factor(Metric, levels = rev(unique(Metric))),
    Tool   = factor(Tool,   levels = tools)
  )

n_metrics_tools <- length(unique(df_tools_long$Metric))

# Tile dimensions for first plot (fixed)
tile_width  <- 0.8          # cm (same for both plots)
tile_height <- 0.3          # cm (fixed for first plot)

panel_aspect_tools <- (n_metrics_tools * tile_height) / (n_tools * tile_width)

puBu_colors <- brewer.pal(7, "PuBu")        # original: light → dark (1→7)
puBu_colors_rev <- rev(puBu_colors)         # reversed: dark → light (1→7)

# FIXED: use df_tools_long, not df_long, and panel_aspect_tools
p1 <- ggplot(df_tools_long, aes(x = Tool, y = Metric, fill = Rank)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = puBu_colors_rev,                # now 1 = darkest, 7 = lightest
    na.value = "white",
    name = "Rank",     # clear legend title
    breaks = 1:7,
    limits = c(1, 7),
    guide = guide_colorbar(reverse = TRUE)   # puts rank 1 at top, rank 7 at bottom
  ) +
  labs(
    title = "Performance ranking of peak callers across metrics",
    x = NULL, y = NULL
  ) +
  nature_heatmap_theme() +
  theme(aspect.ratio = panel_aspect_tools)   # CORRECTED

extra_width  <- 3.5   # cm (shared with second plot)
extra_height <- 2.5   # cm (shared with second plot)

plot_width_tools  <- n_tools * tile_width + extra_width
plot_height_tools <- n_metrics_tools * tile_height + extra_height

# Save first plot
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/rank_figure/heatmap_binary_tool_updated.pdf",
       plot = p1, width = plot_width_tools, height = plot_height_tools,
       device = cairo_pdf, bg = "white")

cat("First plot height:", plot_height_tools, "cm\n")

# ----------------------------------------------------------------------
# 2. SECOND HEATMAP: binary input (rank_input.csv)
# ----------------------------------------------------------------------
df_input <- read_csv("/home/wahid/project_scHMTF/GSE195725_processed_data/rank_figure/rank_input.csv",
                     col_names = TRUE)

conditions <- colnames(df_input)[2:3]
n_conditions <- length(conditions)

df_input_long <- df_input %>%
  pivot_longer(cols = all_of(conditions), names_to = "Condition", values_to = "Value") %>%
  rename(Metric = 1) %>%
  mutate(
    Value = as.numeric(Value),
    Metric = factor(Metric, levels = rev(unique(Metric))),
    Condition = factor(Condition, levels = conditions)
  )

n_metrics_input <- length(unique(df_input_long$Metric))

# ----------------------------------------------------------------------
# Ensure total height matches first plot
# ----------------------------------------------------------------------
# If number of metrics is the same, use same tile_height
if (n_metrics_input == n_metrics_tools) {
  tile_height_input <- tile_height
  message("Number of metrics matches. Using tile_height = ", tile_height, " cm")
} else {
  # Adjust tile_height so that total height equals first plot's height
  # plot_height = n_metrics * tile_height + extra_height
  # => tile_height = (plot_height_tools - extra_height) / n_metrics_input
  tile_height_input <- (plot_height_tools - extra_height) / n_metrics_input
  message("Adjusting tile_height for second plot to ", round(tile_height_input, 3), " cm")
}

panel_aspect_input <- (n_metrics_input * tile_height_input) / (n_conditions * tile_width)

p2 <- ggplot(df_input_long, aes(x = Condition, y = Metric, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "#F0F0F0", high = "darkgreen",
    na.value = "white",
    name = "Valu",
    breaks = c(0, 1),
    limits = c(0, 1)
  ) +
  labs(title = "Binary metrics: noInput vs withInput") +
  nature_heatmap_theme() +
  theme(aspect.ratio = panel_aspect_input)

plot_width_input  <- n_conditions * tile_width + extra_width
plot_height_input <- n_metrics_input * tile_height_input + extra_height

ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/rank_figure/heatmap_binary_input_updated.pdf",
       plot = p2, width = plot_width_input, height = plot_height_input,
       device = cairo_pdf, bg = "white")

cat("Second plot height:", plot_height_input, "cm\n")
