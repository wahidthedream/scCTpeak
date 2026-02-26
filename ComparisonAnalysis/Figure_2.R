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


