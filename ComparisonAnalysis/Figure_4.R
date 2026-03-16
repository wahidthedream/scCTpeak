#############################################################################################
### Figure 4 A, B
#############################################################################################

###========================================================================================###
#### Comparison_PeakCalling_tools_usingFRiPs for HumanPBMC
###========================================================================================###


# Load required packages
library(readr)      # For reading CSV
library(dplyr)      # For data manipulation
library(knitr)      # For kable
library(kableExtra) # For styling table

# Read the CSV file
data <- read_csv("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_frip/HumanPBMC_combined_peak_frip_summary_with_input.csv")

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

  # Display the table in Nature Methods style
data %>%
  kable("html", caption = "Summary of scCUT&Tag with Input Control", align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                full_width = FALSE, 
                position = "center", 
                font_size = 14) %>%
  row_spec(0, bold = TRUE, background = "#f7f7f7")



library(tidyverse)
library(forcats)
library(scales)

# Verify columns exist
stopifnot(c("Histone", "Celltype", "Method", "FRiP") %in% names(data))

# Reorder samples histone-wise
df_plot <- data %>%
  filter(Method %in% c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")) %>%
  group_by(Histone) %>%
  mutate(
    Sample = paste(Histone, Celltype, sep = "_") %>%
      fct_reorder(FRiP, .fun = median) %>%
      fct_rev()
  ) %>%
  ungroup() %>%
  mutate(
    FRiP_percent = FRiP * 100,
    Method = factor(Method, levels = c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", 
                                      "SICER2"))
  ) %>%
  dplyr::select(Sample, Method, FRiP_percent, Histone)

library(ggplot2)
library(scales)
library(extrafont)  # for extra fonts

# Plot 
p <- ggplot(df_plot, aes(x = FRiP_percent, y = Sample, fill = Method)) + 
  geom_col(width = 0.7, position = position_dodge(width = 0.8), color = "black", linewidth = 0.2) +
  facet_wrap(~ Method, nrow = 1) +
  scale_fill_manual(values = hue_pal()(7)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "% of fragments in peaks", 
    y = NULL,
    title = "Fraction of reads in peaks by peak calling method",
    subtitle = "Comparison across different computational approaches"
  ) +
  theme_classic(base_size = 8) +
  theme(
    text = element_text(family = "Helvetica"),  
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7, margin = margin(t = 3)),
    strip.text = element_text(size = 7, face = "bold", margin = margin(b = 2)),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(3, "mm"),
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
    plot.subtitle = element_text(size = 7, hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )

# PDF save 
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_frip/Figure/HumanPBMC_Comparison_PeakCalling_tools_usingFRiP_with_input.pdf",
       plot = p,
       device = "pdf",
       width = 7, 
       height = 8,
       units = "in")
print(p)


###========================================================================================###
#### FRiP SCATTER PLOT for HumanPBMC
#### Without input case 
###========================================================================================###

# Load required packages
library(readr)      # For reading CSV
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting
library(scales)     # For percentage formatting
library(ggsci)      # For color palettes

# Read the CSV file
data <- read_csv("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_frip/HumanPBMC_combined_peak_frip_summary_without_input.csv")

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

# Define Nature Methods color palette for methods
all_methods <- unique(data$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

# Create FRiP summary data
frip_summary <- data %>%
  group_by(Method, Histone) %>%
  summarize(
    median_FRiP = median(FRiP, na.rm = TRUE),
    mean_FRiP = mean(FRiP, na.rm = TRUE),
    sd_FRiP = sd(FRiP, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(Method),
    Histone = factor(Histone)
  )


###========================================================================================###
# SCATTER PLOT 2: With error bars
###========================================================================================###

frip_scatter2 <- ggplot(frip_summary, 
                       aes(x = Histone, 
                           y = median_FRiP,
                           color = Method)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pmax(median_FRiP - sd_FRiP, 0),  # Cap at 0
                    ymax = pmin(median_FRiP + sd_FRiP, 1)), # Cap at 1
                width = 0.3, 
                position = position_dodge(width = 0.5),
                alpha = 0.7) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1),  # Explicitly set 0-100% range
                     breaks = seq(0, 1, by = 0.2)) +  # 0%, 20%, 40%, etc.
  labs(
    title = "FRiP Values with Variability by Histone and Method",
    x = "Histone Modification",
    y = "Median FRiP (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )
print(frip_scatter2)

ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_frip_ex/Figure/HumanPBMC_FRiP_Scatter_with_ErrorBars_Without_Input.pdf",
       plot = frip_scatter2,
       width = 14,
       height = 8,
       device = "pdf",
       bg = "white")


#############################################################################################
### Figure 4 C, D
#############################################################################################


##################
### For Human PBMC 
###################

#!/bin/bash

# Set up logging
LOG_FILE="/home/wahid/project_scHMTF/GSE195725_processed_data/BAM/bw_deeptools/result_deeptools/with_input/analysis.log"
mkdir -p "$(dirname "$LOG_FILE")"
echo "=== Analysis started: $(date) ===" | tee "$LOG_FILE"

# Analyze histone mark signals around peak centers
for histone in H3K27me3 H3K4me1 H3K4me2 H3K4me3 H3K9me3; do 
    for celltype in B CD4T CD8T DC Mono NK other otherT; do
        for caller in DROMPAplus Genrich GoPeaks HOMER MACS2 SEACR SICER2; do
            bw_file="/home/wahid/project_scHMTF/GSE195725_processed_data/BAM/bw_deeptools/${histone}/bam/sorted_bams/bigwig_files/${histone}_${celltype}.bw"
            bed_file="/home/wahid/project_scHMTF/GSE195725_processed_data/BAM/bw_deeptools/peakbed_with_input/${caller}_${histone}_${celltype}.bed"
            output_dir="/home/wahid/project_scHMTF/GSE195725_processed_data/BAM/bw_deeptools/result_deeptools/with_input/"
            output_file="${output_dir}/matrix_${caller}_${histone}_${celltype}.gz"  # Fixed: use ${histone} not H3K27ac
            
            # Create output directory if it doesn't exist
            mkdir -p "$output_dir"
            
            # Check if both files exist and BED file is not empty
            if [[ -f "$bw_file" && -f "$bed_file" ]]; then
                bed_line_count=$(wc -l < "$bed_file" 2>/dev/null)
                if [[ $bed_line_count -gt 0 ]]; then
                    echo "Processing ${histone} ${caller} ${celltype}..." | tee -a "$LOG_FILE"
                    echo "  BED file: $bed_file ($bed_line_count peaks)" | tee -a "$LOG_FILE"
                    
                    if computeMatrix reference-point \
                        -S "$bw_file" \
                        -R "$bed_file" \
                        --referencePoint center \
                        -a 2000 -b 2000 \
                        --binSize 50 \
                        --missingDataAsZero \
                        -p 32 \
                        -o "$output_file" 2>> "$LOG_FILE"; then
                        
                        # Verify the output file was created
                        if [[ -f "$output_file" ]]; then
                            echo "✓ Successfully created matrix for ${histone} ${caller} ${celltype}" | tee -a "$LOG_FILE"
                        else
                            echo "✗ Output file not created for ${histone} ${caller} ${celltype}" | tee -a "$LOG_FILE"
                        fi
                    else
                        echo "✗ computeMatrix failed for ${histone} ${caller} ${celltype}" | tee -a "$LOG_FILE"
                        # Remove failed output file if it exists
                        [[ -f "$output_file" ]] && rm "$output_file"
                    fi
                else
                    echo "⚠ Empty BED file for ${histone} ${caller} ${celltype}: $bed_file" | tee -a "$LOG_FILE"
                fi
            else
                missing_files=""
                [[ ! -f "$bw_file" ]] && missing_files="BigWig"
                [[ ! -f "$bed_file" ]] && missing_files="${missing_files}${missing_files:+, }BED"
                echo "⚠ Missing files for ${histone} ${caller} ${celltype}: $missing_files" | tee -a "$LOG_FILE"
            fi
        done
    done
done 

echo "=== Analysis completed: $(date) ===" | tee -a "$LOG_FILE"


#####==================================================================================================
### profile-plots and heatmap plot for Human PBMC and Mouse Brain (with and without input)
#####==================================================================================================

#!/bin/bash
cd /home/wahid/project_scHMTF/GSE195725_processed_data/BAM/bw_deeptools/result_deeptools/with_input
# Create output directories
mkdir -p profile_plots
mkdir -p heatmaps

# Function to generate plots for one matrix
generate_plots() {
    local matrix_file="$1"
    local base_name=$(basename "$matrix_file" .gz)
    
    # Extract components from filename
    IFS='_' read -ra parts <<< "$base_name"
    local caller="${parts[1]}"
    local histone="${parts[2]}"
    local celltype="${parts[3]}"
    
    echo "Generating plots for: $caller $histone $celltype"
    
    # Generate profile plot
    plotProfile \
        -m "$matrix_file" \
        -o "profile_plots/${base_name}_profile.pdf" \
        --perGroup \
        --plotTitle "${histone} Signal at ${celltype} Peaks (${caller})" \
        --yAxisLabel "${histone} Signal" \
        --legendLocation best \
        --regionsLabel "Peaks" \
        --plotWidth 10 \
        --plotHeight 8 \
        --dpi 300
    
    # Generate heatmap
    plotHeatmap \
        -m "$matrix_file" \
        -o "heatmaps/${base_name}_heatmap.pdf" \
        --colorMap RdBu \
        --whatToShow 'plot, heatmap and colorbar' \
        --plotTitle "${histone} Signal at ${celltype} Peaks (${caller})" \
        --refPointLabel "Peak Center" \
        --xAxisLabel "Distance from peak center (bp)" \
        --regionsLabel "${celltype} Peaks" \
        --sortRegions descend \
        --zMin 0 \
        --legendLocation none \
        --dpi 300
    
    echo "✓ Created plots for $caller $celltype"
}

# Export function for parallel processing
export -f generate_plots

# Generate plots for all matrix files in parallel
find . -maxdepth 1 -name "matrix_*.gz" -type f | \
    parallel -j 8 --progress generate_plots {}

echo "=== All plots generated ==="
echo "Profile plots: $(find profile_plots -name "*.pdf" | wc -l) files"
echo "Heatmaps: $(find heatmaps -name "*.pdf" | wc -l) files"


#####==================================================================================================
### Figure 4 D
#####==================================================================================================


# calculate_overall_means_parsed.R
library(data.table)
library(dplyr)

# Set working directory
setwd("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_deeptools/without_input/middle_values")

# Get all middle_8_values.txt files
files <- list.files(pattern = "_middle_8_values\\.txt$")
cat("Found", length(files), "files to process\n")

# Create empty results table
results <- data.table(
  filename = character(),
  histone = character(),
  tool = character(),
  celltype = character(),
  overall_mean = numeric(),
  regions = integer(),
  total_values = integer()
)

# Process each file
for (i in seq_along(files)) {
  file <- files[i]
  cat("Processing", i, "of", length(files), ":", file, "\n")
  
  tryCatch({
    # Read only the value columns (columns 5-12 for 8 values)
    data <- fread(file, 
                  sep = "\t", 
                  header = FALSE,
                  select = 5:12)  # Only read the 8 value columns
    
    # Calculate overall mean
    all_values <- unlist(data)
    overall_mean <- mean(all_values)
    
    # Count regions (rows) and total values
    regions <- nrow(data)
    total_values <- length(all_values)
    
    # Parse filename into components
    # Example: "matrix_DROMPAplus_H3K27ac-b_B_middle_8_values.txt"
    # Remove "matrix_" prefix and "_middle_8_values.txt" suffix
    base_name <- gsub("^matrix_", "", file)
    base_name <- gsub("_middle_8_values\\.txt$", "", base_name)
    
    # Split by underscore
    parts <- unlist(strsplit(base_name, "_"))
    
    if (length(parts) >= 3) {
      tool <- parts[1]
      histone <- parts[2]
      celltype <- paste(parts[3:length(parts)], collapse = "_")  # In case celltype has underscores
    } else if (length(parts) == 2) {
      tool <- parts[1]
      histone <- parts[2]
      celltype <- "unknown"
    } else {
      tool <- "unknown"
      histone <- "unknown"
      celltype <- "unknown"
    }
    
    # Add to results
    results <- rbind(results, 
                     data.table(filename = file,
                                histone = histone,
                                tool = tool,
                                celltype = celltype,
                                overall_mean = overall_mean,
                                regions = regions,
                                total_values = total_values))
    
  }, error = function(e) {
    cat("Error processing", file, ":", e$message, "\n")
    results <- rbind(results, 
                     data.table(filename = file,
                                histone = "error",
                                tool = "error",
                                celltype = "error",
                                overall_mean = NA,
                                regions = NA,
                                total_values = NA))
  })
}

# Save results
fwrite(results, "overall_means_parsed_without_input.tsv", sep = "\t")

# Print summary
cat("\n=== Processing Complete ===\n")
cat("Processed", nrow(results), "files\n")
cat("Unique histone marks:", n_distinct(results$histone), "\n")
cat("Unique tools:", n_distinct(results$tool), "\n")
cat("Unique celltypes:", n_distinct(results$celltype), "\n")

# Show first few results
print(head(results, 10))

pbmc_with_input<-results

#####==================================================================================================
## Plot - Corrected for Human PBMC data (APA Peak Height)
## Two versions: original scale and log-transformed - SAVED SEPARATELY
#####==================================================================================================

library(ggplot2)
library(dplyr)
library(scales)

# First, check the actual range of your data
cat("=== DATA RANGE CHECK ===\n")
cat("Min overall_mean:", format(min(pbmc_with_input$overall_mean), scientific = TRUE, digits = 3), "\n")
cat("Max overall_mean:", format(max(pbmc_with_input$overall_mean), scientific = TRUE, digits = 3), "\n")
cat("Mean overall_mean:", format(mean(pbmc_with_input$overall_mean), scientific = TRUE, digits = 3), "\n")
cat("SD overall_mean:", format(sd(pbmc_with_input$overall_mean), scientific = TRUE, digits = 3), "\n\n")

# Add log-transformed column (add small constant to avoid log(0))
pbmc_with_input <- pbmc_with_input %>%
  mutate(log_overall_mean = log10(overall_mean + 1e-10))  # Add small constant

cat("=== LOG-TRANSFORMED DATA CHECK ===\n")
cat("Min log_overall_mean:", min(pbmc_with_input$log_overall_mean, na.rm = TRUE), "\n")
cat("Max log_overall_mean:", max(pbmc_with_input$log_overall_mean, na.rm = TRUE), "\n")
cat("Mean log_overall_mean:", mean(pbmc_with_input$log_overall_mean, na.rm = TRUE), "\n\n")

# Create summary statistics for original scale
score_summary_original <- pbmc_with_input %>%
  group_by(tool, histone) %>%
  summarize(
    median_score = median(overall_mean, na.rm = TRUE),
    mean_score = mean(overall_mean, na.rm = TRUE),
    sd_score = sd(overall_mean, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(tool),
    Histone = factor(histone),
    # Calculate error bars
    se_score = sd_score / sqrt(n_samples),
    ci_lower = median_score - 1.96 * se_score,
    ci_upper = median_score + 1.96 * se_score
  )

# Create summary statistics for log-transformed scale
score_summary_log <- pbmc_with_input %>%
  group_by(tool, histone) %>%
  summarize(
    median_log_score = median(log_overall_mean, na.rm = TRUE),
    mean_log_score = mean(log_overall_mean, na.rm = TRUE),
    sd_log_score = sd(log_overall_mean, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(tool),
    Histone = factor(histone),
    # Calculate error bars for log scale
    se_log_score = sd_log_score / sqrt(n_samples),
    ci_lower_log = median_log_score - 1.96 * se_log_score,
    ci_upper_log = median_log_score + 1.96 * se_log_score
  )

# Define colors
all_methods <- unique(score_summary_original$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

###############################################################################
# PLOT 1: Original Scale - SAVE SEPARATELY
###############################################################################

# Calculate y-axis limits for original scale
y_min_orig <- min(score_summary_original$ci_lower, na.rm = TRUE)
y_max_orig <- max(score_summary_original$ci_upper, na.rm = TRUE)
padding_orig <- (y_max_orig - y_min_orig) * 0.1
y_min_plot_orig <- y_min_orig - padding_orig
y_max_plot_orig <- y_max_orig + padding_orig

cat(sprintf("\nY-axis limits (Original): %.6f to %.6f\n", y_min_plot_orig, y_max_plot_orig))

plot1_original <- ggplot(score_summary_original, 
                         aes(x = Histone, 
                             y = median_score,
                             color = Method)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = ci_lower,
        ymax = ci_upper),
    width = 0.2, 
    position = position_dodge(width = 0.7),
    alpha = 0.7,
    size = 0.7
  ) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    limits = c(-2,10),
    breaks = seq(-2,10, by=2)
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "red", 
    alpha = 0.5,
    size = 0.5
  ) +
  labs(
    title = "Human PBMC without Input - Original Scale",
    x = "Histone Modification",
    y = "Mean Peak Signal"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

# Print and save PLOT 1
print(plot1_original)
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_deeptools/signal_intensity_figure/Signal_intensity_Original_hPBMC_withoutinput.pdf", 
       plot1_original, width = 8, height = 6, dpi = 300)



#############################################################################################
### Figure 4 E, G
#############################################################################################

####==================================================================
#### Grouth truth bulk ChiPseq from encode data for Human PBMC
####==================================================================

# === REQUIRED PACKAGES ===
required_packages <- c(
  "rtracklayer", "GenomicRanges", "ggplot2", "viridis",
  "dplyr", "ggpubr", "lme4", "emmeans", "multcomp", "patchwork", "tidyr"
)
suppressWarnings({
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))
  if(length(missing_packages) > 0) install.packages(missing_packages, quiet = TRUE)
  invisible(lapply(required_packages, library, character.only = TRUE))
})

set.seed(2023)

# === PARAMETERS ===
cell_types <- c("B", "CD4T", "CD8T", "Mono", "NK")
methods <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
color_palette <- "viridis"
plot_width <- 20  # cm
plot_height <- 20  # cm
min_overlap <- 100  # bp
significance_level <- 0.05

histones<-c("H3K4me3")
for(hist in histones){
# === DIRECTORIES ===
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkchipseq_encode_correction_label"
sc_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_with_input_peakbed_corrected/bulkchipseq_peakbed/", hist, "_peakbed/")
bulk_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/bulk_ChIPsep_encode/bulk", hist, "_peakbed/")

# === DATA IMPORT & PROCESSING ===

# Function to read BED files with variable columns
read_bed_file <- function(file_path) {
  tryCatch({
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      return(NULL)
    }
    
    # Read first line to determine column count
    first_line <- readLines(file_path, n = 1)
    n_cols <- length(strsplit(first_line, "\t")[[1]])
    
    # Read full file with appropriate columns
    bed <- read.table(file_path, sep = "\t", header = FALSE, 
                     col.names = paste0("V", seq_len(n_cols)),
                     fill = TRUE, stringsAsFactors = FALSE)
    
    # Standardize to 3 columns (chr, start, end)
    if (ncol(bed) >= 3) {
      return(bed[, 1:3])
    } else {
      warning("File has fewer than 3 columns: ", file_path)
      return(NULL)
    }
  }, error = function(e) {
    warning("Error reading file: ", file_path, "\n", e$message)
    return(NULL)
  })
}

# Import bulk peaks
bulk_peaks <- list()
for (cell_type in cell_types) {
  bulk_file <- file.path(bulk_peak_dir, paste0(hist, "_", cell_type, ".bed"))
  bed <- read_bed_file(bulk_file)
  if (!is.null(bed)) {
    suppressWarnings({
      bulk_peaks[[cell_type]] <- GRanges(
        seqnames = bed[,1],
        ranges = IRanges(start = bed[,2] + 1, end = bed[,3])
      )
    })
  } else {
    warning("Failed to import bulk peaks for: ", cell_type)
  }
}

# Import single-cell peaks
import_sc_peaks <- function() {
  sc_peaks <- list()
  for (method in methods) {
    for (cell_type in cell_types) {
      sc_file <- file.path(sc_peak_dir, paste0(method, "_", hist, "_", cell_type, ".bed"))
      bed <- read_bed_file(sc_file)
      if (!is.null(bed)) {
        suppressWarnings({
          sc_peaks[[paste(method, cell_type, sep = "_")]] <- GRanges(
            seqnames = bed[,1],
            ranges = IRanges(start = bed[,2] + 1, end = bed[,3]),
            method = method,
            cell_type = cell_type
          )
        })
      } else {
        warning("Failed to import sc peaks for: ", method, " - ", cell_type)
      }
    }
  }
  return(sc_peaks)
}

sc_peaks <- import_sc_peaks()

# Check if sc_peaks was created successfully
if(length(sc_peaks) == 0) {
  stop("No single-cell peaks were imported. Check your input files and paths.")
}

# Check which methods/cell types were successfully imported
imported_combinations <- sapply(sc_peaks, function(x) 
  paste0(unique(x$method), "_", unique(x$cell_type)))
print("Successfully imported combinations:")
print(imported_combinations)

# === CALCULATE METRICS ===
calculate_enhanced_metrics <- function(sc_peaks, bulk_peaks) {
  # Initialize empty results
  results <- data.frame()
  
  # Get all unique method-cell_type combinations
  existing_combos <- unique(data.frame(
    algorithm = unlist(lapply(sc_peaks, function(x) unique(x$method))),
    cell_type = unlist(lapply(sc_peaks, function(x) unique(x$cell_type))),
    stringsAsFactors = FALSE
  ))
  
  if(nrow(existing_combos) == 0) {
    warning("No valid method-cell_type combinations found")
    return(results)
  }
  
  for(i in 1:nrow(existing_combos)) {
    alg <- existing_combos$algorithm[i]
    ct <- existing_combos$cell_type[i]
    
    gr_name <- paste0(alg, "_", ct)
    if(!gr_name %in% names(sc_peaks)) next
    
    query <- sc_peaks[[gr_name]]
    bulk_gr <- bulk_peaks[[ct]]
    
    if(is.null(bulk_gr)) {
      warning("No bulk peaks for: ", ct)
      next
    }
    
    # Calculate overlaps
    ov <- findOverlaps(query, bulk_gr, minoverlap = min_overlap)
    overlapping_query <- length(unique(queryHits(ov)))
    overlapping_bulk <- length(unique(subjectHits(ov)))
    total_query <- length(query)
    total_bulk <- length(bulk_gr)
    
    # Calculate metrics
    metrics <- list(
      precision = ifelse(total_query > 0, overlapping_query/total_query, NA),
      recall = ifelse(total_bulk > 0, overlapping_bulk/total_bulk, NA),
      f1_score = NA,
      simpson_idx = ifelse(min(total_query, total_bulk) > 0, 
                         overlapping_bulk/min(total_query, total_bulk), NA)
    )
    
    # Calculate F1 if possible
    if(!any(is.na(c(metrics$precision, metrics$recall)))) {
      metrics$f1_score <- 2 * (metrics$precision * metrics$recall) / 
        (metrics$precision + metrics$recall)
    }
    
    # Add to results
    results <- rbind(results, data.frame(
      cell_type = ct,
      algorithm = alg,
      total_peaks = total_query,
      overlapping_peaks = overlapping_query,
      overlap_pct = ifelse(total_query > 0, overlapping_query/total_query * 100, NA),
      simpson_idx = metrics$simpson_idx,
      Precision = metrics$precision,
      Recall = metrics$recall,
      F1_score = metrics$f1_score,
      mean_width = mean(width(query), na.rm = TRUE),
      median_width = median(width(query), na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

results <- calculate_enhanced_metrics(sc_peaks, bulk_peaks)

# Check results
if(nrow(results) == 0) {
  stop("No results calculated. Check input data.")
}

# Process results
results <- results %>%
  mutate(
    algorithm = factor(algorithm, levels = methods),
    cell_type = factor(cell_type, levels = cell_types),
    label = paste0(algorithm, " (", cell_type, ")")
  ) %>%
  filter(!is.na(algorithm), !is.na(cell_type), is.finite(overlap_pct))

# === STATISTICAL ANALYSIS ===
perform_comprehensive_stats <- function(results) {
  cat("\n=== DATA QUALITY CHECK ===\n")
  print(summary(results))
  
  if(length(unique(results$algorithm)) < 2 || length(unique(results$cell_type)) < 2) {
    warning("Insufficient groups for statistical tests")
    return(NULL)
  }
  
  if(all(results$overlap_pct == 0)) {
    warning("All overlap percentages are zero")
    return(NULL)
  }
  
  cat("\n=== NORMALITY TEST ===\n")
  print(shapiro.test(sample(results$overlap_pct, min(500, nrow(results)))))
  
  cat("\n=== MIXED-EFFECTS MODEL ===\n")
  tryCatch({
    fit <- lmer(overlap_pct ~ algorithm + (1|cell_type), data = results)
    print(summary(fit))
    print(anova(fit))
    
    if(requireNamespace("emmeans", quietly = TRUE)) {
      cat("\n=== POST-HOC COMPARISONS ===\n")
      emm <- emmeans(fit, specs = pairwise ~ algorithm)
      print(emm$contrasts)
      return(list(model = fit, contrasts = emm$contrasts))
    }
    return(list(model = fit))
  }, error = function(e) {
    cat("Mixed model failed, using Kruskal-Wallis:\n")
    print(kruskal.test(overlap_pct ~ algorithm, data = results))
    return(NULL)
  })
}

stats_results <- perform_comprehensive_stats(results)

# === FIGURE GENERATION ===
# === FIGURE GENERATION ===
create_benchmark_plot <- function(results) {
  # Create labels in the format "H3K27ac_celltype (method)"
  results <- results %>%
    mutate(label = paste0(cell_type, "_", algorithm),
           # Create a sorting factor that groups by method first, then cell type
           sort_var = interaction(factor(algorithm, levels = methods), 
                                   factor(cell_type, levels = cell_types)))
  
  ggplot(results, aes(x = reorder(label, as.numeric(sort_var)), 
                      y = overlap_pct, fill = algorithm)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", overlap_pct)),
              position = position_stack(vjust = 1.02),
              hjust = 0, size = 3) +
    coord_flip() +
    # y-axis limit set 
    scale_y_continuous(
      limits = c(0, 100),           
      breaks = seq(0, 100, by = 10), 
      #labels = function(x) paste0(x, "%")#, # label
      expand = expansion(mult = c(0, 0.05)) # 
    ) +
    scale_fill_viridis_d(option = color_palette) +
    labs(x = NULL, y = "Overlap with Bulk Peaks (%)",
         title = paste0("Peak Caller Performance on - ", hist),
         fill = "Peak Calling Method") +
    theme_pubr() +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 10))
}
##sort_var = interaction(factor(cell_type, levels = cell_types),
##                       factor(algorithm, levels = methods))
create_prf1_plot <- function(results) {
  prf_data <- results %>%
    dplyr::select(algorithm, Precision, Recall, F1_score) %>%
    tidyr::pivot_longer(-algorithm, names_to = "metric", values_to = "value") %>%
    dplyr::group_by(algorithm, metric) %>%
    dplyr::summarise(mean = mean(value, na.rm = TRUE),
                     sd = sd(value, na.rm = TRUE), .groups = "drop")
  
  ggplot(prf_data, aes(x = algorithm, y = mean, fill = algorithm)) +
    geom_col(position = position_dodge(0.9), width = 0.7) +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
    geom_text(aes(label = sprintf("%.2f", mean)), vjust = -0.5) +
    facet_wrap(~metric, nrow = 1) +
    # y-axis limit 
    scale_y_continuous(
      limits = c(0, 1),               
      breaks = seq(0, 1, by = 0.2)#,    
      #expand = expansion(mult = c(0, 0.05))  
    ) +
    scale_fill_viridis_d(option = color_palette) +
    labs(x = "Peak Calling Algorithm", y = "Metric Value") +
    theme_pubr() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}
# Generate and save plots
plot_benchmark <- create_benchmark_plot(results)
plot_prf1 <- create_prf1_plot(results)

combined_plot <- (plot_benchmark / plot_prf1) + 
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = 'A')

# Save outputs
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, paste0("Figure4E_HumanPBMC_", hist, "_with_input.pdf")), plot_benchmark,
       width = plot_width, height = plot_height, units = "cm", dpi = 600)

ggsave(file.path(out_dir, paste0("Figure4F_HumanPBMC_", hist, "_prf1_plot_with_input.pdf")), plot_prf1,
       width = plot_width, height = plot_height/2, units = "cm", dpi = 600)

ggsave(file.path(out_dir, paste0("HumanPBMC_", hist, "_combined_plot_with_input.pdf")), combined_plot,
       width = plot_width, height = plot_height*1.5, units = "cm", dpi = 600)

# Save results
write.csv(results, file.path(out_dir, paste0("HumanPBMC_", hist, "_results_with_input.csv")), row.names = FALSE)

cat("\n=== ANALYSIS COMPLETED ===\n")
print(Sys.time())
}




#############################################################################################
### Figure 4 F, H
#############################################################################################


library(dplyr)
library(readr)
library(stringr)

# Define the directory where your files are located
input_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkchipseq_encode_correction_label"

# List of files
files <- c(
  "HumanPBMC_H3K27ac-b_results_without_input.csv",
  "HumanPBMC_H3K27ac-s_results_without_input.csv",
  "HumanPBMC_H3K27me3_results_without_input.csv",
  "HumanPBMC_H3K4me1_results_without_input.csv",
  "HumanPBMC_H3K4me3_results_without_input.csv",
  "HumanPBMC_H3K9me3_results_without_input.csv"
)

# Create an empty list to store dataframes
all_data <- list()

# Read and process each file
for (file in files) {
  # Extract histone mark from filename
  # Pattern: HumanPBMC_(histone-mark)_results_with_input.csv
  histone_mark <- str_extract(file, "(?<=HumanPBMC_)[^_]+(?=_results_without_input\\.csv)")
  
  # Full path to file
  file_path <- file.path(input_dir, file)
  
  # Read the CSV file
  if (file.exists(file_path)) {
    cat("Processing:", file, "\n")
    
    # Read the data
    data <- read_csv(file_path, col_types = cols())
    
    # Add histone mark column
    data <- data %>%
      mutate(histones = histone_mark) %>%
      select(histones, everything())  # Move histones to first column
    
    # Add to list
    all_data[[file]] <- data
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# Combine all dataframes
combined_data <- bind_rows(all_data)

# View the structure
cat("\nCombined data structure:\n")
glimpse(combined_data)

# View first few rows
cat("\nFirst few rows of combined data:\n")
print(head(combined_data))

# Save the combined data
output_file <- file.path(input_dir, "HumanPBMC_ALL_histones_combined_results_without_input.csv")
write_csv(combined_data, output_file)

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)

# =============================================================================
# PART 1: Direct Precision vs Recall with 1:1 aspect ratio
# =============================================================================

# Direct Precision vs Recall using individual data points
scatter_direct_pr <- combined_data %>% 
  group_by(algorithm, histones) %>%
  summarize(
    median_Precision = median(Precision, na.rm = TRUE),
    sd_Precision = sd(Precision, na.rm = TRUE),
    median_Recall = median(Recall, na.rm = TRUE),
    sd_Recall = sd(Recall, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(algorithm),
    Histone = factor(histones)
  )

# Get all unique methods and create a consistent color mapping
all_methods <- unique(scatter_direct_pr$Method)

# Nature methods color palette
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)

# Assign colors to methods
if (length(all_methods) <= length(nature_methods_colors)) {
  method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)
} else {
  # If more methods than colors, generate additional colors
  extra_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(all_methods) - length(nature_methods_colors))
  method_colors <- setNames(c(nature_methods_colors, extra_colors)[1:length(all_methods)], all_methods)
}

# Define shapes for histones
histone_names <- unique(scatter_direct_pr$Histone)
histone_shape_map <- setNames(c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10), histone_names)

# Enhanced Nature Methods Theme for 1:1 aspect ratio
nature_methods_theme_1x1 <- function() {
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
      legend.position = "right",
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

# Check the range of direct values
cat("Summary of direct Precision values:\n")
print(summary(scatter_direct_pr$median_Precision))
cat("\nSummary of direct Recall values:\n")
print(summary(scatter_direct_pr$median_Recall))

x_range <- range(scatter_direct_pr$median_Precision, na.rm = TRUE)
y_range <- range(scatter_direct_pr$median_Recall, na.rm = TRUE)

cat(sprintf("\nX-axis (Precision) range: %.3f to %.3f\n", x_range[1], x_range[2]))
cat(sprintf("Y-axis (Recall) range: %.3f to %.3f\n", y_range[1], y_range[2]))

# Main plot - Direct Precision vs Recall (1:1 aspect ratio)
scatter_direct_plot <- ggplot(scatter_direct_pr, 
                             aes(x = median_Precision,
                                 y = median_Recall,
                                 color = Method,
                                 shape = Histone)) +
  geom_point(size = 10, alpha = 0.8, stroke = 1) +
  scale_color_manual(values = method_colors, 
                     name = "Peak-Calling Method") +
  scale_shape_manual(values = histone_shape_map,
                     name = "Histone Modification") +
  scale_x_continuous(
    labels = scales::comma_format(accuracy = 0.01),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 0.01),
    limits = c(0, .8),
    breaks = seq(0, .8, by = 0.1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Median Precision vs Median Recall Across Peak Callers",
    x = "median_Precision",
    y = "median_Recall"
  ) +
  nature_methods_theme_1x1()

print(scatter_direct_plot)

# Save the plot with 1:1 aspect ratio
output_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkchipseq_encode_correction_label/result_GT_bulkchipseq_encode_median_corrected"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save with square dimensions (8x8 inches)
ggsave(file.path(output_dir, "HumanPBMC_Precision_vs_Recall_direct_4_3_without_input.pdf"),
       plot = scatter_direct_plot,
       width = 12,
       height = 9,
       device = "pdf",
       bg = "white")

cat("\n✓ Precision vs Recall plot (1:1 aspect ratio) saved to:", 
    file.path(output_dir, "HumanPBMC_Precision_vs_Recall_direct_4_3_without_input.pdf"), "\n")

# =============================================================================
# PART 2: Direct Overlap vs F1_score with 4:3 aspect ratio
# =============================================================================

# Direct Overlap vs F1_score using individual data points
scatter_direct_of  <- combined_data %>% 
  group_by(algorithm, histones) %>%
  summarize(
    median_overlap_pct = median(overlap_pct, na.rm = TRUE),
    sd_overlap_pct = sd(overlap_pct, na.rm = TRUE),
    median_F1_score = median(F1_score, na.rm = TRUE),
    sd_F1_score = sd(F1_score, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(algorithm),
    Histone = factor(histones)
  )

# Enhanced Nature Methods Theme for 4:3 aspect ratio
nature_methods_theme_4x3 <- function() {
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
      legend.position = "right",
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
# Check the range of direct values
cat("\nSummary of direct Overlap Percentage values:\n")
print(summary(scatter_direct_of$median_overlap_pct))
cat("\nSummary of direct F1 Score values:\n")
print(summary(scatter_direct_of$median_F1_score))

x_range_of <- range(scatter_direct_of$median_overlap_pct, na.rm = TRUE)
y_range_of <- range(scatter_direct_of$F1_smedian_F1_scorecore, na.rm = TRUE)

cat(sprintf("\nX-axis (Overlap %%) range: %.1f to %.1f\n", x_range_of[1], x_range_of[2]))
cat(sprintf("Y-axis (F1 Score) range: %.3f to %.3f\n", y_range_of[1], y_range_of[2]))

# Plot: Direct Overlap % vs F1 Score with 4:3 aspect ratio
scatter_direct_overlap_f1 <- ggplot(scatter_direct_of, 
                                   aes(x = median_overlap_pct,
                                       y = median_F1_score,
                                       color = Method,
                                       shape = Histone)) +
  geom_point(size = 10, alpha = 0.8, stroke = 1) +
  scale_color_manual(values = method_colors, 
                     name = "Peak-Calling Method") +
  scale_shape_manual(values = histone_shape_map,
                     name = "Histone Modification") +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 0.01),
    limits = c(0, .6),
    breaks = seq(0, .6, by = 0.1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Median Overlap Percentage vs Median F1 Score",
    x = "Median Overlap Percentage (%)",
    y = "Median F1 Score"
  ) +
  nature_methods_theme_4x3()

print(scatter_direct_overlap_f1)

# Save with 4:3 dimensions (12x9 inches)
ggsave(file.path(output_dir, "HumanPBMC_Overlap_vs_F1_direct_4_3_without_input.pdf"),
       plot = scatter_direct_overlap_f1,
       width = 12,
       height = 9,  # 12:9 = 4:3 ratio
       device = "pdf",
       bg = "white")

cat("\n✓ Overlap vs F1 plot (4:3 aspect ratio) saved to:", 
    file.path(output_dir, "HumanPBMC_Overlap_vs_F1_direct_4_3_without_input.pdf"), "\n")


