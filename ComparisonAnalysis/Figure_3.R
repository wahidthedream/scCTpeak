#########################################################################################################
## Figure 3 A, B
## Simpson Index / Jaccard Index for all HMs of HumanPBMC (with & without input)
## Author: Wahid (modified)
#########################################################################################################

## Required libraries
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(extrafont)
library(patchwork)
library(grid)        # for unit()
library(purrr)       # for map functions

loadfonts(device = "postscript")

#########################################################################################################
## Universal SEACR Reader (Works For ALL SEACR Versions)
#########################################################################################################
read_seacr_bed <- function(file) {
  df <- read.table(file, header = FALSE, comment.char = "", stringsAsFactors = FALSE)
  
  # Keep only first 3 columns if more exist
  df <- df[, 1:min(3, ncol(df))]
  colnames(df) <- c("chr", "start", "end")[1:ncol(df)]
  
  df <- df[complete.cases(df), ]
  
  # Convert numeric and ensure valid regions
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  df <- df[df$start < df$end, ]
  
  if (nrow(df) == 0) stop("SEACR file has no valid peaks: ", file)
  
  return(df)
}

#########################################################################################################
## Standard BED reader
#########################################################################################################
read_bed <- function(file) {
  df <- read.table(file, header = FALSE, comment.char = "", stringsAsFactors = FALSE)
  df <- df[complete.cases(df), ]
  df <- df[, 1:3]
  colnames(df) <- c("chr", "start", "end")
  
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  df <- df[df$start < df$end, ]
  return(df)
}

#########################################################################################################
## Combine multiple BED files into one GRanges object
#########################################################################################################
combine_bed_files <- function(file_list, method = NULL) {
  gr_list <- lapply(file_list, function(f) {
    if (file.size(f) == 0) return(NULL)
    
    if (!is.null(method) && method == "SEACR") {
      df <- read_seacr_bed(f)
    } else {
      df <- read_bed(f)
    }
    
    makeGRangesFromDataFrame(df, seqnames.field = "chr", start.field = "start", end.field = "end")
  })
  
  gr_list <- gr_list[!sapply(gr_list, is.null)]
  
  if (length(gr_list) == 0) return(NULL)
  
  # Combine all GRanges
  combined_gr <- Reduce(c, gr_list)
  return(combined_gr)
}

#########################################################################################################
## Similarity Matrix Calculation
#########################################################################################################
calculate_similarity_matrices <- function(gr_list) {
  cell_types <- names(gr_list)
  n <- length(cell_types)
  
  jaccard_mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(cell_types, cell_types))
  simpson_mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(cell_types, cell_types))
  
  for (i in 1:n) {
    for (j in i:n) {
      gr1 <- gr_list[[i]]
      gr2 <- gr_list[[j]]
      
      if (is.null(gr1) || is.null(gr2) || length(gr1) == 0 || length(gr2) == 0) {
        jaccard_mat[i, j] <- jaccard_mat[j, i] <- 0
        simpson_mat[i, j] <- simpson_mat[j, i] <- 0
        next
      }
      
      overlaps <- findOverlaps(gr1, gr2)
      unique_gr1_overlaps <- length(unique(queryHits(overlaps)))
      unique_gr2_overlaps <- length(unique(subjectHits(overlaps)))
      
      union_size <- length(unique(c(gr1, gr2)))
      intersection_size <- length(unique(c(unique_gr1_overlaps, unique_gr2_overlaps)))
      
      jaccard <- ifelse(union_size > 0, intersection_size / union_size, 0)
      jaccard_mat[i, j] <- jaccard_mat[j, i] <- jaccard
      
      min_size <- min(length(gr1), length(gr2))
      if (min_size == 0) {
        simpson <- 0
      } else {
        if (length(gr1) <= length(gr2)) {
          simpson <- unique_gr1_overlaps / min_size
        } else {
          simpson <- unique_gr2_overlaps / min_size
        }
      }
      
      simpson_mat[i, j] <- simpson_mat[j, i] <- min(simpson, 1)
    }
  }
  
  diag(jaccard_mat) <- 1
  diag(simpson_mat) <- 1
  return(list(jaccard = jaccard_mat, simpson = simpson_mat))
}

#########################################################################################################
## Nature Methods Style Theme
#########################################################################################################
nature_methods_theme <- function() {
  theme_minimal(base_size = 10) +
    theme(
      text = element_text(family = "Helvetica", color = "black", size = 10),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, 
                                margin = margin(b = 10), color = "#2C3E50"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#566573",
                                  margin = margin(b = 15)),
      axis.title = element_text(face = "bold", size = 10, color = "#2C3E50"),
      axis.title.x = element_text(margin = margin(t = 10), vjust = -0.5),
      axis.title.y = element_text(margin = margin(r = 10), vjust = 2),
      axis.text = element_text(size = 9, color = "#34495E"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.line = element_line(color = "#2C3E50", linewidth = 0.5),
      axis.ticks = element_line(color = "#2C3E50", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      legend.title = element_text(face = "bold", size = 9, color = "#2C3E50",
                                 margin = margin(b = 5)),
      legend.text = element_text(size = 8, color = "#34495E"),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(fill = "white", color = NA),
      plot.caption = element_text(size = 8, color = "#7F8C8D", 
                                 margin = margin(t = 10), hjust = 1)
    )
}

#########################################################################################################
## Heatmap function (Nature Methods Style - NO text labels)
#########################################################################################################
compact_theme <- function() {
  theme_minimal(base_size = 9) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
      axis.text.y = element_text(size = 8, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10, 
                               margin = margin(b = 10)),
      legend.position = "right",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.width = unit(0.6, "cm"),
      legend.key.height = unit(1.2, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Nature Methods Style Sequential Color Palette
nature_sequential_palette <- colorRampPalette(c("#FDE725", "#21918C", "#440154"))(100)

create_nature_heatmap <- function(sim_mat, method, metric, target) {
  sim_mat[sim_mat > 1] <- 1
  
  plot_data <- sim_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Row") %>%
    pivot_longer(-Row, names_to = "Column", values_to = "Similarity") %>%
    mutate(
      Row = factor(Row, levels = rownames(sim_mat)),
      Column = factor(Column, levels = colnames(sim_mat))
    )
  
  # Create the heatmap WITHOUT text labels
  p <- ggplot(plot_data, aes(x = Column, y = Row, fill = Similarity)) +
    geom_tile(color = "white", linewidth = 0.15) +
    # NO TEXT LABELS - Clean Nature Methods style
    
    # Nature Methods style color scale
    scale_fill_gradientn(
      colors = nature_sequential_palette,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      na.value = "grey90",
      guide = guide_colorbar(
        title = paste(metric, "Index"),
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(0.6, "cm"),
        barheight = unit(3, "cm"),
        frame.colour = "black",
        ticks.colour = "black",
        label.position = "right"
      )
    ) +
    labs(
      title = paste(target, "-", method),
      subtitle = paste(metric, "Similarity Matrix"),
      x = "",
      y = ""
    ) +
    compact_theme() +
    coord_fixed()
  
  return(p)
}

#########################################################################################################
## Main analysis function (UPDATED with Nature Methods style) - FIXED VERSION
#########################################################################################################
analyze_peaks <- function(target, analysis_type = "with_control") {
  
  if (analysis_type == "with_control") {
    input_dir  <- paste0("~/HumanPBMC_peakbed/all_with_input_peakbed_corrected/jaccard_simpson_peakbed/", target, "_peakbed/")
    output_dir <- paste0("~/result_jaccard_simpson_additional/", target, "_Jaccard_and_Simpson_with_input/")
  } else {
    input_dir  <- paste0("~/HumanPBMC_peakbed/all_without_input_peakbed_corrected/jaccard_simpson_peakbed/", target, "_peakbed/")
    output_dir <- paste0("~/result_jaccard_simpson_additional/", target, "_Jaccard_and_Simpson_without_input/")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(input_dir)
  
  methods <- c("DROMPAplus", "Genrich", "GoPeaks",  "HOMER", "MACS2", "SEACR" ,"SICER2")
  
  all_jaccard <- list()
  all_simpson <- list()
  
  for (method in methods) {
    message("Processing: ", method)
    
    # Accept method-specific filename patterns (SEACR already uses same naming in your data)
    bed_files <- list.files(pattern = paste0(method, "_", target, "_.*\\.bed$"), full.names = TRUE)
    
    if (length(bed_files) == 0) next
    
    gr_list <- list()
    file_names <- character(0)
    
    for (f in bed_files) {
      if (file.size(f) == 0) next
      
      tryCatch({
        if (method == "SEACR") {
          df <- read_seacr_bed(f)
        } else {
          df <- read_bed(f)
        }
        
        if (nrow(df) == 0) next
        
        gr <- makeGRangesFromDataFrame(df, seqnames.field = "chr", start.field = "start", end.field = "end")
        
        # Create a simple name for the file
        sample_name <- gsub(paste0("^.*", method, "_", target, "_|\\.bed$"), "", basename(f))
        sample_name <- gsub("_$", "", sample_name)  # Remove trailing underscores
        sample_name <- ifelse(nchar(sample_name) == 0, basename(f), sample_name)
        
        gr_list[[sample_name]] <- gr
        file_names <- c(file_names, sample_name)
        
      }, error = function(e) {
        message("  Warning: Could not process file ", basename(f), ": ", e$message)
      })
    }
    
    if (length(gr_list) < 2) {
      message("  Skipping ", method, ": Insufficient valid files (", length(gr_list), " files)")
      next
    }
    
    # Ensure all GRanges have names
    if (length(gr_list) != length(names(gr_list))) {
      names(gr_list) <- paste0("Sample_", seq_along(gr_list))
    }
    
    tryCatch({
      sim <- calculate_similarity_matrices(gr_list)
      
      # Save matrices (CSV files contain all numbers)
      write.csv(sim$jaccard, file.path(output_dir, paste0("Jaccard_Matrix_", method, ".csv")))
      write.csv(sim$simpson, file.path(output_dir, paste0("Simpson_Matrix_", method, ".csv")))
      
      # Create heatmaps WITHOUT text labels
      all_jaccard[[method]] <- create_nature_heatmap(sim$jaccard, method, "Jaccard", target)
      all_simpson[[method]] <- create_nature_heatmap(sim$simpson, method, "Simpson", target)
      
      message("  Completed ", method, ": ", length(gr_list), " samples")
      
    }, error = function(e) {
      message("  Error calculating similarities for ", method, ": ", e$message)
    })
  }
  
  if (length(all_jaccard) > 0) {
    tryCatch({
      # Combine plots
      combined_jaccard <- wrap_plots(all_jaccard, ncol = 3) + 
        plot_layout(guides = "collect") &
        theme(legend.position = "right")
      
      combined_simpson <- wrap_plots(all_simpson, ncol = 3) + 
        plot_layout(guides = "collect") &
        theme(legend.position = "right")
      
      # Save combined plots
      ggsave(file.path(output_dir, paste0(target, "_All_Methods_Jaccard_Heatmaps.pdf")),
             combined_jaccard,
             width = 297, height = 210, units = "mm", device = cairo_pdf, dpi = 1200)
      
      ggsave(file.path(output_dir, paste0(target, "_All_Methods_Simpson_Heatmaps.pdf")),
             combined_simpson,
             width = 297, height = 210, units = "mm", device = cairo_pdf, dpi = 1200)
      
      message("  Saved combined plots for ", target)
      
    }, error = function(e) {
      message("  Error creating combined plots for ", target, ": ", e$message)
    })
  } else {
    message("  No valid data for ", target)
  }
}

#########################################################################################################
## Also fix the combine_bed_files function to be more robust
#########################################################################################################
combine_bed_files <- function(file_list, method = NULL) {
  gr_list <- lapply(seq_along(file_list), function(i) {
    f <- file_list[i]
    
    if (!file.exists(f) || file.size(f) == 0) {
      message("  Warning: File ", basename(f), " does not exist or is empty")
      return(NULL)
    }
    
    tryCatch({
      if (!is.null(method) && method == "SEACR") {
        df <- read_seacr_bed(f)
      } else {
        df <- read_bed(f)
      }
      
      if (nrow(df) == 0) {
        message("  Warning: File ", basename(f), " has no valid peaks")
        return(NULL)
      }
      
      gr <- makeGRangesFromDataFrame(df, seqnames.field = "chr", start.field = "start", end.field = "end")
      return(gr)
      
    }, error = function(e) {
      message("  Warning: Could not process file ", basename(f), ": ", e$message)
      return(NULL)
    })
  })
  
  # Remove NULL elements
  valid_gr <- gr_list[!sapply(gr_list, is.null)]
  
  if (length(valid_gr) == 0) {
    message("  No valid GRanges objects created")
    return(NULL)
  }
  
  # Combine all GRanges
  tryCatch({
    combined_gr <- Reduce(c, valid_gr)
    return(combined_gr)
  }, error = function(e) {
    message("  Warning: Could not combine GRanges: ", e$message)
    return(NULL)
  })
}
#########################################################################################################
## Run ALL analyses
#########################################################################################################
run_all_analyses <- function() {
  
  targets <- c("H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
  methods <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
  
  analysis_types <- c("with_control", "without_control")
  
  for (analysis_type in analysis_types) {
    message(paste0("\n", strrep("=", 80)))
    message("RUNNING ANALYSIS FOR: ", toupper(analysis_type))
    message(strrep("=", 80))
    
    # 1. Original analysis (within each histone mark)
    message("\n1. Running original analysis (within each histone mark)...")
    for (t in targets) {
      message("  Processing: ", t)
      analyze_peaks(t, analysis_type)
    }
    
    # 2. NEW: Among tools analysis
    message("\n2. Running NEW analysis: Among Tools...")
    for (histone in targets) {
      message("  Processing: ", histone)
      analyze_among_tools(histone, analysis_type)
    }
    
    # 3. NEW: Among histones analysis
    message("\n3. Running NEW analysis: Among Histones...")
    for (method in methods) {
      message("  Processing: ", method)
      analyze_among_histones(method, analysis_type)
    }
    
    # 4. Create summary comparisons
    message("\n4. Creating summary comparisons...")
    create_summary_comparisons(analysis_type)
    
    message(paste0("\n", strrep("-", 80)))
    message("COMPLETED: ", toupper(analysis_type), " analysis")
    message(strrep("-", 80))
  }
  
  message(paste0("\n", strrep("=", 80)))
  message("ALL ANALYSES COMPLETED SUCCESSFULLY!")
  message(strrep("=", 80))
  
  # Print summary of outputs
  message("\n📊 OUTPUT SUMMARY:")
  message("├── CSV Files: All similarity matrices with numbers")
  message("├── PDF Files: Nature Methods style heatmaps (NO text labels)")
  message("├── Summary Reports: Overall statistics and comparisons")
  message("└── RData Files: Combined results for easy loading")
}

#########################################################################################################
## EXECUTE THE FULL ANALYSIS
#########################################################################################################

# Run everything with one command
run_all_analyses()




#########################################################################################################
##
## Figure 3 C, D

#########################################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# =============================================================================
# CONFIGURATION
# =============================================================================

base_dir <- "~/result_jaccard_simpson_additional"
analysis_type <- "with_input"

histones <- c("H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K36me3", "H3K4me3", "Olig2", "Rad21")
methods <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")

# Clean color palette
nature_palette <- c(
  "#E64B35",  # Red (DROMPAplus)
  "#4DBBD5",  # Cyan (Genrich)
  "#00A087",  # Teal (GoPeaks)
  "#3C5488",  # Blue (HOMER)
  "#F39B7F",  # Orange (MACS2)
  "#8491B4",  # Lavender (SEACR)
  "#91D1C2"   # Mint (SICER2)
)

names(nature_palette) <- methods

# =============================================================================
# CLEAN NATURE THEME
# =============================================================================

theme_nature_clean <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      axis.title = element_text(face = "plain"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "plain"),
      legend.key = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

# =============================================================================
# DATA EXTRACTION
# =============================================================================

extract_simpson_values <- function() {
  all_data <- data.frame()
  
  for (histone in histones) {
    sim_dir <- file.path(base_dir, paste0(histone, "_Jaccard_and_Simpson_", analysis_type))
    
    if (!dir.exists(sim_dir)) next
    
    for (method in methods) {
      simpson_file <- file.path(sim_dir, paste0("Simpson_Matrix_", method, ".csv"))
      
      if (file.exists(simpson_file)) {
        tryCatch({
          simpson_mat <- read.csv(simpson_file, row.names = 1, check.names = FALSE)
          
          if (nrow(simpson_mat) > 1 && ncol(simpson_mat) > 1 && 
              nrow(simpson_mat) == ncol(simpson_mat)) {
            
            n <- nrow(simpson_mat)
            simpson_values <- numeric()
            
            for (i in 1:(n-1)) {
              for (j in (i+1):n) {
                val <- as.numeric(simpson_mat[i, j])
                if (!is.na(val) && is.numeric(val) && val >= 0 && val <= 1) {
                  simpson_values <- c(simpson_values, val)
                }
              }
            }
            
            if (length(simpson_values) > 0) {
              all_data <- rbind(all_data, data.frame(
                Histone = histone,
                Method = method,
                Simpson_Index = simpson_values,
                stringsAsFactors = FALSE
              ))
            }
          }
        }, error = function(e) NULL)
      }
    }
  }
  
  if (nrow(all_data) == 0) {
    stop("No data extracted.")
  }
  
  all_data$Histone <- factor(all_data$Histone, levels = histones)
  all_data$Method <- factor(all_data$Method, levels = methods)
  
  cat("Data extracted:", nrow(all_data), "points\n")
  return(all_data)
}

# =============================================================================
# FIGURE 1: BOX PLOT - CLEAN VERSION 
# =============================================================================

create_figure1_boxplot <- function(simpson_data) {
  cat("Creating Figure 1...\n")
  
  # Order methods by mean performance
  method_ranking <- simpson_data %>%
    group_by(Method) %>%
    summarise(Mean = mean(Simpson_Index), .groups = "drop") %>%
    arrange(desc(Mean))
  
  method_order <- as.character(method_ranking$Method)
  simpson_data$Method <- factor(simpson_data$Method, levels = method_order)
  plot_colors <- nature_palette[method_order]
  
  # Create box plot - CLEAN VERSION
  p <- ggplot(simpson_data, aes(x = Histone, y = Simpson_Index, fill = Method)) +
    
    # Box plot with minimal outliers
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.7,
      outlier.size = 0.3,      # Very small outliers
      outlier.shape = 16,      # Simple filled circle
      outlier.color = "gray50", # Gray color (not white)
      outlier.alpha = 0.2,     # Very transparent
      lwd = 0.3,
      alpha = 0.9
    ) +
    
    # NO stat_summary - removes white balls completely
    # If you want mean indicators, use a different approach:
    
    # OPTION A: Add small black mean lines (recommended)
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 4,              # X shape instead of circle
      size = 1.0,
      color = "black",        # Black color
      position = position_dodge(width = 0.8),
      show.legend = FALSE
    ) +
    
    # Color scheme
    scale_fill_manual(values = plot_colors, name = "Peak caller") +
    
    # Y-axis
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    # Labels
    labs(
      x = "Histone modification",
      y = "Simpson similarity index"
    ) +
    
    # Theme
    theme_nature_clean() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    )
  
  return(p)
}

# =============================================================================
# ALTERNATIVE: BOX PLOT WITHOUT ANY MEAN POINTS (CLEANEST)
# =============================================================================

create_figure1_boxplot_clean <- function(simpson_data) {
  cat("Creating Figure 1 (clean version)...\n")
  
  # Order methods by mean performance
  method_ranking <- simpson_data %>%
    group_by(Method) %>%
    summarise(Mean = mean(Simpson_Index), .groups = "drop") %>%
    arrange(desc(Mean))
  
  method_order <- as.character(method_ranking$Method)
  simpson_data$Method <- factor(simpson_data$Method, levels = method_order)
  plot_colors <- nature_palette[method_order]
  
  # SIMPLEST VERSION - just box plots, no extra points
  p <- ggplot(simpson_data, aes(x = Histone, y = Simpson_Index, fill = Method)) +
    
    # Clean box plots only
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.7,
      outlier.size = 0.3,
      outlier.shape = 16,
      outlier.color = "black",
      outlier.alpha = 0.3,
      lwd = 0.3,
      alpha = 0.9
    ) +
    
    # Color scheme
    scale_fill_manual(values = plot_colors, name = "Peak caller") +
    
    # Y-axis
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    # Labels
    labs(
      x = "Histone modification",
      y = "Simpson similarity index"
    ) +
    
    # Theme
    theme_nature_clean() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
  
  return(p)
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

run_analysis <- function() {
  cat("Starting analysis...\n")
  
  # Create output directory
  output_dir <- file.path(base_dir, "figures_final")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Extract data
  simpson_data <- extract_simpson_values()
  
  # Create figures (use clean version without white balls)
  fig1 <- create_figure1_boxplot_clean(simpson_data)  # Clean version


  # Save figures
  ggsave(file.path(output_dir, "Figure1_boxplot.pdf"), fig1, 
         width = 180/25.4, height = 120/25.4, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "Figure1_boxplot.png"), fig1,
         width = 180/25.4, height = 120/25.4, dpi = 300, bg = "white")
  
  
  # Save data
  write.csv(simpson_data, file.path(output_dir, "simpson_data.csv"), row.names = FALSE)
  
  cat("Analysis complete. Files saved in:", output_dir, "\n")
  return(list(data = simpson_data, fig1 = fig1, fig2 = fig2))
}

# Run analysis
results <- tryCatch({
  run_analysis()
}, error = function(e) {
  cat("Error:", e$message, "\n")
  stop()
})

cat("Done! No white balls in figures.\n")

# Display the figure
print(results$fig1)
