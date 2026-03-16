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

##====================================================================================================##
## Universal SEACR Reader (Works For ALL SEACR Versions)
##====================================================================================================##
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

##====================================================================================================##
## Standard BED reader
##====================================================================================================##
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

##====================================================================================================##
## Similarity Matrix Calculation
##====================================================================================================##
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

##====================================================================================================##
## Nature Methods Style Theme
##====================================================================================================##
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

##====================================================================================================##
## Heatmap function (Nature Methods Style - NO text labels)
##====================================================================================================##
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

##====================================================================================================##
## Main analysis function (UPDATED with Nature Methods style) - FIXED VERSION
##====================================================================================================##
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

##====================================================================================================##
## Also fix the combine_bed_files function to be more robust
##====================================================================================================##
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
##====================================================================================================##
## Run ALL analyses
##====================================================================================================##
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

##====================================================================================================##
## EXECUTE THE FULL ANALYSIS
##====================================================================================================##

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


#########################################################################################################
##
## Figure 3 E
#########################################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(ggpubr)
library(rstatix)
library(Signac)
library(GenomicRanges)

# Resolve package conflicts
library(conflicted)
conflicts_prefer(base::intersect)
conflicts_prefer(base::union)
conflicts_prefer(base::setdiff)

# Load data
scRNA_seq <- readRDS("~/scRNA_seq/10x_scRNAseq/pbmc_processed.rds")
#scRNA_seq <- readRDS("~/scRNA_seq/10x_scRNAseq/pbmc_processed_with_TPM.rds")
# Check cell type distribution
cat("Cell type distribution:\n")
print(table(scRNA_seq$predicted_celltype))

# Check all available metadata columns
cat("\nAvailable metadata columns:\n")
print(colnames(scRNA_seq@meta.data))

# Check available reductions
cat("\nAvailable reductions:\n")
print(names(scRNA_seq@reductions))

# Option 1: Basic UMAP with default colors - FIXED
p1 <- DimPlot(scRNA_seq, 
              reduction = "umap",
              group.by = "predicted_celltype",
              label = TRUE,
              label.size = 5,           # Increased label size
              repel = TRUE,
              pt.size = .8) +
  ggtitle("Human PBMC UMAP Plot-label by Color") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)



# Save as PDF - Option 2: High quality with specific settings
ggsave("~/scRNA_seq/10x_scRNAseq/HumanPBMC_umap_plot_labeled_by_color.pdf", plot = p1, width = 8, height = 7, dpi = 600, device = cairo_pdf, bg = "white")

cat("\nPDF saved successfully!\n")

#########################################################################################################
## Figure 3 F G
######################################################################################################### 

# Now load other packages
library(Seurat)
library(Matrix)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(Seurat)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(ggrepel)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(purrr)
library(tibble)
library(ggpubr)       
library(viridis)     
library(ggsci)       
library(scales)    
library(conflicted)  

# Resolve namespace conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("setdiff", "dplyr")
# Explicitly prefer dplyr functions over plyr
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::rename)

histones <- c( "H3K27ac-b", "H3K27ac-s",  "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
##====================================================================================================##
# loading the datasets
##====================================================================================================##
# 1. Load the pre-aggregated TPM expression matrix
for (hist in histones){
TPM <- read.delim("/home/wahid/project_scHMTF/GSE195725_processed_data/scRNA_seq/10x_scRNAseq/pbmc_celltype_averaged_TPM.txt.gz", row.names = 1)

# 3. Define parameters
celltypes <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")
peak_callers <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
#sc_peak_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/result/l1_withcontrol/H3K27me3_peakbed"
sc_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_with_input_peakbed_corrected/scrnaseq_peakbed/", hist, "_peakbed")


# 4. Create output directory
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final/result_GT_scrnaseq_with_input/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)




##====================================================================================================##
# Comprehensive BED File Reading Function for each peak caller - FIXED SICER2
##====================================================================================================##
read_peaks <- function(bed_file) {
  tryCatch({
    # Read the BED file
    bed_data <- read.table(bed_file, sep = "\t", header = FALSE, 
                          stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
    
    n_cols <- ncol(bed_data)
    filename <- basename(bed_file)
    peak_caller <- gsub("_.*", "", filename)
    
    cat("  File:", filename, "Columns:", n_cols, "Caller:", peak_caller, "\n")
    
    # Create basic GRanges object (BED is 0-based, GRanges is 1-based)
    peaks <- GRanges(
      seqnames = bed_data[, 1],
      ranges = IRanges(start = bed_data[, 2] + 1, end = bed_data[, 3])
    )
    
    # Handle each peak caller format specifically
    if (peak_caller == "Genrich") {
      # Genrich format: chr, start, end, name, score, strand, signalValue, pValue, qValue, peak
      if (n_cols >= 10) {
        peaks$score <- as.numeric(bed_data[, 7]) # Use signalValue
        mcols(peaks)$name <- bed_data[, 4]
        mcols(peaks)$strand <- bed_data[, 6]
        mcols(peaks)$signalValue <- as.numeric(bed_data[, 7])
        mcols(peaks)$pValue <- as.numeric(bed_data[, 8])
        mcols(peaks)$qValue <- as.numeric(bed_data[, 9])
        mcols(peaks)$peak <- as.numeric(bed_data[, 10])
      }
    } 
    else if (peak_caller == "MACS2") {
      # MACS2 format: chr, start, end, name, score, strand, fold_change, -log10(pvalue), -log10(qvalue), summit
      if (n_cols >= 8) {
        peaks$score <- as.numeric(bed_data[, 7]) # Use fold change
        mcols(peaks)$name <- bed_data[, 4]
        mcols(peaks)$strand <- bed_data[, 6]
        mcols(peaks)$fold_change <- as.numeric(bed_data[, 7])
        mcols(peaks)$neg_log10_pvalue <- as.numeric(bed_data[, 8])
#        mcols(peaks)$neg_log10_qvalue <- as.numeric(bed_data[, 9])
#        mcols(peaks)$summit <- as.numeric(bed_data[, 10])
      }
    }
    else if (peak_caller == "HOMER") {
      # HOMER format: chr, start, end, name, score, strand
      if (n_cols >= 6) {
        peaks$score <- as.numeric(bed_data[, 5]) # Use score column
        mcols(peaks)$name <- bed_data[, 4]
        mcols(peaks)$strand <- bed_data[, 6]
      } else {
        peaks$score <- 1
      }
    }
    else if (peak_caller == "SEACR") {
      # SEACR format: chr, start, end, total_signal, max_signal, summit
      if (n_cols >= 6) {
        peaks$score <- as.numeric(bed_data[, 4]) # Use total_signal
        mcols(peaks)$total_signal <- as.numeric(bed_data[, 4])
        mcols(peaks)$max_signal <- as.numeric(bed_data[, 5])
        mcols(peaks)$summit <- bed_data[, 6]
      } else {
        peaks$score <- 1
      }
    }
    else if (peak_caller == "SICER2") {
      # SICER2 format: chr, start, end, name, (empty), strand - FIXED
      if (n_cols >= 6) {
        peaks$score <- 1 # Default score for SICER2
        mcols(peaks)$name <- bed_data[, 4]
        mcols(peaks)$strand <- bed_data[, 6]  # Strand is in column 6
      } else if (n_cols >= 4) {
        peaks$score <- 1
        mcols(peaks)$name <- bed_data[, 4]
      } else {
        peaks$score <- 1
      }
    }
    else if (peak_caller == "GoPeaks") {
      # GoPeaks format: chr, start, end (only 3 columns)
      peaks$score <- 1 # Default score
    }
    else if (peak_caller == "DROMPAplus") {
      # DROMPAplus format: unknown, use generic approach
      if (n_cols >= 5) {
        peaks$score <- as.numeric(bed_data[, 5])
      } else if (n_cols >= 4) {
        peaks$score <- 1
        mcols(peaks)$name <- bed_data[, 4]
      } else {
        peaks$score <- 1
      }
    }
    else {
      # Default handling for unknown callers
      if (n_cols >= 5) {
        peaks$score <- as.numeric(bed_data[, 5])
      } else if (n_cols >= 4) {
        peaks$score <- 1
        mcols(peaks)$name <- bed_data[, 4]
      } else {
        peaks$score <- 1
      }
    }
    
    # Remove any peaks with invalid scores
    peaks <- peaks[!is.na(peaks$score) & is.finite(peaks$score)]
    
    return(peaks)
    
  }, error = function(e) {
    cat("  ERROR reading", basename(bed_file), ":", conditionMessage(e), "\n")
    
    # Fallback: try basic import
    tryCatch({
      peaks <- import(bed_file, format = "BED")
      if (length(mcols(peaks)) == 0) {
        peaks$score <- 1
      } else if ("score" %in% colnames(mcols(peaks))) {
        peaks$score <- mcols(peaks)$score
      } else {
        peaks$score <- 1
      }
      peaks <- peaks[!is.na(peaks$score) & is.finite(peaks$score)]
      return(peaks)
    }, error = function(e2) {
      cat("  Fallback also failed for", basename(bed_file), "\n")
      return(NULL)
    })
  })
}
#
##====================================================================================================##
# FIXED: Process All Peak Files with Proper RP Model Implementation
# MODIFIED: Ensure same number of gene activities for all peak callers
##====================================================================================================##
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(dplyr)
library(tidyr)

edb <- EnsDb.Hsapiens.v86
genes_gr <- genes(edb)
#seqlevelsStyle(genes_gr) <- "UCSC"
harmonize_seqlevels <- function(gr) {
  sl <- seqlevels(gr)
  
  if (any(grepl("^chr", sl))) {
    return(gr)   #
  } else {
    seqlevels(gr) <- paste0("chr", sl)
    return(gr)
  }
}
genes_gr <- harmonize_seqlevels(genes_gr)
# Select protein-coding genes only
cat("Total genes in EnsDb:", length(genes_gr), "\n")
protein_coding_genes_gr <- genes_gr[genes_gr$gene_biotype == "protein_coding"]
cat("Protein-coding genes:", length(protein_coding_genes_gr), "\n")

# Show distribution of gene biotypes
cat("\nGene biotype distribution:\n")
biotype_counts <- table(genes_gr$gene_biotype)
print(sort(biotype_counts, decreasing = TRUE)[1:10])  # Top 10 biotypes

# Use only protein-coding genes
genes_gr <- protein_coding_genes_gr

# CORRECTION: Create single-base TSS positions first, then define promoter regions
# Step 1: Extract single-base TSS from each gene
tss_gr <- genes_gr

# Get correct TSS positions based on strand for PROMOTER REGION CREATION
tss_positions_for_promoters <- ifelse(strand(genes_gr) == "+", 
                                     start(genes_gr),  # TSS at start for + strand
                                     end(genes_gr))    # TSS at end for - strand

# Set both start and end to TSS position (single-base range)
start(tss_gr) <- tss_positions_for_promoters
end(tss_gr) <- tss_positions_for_promoters

# Step 2: Define the 100KB window around TSS (±100KB)
promoter_region <- promoters(tss_gr, upstream = 100000, downstream = 100000)

# Function to calculate Regulatory Potential (RP) weight
# weight = 2^(-z/10000) where z is absolute distance from TSS
calculate_rp_weight <- function(distance) {
  2^(-abs(distance) / 10000)  # Your specified formula
}

# Function to process peaks for one peak caller with RP model
process_peak_caller <- function(peak_caller) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PROCESSING PEAK CALLER:", peak_caller, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  peak_activity_list <- list()
  processed_celltypes <- 0
  
  # Create a complete template with ALL genes and ALL celltypes
  all_genes_template <- data.frame(
    gene_id = genes_gr$gene_id,
    gene_name = genes_gr$gene_name
  ) %>%
    tidyr::crossing(celltype = celltypes) %>%
    mutate(
      peak_caller = peak_caller,
      peak_count = 0,
      weighted_peak_count = 0,
      mean_distance = 0,
      min_absolute_distance = 0,
      max_weight = 0,
      n_upstream = 0,
      n_downstream = 0
    )
  
  for (ct in celltypes) {
    bed_file <- file.path(sc_peak_dir, paste0(peak_caller, "_", hist, "_", ct, ".bed"))
    
    # Initialize with zeros for all genes for this celltype
    celltype_template <- all_genes_template %>%
      filter(celltype == ct)
    
    if (file.exists(bed_file)) {
      cat("Processing", ct, "...")

      tryCatch({
        # Read peaks
        peaks_gr <- read_peaks(bed_file)
        
        if (is.null(peaks_gr) || length(peaks_gr) == 0) {
          cat(" No peaks found, using zeros for all genes\n")
          peak_activity_list[[ct]] <- celltype_template
          processed_celltypes <- processed_celltypes + 1
          next
        }
        
        seqlevelsStyle(peaks_gr) <- "UCSC"
        
        # Find overlaps between peaks and the 100KB promoter regions
        overlaps <- findOverlaps(peaks_gr, promoter_region)
        
        if (length(overlaps) > 0) {
          # Get matching peaks and genes
          peak_hits <- peaks_gr[queryHits(overlaps)]
          gene_hits <- promoter_region[subjectHits(overlaps)]
          
          # IMPORTANT: Get the original gene indices to calculate correct TSS
          gene_indices <- subjectHits(overlaps)
          original_genes <- genes_gr[gene_indices]
          
          
          # STEP 1: Calculate distance from TSS for each peak-gene pair
          # Get actual TSS positions from original genes (considering gene strand)
          # NOTE: Different variable name to avoid conflict
          gene_tss_positions <- ifelse(strand(original_genes) == "+", 
                                      start(original_genes),  # TSS for positive strand
                                      end(original_genes))    # TSS for negative strand
          
          # Calculate peak centers
          peak_centers <- (start(peak_hits) + end(peak_hits)) / 2
          
          # CORRECTION: Calculate distances (strand-aware)
          # For + strand: upstream = negative, downstream = positive
          # For - strand: upstream = positive, downstream = negative
          distances <- ifelse(strand(original_genes) == "+",
                             peak_centers - gene_tss_positions,  # Positive strand
                             gene_tss_positions - peak_centers)  # Negative strand
          
          # STEP 2: Calculate RP weights using your formula
          rp_weights <- calculate_rp_weight(distances)
          
          # STEP 3: Handle any potential NaN/Inf values
          rp_weights[is.na(rp_weights) | is.infinite(rp_weights)] <- 0
          
          # Create data frame with all information
          gene_activity <- data.frame(
            gene_id = original_genes$gene_id,
            gene_name = original_genes$gene_name,
            celltype = ct,
            peak_weight = rp_weights,
            distance_from_tss = distances,
            strand = as.character(strand(original_genes))
          )
          
          # STEP 5: Summarize both peak count and weighted peak count per gene
          gene_activity_summary <- gene_activity %>%
            group_by(gene_id, gene_name, celltype) %>%
            summarise(
              peak_count = n(),                    # Number of peaks in 100KB window
              weighted_peak_count = sum(peak_weight, na.rm = TRUE),  # Sum of RP weights
              mean_distance = mean(distance_from_tss, na.rm = TRUE), # Average distance from TSS
              min_absolute_distance = min(abs(distance_from_tss), na.rm = TRUE), # Distance of closest peak
              max_weight = max(peak_weight, na.rm = TRUE),       # Highest weight among peaks
              n_upstream = sum(distance_from_tss < 0, na.rm = TRUE), # Peaks upstream of TSS
              n_downstream = sum(distance_from_tss >= 0, na.rm = TRUE), # Peaks downstream of TSS
              .groups = 'drop'
            )
          
          # Replace any remaining NaN values with 0
          gene_activity_summary <- gene_activity_summary %>%
            mutate(across(where(is.numeric), ~ replace(., is.na(.) | is.infinite(.), 0)))
          
          # Merge with template to ensure ALL genes are included
          final_result <- celltype_template %>%
            select(-c(peak_count:n_downstream)) %>%
            left_join(gene_activity_summary, by = c("gene_id", "gene_name", "celltype")) %>%
            mutate(
              across(c(peak_count, weighted_peak_count, mean_distance, 
                       min_absolute_distance, max_weight, n_upstream, n_downstream),
                     ~ replace_na(., 0))
            )
          
          peak_activity_list[[ct]] <- final_result
          processed_celltypes <- processed_celltypes + 1
          cat(" DONE (", nrow(final_result), "genes,", length(overlaps), "peak-gene overlaps)\n")
          
          # Debug output for first few genes with verification
          if (processed_celltypes == 1) {
            cat("Sample RP weights for first 5 peaks:\n")
            sample_data <- head(data.frame(
              Gene = original_genes$gene_name,
              GeneBiotype = original_genes$gene_biotype,
              Strand = as.character(strand(original_genes)),
              TSS = gene_tss_positions,  # Use the corrected variable name
              PeakCenter = round(peak_centers),
              Distance = round(distances),
              RP_Weight = round(rp_weights, 4)
            ), 5)
            print(sample_data)
            
            # Verify negative strand genes
            if (any(strand(original_genes) == "-")) {
              cat("\nNegative strand gene verification (first 2):\n")
              neg_examples <- which(strand(original_genes) == "-")[1:2]
              if (length(neg_examples) > 0) {
                verification <- data.frame(
                  Gene = original_genes$gene_name[neg_examples],
                  GeneBiotype = original_genes$gene_biotype[neg_examples],
                  GeneStart = start(original_genes[neg_examples]),
                  GeneEnd = end(original_genes[neg_examples]),
                  TSS = gene_tss_positions[neg_examples],  # Use the corrected variable name
                  PeakCenter = round(peak_centers[neg_examples]),
                  Distance = round(distances[neg_examples])
                )
                print(verification)
              }
            }
          }
        } else {
          cat(" No peaks near genes, using zeros for all genes\n")
          peak_activity_list[[ct]] <- celltype_template
          processed_celltypes <- processed_celltypes + 1
        }
      }, error = function(e) {
        cat(" ERROR:", conditionMessage(e), "- using zeros for all genes\n")
        peak_activity_list[[ct]] <- celltype_template
        processed_celltypes <- processed_celltypes + 1
      })
      
    } else {
      cat("  File not found - using zeros for all genes\n")
      peak_activity_list[[ct]] <- celltype_template
      processed_celltypes <- processed_celltypes + 1
    }
  }
  
  # Combine all results for this peak caller
  if (length(peak_activity_list) > 0) {
    all_peak_activity <- bind_rows(peak_activity_list)
    cat("\nSUCCESS: Processed", processed_celltypes, "cell types for", peak_caller, "\n")
    cat("Total genes processed:", nrow(all_peak_activity), "\n")
    return(all_peak_activity)
  } else {
    cat("\nUsing template with zeros for all genes\n")
    return(all_genes_template)
  }
}
# Re-process all peak callers with the fixed implementation
all_peak_activity_list <- list()

for (peak_caller in peak_callers) {
  result <- tryCatch({
    process_peak_caller(peak_caller)
  }, error = function(e) {
    cat("FATAL ERROR with", peak_caller, ":", conditionMessage(e), "\n")
    # Return template with zeros if everything fails
    data.frame(
      gene_id = genes_gr$gene_id,
      gene_name = genes_gr$gene_name
    ) %>%
      tidyr::crossing(celltype = celltypes) %>%
      mutate(
        peak_caller = peak_caller,
        peak_count = 0,
        weighted_peak_count = 0,
        mean_distance = 0,
        min_absolute_distance = 0,
        max_weight = 0,
        n_upstream = 0,
        n_downstream = 0
      )
  })
  
  if (!is.null(result)) {
    all_peak_activity_list[[peak_caller]] <- result
  }
  
  # Save intermediate results after each peak caller
  if (length(all_peak_activity_list) > 0) {
    temp_activity <- bind_rows(all_peak_activity_list)
    write.csv(temp_activity, file.path(out_dir, paste0("temp_peak_activity_data_", hist, "_100KB_protein_coding_with_input.csv")), row.names = FALSE)
  }
}

# Combine all results
if (length(all_peak_activity_list) > 0) {
  all_peak_activity <- bind_rows(all_peak_activity_list)
  cat("\n", paste(rep("=", 56), collapse = ""), "\n")
  cat("SUCCESS: Processed", length(all_peak_activity_list), "peak callers\n")
  cat("Total gene-celltype-caller pairs:", nrow(all_peak_activity), "\n")
  cat(paste(rep("=", 56), collapse = ""), "\n")
  
  # Verify that all peak callers have the same number of rows
  row_counts <- all_peak_activity %>%
    group_by(peak_caller) %>%
    summarise(n_genes = n(), .groups = 'drop')
  
  cat("\nGene count per peak caller:\n")
  print(row_counts)
  
  # Check if all have the same number
  if (length(unique(row_counts$n_genes)) == 1) {
    cat("✓ SUCCESS: All peak callers have the same number of gene activities (", 
        unique(row_counts$n_genes), "protein-coding genes)\n")
  } else {
    cat("✗ WARNING: Peak callers have different numbers of gene activities\n")
  }
  
  # Save final results with corrected name
  final_output_file <- file.path(out_dir, paste0("final_peak_counts_with_RP_model_", hist, "_100KB_protein_coding_with_input.csv"))
  write.csv(all_peak_activity, final_output_file, row.names = FALSE)
  cat("Final results saved to:", final_output_file, "\n")
  
  # Print summary statistics
  cat("\nRP Model Summary Statistics (Protein-coding genes only):\n")
  summary_stats <- all_peak_activity %>%
    group_by(peak_caller) %>%
    summarise(
      mean_weighted_count = mean(weighted_peak_count, na.rm = TRUE),
      median_weighted_count = median(weighted_peak_count, na.rm = TRUE),
      mean_peak_count = mean(peak_count, na.rm = TRUE),
      genes_with_peaks = sum(peak_count > 0),
      percent_with_peaks = round(mean(peak_count > 0) * 100, 1),
      .groups = 'drop'
    )
  print(summary_stats)
  
} else {
  stop("No peak callers were successfully processed")
}


##====================================================================================================##
# Prepare and Merge with Expression Data - FIXED with retry logic
##====================================================================================================##

# First, let's examine the structure of your TPM matrix
print("TPM matrix structure:")
print(dim(TPM))
print("First few rows and columns:")
print(TPM[1:5, 1:5])

# Check if gene names are row names
print("Row names (gene names):")
print(head(rownames(TPM)))

# The gene names are ROW NAMES, not a column
# Convert to data frame with gene names as a column
TPM_with_genes <- TPM %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(gene = sub("\\..*", "", gene))  # Clean gene IDs if needed

print("TPM_with_genes structure:")
str(TPM_with_genes[1:5, 1:5])

# Now pivot to long format
TPM_long <- TPM_with_genes %>%
  pivot_longer(cols = -gene,  # All columns except 'gene'
               names_to = "celltype", 
               values_to = "tpm_expression") %>%
  rename(gene_id = gene) %>%
  dplyr::filter(celltype %in% celltypes)  # Filter to relevant cell types

# Let's check what we have
print("TPM_long structure:")
str(TPM_long)
print("First few rows of TPM_long:")
head(TPM_long)

# Check the structure of the fixed peak activity data
print("Fixed all_peak_activity structure:")
str(all_peak_activity[1:5, ])

# Check if celltypes match between datasets
print("Cell types in TPM data:")
print(unique(TPM_long$celltype))
print("Cell types in peak data:")
print(unique(all_peak_activity$celltype))

# Remove version numbers from peak activity gene IDs
all_peak_activity <- all_peak_activity %>%
  mutate(gene_id = sub("\\..*", "", gene_id))

# GENE ID MAPPING SOLUTION (added with retry logic)
print("TPM gene ID examples:")
print(head(unique(TPM_long$gene_id)))
print("Peak activity gene ID examples:")
print(head(unique(all_peak_activity$gene_id)))

# Use biomaRt to map between Ensembl IDs and gene symbols with retry logic
get_gene_mapping_with_retry <- function(gene_ids, max_attempts = 3) {
  attempt <- 1
  while (attempt <= max_attempts) {
    tryCatch({
      cat(sprintf("Attempt %d to connect to Ensembl...\n", attempt))
      
      # Try different mirrors if main server fails
      mirrors <- c("www.ensembl.org", "useast.ensembl.org", "asia.ensembl.org")
      
      for (mirror in mirrors) {
        tryCatch({
          cat(sprintf("  Trying mirror: %s\n", mirror))
          ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = paste0("https://", mirror))
          
          # Get mapping from Ensembl IDs to gene symbols
          gene_mapping <- getBM(
            attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"), filters = "ensembl_gene_id", values = gene_ids, mart = ensembl, useCache = FALSE 
          )
          
          cat(sprintf("  Successfully retrieved %d mappings from %s\n", nrow(gene_mapping), mirror))
          return(gene_mapping)
          
        }, error = function(e) {
          cat(sprintf("  Failed with mirror %s: %s\n", mirror, conditionMessage(e)))
        })
      }
      
      # If all mirrors fail, wait and try again
      Sys.sleep(5)  # Wait 5 seconds before retrying
      attempt <- attempt + 1
      
    }, error = function(e) {
      cat(sprintf("Attempt %d failed: %s\n", attempt, conditionMessage(e)))
      Sys.sleep(5)
      attempt <- attempt + 1
    })
  }
  
  # If all attempts fail, return empty data frame
  cat("All connection attempts failed. Returning empty mapping.\n")
  return(data.frame(ensembl_gene_id = character(),  hgnc_symbol = character(), external_gene_name = character()))
}

# Get unique gene IDs for mapping
unique_gene_ids <- unique(all_peak_activity$gene_id)

# Try to get gene mapping with retry logic
gene_mapping <- get_gene_mapping_with_retry(unique_gene_ids)

# If mapping failed completely, create a simple mapping from available data
if (nrow(gene_mapping) == 0) {
  cat("Creating fallback mapping from available gene names...\n")
  gene_mapping <- all_peak_activity %>%
    distinct(gene_id, gene_name) %>%
    rename(ensembl_gene_id = gene_id, external_gene_name = gene_name) %>%
    mutate(hgnc_symbol = external_gene_name)
}

print("Gene mapping results:")
head(gene_mapping)
print(sprintf("Total mappings: %d", nrow(gene_mapping)))

# Merge with mapping - this is the key fix
merged_with_mapping <- all_peak_activity %>%
  left_join(gene_mapping, by = c("gene_id" = "ensembl_gene_id")) %>%
  # Use external_gene_name if available, otherwise use gene_name from peaks
  mutate(gene_name_to_use = ifelse(!is.na(external_gene_name) & external_gene_name != "", 
                                 external_gene_name, gene_name)) %>%
  left_join(TPM_long, by = c("gene_name_to_use" = "gene_id", "celltype")) %>%
  dplyr::filter(!is.na(tpm_expression) & !is.na(peak_count))

print(paste("Mapped dataset:", nrow(merged_with_mapping), "pairs"))

# Apply log transformations
final_merged_data <- merged_with_mapping %>%
  mutate(
    log_tpm = log2(tpm_expression + 1), log_peak_count = log10(peak_count + 1), log_weighted_count = log10(weighted_peak_count + 1)
  )

print(paste("Final dataset:", nrow(final_merged_data), "gene-celltype-peakcaller pairs"))
print("Final merged data structure:")
str(final_merged_data[1:5, ])


##====================================================================================================##
## COMPLETE CODE: scRNA-seq analysis using top 100 genes from final_merged_data
## FOLLOWING YOUR SPECIFIED WORKFLOW
##====================================================================================================##

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(cowplot)
library(scales)
library(pheatmap)

# Set seed for reproducibility
set.seed(2024)

# Resolve package conflicts
library(conflicted)
conflicts_prefer(base::intersect)
conflicts_prefer(base::union)
conflicts_prefer(base::setdiff)

# Set output directory
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final/result_GT_scrnaseq_with_input"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

##====================================================================================================##
## STEP 1: Select top 100 genes for each cell-type and tool
##====================================================================================================##

cat("=== STEP 1: SELECTING TOP 100 GENES FOR EACH CELL-TYPE AND TOOL ===\n")

# Check if final_merged_data exists
if (!exists("final_merged_data")) {
  stop("Error: final_merged_data not found in R environment!")
}

cat("✓ final_merged_data found! Dimensions:", dim(final_merged_data), "\n")

# Verify required columns
required_cols <- c("peak_caller", "celltype", "gene_name_to_use", "weighted_peak_count")
missing_cols <- required_cols[!required_cols %in% colnames(final_merged_data)]
if (length(missing_cols) > 0) {
  stop("Error: Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Remove duplicates first (keep highest weighted_peak_count for each gene)
top_100_genes <- final_merged_data %>%
  group_by(peak_caller, celltype, gene_name_to_use) %>%
  summarise(
    weighted_peak_count = max(weighted_peak_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Get top 100 genes for each cell-type and tool
  group_by(peak_caller, celltype) %>%
  arrange(desc(weighted_peak_count), .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup()

cat("\n✓ Top 100 genes selected for each combination\n")
print(table(top_100_genes$peak_caller, top_100_genes$celltype))

# Save the top 100 genes
write.csv(top_100_genes, file.path(out_dir, paste0("HumanPBMC_", hist, "_top_100_genes_all_combinations_with_input.csv")), 
          row.names = FALSE)

##====================================================================================================##
## STEP 2: Load scRNA-seq data
##====================================================================================================##

cat("\n=== STEP 2: LOADING SCRNA-SEQ DATA ===\n")

scRNA_seq <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/scRNA_seq/10x_scRNAseq/pbmc_processed_with_TPM.rds")

cat("✓ scRNA-seq data loaded\n")
cat("  Dimensions:", dim(scRNA_seq), "\n")
cat("  Cell types in scRNA-seq:", paste(unique(scRNA_seq$predicted_celltype), collapse = ", "), "\n")

# Use full scRNA-seq data (DO NOT FILTER GENES YET)
scRNA_main <- scRNA_seq

##====================================================================================================##
## STEP 3: Intersect and create gene lists for each combination
##====================================================================================================##

cat("\n=== STEP 3: INTERSECTING GENES AND CREATING GENE LISTS ===\n")

# Get all unique combinations
combinations <- top_100_genes %>%
  distinct(peak_caller, celltype)

cat("Found", nrow(combinations), "unique combinations\n")

# Create intersected gene lists
intersected_gene_lists <- list()

for (i in 1:nrow(combinations)) {
  tool <- combinations$peak_caller[i]
  cell_type <- combinations$celltype[i]
  
  # Get top 100 genes for this combination
  genes <- top_100_genes %>%
    filter(peak_caller == tool, celltype == cell_type) %>%
    pull(gene_name_to_use)
  
  # Intersect with scRNA-seq genes
  intersected_genes <- intersect(genes, rownames(scRNA_main))
  
  # Only keep if we have at least 10 genes
  if (length(intersected_genes) >= 10) {
    list_name <- paste(tool, cell_type, sep = "_")
    intersected_gene_lists[[list_name]] <- intersected_genes
    
    cat(sprintf("%-25s: %3d / 100 genes intersected\n", 
                list_name, length(intersected_genes)))
  }
}

cat("\n✓ Created", length(intersected_gene_lists), "intersected gene lists\n")

# Save intersected gene lists
saveRDS(intersected_gene_lists, file.path(out_dir, "intersected_gene_lists.rds"))

# Also save as text files
lists_dir <- file.path(out_dir, "gene_lists")
if (!dir.exists(lists_dir)) dir.create(lists_dir)

for (list_name in names(intersected_gene_lists)) {
  writeLines(intersected_gene_lists[[list_name]], 
             file.path(lists_dir, paste0(list_name, ".txt")))
}

##====================================================================================================##
## STEP 4: Calculate module scores using AddModuleScore()
##====================================================================================================##

cat("\n=== STEP 4: CALCULATING MODULE SCORES WITH AddModuleScore() ===\n")

# Filter scRNA-seq to common cell types only
common_celltypes <- unique(top_100_genes$celltype)
scRNA_celltypes <- unique(scRNA_main$predicted_celltype)
celltypes_for_analysis <- intersect(common_celltypes, scRNA_celltypes)

if (length(celltypes_for_analysis) == 0) {
  stop("Error: No common cell types between datasets!")
}

cat("Cell types for analysis:", paste(celltypes_for_analysis, collapse = ", "), "\n")

# Filter scRNA-seq to these cell types
scRNA_main <- subset(scRNA_main, subset = predicted_celltype %in% celltypes_for_analysis)
cat("Cells after filtering:", ncol(scRNA_main), "\n")

# Set default assay
DefaultAssay(scRNA_main) <- "RNA"
Idents(scRNA_main) <- scRNA_main$predicted_celltype

# Ensure UMAP exists
if (!"umap" %in% names(scRNA_main@reductions)) {
  cat("Creating UMAP...\n")
  
  # Run PCA if needed
  if (!"pca" %in% names(scRNA_main@reductions)) {
    scRNA_main <- RunPCA(scRNA_main, 
                        features = VariableFeatures(scRNA_main),
                        npcs = min(30, ncol(scRNA_main)-1))
  }
  
  # Run UMAP
  scRNA_main <- RunUMAP(scRNA_main, 
                       reduction = "pca", 
                       dims = 1:min(30, ncol(scRNA_main[["pca"]])))
}

# Calculate module scores for each intersected gene list
cat("\nCalculating module scores for each combination:\n")

module_score_columns <- c()

for (list_name in names(intersected_gene_lists)) {
  genes <- intersected_gene_lists[[list_name]]
  
  if (length(genes) >= 10) {
    cat("  Calculating for", list_name, "... ")
    
    # Clean list name for column name
    clean_name <- gsub("[^[:alnum:]_]", "_", list_name)
    
    # Calculate module score using AddModuleScore()
    scRNA_main <- AddModuleScore(
      object = scRNA_main,
      features = list(genes),  # Must be a list
      name = paste0("Score_", clean_name, "_"),
      pool = NULL,
      nbin = 24,
      ctrl = 100,
      k = FALSE,
      assay = "RNA",
      seed = 2024,
      search = FALSE,
      slot = "data"
    )
    
    # Store column name
    col_name <- paste0("Score_", clean_name, "_1")
    module_score_columns <- c(module_score_columns, col_name)
    
    cat("DONE (", length(genes), " genes)\n")
  }
}

cat("\n✓ Calculated module scores for", length(module_score_columns), "combinations\n")

##====================================================================================================##
## STEP 5: Plot feature visualizations
##====================================================================================================##

cat("\n=== STEP 5: CREATING FEATURE VISUALIZATIONS ===\n")

# Helper function to darken colors
darken_color <- function(color, factor = 0.3) {
  col <- col2rgb(color)
  col <- col * (1 - factor)
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}

# Create directories for plots
plot_dir <- file.path(out_dir, "visualizations")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

umap_dir <- file.path(plot_dir, "umap_plots")
violin_dir <- file.path(plot_dir, "violin_plots")
heatmap_dir <- file.path(plot_dir, "heatmaps")
if (!dir.exists(umap_dir)) dir.create(umap_dir)
if (!dir.exists(violin_dir)) dir.create(violin_dir)
if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir)

# Define colors for tools
tool_colors <- c(
  "DROMPAplus" = "#0072B2",  # Blue
  "Genrich"    = "#E69F00",  # Orange
  "GoPeaks"    = "#009E73",  # Green
  "HOMER"      = "#CC79A7",  # Pink
  "MACS2"      = "#56B4E9",  # Sky Blue
  "SEACR"      = "#D55E00",  # Red
  "SICER2"     = "#F0E442"   # Yellow
)

# 1. Create UMAP feature plots for each module score with tool-specific colors
cat("\n1. Creating UMAP feature plots with tool-specific colors:\n")

for (i in seq_along(module_score_columns)) {
  score_col <- module_score_columns[i]
  list_name <- names(intersected_gene_lists)[i]
  
  cat("  Plotting", list_name, "... ")
  
  # Extract tool name from list_name
  tool <- strsplit(list_name, "_")[[1]][1]
  
  # Get base color for this tool
  if (tool %in% names(tool_colors)) {
    base_color <- tool_colors[tool]
    
    # Create color gradient from light to dark
    plot_colors <- colorRampPalette(c("lightgrey", base_color, darken_color(base_color, 0.4)))(100)
  } else {
    # Default blue gradient if tool not found
    plot_colors <- colorRampPalette(c("lightgrey", "blue", "darkblue"))(100)
  }
  
  # Create feature plot with tool-specific colors
  p <- FeaturePlot(
    scRNA_main,
    features = score_col,
    reduction = "umap",
    pt.size = 0.4,
    order = TRUE,
    label = FALSE,
    combine = TRUE,
    cols = plot_colors  # Use tool-specific color gradient
  ) +
    ggtitle(paste("Module Score:", list_name)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  # Save plot
  ggsave(
    file.path(umap_dir, paste0(hist, "_umap_", gsub("[^[:alnum:]]", "_", list_name), ".png")),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  cat("✓ (", tool, " color)\n")
}

# 2. Create violin plots for each module score (by cell type)
cat("\n2. Creating violin plots by cell type:\n")

# Extract metadata for plotting
metadata_df <- scRNA_main[[]]
metadata_df$CellType <- scRNA_main$predicted_celltype

# Add module scores to metadata
for (score_col in module_score_columns) {
  list_name <- names(intersected_gene_lists)[which(module_score_columns == score_col)]
  
  # Create long format data for this module
  plot_data <- data.frame(
    CellType = metadata_df$CellType,
    Score = metadata_df[[score_col]],
    Module = list_name
  )
  
  # Create violin plot
  p <- ggplot(plot_data, aes(x = CellType, y = Score, fill = CellType)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_summary(fun = median, geom = "point", size = 2, color = "black") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = paste("Module:", list_name),
      x = "Cell Type",
      y = "Module Score"
    ) +
    scale_fill_brewer(palette = "Set3")
  
  # Save plot
  ggsave(
    file.path(violin_dir, paste0("violin_", hist , gsub("[^[:alnum:]]", "_", list_name), ".png")),
    p,
    width = max(8, length(unique(plot_data$CellType)) * 1.2),
    height = 6,
    dpi = 300
  )
  
  cat("  Created violin plot for", list_name, "✓\n")
}

# 3. Create summary heatmap of average module scores
cat("\n3. Creating summary heatmap:\n")

# Prepare data for heatmap
heatmap_data <- data.frame(
  Cell = colnames(scRNA_main),
  CellType = scRNA_main$predicted_celltype
)

# Add module scores
for (score_col in module_score_columns) {
  list_name <- names(intersected_gene_lists)[which(module_score_columns == score_col)]
  heatmap_data[[list_name]] <- metadata_df[[score_col]]
}

# Calculate average scores per cell type
avg_scores <- heatmap_data %>%
  group_by(CellType) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames("CellType")

# Transpose for better visualization
avg_scores_t <- t(avg_scores)

# Create heatmap
png(file.path(heatmap_dir, paste0("HumanPBMC_", hist, "_average_module_scores_heatmap_with_input.png")), 
    width = 12, height = 10, units = "in", res = 300)
pheatmap(
  avg_scores_t,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Average Module Scores by Cell Type",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  fontsize = 11,
  border_color = NA,
  cellwidth = 20,
  cellheight = 15
)
dev.off()

cat("  ✓ Heatmap created\n")

# DEBUG: Check what cell types and combinations actually exist
cat("\n=== DEBUGGING: Checking actual data structure ===\n")

# Check what combinations we actually have
cat("First 10 combinations in intersected_gene_lists:\n")
print(head(names(intersected_gene_lists), 10))

# Extract all unique cell types from the data
all_combinations <- names(intersected_gene_lists)
extracted_celltypes <- sapply(strsplit(all_combinations, "_"), function(x) {
  paste(x[-1], collapse = "_")
})

cat("\nUnique cell types found in intersected_gene_lists:\n")
unique_celltypes <- sort(unique(extracted_celltypes))
print(unique_celltypes)
cat("Total unique cell types:", length(unique_celltypes), "\n")

# Check if "B" cell type exists
cat("\nChecking for B cell types:\n")
b_combinations <- all_combinations[grep("_B", all_combinations)]
if (length(b_combinations) > 0) {
  cat("Found B cell combinations:\n")
  print(b_combinations)
} else {
  cat("No B cell combinations found. Searching for any 'B' pattern...\n")
  b_pattern <- all_combinations[grep("B", all_combinations, ignore.case = TRUE)]
  if (length(b_pattern) > 0) {
    cat("Found these combinations with 'B':\n")
    print(b_pattern)
  } else {
    cat("No combinations with 'B' found at all.\n")
  }
}

# 4. Create a comprehensive grid plot for publication (Auto-detect layout)
cat("\n4. Creating comprehensive grid plot with auto-detected layout:\n")

# Extract all unique tools and cell types from the data
all_tools <- unique(sapply(strsplit(names(intersected_gene_lists), "_"), function(x) x[1]))
all_celltypes <- unique(sapply(strsplit(names(intersected_gene_lists), "_"), function(x) {
  paste(x[-1], collapse = "_")
}))

cat("Auto-detected tools (", length(all_tools), "): ", paste(sort(all_tools), collapse = ", "), "\n", sep = "")
cat("Auto-detected cell types (", length(all_celltypes), "): ", paste(sort(all_celltypes), collapse = ", "), "\n", sep = "")

# Define tool order (prioritize the 7 main tools)
tool_order <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
# Add any other tools found
other_tools <- setdiff(all_tools, tool_order)
tool_order <- c(tool_order, sort(other_tools))

# Sort cell types for consistent ordering
celltype_order <- sort(all_celltypes)

cat("\nFinal tool order (", length(tool_order), "): ", paste(tool_order, collapse = ", "), "\n", sep = "")
cat("Final cell type order (", length(celltype_order), "): ", paste(celltype_order, collapse = ", "), "\n", sep = "")

# Create a matrix to store plots
plot_matrix <- list()

# Generate all combinations in the desired order
for (celltype in celltype_order) {
  for (tool in tool_order) {
    list_name <- paste(tool, celltype, sep = "_")
    
    # Debug for first row
    if (celltype == celltype_order[1]) {
      cat("  Checking: ", list_name, " ... ")
    }
    
    # Check if this combination exists in our data
    if (list_name %in% names(intersected_gene_lists)) {
      # Find the corresponding score column
      score_idx <- which(names(intersected_gene_lists) == list_name)
      score_col <- module_score_columns[score_idx]
      
      # Get tool color
      if (tool %in% names(tool_colors)) {
        base_color <- tool_colors[tool]
        plot_colors <- colorRampPalette(c("lightgrey", base_color, darken_color(base_color, 0.4)))(100)
      } else {
        plot_colors <- colorRampPalette(c("lightgrey", "blue", "darkblue"))(100)
      }
      
      # Create compact UMAP plot
      p <- FeaturePlot(
        scRNA_main,
        features = score_col,
        reduction = "umap",
        pt.size = 0.1,
        order = TRUE,
        label = FALSE,
        combine = TRUE,
        cols = plot_colors
      ) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
          legend.position = "none",
          plot.margin = margin(1, 1, 1, 1, "mm")
        ) +
        ggtitle(paste0(tool, "\n", celltype))
      
      plot_matrix[[list_name]] <- p
      if (celltype == celltype_order[1]) cat("FOUND ✓\n")
    } else {
      # Create empty plot for missing combination
      plot_matrix[[list_name]] <- ggplot() + 
        theme_void() + 
        ggtitle(paste0(tool, "\n", celltype, "\n(No data)")) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 8, color = "gray50")
        )
      if (celltype == celltype_order[1]) cat("NOT FOUND (empty plot)\n")
    }
  }
}

# Arrange in grid (rows = cell types, columns = tools)
cat("\n  Arranging plots in", length(celltype_order), "x", length(tool_order), "grid...\n")

# Create the grid layout
grid_plots <- list()
for (i in 1:length(celltype_order)) {
  for (j in 1:length(tool_order)) {
    list_name <- paste(tool_order[j], celltype_order[i], sep = "_")
    grid_plots[[paste0("r", i, "c", j)]] <- plot_matrix[[list_name]]
  }
}

# Create the grid using patchwork
grid_plot <- wrap_plots(
  grid_plots,
  nrow = length(celltype_order),  # rows for cell types
  ncol = length(tool_order),      # columns for tools
  byrow = TRUE
) +
  plot_annotation(
    title = "Module Scores: Cell Types vs Peak Callers",
    subtitle = paste("UMAP visualizations (", length(celltype_order), " cell types × ", length(tool_order), " tools)", sep = ""),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50")
    )
  )

# Calculate appropriate dimensions based on grid size
plot_width <- max(20, length(tool_order) * 3)
plot_height <- max(20, length(celltype_order) * 3)

# Save as PDF (better for publication)
cat("\n  Saving as PDF...\n")
pdf(
  file.path(plot_dir, paste0("HumanPBMC_", hist, "_module_scores_grid_with_input.pdf")),
  width = plot_width,
  height = plot_height
)
print(grid_plot)
dev.off()

# Also save as high-resolution PNG
cat("  Saving as PNG...\n")
ggsave(
  file.path(plot_dir, paste0("HumanPBMC_", hist, "_module_scores_grid_with_input.png")),
  grid_plot,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

cat("  ✓ Grid plot created (", length(celltype_order), " rows × ", length(tool_order), " columns)\n", sep = "")
cat("    - PDF saved: module_scores_grid_auto.pdf (", plot_width, "×", plot_height, " inches)\n", sep = "")
cat("    - PNG saved: module_scores_grid_auto.png\n")

##====================================================================================================##
## STEP 6: Statistical analysis
###########################################################################

cat("\n=== STEP 6: PERFORMING STATISTICAL ANALYSIS ===\n")

# Prepare data for statistical analysis
stats_data <- data.frame(
  Cell = colnames(scRNA_main),
  CellType = scRNA_main$predicted_celltype
)

# Add module scores
for (score_col in module_score_columns) {
  list_name <- names(intersected_gene_lists)[which(module_score_columns == score_col)]
  stats_data[[list_name]] <- metadata_df[[score_col]]
}

# Perform statistical tests for each module
stats_results <- data.frame()

for (list_name in names(intersected_gene_lists)) {
  if (list_name %in% colnames(stats_data)) {
    # Extract tool and cell type from list name
    parts <- strsplit(list_name, "_")[[1]]
    tool <- parts[1]
    target_cell_type <- paste(parts[-1], collapse = " ")
    
    # Get scores for this module
    scores <- stats_data[[list_name]]
    cell_types <- stats_data$CellType
    
    # Create groups: target cell type vs others
    is_target <- cell_types == target_cell_type
    
    if (sum(is_target) > 5 && sum(!is_target) > 5) {
      target_scores <- scores[is_target]
      other_scores <- scores[!is_target]
      
      # T-test
      t_test <- t.test(target_scores, other_scores)
      
      # Wilcoxon test (non-parametric)
      wilcox_test <- wilcox.test(target_scores, other_scores)
      
      # Calculate effect sizes
      # Cohen's d
      pooled_sd <- sqrt(((length(target_scores)-1)*var(target_scores) + 
                          (length(other_scores)-1)*var(other_scores)) / 
                         (length(target_scores) + length(other_scores) - 2))
      cohens_d <- (mean(target_scores) - mean(other_scores)) / pooled_sd
      
      # Hedge's g (correction for small sample size)
      n_total <- length(target_scores) + length(other_scores)
      hedge_correction <- 1 - (3 / (4*n_total - 9))
      hedges_g <- cohens_d * hedge_correction
      
      # Store results
      stats_results <- rbind(stats_results, data.frame(
        Tool = tool,
        Target_CellType = target_cell_type,
        Module = list_name,
        N_Target = length(target_scores),
        N_Other = length(other_scores),
        Mean_Target = mean(target_scores),
        Mean_Other = mean(other_scores),
        SD_Target = sd(target_scores),
        SD_Other = sd(other_scores),
        T_statistic = t_test$statistic,
        T_p_value = t_test$p.value,
        Wilcoxon_p_value = wilcox_test$p.value,
        Cohens_d = cohens_d,
        Hedges_g = hedges_g,
        Genes_Used = length(intersected_gene_lists[[list_name]])
      ))
    }
  }
}

# Adjust p-values for multiple testing
if (nrow(stats_results) > 0) {
  stats_results$T_p_adj <- p.adjust(stats_results$T_p_value, method = "BH")
  stats_results$Wilcoxon_p_adj <- p.adjust(stats_results$Wilcoxon_p_value, method = "BH")
  
  # Save statistical results
  write.csv(stats_results, file.path(out_dir, paste0(hist, "_statistical_analysis_results_with_input.csv")), row.names = FALSE)
  
  cat("\n✓ Statistical analysis completed\n")
  cat("  Significant results (FDR < 0.05):\n")
  cat("    T-test:", sum(stats_results$T_p_adj < 0.05), "/", nrow(stats_results), "\n")
  cat("    Wilcoxon:", sum(stats_results$Wilcoxon_p_adj < 0.05), "/", nrow(stats_results), "\n")
  
  # Print summary of best results
  best_results <- stats_results %>%
    arrange(T_p_adj) %>%
    head(5)
  
  cat("\n  Top 5 results by adjusted p-value:\n")
  print(best_results[, c("Tool", "Target_CellType", "T_p_adj", "Hedges_g", "Genes_Used")])
}

}



##===================================================================================================
## Figure 3 H, I
##===================================================================================================

allcell_tpm <- read.delim("/home/wahid/project_scHMTF/GSE195725_processed_data/scRNA_seq/10x_scRNAseq/pbmc_celltype_averaged_TPM_with_AllCell_weighted.txt.gz", row.names = 1)



library(dplyr)
library(readr)
library(purrr)

# Define the directory path
base_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_100genes_celltypeswise_bytools/result_GT_scrnaseq_with_input"

# List all CSV files with pattern
csv_files <- list.files(
  path = base_dir,
  pattern = "HumanPBMC_.*_top_100_genes_all_combinations_with_input\\.csv$",
  full.names = TRUE
)

cat("Found", length(csv_files), "CSV files:\n")
print(csv_files)

# Function to extract histone mark from filename
extract_histone_mark <- function(filename) {
  # Extract part between "HumanPBMC_" and "_top_100"
  basename <- basename(filename)
  histone <- gsub("HumanPBMC_", "", basename)
  histone <- gsub("_top_100_genes_all_combinations_with_input\\.csv", "", histone)
  return(histone)
}

# Read and combine all CSV files with histone mark
combined_data <- map_dfr(csv_files, function(file_path) {
  # Read the CSV file
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Extract histone mark from filename
  histone_mark <- extract_histone_mark(file_path)
  
  # Add histone_mark column as first column
  df <- df %>%
    mutate(histone = histone_mark, .before = 1)
  
  return(df)
})

# Reorder columns if needed
combined_data <- combined_data %>%
  select(histone, peak_caller, celltype, gene_name_to_use, 
         weighted_peak_count, log_tpm, tpm_expression, everything())

# Check the combined data
cat("\n=== COMBINED DATA SUMMARY ===\n")
cat("Total rows:", nrow(combined_data), "\n")
cat("Unique histone marks:", paste(unique(combined_data$histone), collapse = ", "), "\n")
cat("Unique peak callers:", paste(unique(combined_data$peak_caller), collapse = ", "), "\n")
cat("Unique cell types:", paste(unique(combined_data$celltype), collapse = ", "), "\n")




#########################################################################################################################
#### adding information of AllCells Tpm
#########################################################################################################################


# Make sure allcell_tpm has row names properly set
# If not, set row names from the Gene column
if ("Gene" %in% colnames(allcell_tpm)) {
  rownames(allcell_tpm) <- allcell_tpm$Gene
}

# Add AllCell_tpm column to combined_data by matching gene names
combined_data$AllCell_tpm <- allcell_tpm[combined_data$gene_name_to_use, "AllCell"]

# Check the result
cat("New column added: AllCell_tpm\n")
cat("Dimensions of combined_data:", dim(combined_data), "\n")
cat("Columns:", colnames(combined_data), "\n\n")

# Show first few rows
cat("First 6 rows of combined_data with AllCell_tpm:\n")
print(head(combined_data))

# Check for any NA values (genes that didn't match)
missing_genes <- sum(is.na(combined_data$AllCell_tpm))
cat("\nNumber of genes with no match (NA values):", missing_genes, "\n")

if (missing_genes > 0) {
  cat("First few unmatched genes:\n")
  unmatched <- unique(combined_data$gene_name_to_use[is.na(combined_data$AllCell_tpm)])
  print(head(unmatched, 10))
}

# If you want to see examples of the matching
cat("\n=== Example of matching ===\n")
example_genes <- head(unique(combined_data$gene_name_to_use), 5)

for (gene in example_genes) {
  tpm_value <- allcell_tpm[gene, "AllCell"]
  cat(gene, "-> AllCell_tpm =", tpm_value, "\n")
}



# Add log_ratio_tpm column = log(tpm_expression / AllCell_tpm)
combined_data$ratio_tpm <- (combined_data$tpm_expression / combined_data$AllCell_tpm)
combined_data$log_ratio_tpm <- log(combined_data$ratio_tpm)

# Check the result
cat("log_ratio_tpm column added to combined_data\n")
cat("Dimensions:", dim(combined_data), "\n")
cat("Columns:", colnames(combined_data), "\n\n")

# Show first few rows
cat("First 6 rows:\n")
print(head(combined_data))



##====================================================================================================##
## Plot - Corrected for Human PBMC data
##====================================================================================================##

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# CRITICAL FIX: Handle infinite values properly
# Create FRiP summary data 
ratio_tpm_summary <- combined_data %>%
  # Filter out infinite values before calculating statistics
  filter(!is.infinite(log_ratio_tpm)) %>%
  group_by(peak_caller, histone, celltype) %>%
  summarize(
    median_log_ratio_tpm = median(log_ratio_tpm, na.rm = TRUE),
    mean_log_ratio_tpm = mean(log_ratio_tpm, na.rm = TRUE),
    sd_log_ratio_tpm = sd(log_ratio_tpm, na.rm = TRUE),
    n_genes = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(peak_caller),
    Histone = factor(histone),
    # Replace any remaining NaN in sd_ratio_tpm with 0
    sd_log_ratio_tpm = ifelse(is.nan(sd_log_ratio_tpm), 0, sd_log_ratio_tpm)
  )

# Check the summary data
cat("Summary data created successfully!\n")
cat("Dimensions:", dim(ratio_tpm_summary), "\n")
cat("Columns:", colnames(ratio_tpm_summary), "\n")
cat("\nSample of data (showing first 10 rows):\n")
print(head(ratio_tpm_summary, 10))

# Check for negative values (log ratios can be negative)
cat("\n=== Data Statistics ===\n")
cat("Range of median_ratio_tpm:", 
    round(min(ratio_tpm_summary$median_log_ratio_tpm, na.rm = TRUE), 3), "to",
    round(max(ratio_tpm_summary$median_log_ratio_tpm, na.rm = TRUE), 3), "\n")
cat("Number of unique cell types:", length(unique(ratio_tpm_summary$celltype)), "\n")
cat("Number of unique histones:", length(unique(ratio_tpm_summary$histone)), "\n")
cat("Number of unique methods:", length(unique(ratio_tpm_summary$Method)), "\n")

# Get all unique cell types
all_celltypes <- unique(ratio_tpm_summary$celltype)
cat("\nCell types found:", paste(all_celltypes, collapse = ", "), "\n")

# Check which histones are available for each cell type
cat("\n=== Available Histones per Cell Type ===\n")
for (cell in all_celltypes) {
  histones <- unique(ratio_tpm_summary$histone[ratio_tpm_summary$celltype == cell])
  cat(sprintf("%s: %s\n", cell, paste(histones, collapse = ", ")))
}

# Define y-axis limits based on your data
# For log ratios, we need to handle negative values
data_min <- min(ratio_tpm_summary$median_log_ratio_tpm - ratio_tpm_summary$sd_log_ratio_tpm, na.rm = TRUE)
data_max <- max(ratio_tpm_summary$median_log_ratio_tpm + ratio_tpm_summary$sd_log_ratio_tpm, na.rm = TRUE)

# Add some padding
global_ymin <- -2.5
global_ymax <- 2.5

# Ensure y-axis includes 0 and 1 for reference

cat(sprintf("\nGlobal y-axis limits: %.2f to %.2f\n", global_ymin, global_ymax))

# Define colors
all_methods <- unique(ratio_tpm_summary$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

# Create a list to store plots
plot_list <- list()

# Create individual plots for each cell type
for (celltype in all_celltypes) {
  
  # Filter data for this cell type
  celltype_data <- ratio_tpm_summary %>%
    filter(celltype == !!celltype)
  
  # Get available histones for this cell type
  available_histones <- unique(celltype_data$histone)
  
  # Order histones consistently if needed
  histone_order <- sort(available_histones)
  celltype_data <- celltype_data %>%
    mutate(Histone = factor(histone, levels = histone_order))
  
  # Create plot for this cell type
  cell_plot <- ggplot(celltype_data, 
                      aes(x = Histone, 
                          y = median_log_ratio_tpm,
                          color = Method)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
    
    geom_errorbar(
      aes(ymin = median_log_ratio_tpm - sd_log_ratio_tpm,
          ymax = median_log_ratio_tpm + sd_log_ratio_tpm),
      width = 0.25, 
      position = position_dodge(width = 0.6),
      alpha = 0.7,
      size = 0.6
    ) +
    
    scale_color_manual(values = method_colors,
                       name = "Method") +
    
    scale_y_continuous(
      limits = c(global_ymin, global_ymax),
      breaks = seq(round(global_ymin), round(global_ymax), by = 0.5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    # Add reference lines
    geom_hline(
      yintercept = 0, 
      linetype = "dashed", 
      color = "blue", 
      alpha = 0.5,
      size = 0.5
    ) +
    
    geom_hline(
      yintercept = 1, 
      linetype = "dashed", 
      color = "red", 
      alpha = 0.5,
      size = 0.6
    ) +
    
    labs(
      title = paste("Cell Type:", celltype),
      x = "Histone Modification",
      y = "Log Expression Ratio",
      subtitle = paste(length(available_histones), "histone modifications")
    ) +
    
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 10),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.3)
    )
  
  # Store plot in list
  plot_list[[celltype]] <- cell_plot
}

# Determine grid layout based on number of cell types
n_celltypes <- length(all_celltypes)

# Create a combined plot using patchwork
if (n_celltypes <= 4) {
  combined_plot <- wrap_plots(plot_list, ncol = 2, nrow = 2)
} else if (n_celltypes <= 6) {
  combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 2)
} else if (n_celltypes <= 9) {
  combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 3)
} else {
  combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 3)
}

# Add overall title and subtitle
combined_plot <- combined_plot +
  plot_annotation(
    title = "Human PBMC with Input: Cell-type Specific Log Expression Ratios",
    caption = "Blue dashed line: y=0 (no change)\nRed dashed line: y=1 (fold change threshold)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 5)),
      plot.caption = element_text(size = 10, hjust = 0, color = "gray40")
    )
  )

# Print combined plot
print(combined_plot)

# Save combined plot
output_path <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_100genes_celltypeswise_bytools/HumanPBMC_all_celltypes_log_grid_plot_correctedout.pdf"

ggsave(output_path,
       plot = combined_plot,
       width = 14,
       height = 10,
       device = "pdf",
       bg = "white")

# Also save as PNG
png_path <- gsub("\\.pdf$", ".png", output_path)
ggsave(png_path,
       plot = combined_plot,
       width = 14,
       height = 10,
       dpi = 300,
       bg = "white")

cat("\n✓ Grid plot saved for all cell types\n")

# Faceted plot with free y-axis (better for comparison)
faceted_plot <- ggplot(ratio_tpm_summary, 
                       aes(x = Histone, 
                           y = median_log_ratio_tpm,
                           color = Method)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  
  geom_errorbar(
    aes(ymin = median_log_ratio_tpm - sd_log_ratio_tpm,
        ymax = median_log_ratio_tpm + sd_log_ratio_tpm),
    width = 0.2, 
    position = position_dodge(width = 0.6),
    alpha = 0.7,
    size = 0.5
  ) +
  
  facet_wrap(~ celltype, 
             ncol = 4,
             scales = "free") +  # Free scales for both axes
  
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
                     # FIXED y-axis scale: -2.5 to 2.5
  scale_y_continuous(
    limits = c(-2.5, 2.5),
    breaks = seq(-2.5, 2.5, by = 0.5),
    expand = expansion(mult = c(0.02, 0.02))
    ) +
  
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "blue", 
    alpha = 0.3,
    size = 0.4
  ) +
  
  geom_hline(
    yintercept = 1, 
    linetype = "dashed", 
    color = "red", 
    alpha = 0.3,
    size = 0.4
  ) +
  
  labs(
    title = "Human PBMC with Input: Cell-type Specific Log Expression Ratios",
    x = "Histone Modification",
    y = "Log Expression Ratio"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray90", color = "gray70"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray92"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
  )

print(faceted_plot)

# Save faceted plot
faceted_path <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_100genes_celltypeswise_bytools/HumanPBMC_all_celltypes_log_faceted_plot_free_scales_with_input.pdf"

ggsave(faceted_path,
       plot = faceted_plot,
       width = 12,
       height = 9,
       device = "pdf",
       bg = "white")

cat("✓ Faceted plot with free scales saved\n")

# Summary statistics
cat("\n=== Final Summary Statistics ===\n")
summary_stats <- ratio_tpm_summary %>%
  group_by(celltype) %>%
  summarise(
    n_observations = n(),
    min_median = min(median_log_ratio_tpm, na.rm = TRUE),
    max_median = max(median_log_ratio_tpm, na.rm = TRUE),
    avg_sd = mean(sd_log_ratio_tpm, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

#sd_ratio_tpm 




                 
