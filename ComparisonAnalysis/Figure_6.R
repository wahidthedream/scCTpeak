###############################################################################################################################
## Figure 6 A
###############################################################################################################################
### ChromHMM Part B R-script
### Using Roadmap data for all blood cell
### Without Input Case
### R script
###############################################################################################################################
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(tibble)

# =============================================================================
# PACKAGE CONFLICT RESOLUTION
# =============================================================================
# Resolve all package conflicts at the beginning
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(base::unname)
conflicts_prefer(scales::viridis_pal)

# =============================================================================
# GLOBAL VARIABLES
# =============================================================================
celltypes <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")  # 8 cell types as ROWS
histones <- c("H3K27ac-b", "H3K27ac-s","H3K27me3","H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
tools <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")  # 7 tools as COLUMNS

# Chromatin state colors (Roadmap Epigenomics standard)
chromatin_colors <- c(
  "1_TssA" = "#FF0000",        "2_TssAFlnk" = "#FF4500",   
  "3_TxFlnk" = "#32CD32",      "4_Tx" = "#008000",         
  "5_TxWk" = "#006400",        "6_EnhG" = "#ADFF2F",        
  "7_Enh" = "#FFFF00",         "8_ZNF/Rpts" = "#66CDAA",    
  "9_Het" = "#AFEEEE",         "10_TssBiv" = "#CD5C5C",     
  "11_BivFlnk" = "#E9967A",    "12_EnhBiv" = "#BDB76B",     
  "13_ReprPC" = "#808080",     "14_ReprPCWk" = "#DCDCDC",   
  "15_Quies" = "#FFFFFF"
)

# State descriptions
state_descriptions <- c(
  "1_TssA" = "Active TSS", "2_TssAFlnk" = "Flanking Active TSS", 
  "3_TxFlnk" = "Transcr. at gene 5' and 3'", "4_Tx" = "Strong transcription",
  "5_TxWk" = "Weak transcription", "6_EnhG" = "Genic enhancers",
  "7_Enh" = "Enhancers", "8_ZNF/Rpts" = "ZNF genes & repeats",
  "9_Het" = "Heterochromatin", "10_TssBiv" = "Bivalent/Poised TSS",
  "11_BivFlnk" = "Flanking Bivalent TSS/Enh", "12_EnhBiv" = "Bivalent Enhancer",
  "13_ReprPC" = "Repressed PolyComb", "14_ReprPCWk" = "Weak Repressed PolyComb",
  "15_Quies" = "Quiescent/Low"
)

# =============================================================================
# IMPROVED DATA PROCESSING FUNCTION - HANDLES ALL TOOL FORMATS
# =============================================================================
summarize_chromatin_states <- function(file_path) {
  # Check if file exists and is readable
  if (!file.exists(file_path)) {
    warning("File does not exist: ", file_path)
    return(NULL)
  }
  
  # Check if file is empty
  if (file.size(file_path) == 0) {
    cat("Empty file:", basename(file_path), "\n")
    return(NULL)
  }
  
  # Read first few lines to determine format
  first_lines <- tryCatch({
    readLines(file_path, n = 10, warn = FALSE)
  }, error = function(e) {
    warning("Cannot read file: ", basename(file_path), " - ", e$message)
    return(NULL)
  })
  
  if (is.null(first_lines) || length(first_lines) == 0) {
    return(NULL)
  }
  
  # Remove empty lines and header lines
  first_lines <- first_lines[first_lines != ""]
  first_lines <- first_lines[!grepl("^track|^browser|^#", first_lines)]
  
  if (length(first_lines) == 0) {
    cat("File contains only header/empty lines:", basename(file_path), "\n")
    return(NULL)
  }
  
  # Determine the number of columns based on file content
  n_cols <- length(str_split(first_lines[1], "\t")[[1]])
  
  cat("File:", basename(file_path), "Columns:", n_cols, "\n")
  
  # If only 1 column, it's likely empty or malformed
  if (n_cols == 1) {
    cat("⚠️ File has only 1 column, likely empty or malformed:", basename(file_path), "\n")
    return(NULL)
  }
  
  # Read data based on number of columns
  df <- tryCatch({
    if (n_cols == 8) {
      # Format 1: DROMPAplus, GoPeaks format (8 columns)
      readr::read_tsv(
        file_path,
        col_names = c("chrom", "start", "end", "peak_chrom", 
                      "peak_start", "peak_end", "state", "overlap"),
        show_col_types = FALSE,
        comment = "#"
      )
    } else if (n_cols == 15) {
      # Format 2: Genrich format (15 columns)
      readr::read_tsv(
        file_path,
        col_names = c("chrom", "start", "end", "peak_id", "score", "strand",
                      "unknown1", "unknown2", "unknown3", "unknown4", "peak_chrom", 
                      "peak_start", "peak_end", "state", "overlap"),
        show_col_types = FALSE,
        comment = "#"
      )
    } else if (n_cols == 14) {
      # Format 3: MACS2 format (11 columns)
      readr::read_tsv(
        file_path,
        col_names = c("chrom", "start", "end", "peak_id", "score", "strand",
                      "unknown1", "unknown2", "unknown3","peak_chrom", "peak_start", "peak_end", "state", "overlap"),
        show_col_types = FALSE,
        comment = "#"
      )
    } else if (n_cols == 11) {
      # Format 4: HOMER format (10 columns)
      readr::read_tsv(
        file_path,
        col_names = c("chrom", "start", "end", "peak_id", "score", "strand",
                      "peak_chrom", "peak_start", "peak_end", "state", "overlap"),
        show_col_types = FALSE,
        comment = "#"
      ) %>%
        mutate(overlap = pmin(end, peak_end) - pmax(start, peak_start))
    } else if (n_cols == 11) {
      # Format 5: SEACR format (9 columns)
      readr::read_tsv(
        file_path,
        col_names = c("chrom", "start", "end", "score1", "score2", "peak_id",
                      "peak_chrom", "peak_start", "peak_end","state", "overlap" ),
        show_col_types = FALSE,
        comment = "#"
      ) %>%
        mutate(
          state = str_extract(peak_id, "\\d+_\\w+"),
          overlap = pmin(end, peak_end) - pmax(start, peak_start)
        )
    } else if (n_cols == 10) {
      # Format 6: SICER2 format (7 columns)
      readr::read_tsv(
        file_path,
        col_names = c("chrom", "start", "end", "peak_id", "strand",
                      "peak_chrom", "peak_start", "peak_end" ,"state", "overlap"),
        show_col_types = FALSE,
        comment = "#",
        col_names = FALSE
      ) %>%
        select(1:8) %>%
        setNames(c("chrom", "start", "end", "peak_id", "strand",
                   "peak_chrom", "peak_start", "peak_end")) %>%
        mutate(
          state = str_extract(peak_id, "\\d+_\\w+"),
          overlap = pmin(end, peak_end) - pmax(start, peak_start)
        )
    } else {
      # Generic format - try to find state column
      temp_df <- readr::read_tsv(file_path, show_col_types = FALSE, comment = "#", col_names = FALSE)
      
      # Find state column (contains pattern like "1_TssA", "15_Quies")
      state_col <- NULL
      for (i in 1:ncol(temp_df)) {
        if (any(grepl("^\\d+_\\w+", temp_df[[i]]))) {
          state_col <- i
          break
        }
      }
      
      if (!is.null(state_col)) {
        temp_df <- temp_df %>%
          mutate(state = .[[state_col]]) %>%
          select(state) %>%
          mutate(
            chrom = "unknown",
            start = 0,
            end = 1,
            overlap = 1
          )
      } else {
        cat("⚠️ Cannot identify state column in:", basename(file_path), "\n")
        return(NULL)
      }
      temp_df
    }
  }, error = function(e) {
    warning("Error reading file ", basename(file_path), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(df) || nrow(df) == 0) {
    cat("⚠️ No data rows in:", basename(file_path), "\n")
    return(NULL)
  }
  
  # Check if state column exists and has valid data
  if (!"state" %in% names(df) || all(is.na(df$state))) {
    cat("⚠️ No valid state information in:", basename(file_path), "\n")
    return(NULL)
  }
  
  # Filter out invalid states and ensure proper format
  valid_states <- paste0(1:15, "_", c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", 
                                    "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", 
                                    "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"))
  
  df <- df %>% 
    filter(!is.na(state)) %>%
    mutate(state = as.character(state)) %>%
    filter(state %in% valid_states)
  
  if (nrow(df) == 0) {
    cat("⚠️ No valid chromatin states found in:", basename(file_path), "\n")
    return(NULL)
  }
  
  # Handle overlap calculation
  if (!"overlap" %in% names(df)) {
    df$overlap <- 1  # Default overlap if not calculable
  }
  
  # Summarize by state
  result <- df %>%
    dplyr::group_by(state) %>%
    dplyr::summarise(
      peak_count = dplyr::n(),
      total_overlap = sum(overlap, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add file information and fractions
  result <- result %>%
    dplyr::mutate(
      filename = basename(file_path),
      sample_id = gsub("\\.bed$", "", filename),
      fraction = peak_count / sum(peak_count),
      fraction_overlap = total_overlap / sum(total_overlap, na.rm = TRUE)
    )
  
  cat("✓ Processed", basename(file_path), "-", nrow(result), "states found\n")
  
  result %>%
    dplyr::arrange(desc(peak_count))
}

# =============================================================================
# INDIVIDUAL TOOL PROCESSING FUNCTION
# =============================================================================
create_celltype_figures <- function(celltype_data, celltype_name, hist, tool) {
  
  # Validate input data
  if (nrow(celltype_data) == 0) {
    warning("No data available for cell type: ", celltype_name)
    return(NULL)
  }
  
  required_cols <- c("eid", "state", "fraction", "peak_count")
  if (!all(required_cols %in% names(celltype_data))) {
    stop("Missing required columns in data for: ", celltype_name)
  }
  
  # Order states numerically
  state_order <- paste0(1:15, "_", c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", 
                                   "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", 
                                   "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"))
  
  celltype_data <- celltype_data %>%
    mutate(state = factor(state, levels = state_order))
  
  # 1. STACKED BAR PLOT - Fraction by EID
  p1 <- ggplot(celltype_data, aes(x = eid, y = fraction, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    scale_fill_manual(values = chromatin_colors, name = "Chromatin State", 
                     labels = state_descriptions, breaks = names(state_descriptions)) +
    labs(
      title = paste("Chromatin State Distribution -", celltype_name),
      subtitle = paste(paste0(hist, " peaks detected by"), tool),
      x = "Sample (EID)",
      y = "Fraction of Peaks"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm")
    )
  
  # 2. HEATMAP - Peak counts across EIDs (only if multiple samples)
  if (length(unique(celltype_data$eid)) > 1) {
    heatmap_data <- celltype_data %>%
      select(eid, state, peak_count) %>%
      pivot_wider(names_from = state, values_from = peak_count, values_fill = 0) %>%
      column_to_rownames("eid") %>%
      as.matrix()
    
    # Reorder columns and ensure proper ordering
    heatmap_data <- heatmap_data[, state_order[state_order %in% colnames(heatmap_data)]]
    
    # Use viridisLite directly to avoid conflicts
    p2 <- Heatmap(heatmap_data,
                  name = "Peak Count",
                  col = colorRamp2(seq(0, max(heatmap_data), length = 100), 
                                  viridisLite::viridis(100)),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  row_names_side = "left",
                  column_names_rot = 45,
                  rect_gp = gpar(col = "white", lwd = 0.3),
                  heatmap_legend_param = list(title = "Peak Count", 
                                            title_gp = gpar(fontsize = 10),
                                            labels_gp = gpar(fontsize = 8)))
  } else {
    p2 <- NULL
  }
  
  # 3. DOT PLOT - Mean fraction across EIDs with proper statistics
  summary_data <- celltype_data %>%
    group_by(state) %>%
    summarise(
      mean_fraction = mean(fraction, na.rm = TRUE),
      sd_fraction = sd(fraction, na.rm = TRUE),
      n_samples = n(),
      sem = sd_fraction/sqrt(n_samples),
      ci_lower = mean_fraction - 1.96 * sem,
      ci_upper = mean_fraction + 1.96 * sem,
      .groups = 'drop'
    )
  
  p3 <- ggplot(summary_data, aes(x = state, y = mean_fraction, color = state)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_fraction - sem, ymax = mean_fraction + sem), 
                  width = 0.2, linewidth = 0.7) +
    scale_color_manual(values = chromatin_colors, guide = "none") +
    labs(
      title = paste("Average Chromatin State Distribution -", celltype_name),
      x = "Chromatin State",
      y = "Mean Fraction ± SEM"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
    )
  
  # 4. INDIVIDUAL EID BAR PLOTS
  eid_plots <- list()
  eids <- unique(celltype_data$eid)
  
  for(eid in eids) {
    eid_data <- celltype_data %>% filter(eid == !!eid)
    
    p <- ggplot(eid_data, aes(x = state, y = fraction, fill = state)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = chromatin_colors, guide = "none") +
      labs(
        title = paste(celltype_name, "-", eid),
        x = "Chromatin State",
        y = "Fraction"
      ) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(linewidth = 0.2)
      )
    
    eid_plots[[eid]] <- p
  }
  
  # Return all plots with metadata
  return(list(
    stacked_bar = p1,
    heatmap = p2,
    dot_plot = p3,
    individual_plots = eid_plots,
    data = celltype_data,
    summary_stats = summary_data,
    n_samples = length(eids),
    processing_date = Sys.time()
  ))
}

# =============================================================================
# GRID PLOT FUNCTION FOR COMBINED FIGURE
# =============================================================================
create_grid_barplot <- function(celltype_data, celltype_name, tool_name) {
  
  # Order states numerically
  state_order <- paste0(1:15, "_", c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", 
                                   "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", 
                                   "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"))
  
  celltype_data <- celltype_data %>%
    mutate(state = factor(state, levels = state_order))
  
  # Create compact barplot for grid
  p <- ggplot(celltype_data, aes(x = eid, y = fraction, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6, linewidth = 0.1) +
    scale_fill_manual(values = chromatin_colors, name = "Chromatin State") +
    labs(
      title = paste0(celltype_name, " - ", tool_name),
      x = NULL,
      y = "Fraction"
    ) +
    theme_classic(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
      axis.text.y = element_text(size = 4),
      axis.title.y = element_text(size = 5, margin = margin(r = 1)),
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
      axis.line = element_line(linewidth = 0.15),
      axis.ticks = element_line(linewidth = 0.15)
    )
  
  return(p)
}

# =============================================================================
# MAIN PROCESSING LOOP - BOTH INDIVIDUAL AND COMBINED
# =============================================================================

for(hist in histones) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("Processing histone:", hist, "\n")
  cat(strrep("=", 80), "\n", sep = "")
  
  # Create main output directory for this histone
  main_output_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_ChromHMM_roadmapdata/Figure_without_input/", hist, "/combined_grid_figures")
  dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Store data for combined grid
  all_grid_data <- list()
  tools_found <- c()
  
  # Process each tool individually
  for(tool in tools) { 
    # Corrected path construction
    bed_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_ChromHMM_roadmapdata/intersect_ChromHMM_without_input/", hist, "/", tool, "/")
    output_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_ChromHMM_roadmapdata/Figure_without_input/", hist, "/", tool, "_celltype_figures")
    dir.create(output_dir, showWarnings = FALSE)
    
    cat("Processing tool:", tool, "\n")
    cat("Input directory:", bed_dir, "\n")
    cat("Output directory:", output_dir, "\n")
    
    # Check if directory exists
    if (!dir.exists(bed_dir)) {
      cat("✗ Directory does not exist:", bed_dir, "\n")
      next
    }
    
    # Process all files with error handling
    bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
    
    if (length(bed_files) == 0) {
      cat("✗ No BED files found in:", bed_dir, "\n")
      next
    }
    
    cat("Found", length(bed_files), "BED files\n")
    
    state_data <- map_dfr(bed_files, function(file) {
      tryCatch({
        result <- summarize_chromatin_states(file)
        if (!is.null(result)) {
          result
        } else {
          NULL
        }
      }, error = function(e) {
        cat("Error processing", basename(file), ":", e$message, "\n")
        return(NULL)
      })
    })
    
    if (is.null(state_data) || nrow(state_data) == 0) {
      cat("✗ No data processed successfully\n")
      next
    }
    
    # Extract cell type and EID from filename
    state_data <- state_data %>%
      mutate(
        cell_type = case_when(
          str_detect(filename, paste0(tool, "_([A-Za-z0-9]+)_", hist)) ~ 
            str_extract(filename, paste0("(?<=", tool, "_)[A-Za-z0-9]+(?=_", hist, ")")),
          str_detect(filename, paste0("_([A-Za-z]+)_", hist)) ~ 
            str_extract(filename, paste0("(?<=_)[A-Za-z]+(?=_", hist, ")")),
          str_detect(filename, "CD4T") ~ "CD4T",
          str_detect(filename, "CD8T") ~ "CD8T",
          str_detect(filename, "Mono") ~ "Mono", 
          str_detect(filename, "DC") ~ "DC",
          str_detect(filename, "NK") ~ "NK",
          str_detect(filename, "otherT") ~ "otherT",
          str_detect(filename, "other") ~ "other",
          str_detect(filename, "B_") ~ "B",
          TRUE ~ "unknown"
        ),
        eid = str_extract(filename, "E[0-9]+"),
        tool = tool,
        histone = hist
      ) %>%
      filter(!is.na(cell_type), cell_type != "unknown")
    
    if (nrow(state_data) == 0) {
      cat("✗ No valid cell types extracted from filenames\n")
      next
    }
    
    # Store data for combined grid
    all_grid_data[[tool]] <- state_data
    tools_found <- c(tools_found, tool)
    
    # Get all cell types for this tool
    cell_types <- unique(state_data$cell_type)
    cat("Found cell types:", paste(cell_types, collapse = ", "), "\n")
    
    # Process individual cell types
    all_celltype_results <- list()
    
    for(celltype in cell_types) {
      cat("Processing cell type:", celltype, "\n")
      
      # Filter data for this cell type
      celltype_data <- state_data %>% filter(cell_type == celltype)
      
      if(nrow(celltype_data) > 0) {
        # Create figures with error handling
        results <- tryCatch({
          create_celltype_figures(celltype_data, celltype, hist, tool)
        }, error = function(e) {
          warning("Failed to create figures for ", celltype, ": ", e$message)
          return(NULL)
        })
        
        if (!is.null(results)) {
          all_celltype_results[[celltype]] <- results
          
          # Create output directory for this cell type
          celltype_dir <- file.path(output_dir, celltype)
          dir.create(celltype_dir, showWarnings = FALSE, recursive = TRUE)
          
          # Save stacked bar plot
          ggsave(file.path(celltype_dir, paste0(celltype, "_", hist, "_", tool, "_stacked_bar_without_input.pdf")),
                 plot = results$stacked_bar, 
                 width = 12, height = 10, units = "cm", dpi = 1200, device = cairo_pdf)
          
          # Save heatmap if it exists
          if (!is.null(results$heatmap)) {
            pdf(file.path(celltype_dir, paste0(celltype, "_", hist, "_", tool, "_heatmap_without_input.pdf")), 
                width = 10/2.54, height = 8/2.54)
            draw(results$heatmap)
            dev.off()
          }
          
          # Save dot plot
          ggsave(file.path(celltype_dir, paste0(celltype, "_", hist, "_", tool, "_dot_plot_without_input.pdf")),
                 plot = results$dot_plot, 
                 width = 10, height = 6, units = "cm", dpi = 1200, device = cairo_pdf)
          
          # Save individual EID plots
          for(eid in names(results$individual_plots)) {
            ggsave(file.path(celltype_dir, paste0(celltype, "_", eid, "_", hist, "_", tool, "_barplot_without_input.pdf")),
                   plot = results$individual_plots[[eid]], 
                   width = 8, height = 6, units = "cm", dpi = 1200, device = cairo_pdf)
          }
          
          # Save data with metadata
          write_tsv(results$data, file.path(celltype_dir, paste0(celltype, "_", hist, "_", tool, "_data_without_input.tsv")))
          write_tsv(results$summary_stats, file.path(celltype_dir, paste0(celltype, "_", tool, "_summary_stats_without_input.tsv")))
          
          cat("✓ Saved figures for", celltype, "-", results$n_samples, "samples\n")
        }
      } else {
        cat("✗ No data for cell type:", celltype, "\n")
      }
    }
    
    # Save session information for individual tool
    sink(file.path(output_dir, paste0("analysis_session_info_", hist, "_", tool, "_without_input.txt")))
    print(sessionInfo())
    sink()
    
    cat("✓ Individual analysis completed successfully for", tool, "!\n")
  }
  
  # =============================================================================
  # CREATE COMBINED GRID FIGURE (8 cell types × 7 tools)
  # =============================================================================
  
  if (length(tools_found) > 0) {
    cat("\nCreating combined grid figure...\n")
    cat("Tools found:", paste(tools_found, collapse = ", "), "\n")
    
    all_plots <- list()
    
    # Create plots for grid
    for(tool in tools_found) {
      tool_data <- all_grid_data[[tool]]
      
      for(celltype in celltypes) {
        celltype_data <- tool_data %>% filter(cell_type == celltype)
        
        if(nrow(celltype_data) > 0) {
          plot_name <- paste0(celltype, "_", tool)
          all_plots[[plot_name]] <- create_grid_barplot(celltype_data, celltype, tool)
          cat("✓ Created grid plot for", celltype, "-", tool, "\n")
        } else {
          # Create empty plot for missing combinations
          empty_plot <- ggplot() +
            geom_blank() +
            labs(
              title = paste0(celltype, " - ", tool),
              x = NULL,
              y = NULL
            ) +
            theme_classic(base_size = 6) +
            theme(
              plot.title = element_text(size = 7, face = "bold", hjust = 0.5, color = "gray50"),
              panel.background = element_rect(fill = "grey95")
            )
          
          plot_name <- paste0(celltype, "_", tool)
          all_plots[[plot_name]] <- empty_plot
          cat("○ Empty plot for", celltype, "-", tool, "(no data)\n")
        }
      }
    }
    
    # Create the grid layout: 8 cell types (rows) × tools (columns)
    n_celltypes <- length(celltypes)
    n_tools <- length(tools_found)
    
    cat("Creating grid with", n_celltypes, "cell types (rows) and", n_tools, "tools (columns)\n")
    
    # Order plots properly: cell types as rows, tools as columns
    ordered_plots <- list()
    
    for(celltype in celltypes) {
      for(tool in tools_found) {
        plot_name <- paste0(celltype, "_", tool)
        if (plot_name %in% names(all_plots)) {
          ordered_plots[[plot_name]] <- all_plots[[plot_name]]
        }
      }
    }
    
    # Create the grid arrangement
    grid_plot <- wrap_plots(
      ordered_plots,
      nrow = n_celltypes,
      ncol = n_tools,
      byrow = TRUE
    )
    
    # Create comprehensive legend
    state_mapping <- data.frame(
      state = names(state_descriptions),
      description = unname(state_descriptions)
    )
    
    legend_plot <- ggplot(state_mapping, aes(x = state, y = 1, fill = state)) +
      geom_bar(stat = "identity", width = 0.5) +
      scale_fill_manual(
        values = chromatin_colors, 
        name = "Chromatin State",
        labels = state_descriptions,
        breaks = names(state_descriptions)
      ) +
      theme_void() +
      theme(
        legend.position = "top",
        legend.text = element_text(size = 5, margin = margin(r = 6)),
        legend.title = element_text(size = 6, face = "bold", margin = margin(b = 2)),
        legend.key.size = unit(0.2, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, "cm")
      ) +
      guides(fill = guide_legend(
        nrow = 3,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5
      ))
    
    legend <- cowplot::get_legend(legend_plot)
    
    # Final composition with title and legend
    final_grid <- wrap_plots(
      # Title
      wrap_elements(grid::textGrob(
        paste("Chromatin State Distribution -", hist, "Histone"),
        gp = grid::gpar(fontsize = 14, fontface = "bold", margin = margin(b = 10))
      )),
      # Grid of plots
      grid_plot,
      # Legend
      legend,
      ncol = 1,
      heights = c(0.05, 0.85, 0.10)
    )
    
    # Save as A4 PDF
    output_filename <- paste0("COMBINED_GRID_", hist, "_", n_celltypes, "x", n_tools, "_A4_FIGURE.pdf")
    ggsave(file.path(main_output_dir, output_filename),
           plot = final_grid,
           width = 29.7, height = 21, units = "cm",
           dpi = 1200, device = cairo_pdf)
    
    cat("✓ Saved combined grid figure:", output_filename, "\n")
    cat("✓ Grid dimensions:", n_celltypes, "cell types (rows) ×", n_tools, "tools (columns)\n")
    
  } else {
    cat("✗ No tools with data found for", hist, "- skipping grid figure\n")
  }
}

cat("\n", strrep("=", 80), "\n", sep = "")
cat("ALL PROCESSING COMPLETED SUCCESSFULLY!\n")
cat(strrep("=", 80), "\n", sep = "")



###############################################################################################################################
## Figure 6 B-C
###############################################################################################################################


library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# =============================================================================
# PATHS
# =============================================================================
output_base <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_ChromHMM_roadmapdata/Figure_without_input/"

# Load the combined data
combined_file <- paste0(output_base, "COMBINED_DATA/ALL_COMBINED_RAW_DATA_without_input.rds")

if(file.exists(combined_file)) {
  cat("Loading combined data from:", combined_file, "\n")
  combined_data <- readRDS(combined_file)
} else {
  combined_file <- paste0(output_base, "COMBINED_DATA/ALL_COMBINED_RAW_DATA_without_input.tsv")
  if(file.exists(combined_file)) {
    cat("Loading combined data from:", combined_file, "\n")
    combined_data <- read_tsv(combined_file, show_col_types = FALSE)
  } else {
    stop("No combined data file found! Please run the data saving script first.")
  }
}

# =============================================================================
# PREPARE SCATTER PLOT DATA WITH CELL TYPES AS DOTS
# =============================================================================

cat("Preparing scatter plot data with cell types as dots...\n")

scatter_data <- combined_data %>%
  mutate(
    state_group = case_when(
      state %in% c("1_TssA", "2_TssAFlnk") ~ "TSS",
      state %in% c("6_EnhG", "7_Enh") ~ "Enhancers",
      TRUE ~ "Other"
    )
  ) %>%
  filter(state_group %in% c("TSS", "Enhancers")) %>%
  group_by(filename, sample_id, cell_type, eid, tool, histone, state_group) %>%
  summarise(
    total_fraction = sum(fraction, na.rm = TRUE),
    total_peak_count = sum(peak_count, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = state_group,
    values_from = c(total_fraction, total_peak_count),
    values_fill = 0
  ) %>%
  rename(
    TSS_fraction = total_fraction_TSS,
    Enhancers_fraction = total_fraction_Enhancers,
    TSS_peak_count = total_peak_count_TSS,
    Enhancers_peak_count = total_peak_count_Enhancers
  ) %>%
  mutate(
    TSS_log10 = log10(TSS_peak_count + 1),
    Enhancers_log10 = log10(Enhancers_peak_count + 1),
    sample_label = paste(cell_type, eid, tool, histone, sep = "_"),
    # Rename columns for clarity
    Method = tool,
    Histone = histone,
    CellType = cell_type
  )

# Aggregate by cell type, tool, and histone (average across samples if multiple)
scatter_agg <- scatter_data %>%
  group_by(Method, Histone) %>%
  summarise(
    median_TSS = median(TSS_fraction, na.rm = TRUE),
    median_Enhancers = median(Enhancers_fraction, na.rm = TRUE),
    sd_TSS = sd(TSS_fraction, na.rm = TRUE),
    sd_Enhancers = sd(Enhancers_fraction, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    # Create label with sample count
    Method = factor(Method),
    Histone = factor(Histone)
  )

cat("✓ Cell types in data:", paste(unique(scatter_agg$CellType), collapse = ", "), "\n")
cat("✓ Total data points after aggregation:", nrow(scatter_agg), "\n")

# =============================================================================
# NATURE METHODS STYLE SETUP
# =============================================================================

# Get all unique methods (tools) for color mapping
all_methods <- unique(scatter_agg$Method)

# Nature methods color palette
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)

# Assign colors to methods
if (length(all_methods) <= length(nature_methods_colors)) {
  method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)
} else {
  extra_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(all_methods) - length(nature_methods_colors))
  method_colors <- setNames(c(nature_methods_colors, extra_colors)[1:length(all_methods)], all_methods)
}

# Define shapes for histones
histone_names <- unique(scatter_agg$Histone)
histone_shape_map <- setNames(c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10), histone_names)

# Enhanced Nature Methods Theme
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

# =============================================================================
# CREATE SCATTER PLOTS WITH CELL TYPES AS DOTS
# =============================================================================

# Create output directory
plot_dir <- paste0(output_base, "MEDIAN_SCATTER_PLOTS/")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nCreating scatter plots with cell types as dots...\n")

# Check the range of values
cat("Data range:\n")
cat(sprintf("TSS Fraction range: %.3f to %.3f\n", 
            min(scatter_agg$median_TSS), max(scatter_agg$median_TSS)))
cat(sprintf("Enhancers Fraction range: %.3f to %.3f\n", 
            min(scatter_agg$median_Enhancers), max(scatter_agg$median_Enhancers)))

# =============================================================================
# PLOT 1: MAIN SCATTER PLOT - Cell Types as Dots (NO LABELS)
# =============================================================================

main_scatter_plot <- ggplot(scatter_agg, 
                            aes(x = median_TSS,
                                y = median_Enhancers,
                                color = Method,
                                shape = Histone)) +
  geom_point(size = 10, alpha = 0.85, stroke = 1.2) +
  scale_color_manual(values = method_colors, 
                     name = "Peak-Calling Method") +
  scale_shape_manual(values = histone_shape_map,
                     name = "Histone Modification") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, .4),
    breaks = seq(0, .4, by = 0.1),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title = "TSS vs Enhancers Distribution",
    x = "Median TSS Fraction (Active TSS + Flanking Active TSS)",
    y = "Median Enhancers Fraction (Genic enhancers + Enhancers)"
  ) +
  nature_methods_enhanced_theme() +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

# Save plots in both aspect ratios
ggsave(
  file.path(plot_dir, "CellType_Dots_TSS_vs_Enhancers_1_1.pdf"),
  plot = main_scatter_plot,
  width = 10, height = 10, units = "in", dpi = 300  # 1:1 ratio
)
cat("✓ Created: CellType_Dots_TSS_vs_Enhancers_1_1.pdf (1:1 ratio)\n")

ggsave(
  file.path(plot_dir, "CellType_Dots_TSS_vs_Enhancers_4_3.pdf"),
  plot = main_scatter_plot,
  width = 12, height = 9, units = "in", dpi = 300  # 4:3 ratio
)
cat("✓ Created: CellType_Dots_TSS_vs_Enhancers_4_3.pdf (4:3 ratio)\n")
print(main_scatter_plot)

###==========================================================================================
### SCATTER BAR PLOT 1 - TSS Median Fraction
### Histone on X-axis, Method as Color, Points = Median TSS Fraction with Error Bars
###==========================================================================================

cat("\n================================================================================\n")
cat("CREATING SCATTER BAR PLOT: TSS MEDIAN FRACTION\n")
cat("================================================================================\n")

# =============================================================================
# PREPARE TSS SCATTER BAR DATA (Histone on X-axis, Method as color)
# =============================================================================

tss_bar_data <- combined_data %>%
  mutate(
    state_group = case_when(
      state %in% c("1_TssA", "2_TssAFlnk") ~ "TSS",
      TRUE ~ "Other"
    )
  ) %>%
  filter(state_group == "TSS") %>%
  # First aggregate by sample to get TSS fraction per sample
  group_by(filename, sample_id, cell_type, eid, tool, histone) %>%
  summarise(
    tss_fraction = sum(fraction, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Then calculate statistics per tool and histone
  group_by(tool, histone) %>%
  summarise(
    median_tss = median(tss_fraction, na.rm = TRUE),
    mean_tss = mean(tss_fraction, na.rm = TRUE),
    sd_tss = sd(tss_fraction, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    # Calculate standard error and 95% confidence intervals
    se_tss = sd_tss / sqrt(n_samples),
    ci_lower = mean_tss - 1.96 * se_tss,
    ci_upper = mean_tss + 1.96 * se_tss,
    # Ensure CI bounds are within [0,1] for fractions
    ci_lower = pmax(0, ci_lower),
    ci_upper = pmin(1, ci_upper),
    # For error bars using SD (more common in biological sciences)
    sd_lower = pmax(0, mean_tss - sd_tss),
    sd_upper = pmin(1, mean_tss + sd_tss),
    # Formatting for plot
    Method = factor(tool, levels = sort(unique(tool))),
    Histone = factor(histone, 
                     levels = c("H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", 
                               "H3K4me2", "H3K4me3", "H3K9me3"))
  ) %>%
  # Remove any NA or invalid values
  filter(!is.na(median_tss), !is.na(sd_tss))

cat("✓ TSS data prepared with", nrow(tss_bar_data), "data points\n")
cat("  - Methods:", paste(unique(tss_bar_data$Method), collapse = ", "), "\n")
cat("  - Histones:", paste(unique(tss_bar_data$Histone), collapse = ", "), "\n")
cat("  - TSS fraction range:", 
    round(min(tss_bar_data$median_tss), 3), "to", 
    round(max(tss_bar_data$median_tss), 3), "\n\n")

# =============================================================================
# TSS SCATTER BAR PLOT - Using Median with SD error bars
# =============================================================================

# Calculate y-axis limits
tss_y_max <- max(tss_bar_data$sd_upper, na.rm = TRUE) * 1.15
tss_y_min <- max(0, min(tss_bar_data$sd_lower, na.rm = TRUE) * 0.95)

tss_scatter_bar <- ggplot(tss_bar_data, 
                         aes(x = Histone, 
                             y = median_tss,
                             color = Method)) +
  # Error bars (SD - most common in biology)
  geom_errorbar(
    aes(ymin = sd_lower,
        ymax = sd_upper),
    width = 0.25,
    position = position_dodge(width = 0.7),
    size = 0.8,
    alpha = 0.7
  ) +
  # Points (median values)
  geom_point(
    size = 4,
    position = position_dodge(width = 0.7),
    stroke = 1.2
  ) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, tss_y_max),
    expand = expansion(mult = c(0.02, 0.05)),
    breaks = seq(0, round(tss_y_max, 1), by = 0.1)
  ) +
  labs(
    title = "TSS Region Peak Distribution",
    subtitle = "Median fraction of peaks overlapping TSS regions (TssA + TssAFlnk) ± SD",
    x = "Histone Modification",
    y = "Median TSS Fraction"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, 
                              margin = margin(b = 8), color = "#2C3E50"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "#566573",
                                 margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 14, color = "#2C3E50"),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.text = element_text(size = 12, color = "#34495E"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.line = element_line(color = "#2C3E50", linewidth = 0.6),
    axis.ticks = element_line(color = "#2C3E50", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    legend.title = element_text(face = "bold", size = 12, color = "#2C3E50"),
    legend.text = element_text(size = 11, color = "#34495E"),
    legend.position = "right",
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(t = 0, b = 0, l = 10, r = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )


# =============================================================================
# SCATTER BAR PLOT 2 - MEDIAN ENHANCERS FRACTION
# =============================================================================

cat("\n================================================================================\n")
cat("CREATING SCATTER BAR PLOT: MEDIAN ENHANCERS FRACTION\n")
cat("================================================================================\n")

# =============================================================================
# PREPARE ENHANCERS SCATTER BAR DATA
# =============================================================================

enhancers_bar_data <- combined_data %>%
  mutate(
    state_group = case_when(
      state %in% c("6_EnhG", "7_Enh") ~ "Enhancers",
      TRUE ~ "Other"
    )
  ) %>%
  filter(state_group == "Enhancers") %>%
  # First aggregate by sample
  group_by(filename, sample_id, cell_type, eid, tool, histone) %>%
  summarise(
    enh_fraction = sum(fraction, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Then calculate statistics per tool and histone
  group_by(tool, histone) %>%
  summarise(
    median_enh = median(enh_fraction, na.rm = TRUE),
    mean_enh = mean(enh_fraction, na.rm = TRUE),
    sd_enh = sd(enh_fraction, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    # Calculate error bars
    se_enh = sd_enh / sqrt(n_samples),
    ci_lower = mean_enh - 1.96 * se_enh,
    ci_upper = mean_enh + 1.96 * se_enh,
    ci_lower = pmax(0, ci_lower),
    ci_upper = pmin(1, ci_upper),
    sd_lower = pmax(0, mean_enh - sd_enh),
    sd_upper = pmin(1, mean_enh + sd_enh),
    # Formatting
    Method = factor(tool, levels = sort(unique(tool))),
    Histone = factor(histone, 
                     levels = c("H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", 
                               "H3K4me2", "H3K4me3", "H3K9me3"))
  ) %>%
  filter(!is.na(median_enh), !is.na(sd_enh))

cat("✓ Enhancers data prepared with", nrow(enhancers_bar_data), "data points\n")
cat("  - Enhancers fraction range:", 
    round(min(enhancers_bar_data$median_enh), 3), "to", 
    round(max(enhancers_bar_data$median_enh), 3), "\n\n")

# =============================================================================
# ENHANCERS SCATTER BAR PLOT - Using Median with SD error bars
# =============================================================================

enh_y_max <- max(enhancers_bar_data$sd_upper, na.rm = TRUE) * 1.15

enhancers_scatter_bar <- ggplot(enhancers_bar_data, 
                               aes(x = Histone, 
                                   y = median_enh,
                                   color = Method)) +
  geom_errorbar(
    aes(ymin = sd_lower,
        ymax = sd_upper),
    width = 0.25,
    position = position_dodge(width = 0.7),
    size = 0.8,
    alpha = 0.7
  ) +
  geom_point(
    size = 4,
    position = position_dodge(width = 0.7),
    stroke = 1.2
  ) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, enh_y_max),
    expand = expansion(mult = c(0.02, 0.05)),
    breaks = seq(0, round(enh_y_max, 1), by = 0.05)
  ) +
  labs(
    title = "Enhancer Region Peak Distribution",
    x = "Histone Modification",
    y = "Median Enhancers Fraction"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, 
                              margin = margin(b = 8), color = "#2C3E50"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "#566573",
                                 margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 14, color = "#2C3E50"),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.text = element_text(size = 12, color = "#34495E"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.line = element_line(color = "#2C3E50", linewidth = 0.6),
    axis.ticks = element_line(color = "#2C3E50", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    legend.title = element_text(face = "bold", size = 12, color = "#2C3E50"),
    legend.text = element_text(size = 11, color = "#34495E"),
    legend.position = "right",
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(t = 0, b = 0, l = 10, r = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

# =============================================================================
# SAVE ALL SCATTER BAR PLOTS
# =============================================================================

# Create output directory for scatter bar plots
scatter_bar_dir <- paste0(output_base, "SCATTER_BAR_PLOTS/")
dir.create(scatter_bar_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n================================================================================\n")
cat("SAVING SCATTER BAR PLOTS\n")
cat("================================================================================\n")

# 1. TSS Median with SD error bars (RECOMMENDED)
ggsave(
  file.path(scatter_bar_dir, "ScatterBar_TSS_Median_with_SD_without_input.pdf"),
  plot = tss_scatter_bar,
  width = 4, height =2, units = "in", dpi = 300
)
ggsave(
  file.path(scatter_bar_dir, "ScatterBar_TSS_Median_with_SD_without_input.png"),
  plot = tss_scatter_bar,
  width = 12, height = 7, units = "in", dpi = 300, bg = "white"
)
cat("✓ Created: ScatterBar_TSS_Median_with_SD_without_input.pdf/png\n")


# 3. Enhancers Median with SD error bars (RECOMMENDED)
ggsave(
  file.path(scatter_bar_dir, "ScatterBar_Enhancers_Median_with_SD_without_input.pdf"),
  plot = enhancers_scatter_bar,
  width = 12, height = 7, units = "in", dpi = 300
)
ggsave(
  file.path(scatter_bar_dir, "ScatterBar_Enhancers_Median_with_SD_without_input.png"),
  plot = enhancers_scatter_bar,
  width = 12, height = 7, units = "in", dpi = 300, bg = "white"
)
cat("✓ Created: ScatterBar_Enhancers_Median_with_SD.pdf/png\n")



###############################################################################################################################
## Figure 6 H
###############################################################################################################################

#################################################################################################################
## COMPLETE MOTIF ANALYSIS PIPELINE WITH HORIZONTAL BAR PLOT AND ADDITIONAL UNIFIED FIGURE
## November 2024
#################################################################################################################

## =============================================================================
## 1. LOAD REQUIRED PACKAGES
## =============================================================================

# Required packages
required_packages <- c(
  "motifmatchr", "BSgenome.Hsapiens.UCSC.hg38", "TFBSTools", "JASPAR2020",
  "GenomicRanges", "rtracklayer", "dplyr", "ggplot2", "pheatmap", "viridis",
  "ggrepel", "forcats", "tidyr", "stringr", "purrr", "readr", "tibble",
  "RColorBrewer", "Cairo", "matrixStats", "data.table", "gridExtra", "patchwork",
  "scales", "BiocGenerics", "Biostrings", "IRanges", "GenomeInfoDb", "cowplot"
)

# Install and load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

## =============================================================================
## 2. DEFINE PARAMETERS
## =============================================================================

# Analysis parameters
histones <- c("H3K27ac-s")  # Start with one histone
tools <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
cell_types <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")

# Create output directory
output_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Subdirectories
for (subdir in c("individual_results", "comparison_plots", "summary_tables", 
                 "diagnostics", "figures", "pooled_motifs", "pooled_motifs_by_tool")) {
  dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

## =============================================================================
## 3. FIXED BED FILE FINDING FUNCTION
## =============================================================================

find_bed_file <- function(tool, histone, cell_type) {
  # Base directory
  base_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_with_input_peakbed_corrected/motif_peakbed"
  
  # Tool-specific directory
  tool_dir <- file.path(base_dir, paste0(tool, "_peakbed"))
  
  if (!dir.exists(tool_dir)) {
    cat("  ✗ Directory not found:", tool_dir, "\n")
    return(NULL)
  }
  
  # List all BED files in the directory
  all_files <- list.files(tool_dir, pattern = "\\.bed$", full.names = FALSE)
  
  if (length(all_files) == 0) {
    cat("  ✗ No BED files in directory\n")
    return(NULL)
  }
  
  # Try multiple patterns (flexible matching)
  patterns_to_try <- c(
    # Exact match: tool_histone_cell.bed
    paste0("^", tool, "_", histone, "_", cell_type, "\\.bed$"),
    
    # With variations in histone name
    paste0("^", tool, "_", gsub("-", "_", histone), "_", cell_type, "\\.bed$"),
    paste0("^", tool, "_.*", cell_type, "\\.bed$"),
    
    # Just cell type pattern
    paste0(".*", cell_type, "\\.bed$")
  )
  
  for (pattern in patterns_to_try) {
    matched_files <- grep(pattern, all_files, value = TRUE, ignore.case = TRUE)
    if (length(matched_files) > 0) {
      cat("  ✓ Found:", matched_files[1], "\n")
      return(file.path(tool_dir, matched_files[1]))
    }
  }
  
  cat("  ✗ No matching file found for", tool, histone, cell_type, "\n")
  cat("    Available files:", paste(all_files, collapse=", "), "\n")
  return(NULL)
}

## =============================================================================
## 4. FIXED PEAK LOADING FUNCTION
## =============================================================================

load_peaks_from_bed <- function(bed_file, cell_type) {
  if (!file.exists(bed_file)) {
    cat("    ✗ File does not exist:", bed_file, "\n")
    return(GRanges())
  }
  
  tryCatch({
    # Read BED file
    peaks_data <- fread(bed_file, header = FALSE)
    
    if (nrow(peaks_data) == 0) {
      cat("    ✗ Empty file\n")
      return(GRanges())
    }
    
    # Minimum columns check
    if (ncol(peaks_data) < 3) {
      cat("    ✗ Invalid BED format\n")
      return(GRanges())
    }
    
    # Create GRanges (BED is 0-based, GRanges is 1-based)
    peaks <- GRanges(
      seqnames = peaks_data$V1,
      ranges = IRanges(start = peaks_data$V2 + 1, end = peaks_data$V3),
      strand = "*"
    )
    
    # Keep only standard chromosomes
    standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
    peaks <- peaks[seqnames(peaks) %in% standard_chroms]
    
    # Remove duplicates and invalid ranges
    peaks <- unique(peaks)
    peaks <- peaks[width(peaks) > 0]
    
    # Trim to chromosome boundaries
    peaks <- trim(peaks)
    
    cat("    ✓ Successfully loaded", length(peaks), "peaks\n")
    
    # Diagnostic: Show some peaks
    if (length(peaks) > 0) {
      cat("    First peak:", as.character(peaks[1]), "\n")
    }
    
    return(peaks)
    
  }, error = function(e) {
    cat("    ✗ Error loading file:", e$message, "\n")
    return(GRanges())
  })
}

## =============================================================================
## 5. BUILD MOTIF DATABASE
## =============================================================================

build_motif_database <- function() {
  cat("Building motif database from JASPAR2020...\n")
  
  opts <- list(
    species = "Homo sapiens",
    collection = "CORE",
    all_versions = FALSE
  )
  
  motifs <- getMatrixSet(JASPAR2020, opts)
  
  pwm_list <- PWMatrixList()
  motif_names <- character()
  
  for (i in seq_along(motifs)) {
    motif_name <- name(motifs[[i]])
    
    tryCatch({
      pwm <- toPWM(motifs[[i]])
      pwm_list[[motif_name]] <- pwm
      motif_names <- c(motif_names, motif_name)
    }, error = function(e) {
      # Skip problematic motifs
    })
  }
  
  cat("✓ Loaded", length(pwm_list), "motifs\n")
  
  # Save motif list
  motif_info <- data.frame(
    Name = motif_names,
    Index = seq_along(motif_names),
    stringsAsFactors = FALSE
  )
  
  write.csv(motif_info,
            file.path(output_dir, "summary_tables", "JASPAR_motifs_with_input.csv"),
            row.names = FALSE)
  
  return(pwm_list)
}

## =============================================================================
## 6. SEQUENCE EXTRACTION FUNCTION
## =============================================================================

extract_peak_sequences <- function(peaks, genome, width = 200) {
  if (length(peaks) == 0) {
    return(DNAStringSet())
  }
  
  # Center peaks and resize
  peak_centers <- resize(peaks, width = 1, fix = "center")
  peak_regions <- resize(peak_centers, width = width, fix = "center")
  
  # Ensure within chromosome boundaries
  peak_regions <- trim(peak_regions)
  
  # Extract sequences
  sequences <- getSeq(genome, peak_regions)
  
  return(sequences)
}

## ===================================================
## 7. GENERATE BACKGROUND SEQUENCES
## ===================================================
generate_background_sequences <- function(genome, n_regions = 1000, width = 200, 
                                          exclude_regions = GRanges()) {
  cat("  Generating background sequences...\n")
  
  chrom_sizes <- seqlengths(genome)
  chrom_sizes <- chrom_sizes[!is.na(chrom_sizes)]
  
  # Keep only standard chromosomes
  standard_chroms <- paste0("chr", c(1:22, "X", "Y"))
  chrom_sizes <- chrom_sizes[names(chrom_sizes) %in% standard_chroms]
  
  set.seed(123)  # For reproducibility
  
  background_regions <- GRanges()
  attempts <- 0
  max_attempts <- n_regions * 5
  
  while (length(background_regions) < n_regions && attempts < max_attempts) {
    attempts <- attempts + 1
    
    # Random chromosome
    chrom <- sample(names(chrom_sizes), 1)
    chrom_len <- chrom_sizes[chrom]
    
    if (chrom_len > width + 200) {
      # Random start position (leave buffer)
      start_pos <- sample(100:(chrom_len - width - 100), 1)
      
      region <- GRanges(
        seqnames = chrom,
        ranges = IRanges(start = start_pos, width = width),
        strand = "*"
      )
      
      # Check if it overlaps with any peak regions
      overlaps <- countOverlaps(region, exclude_regions)
      
      if (overlaps == 0) {
        background_regions <- c(background_regions, region)
      }
    }
  }
  
  # Extract sequences
  background_sequences <- getSeq(genome, background_regions)
  
  cat("  Generated", length(background_sequences), "background sequences\n")
  
  return(background_sequences)
}

## =================================
## 8. CORE MOTIF ENRICHMENT ANALYSIS
## =================================

analyze_motifs_for_sample <- function(peaks, tool_name, histone_name, cell_type_name,
                                     motif_db, n_foreground = 500, n_background = 1000) {
  
  sample_id <- paste(tool_name, histone_name, cell_type_name, sep = "_")
  cat("\n    Analyzing:", sample_id, "\n")
  
  # Diagnostic information
  cat("    Number of input peaks:", length(peaks), "\n")
  
  if (length(peaks) < 50) {
    cat("    ⚠ Too few peaks (<50). Skipping...\n")
    return(data.frame())
  }
  
  ## STEP 1: Extract foreground sequences
  cat("    Extracting foreground sequences...\n")
  genome <- BSgenome.Hsapiens.UCSC.hg38
  
  # Use all peaks or subset if too many
  if (length(peaks) > n_foreground) {
    set.seed(123)
    peaks <- sample(peaks, n_foreground)
  }
  
  foreground_sequences <- extract_peak_sequences(peaks, genome, width = 200)
  
  if (length(foreground_sequences) == 0) {
    cat("    ✗ Failed to extract foreground sequences\n")
    return(data.frame())
  }
  
  cat("    Extracted", length(foreground_sequences), "foreground sequences\n")
  
  ## STEP 2: Generate background sequences
  background_sequences <- generate_background_sequences(
    genome, 
    n_regions = n_background,
    width = 200,
    exclude_regions = peaks
  )
  
  if (length(background_sequences) == 0) {
    cat("    ✗ Failed to generate background sequences\n")
    return(data.frame())
  }
  
  ## STEP 3: Match motifs in both sets
  cat("    Matching motifs...\n")
  
  # Match in foreground
  foreground_matches <- matchMotifs(
    motif_db, 
    foreground_sequences,
    genome = genome,
    out = "scores",
    p.cutoff = 1e-3
  )
  
  foreground_scores <- motifScores(foreground_matches)
  foreground_hits <- foreground_scores > 0  # Binary matrix
  
  # Match in background
  background_matches <- matchMotifs(
    motif_db,
    background_sequences,
    genome = genome,
    out = "scores",
    p.cutoff = 1e-3
  )
  
  background_scores <- motifScores(background_matches)
  background_hits <- background_scores > 0
  
  ## STEP 4: Calculate enrichment for each motif
  cat("    Calculating enrichment...\n")
  
  results <- data.frame()
  
  for (motif_idx in seq_len(ncol(foreground_hits))) {
    motif_name <- colnames(foreground_hits)[motif_idx]
    
    # Foreground statistics
    fore_hits <- sum(foreground_hits[, motif_idx], na.rm = TRUE)
    fore_total <- nrow(foreground_hits)
    fore_freq <- ifelse(fore_total > 0, fore_hits / fore_total, 0)
    
    # Background statistics
    back_hits <- sum(background_hits[, motif_idx], na.rm = TRUE)
    back_total <- nrow(background_hits)
    back_freq <- ifelse(back_total > 0, back_hits / back_total, 0)
    
    # Skip if no hits in either set
    if (fore_hits == 0 && back_hits == 0) next
    
    # Calculate enrichment
    if (back_freq > 0) {
      enrichment <- fore_freq / back_freq
      log2_enrichment <- log2(enrichment + 0.001)
    } else {
      enrichment <- 100  # Large number if no background hits
      log2_enrichment <- log2(100.001)
    }
    
    # Fisher's exact test
    if (fore_total > 0 && back_total > 0) {
      contingency <- matrix(c(fore_hits, fore_total - fore_hits,
                             back_hits, back_total - back_hits),
                           nrow = 2)
      
      fisher_test <- try(fisher.test(contingency), silent = TRUE)
      if (!inherits(fisher_test, "try-error")) {
        p_value <- fisher_test$p.value
      } else {
        p_value <- 1.0
      }
    } else {
      p_value <- 1.0
    }
    
    # Store results
    result_row <- data.frame(
      tool = tool_name,
      histone = histone_name,
      cell_type = cell_type_name,
      motif = motif_name,
      foreground_hits = fore_hits,
      foreground_total = fore_total,
      foreground_frequency = fore_freq,  # This is the frequency
      background_hits = back_hits,
      background_total = back_total,
      background_frequency = back_freq,
      enrichment_ratio = enrichment,
      log2_enrichment = log2_enrichment,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
    
    results <- rbind(results, result_row)
  }
  
  ## STEP 5: Statistical adjustments and ranking
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(
        p_adjusted = p.adjust(p_value, method = "BH"),
        significant = p_adjusted < 0.05,
        significance_star = case_when(
          p_adjusted < 0.001 ~ "***",
          p_adjusted < 0.01 ~ "**",
          p_adjusted < 0.05 ~ "*",
          TRUE ~ ""
        ),
        rank = rank(-log2_enrichment, ties.method = "first")
      ) %>%
      arrange(desc(log2_enrichment))
    
    cat("    ✓ Analyzed", nrow(results), "motifs\n")
    if (nrow(results) > 0) {
      cat("    Top motif:", results$motif[1], 
          "(log2Enrich =", round(results$log2_enrichment[1], 2), ")\n")
    }
  } else {
    cat("    ✗ No motifs analyzed\n")
  }
  
  return(results)
}

## =============================================================================
## 9. SAVE INDIVIDUAL RESULTS
## =============================================================================

save_individual_results <- function(results, tool, histone, cell_type) {
  if (nrow(results) == 0) return()
  
  # FIXED: Use 'histone' parameter instead of undefined 'hist'
  filename <- paste0(tool, "_", histone, "_", cell_type, "_motif_results_with_input.csv")
  filepath <- file.path(output_dir, "individual_results", filename)
  
  # Save full results
  write.csv(results, filepath, row.names = FALSE)
  
  # Save top 25 motifs
  top_25 <- results %>%
    filter(rank <= 25) %>%
    select(tool, histone, cell_type, motif, log2_enrichment, p_adjusted, rank)
  
  # FIXED: Use 'histone' parameter
  top_filename <- paste0(tool, "_", histone, "_", cell_type, "_top25_with_input.csv")
  top_filepath <- file.path(output_dir, "individual_results", top_filename)
  
  write.csv(top_25, top_filepath, row.names = FALSE)
}
## =============================================================================
## 10. SHOW TOP MOTIFS
## =============================================================================

show_top_motifs <- function(results, tool, cell_type) {
  if (nrow(results) == 0) return()
  
  top_5 <- results %>%
    arrange(desc(log2_enrichment)) %>%
    head(5)
  
  cat("    Top 5 motifs for", tool, "-", cell_type, ":\n")
  for (i in seq_len(nrow(top_5))) {
    cat(sprintf("      %d. %s (%.2f)\n", 
                i, top_5$motif[i], top_5$log2_enrichment[i]))
  }
}

## =============================================================================
## 11. HORIZONTAL BAR GRID PLOT WITH MOTIF NAMES ON Y-AXIS
## =============================================================================

create_horizontal_bar_grid_plot <- function(all_results, output_dir) {
  cat("\n================================================================\n")
  cat("CREATING HORIZONTAL BAR GRID PLOT\n")
  cat("Row: Cell types | Column: Tools\n")
  cat("X-axis: Log2 Enrichment | Y-axis: Motif Names (Rank 1 at top)\n")
  cat("================================================================\n")
  
  # Get top 25 motifs for each tool-cell combination
  top25_data <- all_results %>%
    group_by(cell_type, tool) %>%
    arrange(desc(log2_enrichment)) %>%
    mutate(
      rank_within = row_number()
    ) %>%
    filter(rank_within <= 25) %>%
    ungroup() %>%
    mutate(
      # Order motifs by rank (highest rank = highest position on Y-axis)
      motif = factor(motif, levels = unique(motif[order(-rank_within)]))
    )
  
  # Check data
  cat("\nData summary:\n")
  cat("Unique cell types:", n_distinct(top25_data$cell_type), "\n")
  cat("Unique tools:", n_distinct(top25_data$tool), "\n")
  cat("Total motifs in top 25:", nrow(top25_data), "\n")
  
  if (nrow(top25_data) == 0) {
    cat("✗ No data for grid plot\n")
    return(NULL)
  }
  
  # Create the horizontal bar grid plot with motif names on Y-axis
  horizontal_plot <- ggplot(top25_data, 
                           aes(x = log2_enrichment, 
                               y = reorder(motif, rank_within))) +  # Rank 1 at top
    
    # Create horizontal bars
    geom_bar(stat = "identity", aes(fill = tool), alpha = 0.8, width = 0.7) +
    
    # Add value labels at the end of bars
    geom_text(aes(label = round(log2_enrichment, 2)), 
              hjust = -0.1, size = 2.5, color = "black") +
    
    # Add rank labels on Y-axis
    geom_text(aes(x = 0, y = reorder(motif, rank_within), 
                  label = paste0(rank_within, ". ", motif)),
              hjust = 1.1, size = 2.5, color = "black") +
    
    # Add vertical line at x=0 for reference
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    
    # Create facets: cell types in rows, tools in columns
    facet_grid(
      rows = vars(cell_type), 
      cols = vars(tool),
      scales = "free",
      space = "free"
    ) +
    
    # Color scheme
    scale_fill_brewer(palette = "Set2") +
    
    # Labels and titles
    labs(
      title = "Top 25 Motifs by Cell Type and Tool",
      subtitle = "Horizontal bar plot: X-axis = Log2 Enrichment, Y-axis = Motif Names",
      x = "Log2 Enrichment Score",
      y = "Motif (Ranked 1-25)",
      fill = "Tool"
    ) +
    
    # Theme settings
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 15)),
      
      # Facet labels
      strip.text.x = element_text(face = "bold", size = 8, margin = margin(b = 5)),
      strip.text.y = element_text(face = "bold", size = 8, angle = 0, 
                                 margin = margin(r = 5)),
      strip.background.x = element_rect(fill = "gray95", color = NA),
      strip.background.y = element_rect(fill = "gray95", color = NA),
      
      # Axis settings - hide default Y-axis text since we have custom labels
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      axis.title = element_text(size = 9),
      
      # Grid and panel
      panel.spacing = unit(0.8, "lines"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.major.y = element_blank(),
      
      # Legend
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    ) +
    
    # Adjust x-axis to accommodate labels
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)))
  
  # Calculate dimensions
  n_rows <- length(unique(top25_data$cell_type))
  n_cols <- length(unique(top25_data$tool))
  
  # Save the plot
  plot_file <- file.path(output_dir, "comparison_plots", 
                        "horizontal_bar_grid_motif_names_with_input.pdf")
  
  ggsave(plot_file, horizontal_plot, 
         width = max(14, n_cols * 3), 
         height = max(10, n_rows * 3),
         limitsize = FALSE)
  
  cat("\n✓ Horizontal bar grid plot saved:", plot_file, "\n")
  cat("  Dimensions:", n_rows, "rows ×", n_cols, "columns\n")
  cat("  Each panel shows 25 motifs with names on Y-axis\n")
  
  return(horizontal_plot)
}

## =============================================================================
## 12. ADDITIONAL UNIFIED FIGURE: POOLED TOP MOTIFS WITH BUBBLE PLOT (FIXED)
## =============================================================================
## =============================================================================
## 12. ADDITIONAL UNIFIED FIGURE: POOLED TOP MOTIFS WITH BUBBLE PLOT (FIXED) - LOG2 ENRICHMENT VERSION
## =============================================================================

create_pooled_motif_bubble_plot_fixed <- function(results_df, top_motifs_list, output_dir, hist, cell_types, tools, n_top_motifs = 25) {
  cat("\nCreating fixed unified bubble plot with log2 enrichment...\n")
  
  # Step 1: Pool top motifs from all tools for each cell type
  pooled_motifs_data <- list()
  pooled_top_motifs_list <- list()
  
  for (cell_type in cell_types) {
    if (!is.null(top_motifs_list[[cell_type]]) && nrow(top_motifs_list[[cell_type]]) > 0) {
      # Get top N motifs for this cell type
      n_motifs <- min(n_top_motifs, nrow(top_motifs_list[[cell_type]]))
      cell_top_motifs <- top_motifs_list[[cell_type]]$motif[1:n_motifs]
      
      # Save pooled top motifs
      pooled_top_motifs_list[[cell_type]] <- data.frame(
        cell_type = cell_type,
        motif = cell_top_motifs,
        original_rank = 1:length(cell_top_motifs),
        stringsAsFactors = FALSE
      )
      
      # Get log2 enrichment data for these motifs across all tools
      cell_data <- results_df %>%
        filter(cell_type == !!cell_type, motif %in% cell_top_motifs) %>%
        group_by(tool, motif) %>%
        summarise(
          frequency = mean(log2_enrichment, na.rm = TRUE),  # log2_enrichment ব্যবহার
          .groups = 'drop'
        )
      
      if (nrow(cell_data) > 0) {
        # Calculate overall rank for each motif based on log2 enrichment
        motif_ranks <- cell_data %>%
          group_by(motif) %>%
          summarise(
            mean_enrichment = mean(frequency, na.rm = TRUE),  # নাম পরিবর্তন
            .groups = 'drop'
          ) %>%
          arrange(desc(mean_enrichment)) %>%
          mutate(pooled_rank = row_number())
        
        # Merge with tool-specific data
        cell_data_ranked <- cell_data %>%
          left_join(motif_ranks, by = "motif") %>%
          mutate(cell_type = cell_type)
        
        pooled_motifs_data[[cell_type]] <- cell_data_ranked
        cat("  Processed", cell_type, ":", n_motifs, "motifs\n")
      }
    }
  }
  
  # NEW: Save pooled top motifs by tool for each cell type
  if (length(pooled_motifs_data) > 0) {
    dir.create(file.path(output_dir, "pooled_motifs_by_tool"), showWarnings = FALSE, recursive = TRUE)
    
    # Create a list to store tool-wise motifs for each cell type
    tool_wise_motifs_all <- list()
    
    for (cell_type in names(pooled_motifs_data)) {
      cell_data <- pooled_motifs_data[[cell_type]]
      
      if (nrow(cell_data) > 0) {
        # Create separate data frames for each tool
        tool_wise_motifs <- list()
        
        for (tool_name in tools) {
          tool_data <- cell_data %>%
            filter(tool == tool_name) %>%
            arrange(desc(frequency)) %>%
            head(n_top_motifs) %>%
            select(motif, frequency, pooled_rank) %>%
            mutate(
              rank_in_tool = row_number(),
              cell_type = cell_type,
              tool = tool_name
            )
          
          tool_wise_motifs[[tool_name]] <- tool_data
        }
        
        # Combine all tools for this cell type
        cell_tool_motifs <- bind_rows(tool_wise_motifs)
        tool_wise_motifs_all[[cell_type]] <- cell_tool_motifs
        
        # Save to CSV for this cell type
        cell_file <- file.path(output_dir, "pooled_motifs_by_tool", 
                               paste0(hist, "_", cell_type, "_top_motifs_by_tool_with_input.csv"))
        write.csv(cell_tool_motifs, cell_file, row.names = FALSE)
        cat("  Saved tool-wise motifs for", cell_type, "to:", cell_file, "\n")
        
        # Also save individual tool files
        dir.create(file.path(output_dir, "pooled_motifs_by_tool", cell_type), 
                   showWarnings = FALSE, recursive = TRUE)
        
        for (tool_name in tools) {
          if (tool_name %in% cell_tool_motifs$tool) {
            tool_file <- file.path(output_dir, "pooled_motifs_by_tool", cell_type,
                                   paste0(hist, "_", cell_type, "_", tool_name, "_top_motifs_with_input.csv"))
            tool_data <- cell_tool_motifs %>%
              filter(tool == tool_name) %>%
              select(motif, frequency, rank_in_tool, pooled_rank)
            
            write.csv(tool_data, tool_file, row.names = FALSE)
          }
        }
      }
    }
    
    # Save combined tool-wise motifs for all cell types
    if (length(tool_wise_motifs_all) > 0) {
      all_tool_motifs <- bind_rows(tool_wise_motifs_all)
      combined_file <- file.path(output_dir, "pooled_motifs_by_tool",
                                 paste0(hist, "_ALL_celltypes_top_motifs_by_tool_with_input.csv"))
      write.csv(all_tool_motifs, combined_file, row.names = FALSE)
      cat("\n✓ Saved combined tool-wise motifs to:", combined_file, "\n")
    }
  }
  
  # Save pooled top motifs to CSV (original function)
  if (length(pooled_top_motifs_list) > 0) {
    all_pooled_motifs <- bind_rows(pooled_top_motifs_list)
    
    dir.create(file.path(output_dir, "pooled_motifs"), showWarnings = FALSE, recursive = TRUE)
    
    motifs_file <- file.path(output_dir, "pooled_motifs", 
                             paste0(hist, "_pooled_top_motifs_by_celltype_with_input.csv"))
    write.csv(all_pooled_motifs, motifs_file, row.names = FALSE)
    cat("\n✓ Saved pooled top motifs to:", motifs_file, "\n")
    cat("  Total:", nrow(all_pooled_motifs), "motifs across", 
        length(unique(all_pooled_motifs$cell_type)), "cell types\n")
  }
  
  # Combine all cell type data
  all_pooled_data <- bind_rows(pooled_motifs_data)
  
  if (nrow(all_pooled_data) == 0) {
    cat("No data available for pooled motif bubble plot\n")
    return(NULL)
  }
  
  # Debug: Print log2 enrichment statistics
  cat("\n=== LOG2 ENRICHMENT STATISTICS ===\n")
  freq_summary <- all_pooled_data %>%
    summarise(
      min_freq = min(frequency, na.rm = TRUE),
      max_freq = max(frequency, na.rm = TRUE),
      mean_freq = mean(frequency, na.rm = TRUE),
      sd_freq = sd(frequency, na.rm = TRUE),
      q25 = quantile(frequency, 0.25, na.rm = TRUE),
      q75 = quantile(frequency, 0.75, na.rm = TRUE)
    )
  print(freq_summary)
  
  # Apply min-max scaling
  freq_min <- freq_summary$min_freq
  freq_max <- freq_summary$max_freq
  
  all_pooled_data <- all_pooled_data %>%
    mutate(
      frequency_scaled = scales::rescale(frequency, to = c(0.1, 1))
    )
  
  # Step 2: Create bubble plots for each cell type
  create_cell_type_plot <- function(cell_data, cell_name, tools) {
    if (nrow(cell_data) == 0) {
      # Return empty plot
      return(
        ggplot() + 
          theme_void() +
          labs(title = cell_name) +
          theme(
            plot.title = element_text(size = 7, face = "bold", hjust = 0.5)
          )
      )
    }
    
    # Get top 30 motifs for this cell type
    top_motifs_cell <- cell_data %>%
      group_by(motif) %>%
      summarise(avg_rank = mean(pooled_rank, na.rm = TRUE)) %>%
      arrange(avg_rank) %>%
      head(30) %>%
      pull(motif)
    
    plot_data <- cell_data %>%
      filter(motif %in% top_motifs_cell) %>%
      mutate(
        tool = factor(tool, levels = tools),
        motif = factor(motif, levels = rev(unique(motif[order(pooled_rank)])))
      )
    
    # Create bubble plot
    ggplot(plot_data, aes(x = tool, y = motif)) +
      geom_point(aes(size = frequency, color = frequency), 
                 alpha = 0.8, shape = 16) +
      scale_size_continuous(
        name = "Log2 Enrichment",
        range = c(0.5, 3),
        limits = c(freq_min, freq_max)
      ) +
      scale_color_viridis_c(
        name = "Log2 Enrichment",
        option = "plasma",
        limits = c(freq_min, freq_max)
      ) +
      scale_x_discrete(
        limits = tools,
        expand = expansion(mult = c(0.15, 0.15))
      ) +
      scale_y_discrete(
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      theme_minimal(base_size = 5.5) +
      labs(
        x = "Peak Caller",
        y = "Motif",
        title = cell_name
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5),
        axis.text.y = element_text(size = 3.5, face = "bold"),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.05),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(1, 1, 1, 1, "mm"),
        aspect.ratio = 1.5
      )
  }
  
  # Create plots for all cell types
  bubble_plots <- list()
  valid_cell_types <- character()
  
  for (cell_type in cell_types) {
    cell_plot_data <- all_pooled_data %>% filter(cell_type == !!cell_type)
    if (nrow(cell_plot_data) > 0) {
      valid_cell_types <- c(valid_cell_types, cell_type)
    }
    bubble_plots[[cell_type]] <- create_cell_type_plot(cell_plot_data, cell_type, tools)
  }
  
  # Step 3: Create shared legend
  legend_data <- data.frame(frequency = seq(freq_min, freq_max, length.out = 6))
  legend_plot <- ggplot(legend_data, aes(x = frequency, y = 1, 
                                         color = frequency, size = frequency)) +
    geom_point() +
    scale_color_viridis_c(
      option = "plasma",
      name = "Log2 Enrichment",
      limits = c(freq_min, freq_max)
    ) +
    scale_size_continuous(
      name = "Log2 Enrichment",
      range = c(0.3, 2),
      limits = c(freq_min, freq_max)
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 5),
      legend.key.width = unit(1.5, "cm"),
      legend.key.height = unit(0.2, "cm")
    )
  
  shared_legend <- cowplot::get_legend(legend_plot)
  
  # Step 4: Create grid layout
  if (length(valid_cell_types) > 0) {
    n_plots <- length(valid_cell_types)
    
    # Determine grid dimensions
    grid_dims <- switch(
      as.character(min(n_plots, 12)),
      "1" = c(1, 1),
      "2" = c(1, 2),
      "3" = c(1, 3),
      "4" = c(2, 2),
      "5" = c(2, 3),
      "6" = c(2, 3),
      "7" = c(2, 4),
      "8" = c(2, 4),
      "9" = c(3, 3),
      "10" = c(2, 5),
      "11" = c(3, 4),
      "12" = c(3, 4)
    )
    
    n_rows <- grid_dims[1]
    n_cols <- grid_dims[2]
    
    cat("Creating grid with", n_plots, "plots in", n_rows, "rows x", n_cols, "columns\n")
    
    plot_grid <- wrap_plots(
      bubble_plots[valid_cell_types],
      ncol = n_cols,
      nrow = n_rows,
      guides = "collect"
    )
  } else {
    plot_grid <- ggplot() + 
      theme_void() +
      labs(title = "No data available for bubble plot")
  }
  
  # Step 5: Create final composition
  final_composition <- wrap_plots(
    # Title
    wrap_elements(
      grid::textGrob(
        paste(hist, "- Pooled Motif Log2 Enrichment Across Peak Callers"),
        gp = grid::gpar(fontsize = 10, fontface = "bold"),
        hjust = 0.5
      )
    ),
    # Plot grid
    plot_grid,
    # Legend
    shared_legend,
    ncol = 1,
    heights = c(0.05, 0.90, 0.05)
  )
  
  # Step 6: Save as PDF
  dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
  pdf_file <- file.path(output_dir, "figures", 
                        paste0(hist, "_Additional_Unified_Bubble_FIXED_A4_with_input_log2enrichment_with_input.pdf"))
  
  pdf(pdf_file, width = 11.69, height = 8.27)
  print(final_composition)
  dev.off()
  
  cat("\n✓ Fixed unified bubble plot saved as:", pdf_file, "\n")
  cat("✓ Number of cell types with data:", length(valid_cell_types), "\n")
  cat("✓ Log2 Enrichment range:", round(freq_min, 4), "to", round(freq_max, 4), "\n")
  
  # Return results
  return(list(
    data = all_pooled_data,
    freq_summary = freq_summary,
    plots = bubble_plots,
    valid_cell_types = valid_cell_types,
    pooled_top_motifs = if(exists("all_pooled_motifs")) all_pooled_motifs else NULL,
    tool_wise_motifs = if(exists("all_tool_motifs")) all_tool_motifs else NULL
  ))
}

## =============================================================================
## 13. PREPARE DATA FOR BUBBLE PLOT
## =============================================================================

prepare_top_motifs_list <- function(all_results) {
  cat("\nPreparing top motifs list for bubble plot...\n")
  
  top_motifs_list <- list()
  
  for (cell_type in cell_types) {
    # Get top motifs for this cell type across all tools
    cell_data <- all_results %>%
      filter(cell_type == !!cell_type) %>%
      group_by(motif) %>%
      summarise(
        avg_enrichment = mean(log2_enrichment, na.rm = TRUE),
        avg_frequency = mean(foreground_frequency, na.rm = TRUE),
        n_tools = n_distinct(tool),
        .groups = 'drop'
      ) %>%
      arrange(desc(avg_enrichment)) %>%
      head(25)  # Top 25 motifs for this cell type
    
    if (nrow(cell_data) > 0) {
      top_motifs_list[[cell_type]] <- cell_data
      cat("  ", cell_type, ":", nrow(cell_data), "motifs\n")
    }
  }
  
  return(top_motifs_list)
}

## =============================================================================
## 14. MAIN ANALYSIS FUNCTION
## =============================================================================

run_complete_analysis <- function() {
  cat("========================================================================\n")
  cat("STARTING COMPLETE MOTIF ANALYSIS PIPELINE\n")
  cat("Horizontal Bar Plot + Unified Bubble Plot\n")
  cat("========================================================================\n\n")
  
  start_time <- Sys.time()
  
  # Step 1: Build motif database
  cat("1. Building motif database...\n")
  motif_db <- build_motif_database()
  
  # Step 2: Run analysis for each tool-cell combination
  cat("\n2. Running motif analysis for each tool-cell combination...\n")
  
  all_results <- data.frame()
  
  total_combinations <- length(tools) * length(cell_types)
  current_combo <- 0
  
  for (tool in tools) {
    cat("\n--- Processing Tool:", tool, "---\n")
    
    for (cell_type in cell_types) {
      current_combo <- current_combo + 1
      cat(sprintf("\n[%d/%d] ", current_combo, total_combinations))
      
      # FIXED: Specify which histone to use
      # You need to define which histone you want to analyze
      # For example, use the first histone from your list
      histone_to_analyze <- histones[1]  # Or specify which histone
      
      # Find and load peaks
      bed_file <- find_bed_file(tool, histone_to_analyze, cell_type)
      
      if (is.null(bed_file)) {
        cat("  No BED file found\n")
        next
      }
      
      peaks <- load_peaks_from_bed(bed_file, cell_type)
      
      if (length(peaks) == 0) {
        cat("  No peaks loaded\n")
        next
      }
      
      # Run motif analysis
      results <- analyze_motifs_for_sample(
        peaks = peaks,
        tool_name = tool,
        histone_name = histone_to_analyze,  # Pass histone name
        cell_type_name = cell_type,
        motif_db = motif_db,
        n_foreground = 500,
        n_background = 1000
      )
      
      if (nrow(results) > 0) {
        # FIXED: Pass all three parameters
        save_individual_results(results, tool, histone_to_analyze, cell_type)
        
        # Show top motifs
        show_top_motifs(results, tool, cell_type)
        
        # Add to combined results
        all_results <- rbind(all_results, results)
        
        cat("  ✓ Success\n")
      } else {
        cat("  ✗ No results\n")
      }
    }
  }
  
  # Step 3: Save all results
  if (nrow(all_results) > 0) {
    cat("\n\n3. Saving combined results...\n")
    
    # Save complete results
    write.csv(all_results,
              file.path(output_dir, "summary_tables", "ALL_motif_results_with_input.csv"),
              row.names = FALSE)
    
    # Save top 25 summary
    top_25_summary <- all_results %>%
      group_by(tool, cell_type) %>%
      arrange(desc(log2_enrichment)) %>%
      slice(1:25) %>%
      ungroup()
    
    write.csv(top_25_summary,
              file.path(output_dir, "summary_tables", "ALL_top25_summary_with_input.csv"),
              row.names = FALSE)
    
    cat("✓ Saved combined results\n")
    
    # Step 4: Create horizontal bar grid plot
    cat("\n4. Creating horizontal bar grid plot...\n")
    main_plot <- create_horizontal_bar_grid_plot(all_results, output_dir)
    
    # Step 5: Prepare data and create unified bubble plot
    cat("\n5. Creating unified bubble plot...\n")
    top_motifs_list <- prepare_top_motifs_list(all_results)
    
    if (length(top_motifs_list) > 0) {
      bubble_results <- create_pooled_motif_bubble_plot_fixed(
        results_df = all_results,
        top_motifs_list = top_motifs_list,
        output_dir = output_dir,
        hist = histone_to_analyze,  # FIXED: Pass histone name
        cell_types = cell_types,
        tools = tools,
        n_top_motifs = 25
      )
    }
    
    # Step 6: Generate report
    generate_analysis_report(all_results, start_time)
    
  } else {
    cat("\n✗ No results generated from analysis\n")
  }
  
  end_time <- Sys.time()
  cat("\nTotal runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
  
  return(all_results)
}

## =============================================================================
## 15. ANALYSIS REPORT
## =============================================================================

generate_analysis_report <- function(all_results, start_time) {
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  report_file <- file.path(output_dir, "ANALYSIS_REPORT.txt")
  
  sink(report_file)
  
  cat("MOTIF ANALYSIS REPORT - COMPLETE PIPELINE\n")
  cat("=========================================\n\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Runtime:", round(runtime, 2), "minutes\n\n")
  
  cat("ANALYSIS PARAMETERS\n")
  cat("-------------------\n")
  cat("Histone: H3K27ac-b\n")
  cat("Tools:", paste(tools, collapse = ", "), "\n")
  cat("Cell Types:", paste(cell_types, collapse = ", "), "\n\n")
  
  cat("RESULTS SUMMARY\n")
  cat("---------------\n")
  
  if (nrow(all_results) > 0) {
    cat("Total motif records:", nrow(all_results), "\n")
    cat("Cell types with data:", n_distinct(all_results$cell_type), "\n")
    cat("Tools with data:", n_distinct(all_results$tool), "\n")
    cat("Unique motifs found:", n_distinct(all_results$motif), "\n\n")
  } else {
    cat("No results generated\n")
  }
  
  cat("\nOUTPUT FILES GENERATED\n")
  cat("----------------------\n")
  cat("1. Horizontal Bar Plot:\n")
  cat("   - comparison_plots/horizontal_bar_grid_motif_names.pdf\n")
  cat("   - Shows top 25 motifs per cell type-tool combination\n")
  cat("   - X-axis: Log2 enrichment, Y-axis: Motif names\n\n")
  
  cat("2. Unified Bubble Plot:\n")
  cat("   - figures/H3K27ac-b_Additional_Unified_Bubble_FIXED_A4_with_input.pdf\n")
  cat("   - Shows pooled motif frequencies across tools\n")
  cat("   - Bubble size/color: Motif frequency\n\n")
  
  cat("3. Pooled Motifs Data:\n")
  cat("   - pooled_motifs/H3K27ac-b_pooled_top_motifs_by_celltype.csv\n")
  cat("   - pooled_motifs_by_tool/ - Tool-wise motif data\n\n")
  
  cat("4. Individual Results:\n")
  cat("   - individual_results/ - All tool-cell combination results\n")
  cat("   - summary_tables/ - Combined results and summaries\n\n")
  
  sink()
  
  cat("\n✓ Analysis report saved:", report_file, "\n")
}

## =============================================================================
## 16. EXECUTE ANALYSIS
## =============================================================================

# Run the complete analysis
cat("========================================================================\n")
cat("EXECUTING COMPLETE MOTIF ANALYSIS PIPELINE\n")
cat("========================================================================\n\n")

final_results <- run_complete_analysis()

if (!is.null(final_results) && nrow(final_results) > 0) {
  cat("\n========================================================================\n")
  cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("========================================================================\n")
  
  # Show final summary
  cat("\nFINAL OUTPUT SUMMARY:\n")
  cat("--------------------\n")
  cat("1. Main Horizontal Bar Plot:\n")
  cat("   File: comparison_plots/horizontal_bar_grid_motif_names.pdf\n")
  cat("   Features: 8 cell types × 7 tools = 56 panels\n")
  cat("   Each panel: Top 25 motifs with names on Y-axis\n")
  cat("   X-axis: Log2 enrichment values\n\n")
  
  cat("2. Unified Bubble Plot:\n")
  cat("   File: figures/H3K27ac-b_Additional_Unified_Bubble_FIXED_A4_with_input.pdf\n")
  cat("   Features: Pooled motif frequencies across all tools\n")
  cat("   Bubble size/color: Represents motif frequency\n")
  cat("   Grid: Cell types arranged in rows\n\n")
  
  cat("3. Data Files:\n")
  cat("   - Pooled motifs: pooled_motifs/\n")
  cat("   - Tool-wise motifs: pooled_motifs_by_tool/\n")
  cat("   - Individual results: individual_results/\n")
  cat("   - Summary tables: summary_tables/\n\n")
  
  cat("All output files are in:", output_dir, "\n")
  
} else {
  cat("\n========================================================================\n")
  cat("ANALYSIS COMPLETED BUT NO RESULTS GENERATED\n")
  cat("========================================================================\n")
  cat("\nCheck console output above for specific errors.\n")
}

cat("\nTo view plots, check:\n")
cat("1.", file.path(output_dir, "comparison_plots"), "\n")
cat("2.", file.path(output_dir, "figures"), "\n")



# List all the files
file_paths <- c(
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/individual_results/DROMPAplus_H3K27ac-s_B_top25_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/individual_results/Genrich_H3K27ac-s_B_top25_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/individual_results/GoPeaks_H3K27ac-s_B_top25_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/individual_results/HOMER_H3K27ac-s_B_top25_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/individual_results/MACS2_H3K27ac-s_B_top25_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/individual_results/SEACR_H3K27ac-s_B_top25_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/individual_results/SICER2_H3K27ac-s_B_top25_with_input.csv"
)

# Read and combine all files
library(dplyr)
library(readr)

merged_data <- bind_rows(
  lapply(file_paths, read_csv)
)

# Save the merged result
write_csv(merged_data, "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/H3K27ac-s_B_merged_top25_with_input.csv")

# View the structure
glimpse(merged_data)





# Load required libraries
library(ggplot2)
library(dplyr)

# Read your data
data <- read.csv("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/H3K27ac-s_B_merged_top25_without_input.csv")

# Create pooled data ranked by MAXIMUM enrichment
pooled_data <- data %>%
  group_by(motif) %>%
  summarise(
    max_log2_enrichment = max(log2_enrichment, na.rm = TRUE),
    avg_log2_enrichment = mean(log2_enrichment, na.rm = TRUE),
    n_tools = n(),
    max_tool = tool[which.max(log2_enrichment)],
    max_value = max(log2_enrichment, na.rm = TRUE)
  ) %>%
  arrange(desc(max_log2_enrichment)) %>%
  mutate(
    pooled_rank = row_number(),
    cell_type = "B"
  )

# Select top 10 motifs based on maximum log2_enrichment
top_10_motifs <- pooled_data %>%
  arrange(desc(max_log2_enrichment)) %>%
  slice_head(n = 15) %>%
  pull(motif)

# Filter original data for these top 10 motifs
top_10_data <- data %>%
  filter(motif %in% top_10_motifs) %>%
  mutate(
    # Order motifs by maximum enrichment (rank 1 at top)
    motif = factor(motif, levels = top_10_motifs)
  )

# COMPACT VERSION 1: Minimal design
compact_bubble <- ggplot(top_10_data, aes(x = tool, y = motif)) +
  geom_point(aes(size = log2_enrichment), 
             color = "black",
             alpha = 0.7,
             shape = 16) +
  scale_size_continuous(
    range = c(2, 8),  # Smaller size range
    name = "log2 Enrichment"
  ) +
  labs(
    x = "Tool",
    y = "Motif"
  ) +
  theme_minimal(base_size = 9) +  # Smaller base size
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7, face = "bold"),
    axis.title = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(5, 5, 5, 5, "pt")  # Smaller margins
  ) +
  scale_y_discrete(limits = rev)

# COMPACT VERSION 2: Even more minimal
ultra_compact <- ggplot(top_10_data, aes(x = tool, y = motif)) +
  geom_point(aes(size = log2_enrichment), 
             color = "black",
             alpha = 0.8,
             shape = 16) +
  scale_size_continuous(
    range = c(1.5, 6),
    guide = "none"  # Remove legend for ultra compact
  ) +
  labs(
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "bold"),
    axis.text.y = element_text(size = 6, face = "bold"),
    panel.grid.major = element_line(color = "grey97", linewidth = 0.1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  scale_y_discrete(limits = rev)

# COMPACT VERSION 3: Square grid format
square_compact <- ggplot(top_10_data, aes(x = tool, y = motif)) +
  geom_point(aes(size = log2_enrichment), 
             color = "black",
             alpha = 0.8,
             shape = 16) +
  scale_size_continuous(
    range = c(2, 7),
    name = ""
  ) +
  labs(
    x = "Tool",
    y = "Motif"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.text = element_text(size = 6),
    plot.margin = margin(3, 3, 3, 3, "mm"),
    aspect.ratio = 1.2  # Control aspect ratio
  ) +
  scale_y_discrete(limits = rev)

# COMPACT VERSION 4: Professional publication style
pub_style <- ggplot(top_10_data, aes(x = tool, y = motif)) +
  geom_point(aes(size = log2_enrichment), 
             color = "black",
             fill = "black",
             alpha = 0.9,
             shape = 21,
             stroke = 0.3) +
  scale_size_continuous(
    range = c(1.5, 6),
    name = expression(paste(log[2], " Enrichment")),
    breaks = c(1, 3, 5, 7)
  ) +
  labs(
    x = "Peak Calling Tool",
    y = "Transcription Factor Motif"
  ) +
  theme_bw(base_size = 7) +  # Use bw theme for publication
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 5.5, color = "black", face = "italic"),
    axis.title.x = element_text(size = 7, margin = margin(t = 3)),
    axis.title.y = element_text(size = 7, margin = margin(r = 3)),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.15),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "white"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  scale_y_discrete(limits = rev)

# Save as PDF
output_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_motif_final_trail_V8/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save all compact versions
ggsave(
  filename = file.path(output_dir, "H3K27ac-s_compact_bubble_1_without_input_B.pdf"),
  plot = compact_bubble,
  width = 6,  # Smaller width
  height = 5, # Smaller height
  device = "pdf"
)

ggsave(
  filename = file.path(output_dir, "H3K27ac-s_compact_bubble_2_ultra_without_input_B.pdf"),
  plot = ultra_compact,
  width = 5,
  height = 4,
  device = "pdf"
)

ggsave(
  filename = file.path(output_dir, "H3K27ac-s_compact_bubble_3_square_without_input_B.pdf"),
  plot = square_compact,
  width = 5.5,
  height = 5,
  device = "pdf"
)

ggsave(
  filename = file.path(output_dir, "H3K27ac-s_compact_bubble_4_publication_without_input_B.pdf"),
  plot = pub_style,
  width = 6,
  height = 5,
  device = "pdf"
)

# Print summary
cat("\n✓ Compact bubble plots saved:\n")
cat("1. compact_bubble_1.pdf (6x5 inches)\n")
cat("2. compact_bubble_2_ultra.pdf (5x4 inches, no legend)\n")
cat("3. compact_bubble_3_square.pdf (5.5x5 inches)\n")
cat("4. compact_bubble_4_publication.pdf (6x5 inches, publication style)\n")
cat("\nTop 20 motifs (ranked by maximum enrichment):\n")
for(i in 1:min(20, length(top_10_motifs))) {
  cat(sprintf("%2d. %s\n", i, top_10_motifs[i]))
}
