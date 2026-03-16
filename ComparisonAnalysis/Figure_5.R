
###########################################################################################################
## Figure 5 A 
###########################################################################################################
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
library(ggpubr)        # For publication-ready themes
library(viridis)       # For better color scales
library(ggsci)         # For Nature-style color palettes
library(scales)        # For better axis formatting
library(conflicted)    # For namespace management

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
# Define histones to process
#histones <- c( "H3K27ac-b", "H3K27ac-s")

histones <- c( "H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
##=====================================================================
# loading the datasets
##=====================================================================
# 1. Load the pre-aggregated TPM expression matrix
for (hist in histones){
TPM <- read.delim("/home/wahid/project_scHMTF/GSE195725_processed_data/scRNA_seq/10x_scRNAseq/pbmc_celltype_averaged_TPM.txt.gz", row.names = 1)

# 3. Define parameters
celltypes <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")
peak_callers <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
#sc_peak_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/result/l1_withcontrol/H3K27me3_peakbed"
sc_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_without_input_peakbed_corrected/scrnaseq_peakbed/", hist, "_peakbed")


# 4. Create output directory
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_ex/result_GT_scrnaseq_withhout_input/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)




##=====================================================================
# Comprehensive BED File Reading Function for each peak caller - FIXED SICER2
##=====================================================================
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
##=====================================================================
# FIXED: Process All Peak Files with Proper RP Model Implementation
# MODIFIED: Ensure same number of gene activities for all peak callers
##=====================================================================
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
    write.csv(temp_activity, file.path(out_dir, paste0("temp_peak_activity_data_", hist, "_100KB_protein_coding_withhout_input.csv")), row.names = FALSE)
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
  final_output_file <- file.path(out_dir, paste0("final_peak_counts_with_RP_model_", hist, "_100KB_protein_coding_withhout_input.csv"))
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



####################################################################
# Prepare and Merge with Expression Data - FIXED with retry logic
####################################################################

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


################TOP1000genes_final_merged_data by peak caller and celltype [Top gene selected by higher peaks]

# Get top 1000 genes for each peak_caller method based on TPM expression
TOP_final_merged_data <- final_merged_data %>%
  group_by(peak_caller,celltype) %>%
  arrange(desc(weighted_peak_count), .by_group = TRUE) %>%
  slice_head(n = 1000) %>%
  ungroup()



####################################################################
# Correlation analysis
####################################################################

# Function to calculate correlations for each cell type and peak caller
correlation_peak_count <- function(data, dataset_name = "") {
  
  # Group by both peak_caller and celltype
  result <- data %>%
    group_by(peak_caller, celltype) %>%
    summarise(
      spearman_cor = suppressWarnings(cor(peak_count, tpm_expression, method = "spearman", use = "pairwise.complete.obs")),
      spearman_pvalue = suppressWarnings(cor.test(peak_count, tpm_expression, method = "spearman", exact = FALSE)$p.value),
      n_genes = n(), .groups = "drop"
    ) %>%
    
    # Add the crucial p-value adjustment step here
    group_by(celltype) %>%
    mutate(
      spearman_pvalue_adj = p.adjust(spearman_pvalue, method = "BH"), spearman_significant = spearman_pvalue_adj < 0.05 & !is.na(spearman_pvalue_adj),
      dataset = dataset_name
    ) %>%
    ungroup()
  
  return(result)
}

# Apply the function to all datasets


# Function to calculate correlations for each cell type and peak caller
correlation_weighted <- function(data, dataset_name = "") {
  
  # Group by both peak_caller and celltype
  result <- data %>%
    group_by(peak_caller, celltype) %>%
    summarise(
      spearman_cor = suppressWarnings(cor(weighted_peak_count, tpm_expression, method = "spearman", use = "pairwise.complete.obs")),
      spearman_pvalue = suppressWarnings(cor.test(weighted_peak_count, tpm_expression, method = "spearman", exact = FALSE)$p.value),
      n_genes = n(), .groups = "drop"
    ) %>%
    
    # Add the crucial p-value adjustment step here
    group_by(celltype) %>%
    mutate(
      spearman_pvalue_adj = p.adjust(spearman_pvalue, method = "BH"), spearman_significant = spearman_pvalue_adj < 0.05 & !is.na(spearman_pvalue_adj),
      dataset = dataset_name
    ) %>%
    ungroup()
  
  return(result)
}

# Apply the function to all datasets
correlation_peak_count_all <- correlation_peak_count(final_merged_data, "Unified_all_genes")
correlation_peak_count_TOP <- correlation_peak_count(TOP_final_merged_data, "TOP1KGENE")

correlation_weighted_all <- correlation_weighted(final_merged_data, "Unified_all_genes")
correlation_weighted_TOP <- correlation_weighted(TOP_final_merged_data, "TOP1KGENE")


cor_matrix_peak_count_all <- correlation_peak_count_all %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_peak_count_all) <- cor_matrix_peak_count_all$celltype
cor_matrix_peak_count_all <- as.matrix(cor_matrix_peak_count_all[, -1])
cor_matrix_peak_count_all



cor_matrix_peak_count_TOP <- correlation_peak_count_TOP %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_peak_count_TOP) <- cor_matrix_peak_count_TOP$celltype
cor_matrix_peak_count_TOP <- as.matrix(cor_matrix_peak_count_TOP[, -1])
cor_matrix_peak_count_TOP





#correlation_weighted_TOP1KGENE
cor_matrix_weighted_all <- correlation_weighted_all %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_weighted_all) <- cor_matrix_weighted_all$celltype
cor_matrix_weighted_all <- as.matrix(cor_matrix_weighted_all[, -1])
cor_matrix_weighted_all

cor_matrix_weighted_TOP <- correlation_weighted_TOP %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_weighted_TOP) <- cor_matrix_weighted_TOP$celltype
cor_matrix_weighted_TOP <- as.matrix(cor_matrix_weighted_TOP[, -1])
cor_matrix_weighted_TOP




####################################################################
# Enhanced Nature Genetics Style Theme
####################################################################

theme_nature_genetics <- function(base_size = 7, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black", size = base_size),
      axis.line = element_line(color = "black", linewidth = 0.25),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = unit(0.5, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1, 
                               hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, 
                                  color = "grey40", margin = margin(b = 4)),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(face = "bold", size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      legend.key.size = unit(3, "mm"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = unit(1, "mm"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

theme_set(theme_nature_genetics())

nature_genetics_palette <- c(
  "#E64B35", "#8bb5bfff", "#0cc3a7ff", "#4e9a6aff", 
  "#274ab1ff", "#4e77e7ff", "#3e4f4bff", "#13c758ff"
)




####################################################################
# NEW FIGURE 4: Performance Summary - Boxplot Style
####################################################################



# Clean data
correlation_weighted_all_clean <- correlation_weighted_all %>%
  filter(!is.na(spearman_cor) & is.finite(spearman_cor))

# Define FIXED y-axis order (vertical axis)
method_order <- c( "SICER2", "SEACR", "MACS2",  "HOMER", "GoPeaks", "Genrich", "DROMPAplus")

# Calculate mean/median for ordering by performance
method_stats <- correlation_weighted_all_clean %>%
  group_by(peak_caller) %>%
  summarise(
    mean_cor = mean(spearman_cor, na.rm = TRUE),
    median_cor = median(spearman_cor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cor))  # Best to worst

cat("Method ranking by mean correlation:\n")
print(method_stats)

# OPTION A: Order by performance (best to worst) - RECOMMENDED
performance_order <- as.character(method_stats$peak_caller)

# OPTION B: Use fixed order
fixed_order <- method_order

# Choose which order you want
selected_order <- fixed_order  # Change to performance_order if you want ranking

# Create plot with HORIZONTAL boxplots (no coord_flip)
performance_boxplot_weighted_all <- ggplot(
  correlation_weighted_all_clean, 
  aes(x = spearman_cor, 
      y = factor(peak_caller, levels = selected_order))  # Y-axis: fixed method order
) +
  
  geom_boxplot(
    fill = "grey90", 
    color = "black", 
    linewidth = 0.3, 
    outlier.size = 0.5,
    orientation = "y"  # Horizontal boxplot
  ) +
  
  geom_point(
    aes(color = celltype), 
    size = 1, 
    alpha = 0.7, 
    position = position_jitter(height = 0.2)
  ) +
  
  scale_color_manual(values = nature_genetics_palette) +
  scale_y_discrete(limits = selected_order) +  # Fix y-axis order
  
  labs(
    title = "Peak Caller Performance Distribution", 
    x = "Spearman Correlation",  # X-axis: correlation value
    y = "Peak Caller",           # Y-axis: methods in fixed order
    color = "Cell Type"
  ) +
  
  theme_nature_genetics() +  # Use your theme
  theme(
    legend.position = "bottom"
  )

# Save plot
ggsave(
  file.path(out_dir, paste0("Figure5A_Performance_Boxplot_", hist, "_byweighted_ALL_withhout_input.pdf")), 
  plot = performance_boxplot_weighted_all, 
  width = 10,    # Wider for horizontal
  height = 8,  
  units = "cm", 
  dpi = 1200
)

}



###########################################################################################################
## Figure 5 B
###########################################################################################################



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
library(ggpubr)        # For publication-ready themes
library(viridis)       # For better color scales
library(ggsci)         # For Nature-style color palettes
library(scales)        # For better axis formatting
library(conflicted)    # For namespace management
library(gridExtra)     # For arranging multiple plots

# Resolve namespace conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("slice_head", "dplyr")

# Function to get gene mapping with retry logic
get_gene_mapping_with_retry <- function(ensembl_ids, max_retries = 3) {
  for (attempt in 1:max_retries) {
    tryCatch({
      # Get gene information from EnsDb
      gene_info <- genes(EnsDb.Hsapiens.v86, 
                        filter = GeneIdFilter(ensembl_ids),
                        columns = c("gene_id", "gene_name", "symbol"))
      
      # Convert to data frame and clean up
      gene_df <- as.data.frame(gene_info) %>%
        dplyr::select(ensembl_gene_id = gene_id, 
                     hgnc_symbol = symbol, 
                     external_gene_name = gene_name) %>%
        distinct()
      
      return(gene_df)
      
    }, error = function(e) {
      if (attempt == max_retries) {
        cat("Failed after", max_retries, "attempts:", conditionMessage(e), "\n")
        
        # Fallback: return basic mapping
        return(data.frame(
          ensembl_gene_id = ensembl_ids,
          hgnc_symbol = NA_character_,
          external_gene_name = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
      cat("Attempt", attempt, "failed, retrying...\n")
      Sys.sleep(2) # Wait before retrying
    })
  }
}

 # Add this line
#histones <- c( "H3K27ac-b")
histones <- c( "H3K27ac-b", "H3K27ac-s", "H3K4me1", "H3K4me2", "H3K4me3")
##=====================================================================
# loading the datasets
##=====================================================================
# 1. Load the pre-aggregated TPM expression matrix
for (hist in histones){
TPM <- read.delim("/home/wahid/project_scHMTF/GSE195725_processed_data/bulk_RNA_seq/GSE107011_Processed_data_AggregateTPM.txt.gz", row.names = 1)
print("TPM matrix dimensions:")
dim(TPM)
print("Column names (cell types):")
colnames(TPM)


# 3. Define parameters
celltypes <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")
peak_callers <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
sc_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_with_input_peakbed_corrected/bulkrnaseq_peakbed/", hist, "_peakbed")

# 4. Create output directory
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_correlation_range/result_GT_bulkrnaseq_with_input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##=====================================================================
# Comprehensive BED File Reading Function for each peak caller - FIXED SICER2
##=====================================================================
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
      if (n_cols >= 9) {
        peaks$score <- as.numeric(bed_data[, 7]) # Use fold change
        mcols(peaks)$name <- bed_data[, 4]
        mcols(peaks)$strand <- bed_data[, 6]
        mcols(peaks)$fold_change <- as.numeric(bed_data[, 7])
        mcols(peaks)$neg_log10_pvalue <- as.numeric(bed_data[, 8])
        mcols(peaks)$neg_log10_qvalue <- as.numeric(bed_data[, 9])
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
##=====================================================================
# FIXED: Process All Peak Files with Proper RP Model Implementation
# MODIFIED: Ensure same number of gene activities for all peak callers
##=====================================================================
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



####################################################################
# Prepare and Merge with Expression Data - CORRECTED
####################################################################

# Resolve function conflicts at the beginning
conflicts_prefer(dplyr::slice)
conflicts_prefer(dplyr::slice_head)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::group_by)
conflicts_prefer(dplyr::ungroup)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::arrange)

# First, let's fix the TPM data structure
print("Original TPM matrix structure:")
print(dim(TPM))
print("First few rows:")
print(TPM[1:5, 1:5])

# The issue: your TPM data has a 'gene' column AND row names
# Let's extract the proper gene names from the 'gene' column
TPM_with_genes <- TPM %>%
  as.data.frame() %>%
  # Use the 'gene' column as the gene identifier
  mutate(ensembl_gene_id = sub("\\..*", "", gene)) %>%
  dplyr::select(-gene) %>%  # Remove the original gene column
  # Convert all sample columns to numeric
  mutate(across(all_of(celltypes), as.numeric))

print("Fixed TPM_with_genes structure:")
print(TPM_with_genes[1:5, 1:5])

# Check for duplicates in TPM data
duplicate_genes <- TPM_with_genes %>%
  dplyr::count(ensembl_gene_id) %>%
  dplyr::filter(n > 1) %>%
  nrow()

print(paste("Number of duplicate gene entries in TPM data:", duplicate_genes))

# Remove duplicates from TPM data by taking the mean expression
TPM_with_genes_unique <- TPM_with_genes %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(across(all_of(celltypes), mean, na.rm = TRUE)) %>%
  dplyr::ungroup()

print("TPM data after removing duplicates:")
print(dim(TPM_with_genes_unique))

# Now pivot to long format
TPM_long <- TPM_with_genes_unique %>%
  pivot_longer(
    cols = -ensembl_gene_id,
    names_to = "celltype",
    values_to = "tpm_expression"
  ) %>%
  dplyr::filter(celltype %in% celltypes)

print("TPM_long structure:")
print(head(TPM_long))

# Check for duplicates in TPM_long
duplicate_pairs <- TPM_long %>%
  dplyr::count(ensembl_gene_id, celltype) %>%
  dplyr::filter(n > 1) %>%
  nrow()

print(paste("Number of duplicate gene-celltype pairs in TPM_long:", duplicate_pairs))

# Remove version numbers from peak activity gene IDs
all_peak_activity <- all_peak_activity %>%
  dplyr::mutate(ensembl_gene_id = sub("\\..*", "", gene_id))

print("Peak activity gene ID examples:")
print(head(unique(all_peak_activity$ensembl_gene_id)))
print("TPM gene ID examples:")
print(head(unique(TPM_long$ensembl_gene_id)))

# Check for common genes
common_genes <- intersect(unique(all_peak_activity$ensembl_gene_id), 
                         unique(TPM_long$ensembl_gene_id))
print(paste("Number of common genes:", length(common_genes)))

# Get gene mapping for the common genes only
gene_mapping <- get_gene_mapping_with_retry(common_genes)

# Simplify the mapping to avoid many-to-many relationships
# Keep only one mapping per Ensembl ID (prefer HGNC symbol)
gene_mapping_unique <- gene_mapping %>%
  dplyr::group_by(ensembl_gene_id) %>%
  # Prefer non-empty HGNC symbols, then external names
  dplyr::arrange(desc(hgnc_symbol != ""), desc(external_gene_name != "")) %>%
  dplyr::slice(1) %>%  # Explicitly use dplyr::slice
  dplyr::ungroup()

print("Unique gene mapping results:")
print(head(gene_mapping_unique))

# Check for duplicates in peak data
peak_duplicates <- all_peak_activity %>%
  dplyr::count(ensembl_gene_id, celltype, peak_caller) %>%
  dplyr::filter(n > 1) %>%
  nrow()

print(paste("Number of duplicate gene-celltype-peak_caller entries in peak data:", peak_duplicates))

# If there are duplicates in peak data, aggregate them
if (peak_duplicates > 0) {
  all_peak_activity <- all_peak_activity %>%
    dplyr::group_by(ensembl_gene_id, celltype, peak_caller, gene_id, gene_name) %>%
    dplyr::summarise(
      peak_count = sum(peak_count, na.rm = TRUE),
      weighted_peak_count = sum(weighted_peak_count, na.rm = TRUE),
      .groups = "drop"
    )
  print("Aggregated duplicate peak entries")
}

# Merge the data step by step
# First merge peak activity with gene mapping
peak_with_names <- all_peak_activity %>%
  dplyr::inner_join(gene_mapping_unique, by = "ensembl_gene_id", relationship = "many-to-one") %>%
  # Use the best available gene name
  dplyr::mutate(gene_name_final = coalesce(hgnc_symbol, external_gene_name, gene_name))

print("Peak data with gene names:")
print(head(peak_with_names))

# Now merge with TPM data - use left_join instead of inner_join to see what's happening
merged_data <- peak_with_names %>%
  dplyr::left_join(TPM_long, by = c("ensembl_gene_id", "celltype"), relationship = "many-to-one")

print(paste("After left join:", nrow(merged_data), "rows"))
print(paste("Rows with TPM data:", sum(!is.na(merged_data$tpm_expression))))

# Filter and transform
final_merged_data <- merged_data %>%
  dplyr::filter(!is.na(tpm_expression) & !is.na(peak_count)) %>%
  dplyr::mutate(
    log_tpm = log2(tpm_expression + 1),
    log_peak_count = log10(peak_count + 1),
    log_weighted_count = log10(weighted_peak_count + 1)
  )

print(paste("Final dataset:", nrow(final_merged_data), "gene-celltype-peakcaller pairs"))
print("Final merged data structure:")
print(head(final_merged_data))

# Get top 1000 genes for each peak_caller method based on weighted peak count
TOP_final_merged_data <- final_merged_data %>%
  dplyr::group_by(peak_caller, celltype) %>%
  dplyr::arrange(desc(weighted_peak_count), .by_group = TRUE) %>%
  dplyr::slice_head(n = 1000) %>%  # Explicitly use dplyr::slice_head
  dplyr::ungroup()

print(paste("Top dataset:", nrow(TOP_final_merged_data), "rows"))


####################################################################
# Correlation analysis
####################################################################

# Function to calculate correlations for each cell type and peak caller
correlation_peak_count <- function(data, dataset_name = "") {
  
  # Group by both peak_caller and celltype
  result <- data %>%
    group_by(peak_caller, celltype) %>%
    summarise(
      spearman_cor = suppressWarnings(cor(peak_count, tpm_expression, method = "spearman", use = "pairwise.complete.obs")),
      spearman_pvalue = suppressWarnings(cor.test(peak_count, tpm_expression, method = "spearman", exact = FALSE)$p.value),
      n_genes = n(), .groups = "drop"
    ) %>%
    
    # Add the crucial p-value adjustment step here
    group_by(celltype) %>%
    mutate(
      spearman_pvalue_adj = p.adjust(spearman_pvalue, method = "BH"), spearman_significant = spearman_pvalue_adj < 0.05 & !is.na(spearman_pvalue_adj),
      dataset = dataset_name
    ) %>%
    ungroup()
  
  return(result)
}

# Apply the function to all datasets


# Function to calculate correlations for each cell type and peak caller
correlation_weighted <- function(data, dataset_name = "") {
  
  # Group by both peak_caller and celltype
  result <- data %>%
    group_by(peak_caller, celltype) %>%
    summarise(
      spearman_cor = suppressWarnings(cor(weighted_peak_count, tpm_expression, method = "spearman", use = "pairwise.complete.obs")),
      spearman_pvalue = suppressWarnings(cor.test(weighted_peak_count, tpm_expression, method = "spearman", exact = FALSE)$p.value),
      n_genes = n(), .groups = "drop"
    ) %>%
    
    # Add the crucial p-value adjustment step here
    group_by(celltype) %>%
    mutate(
      spearman_pvalue_adj = p.adjust(spearman_pvalue, method = "BH"), spearman_significant = spearman_pvalue_adj < 0.05 & !is.na(spearman_pvalue_adj),
      dataset = dataset_name
    ) %>%
    ungroup()
  
  return(result)
}

# Apply the function to all datasets
correlation_peak_count_all <- correlation_peak_count(final_merged_data, "Unified_all_genes")
correlation_peak_count_TOP <- correlation_peak_count(TOP_final_merged_data, "TOP1KGENE")

correlation_weighted_all <- correlation_weighted(final_merged_data, "Unified_all_genes")
correlation_weighted_TOP <- correlation_weighted(TOP_final_merged_data, "TOP1KGENE")


cor_matrix_peak_count_all <- correlation_peak_count_all %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_peak_count_all) <- cor_matrix_peak_count_all$celltype
cor_matrix_peak_count_all <- as.matrix(cor_matrix_peak_count_all[, -1])
cor_matrix_peak_count_all



cor_matrix_peak_count_TOP <- correlation_peak_count_TOP %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_peak_count_TOP) <- cor_matrix_peak_count_TOP$celltype
cor_matrix_peak_count_TOP <- as.matrix(cor_matrix_peak_count_TOP[, -1])
cor_matrix_peak_count_TOP





#correlation_weighted_TOP1KGENE
cor_matrix_weighted_all <- correlation_weighted_all %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_weighted_all) <- cor_matrix_weighted_all$celltype
cor_matrix_weighted_all <- as.matrix(cor_matrix_weighted_all[, -1])
cor_matrix_weighted_all

cor_matrix_weighted_TOP <- correlation_weighted_TOP %>%
  dplyr::select(celltype, peak_caller, spearman_cor) %>%
  tidyr::pivot_wider(names_from = peak_caller, values_from = spearman_cor) %>%
  as.data.frame()
rownames(cor_matrix_weighted_TOP) <- cor_matrix_weighted_TOP$celltype
cor_matrix_weighted_TOP <- as.matrix(cor_matrix_weighted_TOP[, -1])
cor_matrix_weighted_TOP




####################################################################
# Enhanced Nature Genetics Style Theme
####################################################################

theme_nature_genetics <- function(base_size = 7, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black", size = base_size),
      axis.line = element_line(color = "black", linewidth = 0.25),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = unit(0.5, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1, 
                               hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, 
                                  color = "grey40", margin = margin(b = 4)),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(face = "bold", size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      legend.key.size = unit(3, "mm"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = unit(1, "mm"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

theme_set(theme_nature_genetics())

nature_genetics_palette <- c(
  "#E64B35", "#8bb5bfff", "#0cc3a7ff", "#4e9a6aff", 
  "#274ab1ff", "#4e77e7ff", "#3e4f4bff", "#13c758ff"
)



##=====================================================================
# NEW FIGURE 4: Performance Summary - Boxplot Style
##=====================================================================
correlation_peak_count_all <- correlation_peak_count_all %>% filter(is.finite(spearman_cor))
correlation_peak_count_TOP<- correlation_peak_count_TOP %>% filter(is.finite(spearman_cor))

correlation_weighted_all <- correlation_weighted_all %>% filter(is.finite(spearman_cor))
correlation_weighted_TOP<- correlation_weighted_TOP %>% filter(is.finite(spearman_cor))



# Clean data
correlation_weighted_all_clean <- correlation_weighted_all %>%
  filter(!is.na(spearman_cor) & is.finite(spearman_cor))

# Define FIXED y-axis order (vertical axis)
method_order <- c( "SICER2", "SEACR", "MACS2",  "HOMER", "GoPeaks", "Genrich", "DROMPAplus")

# Calculate mean/median for ordering by performance
method_stats <- correlation_weighted_all_clean %>%
  group_by(peak_caller) %>%
  summarise(
    mean_cor = mean(spearman_cor, na.rm = TRUE),
    median_cor = median(spearman_cor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cor))  # Best to worst

cat("Method ranking by mean correlation:\n")
print(method_stats)

# OPTION A: Order by performance (best to worst) - RECOMMENDED
performance_order <- as.character(method_stats$peak_caller)

# OPTION B: Use fixed order
fixed_order <- method_order

# Choose which order you want
selected_order <- fixed_order  # Change to performance_order if you want ranking

# Create plot with HORIZONTAL boxplots (no coord_flip)
performance_boxplot_weighted_all <- ggplot(
  correlation_weighted_all_clean, 
  aes(x = spearman_cor, 
      y = factor(peak_caller, levels = selected_order))  # Y-axis: fixed method order
) +
  
  geom_boxplot(
    fill = "grey90", 
    color = "black", 
    linewidth = 0.3, 
    outlier.size = 0.5,
    orientation = "y"  # Horizontal boxplot
  ) +
  
  geom_point(
    aes(color = celltype), 
    size = 1, 
    alpha = 0.7, 
    position = position_jitter(height = 0.2)
  ) +
  
  scale_color_manual(values = nature_genetics_palette) +
  scale_y_discrete(limits = selected_order) +  # Fix y-axis order
  scale_x_continuous(
    limits = c(0, 0.7),  # সীমা ঠিক করে দিলাম
    breaks = seq(0, 0.7, by = 0.1)  # 0, 0.1, 0.2, ..., 0.7
  ) +
  labs(
    title = "Peak Caller Performance Distribution", 
    x = "Spearman Correlation",  # X-axis: correlation value
    y = "Peak Caller",           # Y-axis: methods in fixed order
    color = "Cell Type"
  ) +
  
  theme_nature_genetics() +  # Use your theme
  theme(
    legend.position = "bottom"
  )

# Save plot
ggsave(
  file.path(out_dir, paste0("Fgure5B_Performance_Boxplot_", hist, "_byweighted_ALL_with_input.pdf")), 
  plot = performance_boxplot_weighted_all, 
  width = 10,    # Wider for horizontal
  height = 8,  
  units = "cm", 
  dpi = 1200
)

}


###########################################################################################################
## Figure 5 C-E
###########################################################################################################

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
library(ggpubr)        # For publication-ready themes
library(viridis)       # For better color scales
library(ggsci)         # For Nature-style color palettes
library(scales)        # For better axis formatting
library(conflicted)    # For namespace management

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
# Define histones to process
#histones <- c( "H3K27ac-b", "H3K27ac-s")

histones <- c( "H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
###===========================================================================
# loading the datasets
###===========================================================================
# 1. Load the pre-aggregated TPM expression matrix
for (hist in histones){
TPM <- read.delim("/home/wahid/project_scHMTF/GSE195725_processed_data/scRNA_seq/10x_scRNAseq/pbmc_celltype_averaged_TPM.txt.gz", row.names = 1)

# 3. Define parameters
celltypes <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")
peak_callers <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
#sc_peak_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/result/l1_withcontrol/H3K27me3_peakbed"
sc_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_with_input_peakbed_corrected/scrnaseq_peakbed/", hist, "_peakbed")


# 4. Create output directory
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)




###===========================================================================
# Comprehensive BED File Reading Function for each peak caller - FIXED SICER2
###===========================================================================
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
##=====================================================================
# FIXED: Process All Peak Files with Proper RP Model Implementation
# MODIFIED: Ensure same number of gene activities for all peak callers
##=====================================================================
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



##=====================================================================
# Prepare and Merge with Expression Data - FIXED with retry logic
##=====================================================================

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


################TOP1000genes_final_merged_data by peak caller and celltype [Top gene selected by higher peaks]

# Get top 1000 genes for each peak_caller method based on TPM expression
TOP_final_merged_data <- final_merged_data %>%
  group_by(peak_caller,celltype) %>%
  arrange(desc(weighted_peak_count), .by_group = TRUE) %>%
  slice_head(n = 1000) %>%
  ungroup()




##=====================================================================
# FIGURE 07: Cell Type Specificity Analysis for TOP1000 Genes
# Enhanced with professional styling and optimal visualization
##=====================================================================

# Calculate cell type specificity scores (Shannon Entropy-based)
specificity_analysis2 <- TOP_final_merged_data %>%
  group_by(gene_name, peak_caller) %>%
  mutate(
    # Add small epsilon to avoid log(0)
    tpm_normalized = tpm_expression + 1e-10,
    total_expression = sum(tpm_normalized, na.rm = TRUE),
    expression_fraction = tpm_normalized / total_expression,
    # Calculate Shannon entropy (measure of distribution uniformity)
    entropy = -sum(expression_fraction * log(expression_fraction), na.rm = TRUE),
    # Convert entropy to specificity (1 - normalized entropy)
    max_entropy = log(n()),  # Maximum possible entropy for n cell types
    normalized_entropy = entropy / max_entropy,
    specificity_index = 1 - normalized_entropy  # Normalized
  ) %>%
  ungroup() %>%
  group_by(celltype, peak_caller) %>%
  summarise(
    mean_specificity = mean(specificity_index, na.rm = TRUE),
    median_specificity = median(specificity_index, na.rm = TRUE),
    n_genes = n(),
    .groups = 'drop'
  )

# Calculate the actual min and max values for the color scale
#min_specificity <- min(specificity_analysis2$mean_specificity, na.rm = TRUE)
#max_specificity <- max(specificity_analysis2$mean_specificity, na.rm = TRUE)
min_specificity <- 0.03
max_specificity <- 0.43
# Create custom breaks for the color scale
specificity_breaks <- pretty(c(min_specificity, max_specificity), n = 5)
midpoint <- mean(c(min_specificity, max_specificity))

# Reorder factors for better organization (if needed)
# specificity_analysis2$peak_caller <- factor(specificity_analysis2$peak_caller, 
#                                            levels = c("Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2"))

# Create a custom color palette for better contrast
custom_magma <- colorRampPalette(cividis(12))(256)

# Specificity heatmap with professional aesthetics
specificity_heatmap2 <- ggplot(specificity_analysis2, 
                              aes(x = peak_caller, y = celltype, fill = mean_specificity)) +
  # Background tile for subtle grid effect
  geom_tile(aes(x = peak_caller, y = celltype), 
            width = 1, height = 1, fill = "grey95", alpha = 0.3) +
  # Main heatmap tiles with elegant border
  geom_tile(color = "white", size = 0.8, width = 0.92, height = 0.92) +
  
  # Value labels with optimal contrast
  geom_text(aes(label = sprintf("%.3f", mean_specificity),
                color = ifelse(mean_specificity > midpoint, "white", "black")), 
            size = 3.8, fontface = "bold", family = "Helvetica") +
  scale_color_identity() +
  
  # Enhanced color scale with professional legend
  scale_fill_gradientn(
    name = "Specificity Index",
    colors = custom_magma,
    limits = c(min_specificity, max_specificity),
    breaks = specificity_breaks,
    labels = function(x) sprintf("%.3f", x),
    guide = guide_colorbar(
      barwidth = 18, 
      barheight = 0.8, 
      title.position = "top", 
      title.hjust = 0.5, 
      direction = "horizontal",
      frame.colour = "black", 
      ticks.colour = "black",
      frame.linewidth = 0.8,
      ticks.linewidth = 0.8,
      title.theme = element_text(
        size = 11, face = "bold", family = "Helvetica", 
        margin = margin(b = 4), color = "#2c3e50"
      ),
      label.theme = element_text(size = 9, family = "Helvetica", color = "#2c3e50")
    )
  ) +
  
  # Clear, informative labels
  labs(
    x = "Peak Caller Method", 
    y = "Cell Type", 
    title = paste0("Cell-Type Specificity of TOP1000 Genes with ", hist, " Peaks"),
    subtitle = "Higher values = more cell-type-specific expression"
  ) +
  
  # Professional theme with refined styling
  theme_minimal(base_size = 13, base_family = "Helvetica") +
  theme(
    # Axis styling
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = 11, 
      face = "bold", color = "#2c3e50", 
      margin = margin(t = 3)
    ),
    axis.text.y = element_text(
      size = 11, face = "bold", hjust = 1, 
      color = "#2c3e50", margin = margin(r = 3)
    ),
    axis.title.x = element_text(
      face = "bold", size = 12, color = "#2c3e50",
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      face = "bold", size = 12, color = "#2c3e50",
      margin = margin(r = 8)
    ),
    
    # Panel and grid styling
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Legend positioning and styling
    legend.position = "top",
    legend.margin = margin(b = 8),
    legend.box.margin = margin(b = 5),
    
    # Title and subtitle styling
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 16, 
      margin = margin(b = 6), color = "#2c3e50"
    ),
    plot.subtitle = element_text(
      hjust = 0.5, size = 12, color = "#7f8c8d",
      lineheight = 1.2, margin = margin(b = 10)
    ),
    plot.caption = element_text(
      hjust = 0.5, size = 10, color = "#95a5a6",
      margin = margin(t = 12), face = "italic"
    ),
    plot.margin = margin(20, 25, 20, 25)
  )

# Save with optimal dimensions for publication
ggsave(file.path(out_dir, paste0("Figure7_CellType_Specificity_TOP_", hist, "_with_input.pdf")), plot = specificity_heatmap2, 
       width = 19,  # Optimal width for clear labels
       height = 15, # Optimal height for good proportions
       units = "cm", 
       dpi = 1200, 
       device = cairo_pdf)


##=====================================================================
# SAVE SPECIFICITY_ANALYSIS2 DATAFRAME
##=====================================================================
cat("\n=== SAVING SPECIFICITY_ANALYSIS2 DATAFRAME ===\n")

# Save specificity_analysis2 as CSV
specificity_csv_file <- file.path(out_dir, paste0("HumanPBMC_specificity_", hist, "_with_input.csv"))
write.csv(specificity_analysis2, specificity_csv_file, row.names = FALSE)
cat("✓ specificity_analysis2 saved as CSV:", specificity_csv_file, "\n")

# Save as RDS (preserves data types and structure)
specificity_rds_file <- file.path(out_dir, paste0("HumanPBMC_specificity_", hist, "_with_input.rds"))
saveRDS(specificity_analysis2, specificity_rds_file)
cat("✓ specificity_analysis2 saved as RDS:", specificity_rds_file, "\n")

# Save with metadata (most complete version)
specificity_results_complete <- list(
  specificity_data = specificity_analysis2,
  metadata = list(
    histone = hist,
    cell_types = celltypes,
    peak_callers = peak_callers,
    calculation_date = Sys.time(),
    n_celltypes = length(unique(specificity_analysis2$celltype)),
    n_peak_callers = length(unique(specificity_analysis2$peak_caller)),
    n_observations = nrow(specificity_analysis2),
    specificity_range = c(
      min = min(specificity_analysis2$mean_specificity, na.rm = TRUE),
      max = max(specificity_analysis2$mean_specificity, na.rm = TRUE)
    )
  )
)

saveRDS(specificity_results_complete,
        file.path(out_dir, 
                 paste0("specificity_complete_", hist, "_with_input.rds")))
cat("✓ Complete specificity results saved with metadata\n")

}

##=====================================================================
# Error bar plot 
##=====================================================================
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Define the directory where your files are located
input_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input"

# List of files with full paths
files <- c(
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/HumanPBMC_specificity_H3K27ac-b_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/HumanPBMC_specificity_H3K27ac-s_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/HumanPBMC_specificity_H3K27me3_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/HumanPBMC_specificity_H3K4me1_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/HumanPBMC_specificity_H3K4me2_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/HumanPBMC_specificity_H3K4me3_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/result_GT_scrnaseq_with_input/HumanPBMC_specificity_H3K9me3_with_input.csv"
)

# Create an empty list to store dataframes
all_data <- list()

cat("Starting to merge specificity files...\n")
cat("==========================================\n")

# Read and process each file
for (file in files) {
  if (file.exists(file)) {
    # Extract histone mark from filename
    # Pattern: HumanPBMC_specificity_(histone-mark)_with_input.csv
    histone_mark <- str_extract(file, "(?<=HumanPBMC_specificity_)[^_]+(?=_with_input\\.csv)")
    
    if (is.na(histone_mark)) {
      # Try alternative pattern if needed
      histone_mark <- str_extract(file, "H3K[0-9a-zA-Z-]+(?=_with_input)")
    }
    
    cat("Processing:", basename(file), "\n")
    cat("Histone mark extracted:", histone_mark, "\n")
    
    # Read the CSV file
    data <- read_csv(file, col_types = cols())
    
    # Check the structure of the data
    cat("  Columns found:", paste(colnames(data), collapse = ", "), "\n")
    cat("  Number of rows:", nrow(data), "\n")
    
    # Add histone mark column as first column
    data <- data %>%
      mutate(histone = histone_mark) %>%
      select(histone, everything())  # Move histone to first column
    
    # Add to list
    all_data[[file]] <- data
    
    # Show first few rows
    cat("  First few rows:\n")
    print(head(data, 3))
    cat("\n")
    
  } else {
    cat("ERROR: File not found:", file, "\n")
  }
}

cat("==========================================\n")
cat("All files processed. Now combining...\n")

# Check if we have any data
if (length(all_data) == 0) {
  stop("No files were successfully read. Please check file paths.")
}

# Combine all dataframes
combined_specificity_data <- bind_rows(all_data)

# View the structure of combined data
cat("\nCombined data structure:\n")
cat("Number of rows:", nrow(combined_specificity_data), "\n")
cat("Number of columns:", ncol(combined_specificity_data), "\n")
cat("Column names:", paste(colnames(combined_specificity_data), collapse = ", "), "\n")

# Count unique values
cat("\nUnique values:\n")
cat("Histones:", paste(unique(combined_specificity_data$histone), collapse = ", "), "\n")
cat("Cell types:", paste(unique(combined_specificity_data$celltype), collapse = ", "), "\n")
cat("Peak callers:", paste(unique(combined_specificity_data$peak_caller), collapse = ", "), "\n")

# View first few rows of combined data
cat("\nFirst 10 rows of combined data:\n")
print(head(combined_specificity_data, 10))

# Save the combined data
output_file <- file.path(input_dir, "HumanPBMC_ALL_histones_specificity_combined_with_input.csv")
write_csv(combined_specificity_data, output_file)
cat("\n✓ Combined data saved to:", output_file, "\n")


# Create specificity summary data (corrected)
specificity_summary <- combined_specificity_data %>%
  group_by(peak_caller, histone) %>%
  summarize(
    mean_specificity_overall = mean(mean_specificity, na.rm = TRUE),
    median_specificity_overall = median(median_specificity, na.rm = TRUE),
    sd_specificity = sd(mean_specificity, na.rm = TRUE),
    se_specificity = sd_specificity / sqrt(n()),
    n_celltypes = n(),
    total_genes = sum(n_genes),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(peak_caller),
    Histone = factor(histone)
  )

# Check the summary
cat("\nSpecificity summary structure:\n")
glimpse(specificity_summary)

cat("\nFirst few rows of specificity summary:\n")
print(head(specificity_summary))



all_methods <- unique(data$peak_caller)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)


library(readr)      # For reading CSV
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting
library(scales)     # For percentage formatting
library(ggsci)    

###############################################################################
# SCATTER PLOT 2: With error bars
###############################################################################

specificity_scatter2 <- ggplot(specificity_summary, 
                       aes(x = Histone, 
                           y = mean_specificity_overall,
                           color = Method)) +
  geom_point(size = 10, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pmax(mean_specificity_overall - sd_specificity, 0),  # Cap at 0
                    ymax = pmin(mean_specificity_overall + sd_specificity, 1)), # Cap at 1
                width = 0.3, 
                position = position_dodge(width = 0.5),
                alpha = 0.7) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.03, .43),  # Explicitly set 0-100% range
                     breaks = seq(0.03, .43, by = 0.1)) +  # 0%, 20%, 40%, etc.
  labs(
    title = "Specificity Values with Variability by Histone and Method",
    x = "Histone Modification",
    y = "Specificity (%)"
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
print(specificity_scatter2)


ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_specificity/HumanPBMC_scRNA_Specificity_Scatter_Plot_With_Input.pdf",
       plot = specificity_scatter2,
       width = 12,
       height = 9,
       device = "pdf",
       bg = "white")








###########################################################################################################
## Figure 5 F-H 
###########################################################################################################

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
library(ggpubr)        # For publication-ready themes
library(viridis)       # For better color scales
library(ggsci)         # For Nature-style color palettes
library(scales)        # For better axis formatting
library(conflicted)    # For namespace management
library(gridExtra)     # For arranging multiple plots

# Resolve namespace conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("slice_head", "dplyr")

# Function to get gene mapping with retry logic
get_gene_mapping_with_retry <- function(ensembl_ids, max_retries = 3) {
  for (attempt in 1:max_retries) {
    tryCatch({
      # Get gene information from EnsDb
      gene_info <- genes(EnsDb.Hsapiens.v86, 
                        filter = GeneIdFilter(ensembl_ids),
                        columns = c("gene_id", "gene_name", "symbol"))
      
      # Convert to data frame and clean up
      gene_df <- as.data.frame(gene_info) %>%
        dplyr::select(ensembl_gene_id = gene_id, 
                     hgnc_symbol = symbol, 
                     external_gene_name = gene_name) %>%
        distinct()
      
      return(gene_df)
      
    }, error = function(e) {
      if (attempt == max_retries) {
        cat("Failed after", max_retries, "attempts:", conditionMessage(e), "\n")
        
        # Fallback: return basic mapping
        return(data.frame(
          ensembl_gene_id = ensembl_ids,
          hgnc_symbol = NA_character_,
          external_gene_name = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
      cat("Attempt", attempt, "failed, retrying...\n")
      Sys.sleep(2) # Wait before retrying
    })
  }
}

 # Add this line
#histones <- c( "H3K27ac-b")
histones <- c( "H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
###==========================================================================
# loading the datasets
###==========================================================================
# 1. Load the pre-aggregated TPM expression matrix
for (hist in histones){
TPM <- read.delim("/home/wahid/project_scHMTF/GSE195725_processed_data/bulk_RNA_seq/GSE107011_Processed_data_AggregateTPM.txt.gz", row.names = 1)
print("TPM matrix dimensions:")
dim(TPM)
print("Column names (cell types):")
colnames(TPM)


# 3. Define parameters
celltypes <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")
peak_callers <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
sc_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_with_input_peakbed_corrected/bulkrnaseq_peakbed/", hist, "_peakbed")

# 4. Create output directory
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###==========================================================================
# Comprehensive BED File Reading Function for each peak caller - FIXED SICER2
###==========================================================================
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
      if (n_cols >= 9) {
        peaks$score <- as.numeric(bed_data[, 7]) # Use fold change
        mcols(peaks)$name <- bed_data[, 4]
        mcols(peaks)$strand <- bed_data[, 6]
        mcols(peaks)$fold_change <- as.numeric(bed_data[, 7])
        mcols(peaks)$neg_log10_pvalue <- as.numeric(bed_data[, 8])
        mcols(peaks)$neg_log10_qvalue <- as.numeric(bed_data[, 9])
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
###==========================================================================
# FIXED: Process All Peak Files with Proper RP Model Implementation
# MODIFIED: Ensure same number of gene activities for all peak callers
###==========================================================================
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



###==========================================================================
# Prepare and Merge with Expression Data - CORRECTED
###==========================================================================

# Resolve function conflicts at the beginning
conflicts_prefer(dplyr::slice)
conflicts_prefer(dplyr::slice_head)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::group_by)
conflicts_prefer(dplyr::ungroup)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::arrange)

# First, let's fix the TPM data structure
print("Original TPM matrix structure:")
print(dim(TPM))
print("First few rows:")
print(TPM[1:5, 1:5])

# The issue: your TPM data has a 'gene' column AND row names
# Let's extract the proper gene names from the 'gene' column
TPM_with_genes <- TPM %>%
  as.data.frame() %>%
  # Use the 'gene' column as the gene identifier
  mutate(ensembl_gene_id = sub("\\..*", "", gene)) %>%
  dplyr::select(-gene) %>%  # Remove the original gene column
  # Convert all sample columns to numeric
  mutate(across(all_of(celltypes), as.numeric))

print("Fixed TPM_with_genes structure:")
print(TPM_with_genes[1:5, 1:5])

# Check for duplicates in TPM data
duplicate_genes <- TPM_with_genes %>%
  dplyr::count(ensembl_gene_id) %>%
  dplyr::filter(n > 1) %>%
  nrow()

print(paste("Number of duplicate gene entries in TPM data:", duplicate_genes))

# Remove duplicates from TPM data by taking the mean expression
TPM_with_genes_unique <- TPM_with_genes %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(across(all_of(celltypes), mean, na.rm = TRUE)) %>%
  dplyr::ungroup()

print("TPM data after removing duplicates:")
print(dim(TPM_with_genes_unique))

# Now pivot to long format
TPM_long <- TPM_with_genes_unique %>%
  pivot_longer(
    cols = -ensembl_gene_id,
    names_to = "celltype",
    values_to = "tpm_expression"
  ) %>%
  dplyr::filter(celltype %in% celltypes)

print("TPM_long structure:")
print(head(TPM_long))

# Check for duplicates in TPM_long
duplicate_pairs <- TPM_long %>%
  dplyr::count(ensembl_gene_id, celltype) %>%
  dplyr::filter(n > 1) %>%
  nrow()

print(paste("Number of duplicate gene-celltype pairs in TPM_long:", duplicate_pairs))

# Remove version numbers from peak activity gene IDs
all_peak_activity <- all_peak_activity %>%
  dplyr::mutate(ensembl_gene_id = sub("\\..*", "", gene_id))

print("Peak activity gene ID examples:")
print(head(unique(all_peak_activity$ensembl_gene_id)))
print("TPM gene ID examples:")
print(head(unique(TPM_long$ensembl_gene_id)))

# Check for common genes
common_genes <- intersect(unique(all_peak_activity$ensembl_gene_id), 
                         unique(TPM_long$ensembl_gene_id))
print(paste("Number of common genes:", length(common_genes)))

# Get gene mapping for the common genes only
gene_mapping <- get_gene_mapping_with_retry(common_genes)

# Simplify the mapping to avoid many-to-many relationships
# Keep only one mapping per Ensembl ID (prefer HGNC symbol)
gene_mapping_unique <- gene_mapping %>%
  dplyr::group_by(ensembl_gene_id) %>%
  # Prefer non-empty HGNC symbols, then external names
  dplyr::arrange(desc(hgnc_symbol != ""), desc(external_gene_name != "")) %>%
  dplyr::slice(1) %>%  # Explicitly use dplyr::slice
  dplyr::ungroup()

print("Unique gene mapping results:")
print(head(gene_mapping_unique))

# Check for duplicates in peak data
peak_duplicates <- all_peak_activity %>%
  dplyr::count(ensembl_gene_id, celltype, peak_caller) %>%
  dplyr::filter(n > 1) %>%
  nrow()

print(paste("Number of duplicate gene-celltype-peak_caller entries in peak data:", peak_duplicates))

# If there are duplicates in peak data, aggregate them
if (peak_duplicates > 0) {
  all_peak_activity <- all_peak_activity %>%
    dplyr::group_by(ensembl_gene_id, celltype, peak_caller, gene_id, gene_name) %>%
    dplyr::summarise(
      peak_count = sum(peak_count, na.rm = TRUE),
      weighted_peak_count = sum(weighted_peak_count, na.rm = TRUE),
      .groups = "drop"
    )
  print("Aggregated duplicate peak entries")
}

# Merge the data step by step
# First merge peak activity with gene mapping
peak_with_names <- all_peak_activity %>%
  dplyr::inner_join(gene_mapping_unique, by = "ensembl_gene_id", relationship = "many-to-one") %>%
  # Use the best available gene name
  dplyr::mutate(gene_name_final = coalesce(hgnc_symbol, external_gene_name, gene_name))

print("Peak data with gene names:")
print(head(peak_with_names))

# Now merge with TPM data - use left_join instead of inner_join to see what's happening
merged_data <- peak_with_names %>%
  dplyr::left_join(TPM_long, by = c("ensembl_gene_id", "celltype"), relationship = "many-to-one")

print(paste("After left join:", nrow(merged_data), "rows"))
print(paste("Rows with TPM data:", sum(!is.na(merged_data$tpm_expression))))

# Filter and transform
final_merged_data <- merged_data %>%
  dplyr::filter(!is.na(tpm_expression) & !is.na(peak_count)) %>%
  dplyr::mutate(
    log_tpm = log2(tpm_expression + 1),
    log_peak_count = log10(peak_count + 1),
    log_weighted_count = log10(weighted_peak_count + 1)
  )

print(paste("Final dataset:", nrow(final_merged_data), "gene-celltype-peakcaller pairs"))
print("Final merged data structure:")
print(head(final_merged_data))

# Get top 1000 genes for each peak_caller method based on weighted peak count
TOP_final_merged_data <- final_merged_data %>%
  dplyr::group_by(peak_caller, celltype) %>%
  dplyr::arrange(desc(weighted_peak_count), .by_group = TRUE) %>%
  dplyr::slice_head(n = 1000) %>%  # Explicitly use dplyr::slice_head
  dplyr::ungroup()

print(paste("Top dataset:", nrow(TOP_final_merged_data), "rows"))


###==========================================================================
#  Cell Type Specificity Analysis for TOP1000 Genes
# Enhanced with professional styling and optimal visualization
###==========================================================================

# Calculate cell type specificity scores (Shannon Entropy-based)
specificity_analysis2 <- TOP_final_merged_data %>%
  group_by(gene_name, peak_caller) %>%
  mutate(
    # Add small epsilon to avoid log(0)
    tpm_normalized = tpm_expression + 1e-10,
    total_expression = sum(tpm_normalized, na.rm = TRUE),
    expression_fraction = tpm_normalized / total_expression,
    # Calculate Shannon entropy (measure of distribution uniformity)
    entropy = -sum(expression_fraction * log(expression_fraction), na.rm = TRUE),
    # Convert entropy to specificity (1 - normalized entropy)
    max_entropy = log(n()),  # Maximum possible entropy for n cell types
    normalized_entropy = entropy / max_entropy,
    specificity_index = 1 - normalized_entropy  # Now ranges 0-1
  ) %>%
  ungroup() %>%
  group_by(celltype, peak_caller) %>%
  summarise(
    mean_specificity = mean(specificity_index, na.rm = TRUE),
    median_specificity = median(specificity_index, na.rm = TRUE),
    n_genes = n(),
    .groups = 'drop'
  )

# Calculate the actual min and max values for the color scale
min_specificity <- 0.035
max_specificity <- 0.320
#min_specificity <- min(specificity_analysis2$mean_specificity, na.rm = TRUE)
#max_specificity <- max(specificity_analysis2$mean_specificity, na.rm = TRUE)
# Create custom breaks for the color scale
specificity_breaks <- pretty(c(min_specificity, max_specificity), n = 5)
midpoint <- mean(c(min_specificity, max_specificity))

# Reorder factors for better organization (if needed)
# specificity_analysis2$peak_caller <- factor(specificity_analysis2$peak_caller, 
#                                            levels = c("Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2"))

# Create a custom color palette for better contrast
custom_magma <- colorRampPalette(viridisLite::magma(12))(256)

# Specificity heatmap with professional aesthetics
specificity_heatmap2 <- ggplot(specificity_analysis2, 
                              aes(x = peak_caller, y = celltype, fill = mean_specificity)) +
  # Background tile for subtle grid effect
  geom_tile(aes(x = peak_caller, y = celltype), 
            width = 1, height = 1, fill = "grey95", alpha = 0.3) +
  # Main heatmap tiles with elegant border
  geom_tile(color = "white", size = 0.8, width = 0.92, height = 0.92) +
  
  # Value labels with optimal contrast
  geom_text(aes(label = sprintf("%.3f", mean_specificity),
                color = ifelse(mean_specificity > midpoint, "white", "black")), 
            size = 3.8, fontface = "bold", family = "Helvetica") +
  scale_color_identity() +
  
  # Enhanced color scale with professional legend
  scale_fill_gradientn(
    name = "Specificity Index",
    colors = custom_magma,
    limits = c(min_specificity, max_specificity),
    breaks = specificity_breaks,
    labels = function(x) sprintf("%.3f", x),
    guide = guide_colorbar(
      barwidth = 18, 
      barheight = 0.8, 
      title.position = "top", 
      title.hjust = 0.5, 
      direction = "horizontal",
      frame.colour = "black", 
      ticks.colour = "black",
      frame.linewidth = 0.8,
      ticks.linewidth = 0.8,
      title.theme = element_text(
        size = 11, face = "bold", family = "Helvetica", 
        margin = margin(b = 4), color = "#2c3e50"
      ),
      label.theme = element_text(size = 9, family = "Helvetica", color = "#2c3e50")
    )
  ) +
  
  # Clear, informative labels
  labs(
    x = "Peak Caller Method", 
    y = "Cell Type", 
    title = paste0("Cell-Type Specificity of TOP1000 Genes with ", hist, " Peaks")
  ) +
  
  # Professional theme with refined styling
  theme_minimal(base_size = 13, base_family = "Helvetica") +
  theme(
    # Axis styling
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = 11, 
      face = "bold", color = "#2c3e50", 
      margin = margin(t = 3)
    ),
    axis.text.y = element_text(
      size = 11, face = "bold", hjust = 1, 
      color = "#2c3e50", margin = margin(r = 3)
    ),
    axis.title.x = element_text(
      face = "bold", size = 12, color = "#2c3e50",
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      face = "bold", size = 12, color = "#2c3e50",
      margin = margin(r = 8)
    ),
    
    # Panel and grid styling
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Legend positioning and styling
    legend.position = "top",
    legend.margin = margin(b = 8),
    legend.box.margin = margin(b = 5),
    
    # Title and subtitle styling
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 16, 
      margin = margin(b = 6), color = "#2c3e50"
    ),
    plot.subtitle = element_text(
      hjust = 0.5, size = 12, color = "#7f8c8d",
      lineheight = 1.2, margin = margin(b = 10)
    ),
    plot.caption = element_text(
      hjust = 0.5, size = 10, color = "#95a5a6",
      margin = margin(t = 12), face = "italic"
    ),
    plot.margin = margin(20, 25, 20, 25)
  )

# Save with optimal dimensions for publication
ggsave(file.path(out_dir, paste0("Figure7_CellType_Specificity_TOP_", hist, "_with_input.pdf")), plot = specificity_heatmap2, 
       width = 19,  # Optimal width for clear labels
       height = 15, # Optimal height for good proportions
       units = "cm", 
       dpi = 1200, 
       device = cairo_pdf)

###==========================================================================
# SAVE SPECIFICITY_ANALYSIS2 DATAFRAME
###==========================================================================
cat("\n=== SAVING SPECIFICITY_ANALYSIS2 DATAFRAME ===\n")

# Save specificity_analysis2 as CSV
specificity_csv_file <- file.path(out_dir, paste0("HumanPBMC_specificity_", hist, "_with_input.csv"))
write.csv(specificity_analysis2, specificity_csv_file, row.names = FALSE)
cat("✓ specificity_analysis2 saved as CSV:", specificity_csv_file, "\n")

# Save as RDS (preserves data types and structure)
specificity_rds_file <- file.path(out_dir, paste0("HumanPBMC_specificity_", hist, "_with_input.rds"))
saveRDS(specificity_analysis2, specificity_rds_file)
cat("✓ specificity_analysis2 saved as RDS:", specificity_rds_file, "\n")

# Save with metadata (most complete version)
specificity_results_complete <- list(
  specificity_data = specificity_analysis2,
  metadata = list(
    histone = hist,
    cell_types = celltypes,
    peak_callers = peak_callers,
    calculation_date = Sys.time(),
    n_celltypes = length(unique(specificity_analysis2$celltype)),
    n_peak_callers = length(unique(specificity_analysis2$peak_caller)),
    n_observations = nrow(specificity_analysis2),
    specificity_range = c(
      min = min(specificity_analysis2$mean_specificity, na.rm = TRUE),
      max = max(specificity_analysis2$mean_specificity, na.rm = TRUE)
    )
  )
)

saveRDS(specificity_results_complete,
        file.path(out_dir, 
                 paste0("specificity_complete_", hist, "_with_input.rds")))
cat("✓ Complete specificity results saved with metadata\n")

}


###==================================================================================
# Figure 5 G
###==================================================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Define the directory where your files are located
input_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input"

# List of files with full paths
files <- c(
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input/HumanPBMC_specificity_H3K27ac-b_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input/HumanPBMC_specificity_H3K27ac-s_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input/HumanPBMC_specificity_H3K27me3_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input/HumanPBMC_specificity_H3K4me1_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input/HumanPBMC_specificity_H3K4me2_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input/HumanPBMC_specificity_H3K4me3_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/result_GT_bulkrnaseq_with_input/HumanPBMC_specificity_H3K9me3_with_input.csv"
)

# Create an empty list to store dataframes
all_data <- list()

cat("Starting to merge specificity files...\n")
cat("==========================================\n")

# Read and process each file
for (file in files) {
  if (file.exists(file)) {
    # Extract histone mark from filename
    # Pattern: HumanPBMC_specificity_(histone-mark)_with_input.csv
    histone_mark <- str_extract(file, "(?<=HumanPBMC_specificity_)[^_]+(?=_with_input\\.csv)")
    
    if (is.na(histone_mark)) {
      # Try alternative pattern if needed
      histone_mark <- str_extract(file, "H3K[0-9a-zA-Z-]+(?=_with_input)")
    }
    
    cat("Processing:", basename(file), "\n")
    cat("Histone mark extracted:", histone_mark, "\n")
    
    # Read the CSV file
    data <- read_csv(file, col_types = cols())
    
    # Check the structure of the data
    cat("  Columns found:", paste(colnames(data), collapse = ", "), "\n")
    cat("  Number of rows:", nrow(data), "\n")
    
    # Add histone mark column as first column
    data <- data %>%
      mutate(histone = histone_mark) %>%
      select(histone, everything())  # Move histone to first column
    
    # Add to list
    all_data[[file]] <- data
    
    # Show first few rows
    cat("  First few rows:\n")
    print(head(data, 3))
    cat("\n")
    
  } else {
    cat("ERROR: File not found:", file, "\n")
  }
}

cat("==========================================\n")
cat("All files processed. Now combining...\n")

# Check if we have any data
if (length(all_data) == 0) {
  stop("No files were successfully read. Please check file paths.")
}

# Combine all dataframes
combined_specificity_data <- bind_rows(all_data)

# View the structure of combined data
cat("\nCombined data structure:\n")
cat("Number of rows:", nrow(combined_specificity_data), "\n")
cat("Number of columns:", ncol(combined_specificity_data), "\n")
cat("Column names:", paste(colnames(combined_specificity_data), collapse = ", "), "\n")

# Count unique values
cat("\nUnique values:\n")
cat("Histones:", paste(unique(combined_specificity_data$histone), collapse = ", "), "\n")
cat("Cell types:", paste(unique(combined_specificity_data$celltype), collapse = ", "), "\n")
cat("Peak callers:", paste(unique(combined_specificity_data$peak_caller), collapse = ", "), "\n")

# View first few rows of combined data
cat("\nFirst 10 rows of combined data:\n")
print(head(combined_specificity_data, 10))

# Save the combined data
output_file <- file.path(input_dir, "HumanPBMC_ALL_histones_specificity_combined_with_input.csv")
write_csv(combined_specificity_data, output_file)
cat("\n✓ Combined data saved to:", output_file, "\n")


# Create specificity summary data (corrected)
specificity_summary <- combined_specificity_data %>%
  group_by(peak_caller, histone) %>%
  summarize(
    mean_specificity_overall = mean(mean_specificity, na.rm = TRUE),
    median_specificity_overall = median(median_specificity, na.rm = TRUE),
    sd_specificity = sd(mean_specificity, na.rm = TRUE),
    se_specificity = sd_specificity / sqrt(n()),
    n_celltypes = n(),
    total_genes = sum(n_genes),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(peak_caller),
    Histone = factor(histone)
  )

# Check the summary
cat("\nSpecificity summary structure:\n")
glimpse(specificity_summary)

cat("\nFirst few rows of specificity summary:\n")
print(head(specificity_summary))



all_methods <- unique(data$peak_caller)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)


library(readr)      # For reading CSV
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting
library(scales)     # For percentage formatting
library(ggsci)    

###==================================================================================
# SCATTER PLOT 2: With error bars
###==================================================================================

specificity_scatter2 <- ggplot(specificity_summary, 
                       aes(x = Histone, 
                           y = mean_specificity_overall,
                           color = Method)) +
  geom_point(size = 10, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pmax(mean_specificity_overall - sd_specificity, 0),  # Cap at 0
                    ymax = pmin(mean_specificity_overall + sd_specificity, 1)), # Cap at 1
                width = 0.3, 
                position = position_dodge(width = 0.5),
                alpha = 0.7) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.03, .33),  # Explicitly set 0-100% range
                     breaks = seq(.03, .33, by = 0.1)) +  # 0%, 20%, 40%, etc.
  labs(
    title = "Specificity Values with Variability by Histone and Method",
    x = "Histone Modification",
    y = "Specificity (%)"
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
print(specificity_scatter2)


ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_bulkrnaseq_corrected_specificity/HumanPBMC_bulkRNA_Specificity_Scatter_Plot_With_Input.pdf",
       plot = specificity_scatter2,
       width = 12,
       height = 9,
       device = "pdf",
       bg = "white")






###########################################################################################################
## Figure 5 I J
###########################################################################################################

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
library(ggpubr)        # For publication-ready themes
library(viridis)       # For better color scales
library(ggsci)         # For Nature-style color palettes
library(scales)        # For better axis formatting
library(conflicted)    # For namespace management

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
# Define histones to process
#histones <- c( "H3K27ac-b", "H3K27ac-s")

histones <- c( "H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
#=======================================================================================
# loading the datasets
#=======================================================================================
# 1. Load the pre-aggregated TPM expression matrix
for (hist in histones){
TPM <- read.delim("/home/wahid/project_scHMTF/GSE195725_processed_data/scRNA_seq/10x_scRNAseq/pbmc_celltype_averaged_TPM.txt.gz", row.names = 1)

# 3. Define parameters
celltypes <- c("B", "CD4T", "CD8T", "DC", "Mono", "NK", "otherT", "other")
peak_callers <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
#sc_peak_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/result/l1_withcontrol/H3K27me3_peakbed"
sc_peak_dir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_with_input_peakbed_corrected/scrnaseq_peakbed/", hist, "_peakbed")


# 4. Create output directory
out_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)




#=======================================================================================
# Comprehensive BED File Reading Function for each peak caller - FIXED SICER2
#=======================================================================================
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
#=======================================================================================
# FIXED: Process All Peak Files with Proper RP Model Implementation
# MODIFIED: Ensure same number of gene activities for all peak callers
#=======================================================================================
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
    return(gr)   # already UCSC
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



#=======================================================================================
# Prepare and Merge with Expression Data - FIXED with retry logic
#=======================================================================================

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


# Get top 1000 genes for each peak_caller method based on TPM expression
TOP_final_merged_data <- final_merged_data %>%
  group_by(peak_caller,celltype) %>%
  arrange(desc(weighted_peak_count), .by_group = TRUE) %>%
  slice_head(n = 1000) %>%
  ungroup()



library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

#=======================================================================================
# REMOVE MT-GENES AND RIBOSOMAL GENES
#=======================================================================================
cat("\n=== REMOVING MT-GENES AND RIBOSOMAL GENES FROM DATASET ===\n")

# Store initial data for comparison
initial_rows <- nrow(final_merged_data)
initial_genes <- unique(final_merged_data$gene_name)

# Identify genes to be removed
genes_to_remove <- final_merged_data %>%
  filter(grepl("^MT-|^RP[SL]", gene_name, ignore.case = TRUE)) %>%
  pull(gene_name) %>%
  unique()

cat(sprintf("Identified %d unique MT/ribosomal genes to remove\n", length(genes_to_remove)))

# Remove MT-genes and ribosomal genes
final_merged_data <- final_merged_data %>%
  filter(!grepl("^MT-|^RP[SL]", gene_name, ignore.case = TRUE))

# Calculate statistics
removed_mt_rp <- initial_rows - nrow(final_merged_data)
cat(sprintf("  Removed %d rows containing MT/ribosomal genes\n", removed_mt_rp))
cat(sprintf("  %d rows remaining\n", nrow(final_merged_data)))
cat(sprintf("  Kept %d unique genes (removed %d unique genes)\n", 
            length(unique(final_merged_data$gene_name)),
            length(genes_to_remove)))

# Store tool order BEFORE any filtering that might remove tools
tool_order <- unique(final_merged_data$peak_caller)
cat("Tool order for plotting:", paste(tool_order, collapse = ", "), "\n")

#=======================================================================================
### Select the DEG from seurat 
#=======================================================================================
cat("\n=== LOADING AND PROCESSING SCRNA-SEQ DATA ===\n")

# Load data and set cell identities
scRNA_seq <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/scRNA_seq/10x_scRNAseq/pbmc_processed_with_TPM.rds")
Idents(scRNA_seq) <- scRNA_seq$predicted_celltype

# Find all positive markers with enhanced parameters
cat("Finding marker genes using Seurat...\n")
scRNA_seq.markers <- FindAllMarkers(
  scRNA_seq,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

cat(sprintf("Found %d total marker entries\n", nrow(scRNA_seq.markers)))
cat(sprintf("Unique marker genes: %d\n", length(unique(scRNA_seq.markers$gene))))

#=======================================================================================
### Filtering DEG
#=======================================================================================
cat("\n=== FILTERING SIGNIFICANT MARKERS ===\n")

significant_markers_all <- scRNA_seq.markers %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

# Count significant markers per cell type
cat("Significant markers per cell type (log2FC > 1, FDR < 0.05):\n")
print(table(significant_markers_all$cluster))

cat(sprintf("\nTotal significant marker entries: %d\n", nrow(significant_markers_all)))
cat(sprintf("Unique significant marker genes: %d\n", length(unique(significant_markers_all$gene))))

#=======================================================================================
### Filter final_merged_data by DEG
#=======================================================================================
cat("\n=== FILTERING FINAL_MERGED_DATA BY MARKER GENES ===\n")

# Create standardized gene name for matching
final_merged_data <- final_merged_data %>%
  mutate(gene_upper = toupper(gene_name))

significant_markers_all <- significant_markers_all %>%
  mutate(gene_upper = toupper(gene))

# Check gene name formats
cat("Checking gene name formats...\n")
cat("Sample from markers:", paste(head(significant_markers_all$gene_upper, 5), collapse = ", "), "\n")
cat("Sample from final_merged_data:", paste(head(unique(final_merged_data$gene_upper), 5), collapse = ", "), "\n")

# Filter using case-insensitive matching
rows_before <- nrow(final_merged_data)
final_merged_data <- final_merged_data %>%
  filter(gene_upper %in% significant_markers_all$gene_upper)
rows_after <- nrow(final_merged_data)

# Calculate matching statistics
matched_genes <- unique(final_merged_data$gene_upper)
expected_genes <- unique(significant_markers_all$gene_upper)

cat(sprintf("\nFiltering results:\n"))
cat(sprintf("  Rows before filtering: %d\n", rows_before))
cat(sprintf("  Rows after filtering: %d\n", rows_after))
cat(sprintf("  Rows removed: %d\n", rows_before - rows_after))
cat(sprintf("  Genes matched: %d/%d (%.1f%%)\n",
            length(matched_genes),
            length(expected_genes),
            length(matched_genes)/length(expected_genes)*100))

# Check celltype distribution
cat("\nCelltype distribution after filtering:\n")
print(table(final_merged_data$celltype))

#=======================================================================================
### VALIDATION CHECKS
#=======================================================================================
cat("\n=== VALIDATION CHECKS ===\n")

# Check if any data remains
if(nrow(final_merged_data) == 0) {
  stop("ERROR: All data filtered out! Check gene name matching.")
}

# Check if any celltypes were completely filtered out
remaining_celltypes <- unique(final_merged_data$celltype)
expected_celltypes <- unique(significant_markers_all$cluster)
missing_celltypes <- setdiff(expected_celltypes, remaining_celltypes)

if(length(missing_celltypes) > 0) {
  cat("Warning: The following celltypes have no data after filtering:\n")
  cat(paste(missing_celltypes, collapse = ", "), "\n")
}

# Check for NA values
na_counts <- colSums(is.na(final_merged_data))
if(any(na_counts > 0)) {
  cat("\nNA counts per column:\n")
  print(na_counts[na_counts > 0])
}

#=======================================================================================
# FIGURE 5 i: Celltype-Wise Top 10 Genes (EXACTLY LIKE FIGURE 3 STYLE)
#=======================================================================================
cat("\n=== CREATING FIGURE 7: CELLTYPE-WISE TOP 10 GENES ===\n")

# Create directory for celltype-wise plots
celltype_wise_dir <- file.path(out_dir, paste0("celltype_wise_plots_", hist))
if (!dir.exists(celltype_wise_dir)) dir.create(celltype_wise_dir, recursive = TRUE)

# Get all unique cell types
all_celltypes <- unique(final_merged_data$celltype)

cat("Creating Figure 7 style plots for", length(all_celltypes), "cell types:\n")

# Create list to store plots for each cell type
plot_list_celltype <- list()

for (celltype in all_celltypes) {
  cat("  Creating plot for", celltype, "... ")
  
  # Filter data for this cell type
  celltype_data <- final_merged_data %>%
    filter(celltype == !!celltype)
  
  if (nrow(celltype_data) == 0) {
    cat("No data found\n")
    next
  }
  
  # Check for and handle duplicates by taking the mean
  celltype_data_clean <- celltype_data %>%
    group_by(gene_name, peak_caller) %>%
    summarise(
      log_tpm = mean(log_tpm, na.rm = TRUE),
      log_weighted_count = mean(log_weighted_count, na.rm = TRUE),
      tpm_expression = mean(tpm_expression, na.rm = TRUE),
      weighted_peak_count = mean(weighted_peak_count, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Get top 10 genes for this cell type (based on TPM)
  top_genes_celltype <- celltype_data_clean %>%
    group_by(gene_name) %>%
    summarise(
      av_log_tpm = mean(log_tpm, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(av_log_tpm)) %>%
    slice_head(n = 10)
  
  if (nrow(top_genes_celltype) == 0) {
    cat("No valid TPM data\n")
    next
  }
  
  # Get data for top 10 genes - use cleaned data
  plot_data <- celltype_data_clean %>%
    filter(gene_name %in% top_genes_celltype$gene_name) %>%
    # Ensure we have all tools for comparison
    complete(gene_name, peak_caller, fill = list(
      weighted_peak_count = 0,
      log_weighted_count = 0,
      tpm_expression = 0,
      log_tpm = 0
    )) %>%
    # Reorder genes by average TPM
    mutate(gene_name = factor(gene_name, levels = rev(top_genes_celltype$gene_name))) %>%
    # Reorder tools
    mutate(peak_caller = factor(peak_caller, levels = tool_order))
  
  # Create plot EXACTLY LIKE FIGURE 3 STYLE
  p <- ggplot(plot_data, 
              aes(x = peak_caller,  # X-axis: Different tools
                  y = gene_name,    # Y-axis: Top 10 genes for this celltype
                  size = log_weighted_count, 
                  color = log_tpm)) +
    geom_point(alpha = 0.8, shape = 16) +  # Use solid circles
    scale_size_continuous(
      name = expression(log[2]("Weighted Peak Count + 1")),
      range = c(1, 4),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_color_viridis(
      name = expression(log[10]("TPM + 1")),
      option = "viridis",
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_x_discrete(drop = FALSE) +  # Show all tools
    labs(
      title = paste("Top Genes -", celltype),  # SAME TITLE FORMAT AS FIGURE 3
      x = "Peak Caller",  # X-axis label
      y = "Gene"          # Y-axis label
    ) +
    theme_classic(base_size = 8) +  # SAME THEME AS FIGURE 3
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, face = "plain"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.box = "vertical",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "mm")
    )
  
  # Store plot in list
  plot_list_celltype[[celltype]] <- p
  
  cat("✓ Created\n")
}

#=======================================================================================
# SAVE INDIVIDUAL CELLTYPE PLOTS
#=======================================================================================

cat("\nSaving individual celltype plots...\n")

for (celltype in names(plot_list_celltype)) {
  # Save individual plot for each cell type
  ggsave(
    file.path(celltype_wise_dir, 
              paste0("Figure4h_",  hist, "_Celltype_", gsub("[^[:alnum:]]", "_", celltype), 
                     "_Top10Genes_with_input.pdf")),
    plot = plot_list_celltype[[celltype]],
    width = 8.5,
    height = 10,
    units = "cm",
    dpi = 600
  )
  
  # Also save as PNG
  ggsave(
    file.path(celltype_wise_dir, 
              paste0("Figure4h_",  hist, "_Celltype_", gsub("[^[:alnum:]]", "_", celltype), 
                     "_Top10Genes_with_input.png")),
    plot = plot_list_celltype[[celltype]],
    width = 8.5,
    height = 10,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
}

cat("✓ Saved individual plots\n")


#=======================================================================================
# Saving the top 10 genes information 
#=======================================================================================

#=======================================================================================
# SAVE DATA FROM FIGURE 3 ANALYSIS
#=======================================================================================
cat("\n=== SAVING DATA FROM FIGURE 3 ANALYSIS ===\n")

# Create directory for saving data
data_save_dir <- file.path(out_dir, paste0("Figure_data_", hist))
dir.create(data_save_dir, recursive = TRUE, showWarnings = FALSE)

# Store all data in a list
Figure_data <- list()

for (celltype in all_celltypes) {
  cat("  Saving data for", celltype, "... ")
  
  # Filter data for this cell type (same as your Figure 3 code)
  celltype_data <- final_merged_data %>%
    filter(celltype == !!celltype)
  
  if (nrow(celltype_data) == 0) {
    cat("No data found\n")
    next
  }
  
  # Check for and handle duplicates by taking the mean (same as your code)
  celltype_data_clean <- celltype_data %>%
    group_by(gene_name, peak_caller) %>%
    summarise(
      log_tpm = mean(log_tpm, na.rm = TRUE),
      log_weighted_count = mean(log_weighted_count, na.rm = TRUE),
      tpm_expression = mean(tpm_expression, na.rm = TRUE),
      weighted_peak_count = mean(weighted_peak_count, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Get top 10 genes for this cell type (based on TPM - same as your code)
  top_genes_celltype <- celltype_data_clean %>%
    group_by(gene_name) %>%
    summarise(
      av_log_tpm = mean(log_tpm, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(av_log_tpm)) %>%
    slice_head(n = 10)
  
  if (nrow(top_genes_celltype) == 0) {
    cat("No valid TPM data\n")
    next
  }
  
  # Get data for top 10 genes (same as your plot_data)
  plot_data <- celltype_data_clean %>%
    filter(gene_name %in% top_genes_celltype$gene_name) %>%
    complete(gene_name, peak_caller, fill = list(
      weighted_peak_count = 0,
      log_weighted_count = 0,
      tpm_expression = 0,
      log_tpm = 0
    )) %>%
    mutate(gene_name = factor(gene_name, levels = rev(top_genes_celltype$gene_name))) %>%
    mutate(peak_caller = factor(peak_caller, levels = tool_order))
  
  # Store in list
  Figure_data[[celltype]] <- list(
    top_genes = top_genes_celltype,
    plot_data = plot_data,
    celltype = celltype,
    n_genes = nrow(top_genes_celltype),
    gene_names = top_genes_celltype$gene_name
  )
  
  # Also save individual celltype data as CSV
  csv_file <- file.path(data_save_dir, paste0("Figure_data_", hist, "_", gsub("[^[:alnum:]]", "_", celltype), "_with_input.csv"))
  write.csv(plot_data, csv_file, row.names = FALSE)
  
  cat("✓ Saved\n")
}

# Save combined data
all_top_genes <- unique(unlist(lapply(Figure_data, function(x) x$gene_names)))

combined_data <- final_merged_data %>%
  filter(gene_name %in% all_top_genes) %>%
  select(gene_name, peak_caller, celltype,
         weighted_peak_count, tpm_expression, log_tpm, log_weighted_count)

combined_csv <- file.path(data_save_dir, paste0("Figure_combined_data_", hist, "_with_input.csv"))
write.csv(combined_data, combined_csv, row.names = FALSE)


####################################################################
# CREATE COMBINED A4 PAGE (ALL CELLTYPES TOGETHER)
####################################################################

cat("\nCreating combined A4 page with all celltypes...\n")

# Arrange all celltype plots on A4 page (landscape orientation)
n_celltypes <- length(plot_list_celltype)

# Determine optimal grid layout (SAME LOGIC AS FIGURE 3)
if (n_celltypes <= 4) {
  n_col <- 2
  n_row <- ceiling(n_celltypes / n_col)
} else if (n_celltypes <= 6) {
  n_col <- 3
  n_row <- ceiling(n_celltypes / n_col)
} else {
  n_col <- ceiling(sqrt(n_celltypes))
  n_row <- ceiling(n_celltypes / n_col)
}

# Combine plots using patchwork (preserves legends)
combined_celltype_plot <- wrap_plots(plot_list_celltype, 
                                     ncol = n_col, 
                                     nrow = n_row) +
  plot_annotation(
    title = paste("Top 10 Genes for Each Cell Type -", hist),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray50")
    )
  )

# Save combined plot on A4 page
ggsave(
  file.path(celltype_wise_dir, 
            paste0("Figure5i_",  hist, "_All_Celltypes_TopGenes_Combined_with_input.pdf")),
  plot = combined_celltype_plot,
  width = 29.7,
  height = 21.0,
  units = "cm",
  dpi = 600
)

# Also save as PNG
ggsave(
  file.path(celltype_wise_dir, 
            paste0("Figure5i_",  hist, "_All_Celltypes_TopGenes_Combined_with_input.png")),
  plot = combined_celltype_plot,
  width = 29.7,
  height = 21.0,
  units = "cm",
  dpi = 300,
  bg = "white"
)

cat("✓ Saved combined A4 page\n")
}

###########################################################################################################
## Figure 5 K-N
###########################################################################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Define the base directory
base_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input"

# List of files with full paths
files <- c(
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/Figure_data_H3K27ac-b/Figure_combined_data_H3K27ac-b_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/Figure_data_H3K27ac-s/Figure_combined_data_H3K27ac-s_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/Figure_data_H3K27me3/Figure_combined_data_H3K27me3_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/Figure_data_H3K4me1/Figure_combined_data_H3K4me1_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/Figure_data_H3K4me2/Figure_combined_data_H3K4me2_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/Figure_data_H3K4me3/Figure_combined_data_H3K4me3_with_input.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter/result_GT_scrnaseq_with_input/Figure_data_H3K9me3/Figure_combined_data_H3K9me3_with_input.csv"
)

# Create an empty list to store dataframes
all_data <- list()

cat("Starting to merge Figure 3 combined data files...\n")
cat("==========================================\n")

# Read and process each file
for (file in files) {
  if (file.exists(file)) {
    # Extract histone mark from file path
    # Pattern: Figure_data_(histone-mark)/figure3_combined_data_(histone-mark)_with_input.csv
    histone_mark <- str_extract(file, "(?<=Figure_data_)[^/]+")
    
    # Alternative: extract from filename
    if (is.na(histone_mark)) {
      histone_mark <- str_extract(file, "(?<=Figure_combined_data_)[^_]+")
    }
    
    cat("Processing:", basename(dirname(file)), "\n")
    cat("Histone mark extracted:", histone_mark, "\n")
    
    # Read the CSV file
    data <- read_csv(file, col_types = cols())
    
    # Check the structure
    cat("  Columns found:", paste(colnames(data), collapse = ", "), "\n")
    cat("  Number of rows:", nrow(data), "\n")
    
    # Add histone mark column as FIRST column
    data <- data %>%
      mutate(histone = histone_mark) %>%
      select(histone, everything())  # Move histone to first position
    
    # Check if all required columns exist
    required_cols <- c("gene_name", "peak_caller", "celltype", 
                      "weighted_peak_count", "tpm_expression", 
                      "log_tpm", "log_weighted_count")
    
    missing_cols <- setdiff(required_cols, colnames(data))
    
    if (length(missing_cols) > 0) {
      cat("  WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
      # Add missing columns with NA values
      for (col in missing_cols) {
        if (!col %in% colnames(data)) {
          data[[col]] <- NA
        }
      }
    }
    
    # Reorder columns to match desired format
    data <- data %>%
      select(histone, gene_name, peak_caller, celltype,
             weighted_peak_count, tpm_expression, log_tpm, log_weighted_count,
             everything())
    
    # Add to list
    all_data[[file]] <- data
    
    # Show first few rows
    cat("  First few rows:\n")
    print(head(data, 3))
    cat("\n")
    
  } else {
    cat("ERROR: File not found:", file, "\n")
  }
}

cat("==========================================\n")
cat("All files processed. Now combining...\n")

# Check if we have any data
if (length(all_data) == 0) {
  stop("No files were successfully read. Please check file paths.")
}

# Combine all dataframes
combined_all_histones_data <- bind_rows(all_data)

# View the structure of combined data
cat("\nCombined data structure:\n")
cat("Number of rows:", nrow(combined_all_histones_data), "\n")
cat("Number of columns:", ncol(combined_all_histones_data), "\n")
cat("Column names:", paste(colnames(combined_all_histones_data), collapse = ", "), "\n")

# Count unique values
cat("\nUnique values:\n")
cat("Histones:", paste(unique(combined_all_histones_data$histone), collapse = ", "), "\n")
cat("Genes:", length(unique(combined_all_histones_data$gene_name)), "\n")
cat("Cell types:", paste(unique(combined_all_histones_data$celltype), collapse = ", "), "\n")
cat("Peak callers:", paste(unique(combined_all_histones_data$peak_caller), collapse = ", "), "\n")

# View first few rows of combined data
cat("\nFirst 10 rows of combined data:\n")
print(head(combined_all_histones_data, 10))

# Save the combined data
output_dir <- file.path(base_dir, "ALL_HISTONES_COMBINED")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(output_dir, "ALL_histones_figure3_combined_data_with_input.csv")
write_csv(combined_all_histones_data, output_file)
cat("\n✓ Combined data saved to:", output_file, "\n")


combined_all_histones_data[is.na(combined_all_histones_data)]<-0
combined_all_histones_data[is.na(combined_all_histones_data)]<-0

top10_summary <- combined_all_histones_data %>%
  group_by(peak_caller, histone, celltype) %>%
  summarize(
    median_log_weighted_count = median(log_weighted_count, na.rm = TRUE),
    mean_log_weighted_count = mean(log_weighted_count, na.rm = TRUE),
    sd_log_weighted_count = sd(log_weighted_count, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(peak_caller),
    Histone = factor(histone)
  )



# Install if needed
# install.packages("data.table")
# install.packages("dplyr")

library(data.table)
library(dplyr)

# Read the CSV file
df <- fread('/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_frip/HumanPBMC_combined_peak_frip_summary_with_input.csv')

# View the data structure
head(df)


# Make sure column names match - adjust as needed
top10_summary_with_peaks <- top10_summary %>%
  left_join(df %>% 
              select(Histone, Method, Celltype, Number_of_Peaks),
            by = c("Histone", "Method", "celltype" = "Celltype"))

# Check the result
head(top10_summary_with_peaks)

# Add log_peaks and norm_ratio columns
top10_summary_with_peaks <- top10_summary_with_peaks %>%
  mutate(
    # Calculate log of Number_of_Peaks (natural log)
    log_peaks = log10(Number_of_Peaks),
    
    # Calculate norm_ratio = median_log_weighted_count / log_peaks
    norm_ratio = median_log_weighted_count / log_peaks,
  )

# Check the result
head(top10_summary_with_peaks)
# Check if there are any NaN values in norm_ratio
nan_count <- sum(is.nan(top10_summary_with_peaks$norm_ratio))
print(paste("Number of NaN values in norm_ratio:", nan_count))

# Simple solution: Replace NaN with 0
top10_summary_with_peaks <- top10_summary_with_peaks %>%
  mutate(
    norm_ratio = ifelse(is.nan(norm_ratio), 0, norm_ratio)
  )

# Check results
print(paste("NaN values after replacement:", sum(is.nan(top10_summary_with_peaks$norm_ratio))))

top10_summary2 <- top10_summary_with_peaks %>%
  group_by(peak_caller, histone) %>%
  summarize(
    median_norm_ratio = median(norm_ratio, na.rm = TRUE),
    mean_norm_ratio = mean(norm_ratio, na.rm = TRUE),
    sd_norm_ratio = sd(norm_ratio, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(peak_caller),
    Histone = factor(histone)
  )

  all_methods <- unique(top10_summary2$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#631879", "#9C9EDE", "#637939", "#8C6D31", "#BD9E39"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)


library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)

# First clean the data
top10_summary_clean <- top10_summary2 %>%
  filter(!is.na(median_norm_ratio) & !is.na(sd_norm_ratio))

# Create scatter plot WITHOUT percentage formatting
top10_scatter2 <- ggplot(top10_summary_clean, 
                       aes(x = Histone, 
                           y = median_norm_ratio,
                           color = Method)) +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = median_norm_ratio - sd_norm_ratio,
                    ymax = median_norm_ratio + sd_norm_ratio),
                width = 0.3, 
                position = position_dodge(width = 0.5),
                alpha = 0.7) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  # WITHOUT percent_format - use default numeric formatting
  scale_y_continuous(
    limits = c(-.05, .25),
    breaks = seq(-.05, .25, by=.05)
  ) +
  labs(
    title = "HumanPBMC with input",
    x = "Histone Modification",
    y = "Norm_ratio"
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

print(top10_scatter2)

ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_GT_scrnaseq_corrected_final_top10gene_scatter_norm/HumanPBMC_top10_Scatter_with_ErrorBars_With_Input.pdf",
       plot = top10_scatter2,
       width = 8,
       height = 6,
       device = "pdf",
       bg = "white")

