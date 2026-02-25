#!/bin/bash
#############################################################################################
### scCTpeak: A Complete Benchmark Suite for scCUT&Tag Peak-Calling
### Unified Function for All Tools, Datasets, and Configurations
### WITH DATA PROCESSING MODULE & UNIFIED OUTPUT FORMAT
### Version: 3.1.0
#############################################################################################

# ===========================================================================================
# GLOBAL CONFIGURATION
# ===========================================================================================
SCCTPEAK_VERSION="3.1.0"
SCCTPEAK_AUTHOR="Wahiduzzaman, M., et al."

# Default paths (UPDATE THESE TO MATCH YOUR SYSTEM)
SCCTPEAK_HOME="/home/wahid/project_scHMTF"
HUMAN_DIR="${SCCTPEAK_HOME}/GSE195725_processed_data"
MOUSE_DIR="${SCCTPEAK_HOME}/GSE157637_processed_data"

# Tool paths (UPDATE THESE TO MATCH YOUR INSTALLATIONS)
SEACR_SCRIPT="/home/wahid/tools/SEACR/SEACR_1.3.sh"
DROMPA_DIR="/home/wahid/DROMPAplus/bin"
SAMTOOLS="/bin/samtools"
RSCRIPT_CMD="/home/wahid/anaconda3/envs/bioconda/bin/Rscript"

# R script paths
SCCTPEAK_R_SCRIPT="/home/wahid/project_scHMTF/scripts/scCTpeak_data_processing.R"

# Threads for parallel processing
THREADS=32

# Split BAM directories (where combined BAM and barcode files are stored)
SPLIT_BAM_HUMAN_DIR="${HUMAN_DIR}/splitbam_realbam"
SPLIT_BAM_MOUSE_DIR="${MOUSE_DIR}/splitbam_realbam"

# ===========================================================================================
# UNIFIED OUTPUT FUNCTIONS
# ===========================================================================================

# Generate unified filename
get_unified_filename() {
    local tool="$1"
    local histone="$2"
    local cell="$3"
    local mode="$4"
    local control="$5"
    
    # Tool name mapping to uppercase
    local tool_name
    case "$tool" in
        "drompa") tool_name="DROMPAplus" ;;
        "macs2") tool_name="MACS2" ;;
        "genrich") tool_name="Genrich" ;;
        "gopeaks") tool_name="GoPeaks" ;;
        "homer") tool_name="HOMER" ;;
        "seacr") tool_name="SEACR" ;;
        "sicer2") tool_name="SICER2" ;;
        *) tool_name="$tool" ;;
    esac
    
    # Mark name with mode suffix for H3K27ac
    local mark_name="$histone"
    if [[ "$histone" == "H3K27ac" ]]; then
        if [[ "$mode" == "broad" ]]; then
            mark_name="${histone}-b"
        elif [[ "$mode" == "narrow" ]]; then
            mark_name="${histone}-s"
        fi
    fi
    
    # Control suffix (optional)
    local control_suffix=""
    [[ -n "$control" ]] && control_suffix="_${control}"
    
    echo "${tool_name}_${mark_name}_${cell}${control_suffix}.bed"
}

# Sort and save unified BED file
create_unified_bed() {
    local input_file="$1"
    local tool="$2"
    local histone="$3"
    local cell="$4"
    local mode="$5"
    local control="$6"
    local output_dir="$7"
    
    if [[ ! -f "$input_file" ]]; then
        echo "  ⚠️ Input file not found: $input_file"
        return 1
    fi
    
    # Get unified filename
    local unified_name=$(get_unified_filename "$tool" "$histone" "$cell" "$mode" "$control")
    local output_file="${output_dir}/${unified_name}"
    
    # Create output directory if needed
    mkdir -p "$output_dir"
    
    echo "  📦 Creating unified BED: $unified_name"
    
    # Process different input formats
    local temp_file="${output_file}.tmp"
    
    # Check file type and convert to BED6
    if [[ "$input_file" == *.narrowPeak || "$input_file" == *.broadPeak ]]; then
        # narrowPeak/broadPeak format (10 columns)
        awk 'BEGIN{OFS="\t"} {
            if ($1 ~ /^chr/) {
                chrom = $1
                start = $2
                end = $3
                name = "peak_"NR
                if ($4 != ".") name = $4
                score = $5
                strand = $6
                signal = $7
                pval = $8
                qval = $9
                summit = $10
                
                # Use score as peak intensity
                print chrom, start, end, name, score, strand
            }
        }' "$input_file" > "$temp_file"
        
    elif [[ "$input_file" == *.bed && $(awk 'NR==1{print NF}' "$input_file") -ge 6 ]]; then
        # Already BED6+ format
        awk 'BEGIN{OFS="\t"} {
            if ($1 ~ /^chr/) {
                chrom = $1
                start = $2
                end = $3
                name = ($4 != "." && $4 != "" && $4 != " ") ? $4 : "peak_"NR
                score = ($5 != "." && $5 != "") ? $5 : 0
                strand = ($6 != "." && $6 != "") ? $6 : "."
                
                # Ensure valid values
                if (score == "") score = 0
                if (score == ".") score = 0
                if (score < 0) score = 0
                if (score > 1000) score = 1000
                if (strand == "") strand = "."
                
                print chrom, start, end, name, score, strand
            }
        }' "$input_file" > "$temp_file"
        
    else
        # Generic format - try to parse
        awk 'BEGIN{OFS="\t"} {
            if ($1 ~ /^chr/) {
                chrom = $1
                start = $2
                end = $3
                name = "peak_"NR
                score = 0
                strand = "."
                
                if (NF >= 4 && $4 != ".") name = $4
                if (NF >= 5 && $5 != ".") score = $5 + 0
                if (NF >= 6 && $6 != ".") strand = $6
                
                print chrom, start, end, name, score, strand
            }
        }' "$input_file" > "$temp_file"
    fi
    
    # Sort and ensure unique peaks
    sort -k1,1 -k2,2n "$temp_file" | awk 'BEGIN{OFS="\t"; prev_chr=""; prev_start=0; prev_end=0} {
        # Remove duplicate peaks (same chr:start-end)
        if (!($1 == prev_chr && $2 == prev_start && $3 == prev_end)) {
            print $0
            prev_chr = $1
            prev_start = $2
            prev_end = $3
        }
    }' > "$output_file"
    
    # Clean up
    rm -f "$temp_file"
    
    local peak_count=$(wc -l < "$output_file")
    echo "  ✅ Created: $unified_name ($peak_count peaks)"
    
    return 0
}

# ===========================================================================================
# DATA PROCESSING MODULE
# ===========================================================================================

print_header() {
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
    echo "║ $1"
    echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
}

print_step() {
    echo "↳ $1"
}

# -------------------------------------------------------------------------------------------
# Create R script for data processing (FIXED VERSION)
# -------------------------------------------------------------------------------------------
create_data_processing_r_script() {
    local script_path="$1"
    
    cat > "$script_path" << 'EOF'
#############################################################################################################################
#### scCTpeak DATA PROCESSING MODULE
#### Making 50 BP Sorted BAM Files
#### Human PBMC (HMs) and Mouse Brain (HMs + TFs)
#### CellType-wise
#### Author: Md Wahiduzzaman
#############################################################################################################################

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(data.table)
  library(parallel)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--dataset_path"), type="character", help="Path to Seurat RDS file"),
  make_option(c("--frag_path"), type="character", help="Path to fragments TSV.GZ file"),
  make_option(c("--assay_name"), type="character", default="peaks", help="Assay name [default: peaks]"),
  make_option(c("--genome_file"), type="character", help="Path to genome file"),
  make_option(c("--outdir"), type="character", help="Output directory for BAM files"),
  make_option(c("--cell_id_col"), type="character", default="celltype", help="Cell ID column name"),
  make_option(c("--threads"), type="integer", default=32, help="Number of threads [default: 32]"),
  make_option(c("--histone"), type="character", help="Histone/TF name (for naming)"),
  make_option(c("--output_base"), type="character", help="Base output directory for scCTpeak")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
required_args <- c("dataset_path", "frag_path", "genome_file", "outdir", "histone")
missing_args <- required_args[!required_args %in% names(opt)]
if(length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
}

# Set global parameters
SAMTOOLS <- "/bin/samtools"
THREADS <- opt$threads

if(!file.exists(SAMTOOLS)) {
  stop("❌ samtools not found at: ", SAMTOOLS)
}

# Core Functions
convert_to_bam_sorted_named <- function(samfile){
  bamfile <- sub("\\.sam$", ".bam", samfile)
  cmd1 <- paste(SAMTOOLS, "view -@", THREADS, "-bS", shQuote(samfile), "|",
                SAMTOOLS, "sort -@", THREADS, "-o", shQuote(bamfile))
  system(cmd1)
  
  cmd2 <- paste(SAMTOOLS, "index", shQuote(bamfile))
  system(cmd2)
  
  file.remove(samfile)
  return(bamfile)
}

write_pe_sam <- function(df, outfile, header){
  sam_lines <- lapply(seq_len(nrow(df)), function(i){
    row <- df[i]
    read_id <- paste0(row$cell, "_", i)
    chr <- row$chr
    start1 <- format(as.integer(row$start)+1, scientific=FALSE)
    start2 <- format(as.integer(row$end)-50, scientific=FALSE)
    if(as.integer(start2)<1) start2 <- 1
    flag1 <- 99
    flag2 <- 147
    r1 <- paste(read_id, flag1, chr, start1, 255, "50M", "=", start2,
                as.integer(start2)-as.integer(start1), "*", "*", sep="\t")
    r2 <- paste(read_id, flag2, chr, start2, 255, "50M", "=", start1,
                as.integer(start1)-as.integer(start2), "*", "*", sep="\t")
    c(r1, r2)
  })
  writeLines(c(header, unlist(sam_lines)), con=outfile)
}


process_dataset <- function(dataset_path, frag_path, assay_name, genome_file, outdir, cell_id_col, histone_name){
  
  cat("🔧 Starting scCTpeak Data Processing\n")
  cat("   Dataset:", basename(dataset_path), "\n")
  cat("   Histone/TF:", histone_name, "\n")
  cat("   Output:", outdir, "\n")
  
  # Load object
  obj <- readRDS(dataset_path)
  obj <- UpdateSeuratObject(obj)
  
  # Fragment object
  fragments <- CreateFragmentObject(path=frag_path, cells=colnames(obj))
  
  # Genome header
  chrom_sizes <- read.table(genome_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(chrom_sizes) <- c("seqnames","seqlengths")
  seqinfo_obj <- Seqinfo(seqnames=chrom_sizes$seqnames,
                         seqlengths=chrom_sizes$seqlengths,
                         genome=ifelse(grepl("mm10", genome_file), "mm10", "hg38"))
  
  chrom_assay <- CreateChromatinAssay(counts=GetAssayData(obj[[assay_name]], slot="counts"),
                                     fragments=fragments,
                                     sep=c("-", "-"),
                                     genome=seqinfo_obj)
  obj[[assay_name]] <- chrom_assay
  DefaultAssay(obj) <- assay_name
  
  # Read fragments as data.table
  frags <- fread(cmd=paste("zcat", Fragments(obj)[[1]]@path),
                 header=FALSE,
                 colClasses=c("character","integer","integer","character","integer"))
  colnames(frags) <- c("chr","start","end","cell","score")
  
  # Output dir - create scCTpeak compatible structure
  bam_dir <- file.path(outdir, "celltypeswisebam_50bp", histone_name)
  dir.create(bam_dir, recursive=TRUE, showWarnings=FALSE)
  
  # Header
  header <- c("@HD\tVN:1.6\tSO:unsorted",
              paste0("@SQ\tSN:", chrom_sizes$seqnames, "\tLN:", chrom_sizes$seqlengths))
    
  # Assign clusters
  Idents(obj) <- obj[[cell_id_col]][,1]
  clusters <- unique(Idents(obj))
  
  # Clean cluster names: remove spaces AND underscores
  clusters_clean <- gsub("[ _]", "", clusters)  # FIXED: removes both space and underscore
  
  cat("📊 Found", length(clusters), "cell types:", paste(clusters, collapse=", "), "\n")
  cat("📊 Cleaned names:", paste(clusters_clean, collapse=", "), "\n")
  
  # Remove old files
  unlink(list.files(bam_dir, pattern="\\.(sam|bam|bai)$", full.names=TRUE))
  
  # Parallel processing
  results <- mclapply(seq_along(clusters), function(i){
    cluster <- clusters[i]
    cluster_clean <- clusters_clean[i]
    
    cat("  Processing cluster:", cluster, "->", cluster_clean, "\n")
    
    treatment_cells <- colnames(obj)[Idents(obj)==cluster]
    control_cells   <- colnames(obj)[Idents(obj)!=cluster]
    
    frag_treat <- frags[cell %in% treatment_cells]
    frag_ctrl  <- frags[cell %in% control_cells]
    
    # Create scCTpeak compatible filenames with CLEANED names
    treatment_bam <- file.path(bam_dir, paste0(histone_name, "_", cluster_clean, ".bam"))
    control_bam <- file.path(bam_dir, paste0("input_", cluster_clean, ".bam"))
    
    # Temporary SAM files
    sam_treat <- sub("\\.bam$", ".sam", treatment_bam)
    sam_ctrl <- sub("\\.bam$", ".sam", control_bam)
    
    # Write SAM files
    if(nrow(frag_treat) > 0) {
      write_pe_sam(frag_treat, sam_treat, header)
      convert_to_bam_sorted_named(sam_treat)
      cat("    ✅ Created treatment BAM:", basename(treatment_bam), "\n")
    } else {
      cat("    ⚠️ No fragments for treatment cluster", cluster_clean, "\n")
    }
    
    if(nrow(frag_ctrl) > 0) {
      write_pe_sam(frag_ctrl, sam_ctrl, header)
      convert_to_bam_sorted_named(sam_ctrl)
      cat("    ✅ Created control BAM:", basename(control_bam), "\n")
    } else {
      cat("    ⚠️ No fragments for control cluster", cluster_clean, "\n")
    }
    
    # Return file info
    list(
      cluster = cluster,
      cluster_clean = cluster_clean,
      treatment_bam = if(file.exists(treatment_bam)) treatment_bam else NULL,
      control_bam = if(file.exists(control_bam)) control_bam else NULL,
      treatment_fragments = nrow(frag_treat),
      control_fragments = nrow(frag_ctrl)
    )
  }, mc.cores = min(opt$threads, 64))

  # Summary
  cat("\n✅ Processing Complete!\n")
  cat("   Output directory:", bam_dir, "\n\n")
  
  # Generate summary
  cat("📋 Generated BAM Files:\n")
  for(result in results) {
    if(!is.null(result$treatment_bam)) {
      cat("  - Treatment:", basename(result$treatment_bam), 
          sprintf("(%d fragments)", result$treatment_fragments), "\n")
    }
    if(!is.null(result$control_bam)) {
      cat("  - Control:", basename(result$control_bam), 
          sprintf("(%d fragments)", result$control_fragments), "\n")
    }
  }
  
  # Create symlinks for scCTpeak compatibility
  if(!is.null(opt$output_base)) {
    scctpeak_bam_dir <- file.path(opt$output_base, "simulation_with_scCTpeak", histone_name)
    dir.create(dirname(scctpeak_bam_dir), recursive=TRUE, showWarnings=FALSE)
    
    cat("🔗 Creating symlinks for scCTpeak:", scctpeak_bam_dir, "\n")
    
    # Create symlinks for all BAM files
    bam_files <- list.files(bam_dir, pattern="\\.bam$", full.names=TRUE)
    for(bam_file in bam_files) {
      link_name <- file.path(scctpeak_bam_dir, basename(bam_file))
      if(!file.exists(link_name)) {
        file.symlink(bam_file, link_name)
      }
    }
    
    # Create BAI symlinks
    bai_files <- list.files(bam_dir, pattern="\\.bai$", full.names=TRUE)
    for(bai_file in bai_files) {
      link_name <- file.path(scctpeak_bam_dir, basename(bai_file))
      if(!file.exists(link_name)) {
        file.symlink(bai_file, link_name)
      }
    }
  }
  
  return(list(bam_dir = bam_dir, results = results))
}

# Main execution
main <- function() {
  cat("================================================================\n")
  cat("                    scCTpeak DATA PROCESSING                    \n")
  cat("================================================================\n")
  
  # Process dataset
  result <- process_dataset(
    dataset_path = opt$dataset_path,
    frag_path = opt$frag_path,
    assay_name = opt$assay_name,
    genome_file = opt$genome_file,
    outdir = opt$outdir,
    cell_id_col = opt$cell_id_col,
    histone_name = opt$histone
  )
  
  # Save processing summary
  summary_file <- file.path(result$bam_dir, "processing_summary.txt")
  sink(summary_file)
  cat("scCTpeak Data Processing Summary\n")
  cat("================================\n")
  cat("Timestamp:", as.character(Sys.time()), "\n")
  cat("Dataset:", basename(opt$dataset_path), "\n")
  cat("Histone/TF:", opt$histone, "\n")
  cat("Output directory:", result$bam_dir, "\n")
  cat("\nGenerated BAM files:\n")
  
  total_bam <- 0
  for(res in result$results) {
    if(!is.null(res$treatment_bam)) {
      cat("  ", basename(res$treatment_bam), 
          sprintf("(%d fragments)", res$treatment_fragments), "\n")
      total_bam <- total_bam + 1
    }
    if(!is.null(res$control_bam)) {
      cat("  ", basename(res$control_bam), 
          sprintf("(%d fragments)", res$control_fragments), "\n")
      total_bam <- total_bam + 1
    }
  }
  
  cat("\nTotal BAM files:", total_bam, "\n")
  cat("Total cell types:", length(result$results), "\n")
  sink()
  
  cat("\n📄 Summary saved to:", summary_file, "\n")
  cat("🎉 scCTpeak data processing completed successfully!\n")
  
  # Return success
  return(0)
}

# Run main function
if(!interactive()) {
  quit(status = main())
}
EOF
    
    chmod +x "$script_path"
    echo "✅ Created R script: $script_path"
}
# -------------------------------------------------------------------------------------------
# Run Data Processing
# -------------------------------------------------------------------------------------------
run_data_processing() {
    local dataset="$1"
    local histone="$2"
    
    print_header "scCTpeak Data Processing | Dataset: $dataset | Histone/TF: $histone"
    
    # Create R script if it doesn't exist
    if [[ ! -f "$SCCTPEAK_R_SCRIPT" ]]; then
        echo "📝 Creating R data processing script..."
        mkdir -p "$(dirname "$SCCTPEAK_R_SCRIPT")"
        create_data_processing_r_script "$SCCTPEAK_R_SCRIPT"
    fi
    
    # Set parameters based on dataset
    local dataset_path frag_path assay_name genome_file cell_id_col outdir_base
    
    if [[ "$dataset" == "human" ]]; then
        outdir_base="${HUMAN_DIR}"
        genome_file="${HUMAN_DIR}/ref/hg38.genome"
        cell_id_col="predicted.celltype.l1"
        
        case "$histone" in
            H3K27ac|H3K27me3|H3K4me1|H3K4me2|H3K4me3|H3K9me3)
                dataset_path="${HUMAN_DIR}/${histone}.rds"
                frag_path="${HUMAN_DIR}/${histone}_fragments.tsv.gz"
                assay_name="tiles"
                ;;
            *)
                echo "❌ Unknown histone for human: $histone"
                echo "   Valid options: H3K27ac, H3K27me3, H3K4me1, H3K4me2, H3K4me3, H3K9me3"
                return 1
                ;;
        esac
    elif [[ "$dataset" == "mouse" ]]; then
        outdir_base="${MOUSE_DIR}"
        genome_file="${MOUSE_DIR}/ref/mm10.genome"
        cell_id_col="cell_type"
        
        case "$histone" in
            H3K27ac|H3K27me3|H3K36me3|H3K4me3)
                dataset_path="${MOUSE_DIR}/${histone}_seurat_object.Rds"
                frag_path="${MOUSE_DIR}/${histone}_fragments.tsv.gz"
                assay_name="peaks"
                ;;
            Olig2|Rad21)
                dataset_path="${MOUSE_DIR}/${histone}_seurat_object.Rds"
                frag_path="${MOUSE_DIR}/${histone}_fragments.tsv.gz"
                assay_name="peaks"
                ;;
            *)
                echo "❌ Unknown histone/TF for mouse: $histone"
                echo "   Valid options: H3K27ac, H3K27me3, H3K36me3, H3K4me3, Olig2, Rad21"
                return 1
                ;;
        esac
    else
        echo "❌ Unknown dataset: $dataset"
        return 1
    fi
    
    # Check if input files exist
    if [[ ! -f "$dataset_path" ]]; then
        echo "❌ Dataset file not found: $dataset_path"
        return 1
    fi
    
    if [[ ! -f "$frag_path" ]]; then
        echo "❌ Fragments file not found: $frag_path"
        return 1
    fi
    
    if [[ ! -f "$genome_file" ]]; then
        echo "❌ Genome file not found: $genome_file"
        return 1
    fi
    
    # Create output directory
    local outdir="${outdir_base}/simulation_with_scCTpeak"
    mkdir -p "$outdir"
    
    print_step "Dataset: $(basename "$dataset_path")"
    print_step "Fragments: $(basename "$frag_path")"
    print_step "Genome: $(basename "$genome_file")"
    print_step "Output: $outdir"
    
    # Run R script
    local r_cmd="$RSCRIPT_CMD '$SCCTPEAK_R_SCRIPT' \
        --dataset_path '$dataset_path' \
        --frag_path '$frag_path' \
        --assay_name '$assay_name' \
        --genome_file '$genome_file' \
        --outdir '$outdir' \
        --cell_id_col '$cell_id_col' \
        --threads $THREADS \
        --histone '$histone' \
        --output_base '$outdir_base'"
    
    echo "Running R script..."
    echo ""
    
    if eval "$r_cmd"; then
        echo "✅ Data processing completed successfully!"
        echo ""
        
        # Show generated files
        local bam_dir="${outdir}/simulation_with_scCTpeak/${histone}"
        if [[ -d "$bam_dir" ]]; then
            echo "📁 Generated BAM files in: $bam_dir"
            local bam_count=$(find "$bam_dir" -name "*.bam" | wc -l)
            echo "   Total BAM files: $bam_count"
            
            # List files
            if [[ $bam_count -gt 0 ]]; then
                print_step "Generated files:"
                find "$bam_dir" -name "*.bam" -exec basename {} \; | sort | head -20 | while read -r file; do
                    echo "  - $file"
                done
                
                if [[ $bam_count -gt 20 ]]; then
                    echo "  ... and $((bam_count - 20)) more"
                fi
            fi
        fi
        
        # Also show scCTpeak compatible directory
        local scctpeak_bam_dir="${outdir_base}/simulation_with_scCTpeak/${histone}"
        if [[ -d "$scctpeak_bam_dir" ]]; then
            echo ""
            echo "🔗 scCTPEAK compatible directory: $scctpeak_bam_dir"
            echo "   (Symlinks created for peak calling tools)"
        fi
    else
        echo "❌ Data processing failed"
        return 1
    fi
    
    echo ""
    return 0
}

# -------------------------------------------------------------------------------------------
# Batch Data Processing
# -------------------------------------------------------------------------------------------
run_batch_data_processing() {
    local dataset="$1"
    
    print_header "Batch Data Processing | Dataset: $dataset"
    
    if [[ "$dataset" == "human" ]]; then
        local histones=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
        
        for histone in "${histones[@]}"; do
            print_step "Processing $histone..."
            if ! run_data_processing "$dataset" "$histone"; then
                echo "⚠️ Failed to process $histone, continuing..."
            fi
            
            # Small delay between runs
            sleep 3
        done
        
    elif [[ "$dataset" == "mouse" ]]; then
        # Process histones
        local histones=("H3K27ac" "H3K27me3" "H3K36me3" "H3K4me3")
        
        for histone in "${histones[@]}"; do
            print_step "Processing histone $histone..."
            if ! run_data_processing "$dataset" "$histone"; then
                echo "⚠️ Failed to process $histone, continuing..."
            fi
            
            # Small delay between runs
            sleep 3
        done
        
        # Process TFs
        local tfs=("Olig2" "Rad21")
        
        for tf in "${tfs[@]}"; do
            print_step "Processing TF $tf..."
            if ! run_data_processing "$dataset" "$tf"; then
                echo "⚠️ Failed to process $tf, continuing..."
            fi
            
            # Small delay between runs
            sleep 3
        done
        
    else
        echo "❌ Unknown dataset: $dataset"
        return 1
    fi
    
    echo "✅ Batch data processing completed for $dataset"
    echo ""
    return 0
}

# -------------------------------------------------------------------------------------------
# Verify BAM Files
# -------------------------------------------------------------------------------------------
verify_processed_data() {
    local dataset="$1"
    local histone="$2"
    local cell="$3"
    
    print_header "Verifying Processed Data | Dataset: $dataset | Histone: $histone | Cell: $cell"
    
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak"
    fi
    
    local treatment_bam="${BASE_DIR}/${histone}/${histone}_${cell}.bam"
    local input_bam="${BASE_DIR}/${histone}/input_${cell}.bam"
    
    print_step "Checking treatment BAM: $(basename "$treatment_bam")"
    if [[ -f "$treatment_bam" ]]; then
        # Check BAM file integrity
        if "$SAMTOOLS" quickcheck "$treatment_bam" 2>/dev/null; then
            local read_count=$("$SAMTOOLS" view -c "$treatment_bam" 2>/dev/null || echo "0")
            echo "  ✅ Valid BAM file with $read_count reads"
        else
            echo "  ❌ Invalid BAM file"
            return 1
        fi
    else
        echo "  ❌ Treatment BAM not found"
        return 1
    fi
    
    print_step "Checking input BAM: $(basename "$input_bam")"
    if [[ -f "$input_bam" ]]; then
        # Check BAM file integrity
        if "$SAMTOOLS" quickcheck "$input_bam" 2>/dev/null; then
            local read_count=$("$SAMTOOLS" view -c "$input_bam" 2>/dev/null || echo "0")
            echo "  ✅ Valid BAM file with $read_count reads"
        else
            echo "  ❌ Invalid BAM file"
            return 1
        fi
    else
        echo "  ❌ Input BAM not found"
        return 1
    fi
    
    echo "✅ All BAM files verified successfully"
    echo ""
    return 0
}

# ===========================================================================================
# DATASET CONFIGURATIONS WITH CELL TYPES
# ===========================================================================================
declare -A DATASET_CONFIGS=(
    # Human PBMC dataset
    ["human_broad_marks"]="H3K27ac H3K27me3 H3K9me3"
    ["human_narrow_marks"]="H3K27ac H3K4me1 H3K4me2 H3K4me3"
    ["human_all_cell_types"]="B CD4T CD8T DC Mono NK otherT other"
    ["human_genome_size"]="2.7e9"
    ["human_genome"]="hg38"
    ["human_chromsizes"]="${HUMAN_DIR}/ref/hg38.chrom.sizes"
    ["human_genome_file"]="${HUMAN_DIR}/ref/genome_file.txt"
    ["human_gene_annot"]="${HUMAN_DIR}/ref/refFlat.dupremoved.txt"
    ["human_data_dir"]="${HUMAN_DIR}/simulation_with_scCTpeak"
    
    # Mouse Brain dataset
    ["mouse_broad_marks"]="H3K27ac H3K27me3 H3K36me3"
    ["mouse_narrow_marks"]="H3K27ac H3K4me3 Olig2 Rad21"
    ["mouse_genome_size"]="1.87e9"
    ["mouse_genome"]="mm10"
    ["mouse_chromsizes"]="${MOUSE_DIR}/ref/mm10.chrom.sizes"
    ["mouse_genome_file"]="${MOUSE_DIR}/ref/mm10.genome_file.txt"
    ["mouse_gene_annot"]="${MOUSE_DIR}/ref/refFlat.dupremoved.txt"
    ["mouse_data_dir"]="${MOUSE_DIR}/simulation_with_scCTpeak"
)

# Cell type mapping for mouse brain
get_mouse_cell_types() {
    local histone="$1"
    case "$histone" in
        H3K27ac) echo "Astrocytes mOL OEC OPC VLMC" ;;
        H3K27me3) echo "Astrocytes Microglia mOL Neurons1 Neurons3 OEC OPC VLMC" ;;
        H3K36me3) echo "Astrocytes mOL OEC OPC" ;;
        H3K4me3) echo "Astrocytes Microglia mOL Neurons1 Neurons2 Neurons3 OEC OPC VLMC" ;;
        Olig2|Rad21) echo "Astrocytes mOL OEC Unknown" ;;
        *) echo "" ;;
    esac
}

# Get specific cell types for mouse brain
get_mouse_cell_type() {
    local histone="$1"
    local cell_type="$2"
    
    IFS=' ' read -r -a ALL_CELLS <<< "$(get_mouse_cell_types "$histone")"
    
    if [[ "$cell_type" == "all" ]]; then
        echo "${ALL_CELLS[@]}"
        return
    fi
    
    for cell in "${ALL_CELLS[@]}"; do
        if [[ "$cell" == "$cell_type" ]]; then
            echo "$cell_type"
            return
        fi
    done
    
    echo ""
}

# Get cell types for any dataset
get_cell_types() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    
    if [[ "$dataset" == "human" ]]; then
        if [[ "$cell_type" == "all" ]]; then
            echo "${DATASET_CONFIGS[human_all_cell_types]}"
        else
            echo "$cell_type"
        fi
    else
        if [[ "$cell_type" == "all" ]]; then
            get_mouse_cell_types "$histone"
        else
            get_mouse_cell_type "$histone" "$cell_type"
        fi
    fi
}

# ===========================================================================================
# UTILITY FUNCTIONS
# ===========================================================================================
validate_bam_files() {
    local base_dir="$1"
    local histone="$2"
    local cell="$3"
    
    local treatment_bam="${base_dir}/${histone}/${histone}_${cell}.bam"
    local input_bam="${base_dir}/${histone}/input_${cell}.bam"
    
    if [[ ! -f "$treatment_bam" ]]; then
        echo "❌ Treatment BAM not found: $treatment_bam"
        return 1
    fi
    
    # Check if input exists (not required for without_input mode)
    if [[ "$4" == "check_input" ]] && [[ ! -f "$input_bam" ]]; then
        echo "❌ Input BAM not found: $input_bam"
        return 1
    fi
    
    return 0
}

# -------------------------------------------------------------------------------------------
# Split BAM Functions
# -------------------------------------------------------------------------------------------

# Get list of cell types from barcode files in a histone directory
get_cell_types_from_barcodes() {
    local histone_dir="$1"
    local barcode_files=("$histone_dir"/*_barcodes.txt)
    local cell_types=()
    for f in "${barcode_files[@]}"; do
        [[ ! -f "$f" ]] && continue
        base=$(basename "$f" _barcodes.txt)
        # Skip temp files
        [[ "$base" == "temp_bam" ]] && continue
        cell_types+=("$base")
    done
    echo "${cell_types[@]}"
}

# Split a combined BAM by cell type
run_split_bam() {
    local dataset="$1"
    local histone="$2"
    
    print_header "Split BAM | Dataset: $dataset | Histone: $histone"
    
    # Set base directories
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${SPLIT_BAM_HUMAN_DIR}"
        OUTPUT_SUBDIR="split_celltype_bams_scCTpeak"   # as in user's script
        REMOVE_BARCODE_SUFFIX=true                      # human needs -1/-2 removal
    else
        BASE_DIR="${SPLIT_BAM_MOUSE_DIR}"
        OUTPUT_SUBDIR="split_celltype_bams_scCTpeak"
        REMOVE_BARCODE_SUFFIX=false
    fi
    
    local histone_dir="${BASE_DIR}/${histone}"
    local combined_bam="${histone_dir}/${histone}.bam"
    local out_dir="${histone_dir}/${OUTPUT_SUBDIR}"
    
    # Check prerequisites
    if [[ ! -d "$histone_dir" ]]; then
        echo "❌ Histone directory not found: $histone_dir"
        return 1
    fi
    if [[ ! -f "$combined_bam" ]]; then
        echo "❌ Combined BAM not found: $combined_bam"
        return 1
    fi
    
    mkdir -p "$out_dir"
    cd "$out_dir" || return 1
    
    # Get cell types from barcode files
    local cell_types=($(get_cell_types_from_barcodes "$histone_dir"))
    if [[ ${#cell_types[@]} -eq 0 ]]; then
        echo "❌ No barcode files found in $histone_dir"
        return 1
    fi
    
    echo "Found cell types: ${cell_types[*]}"
    
    # Step 1: Split BAM by cell type
    for cell in "${cell_types[@]}"; do
        local barcode_file="${histone_dir}/${cell}_barcodes.txt"
        local out_bam="${out_dir}/${histone}_${cell}.bam"
        local temp_sam="${out_dir}/temp_${cell}.sam"
        
        echo "  Splitting $cell ..."
        
        # Build awk script with optional suffix removal
        if [[ "$REMOVE_BARCODE_SUFFIX" == true ]]; then
            # Human: remove trailing -1/-2 from barcodes
            awk_script='
            BEGIN {
                while ((getline line < barcode_file) > 0) {
                    gsub(/-.*$/, "", line)
                    barcodes[line] = 1
                }
                close(barcode_file)
            }
            /^@/ { print; next }
            {
                for (i = 12; i <= NF; i++) {
                    if ($i ~ /^CB:Z:/) {
                        split($i, arr, ":")
                        cb = arr[3]
                        gsub(/-.*$/, "", cb)
                        if (cb in barcodes) {
                            print
                            break
                        }
                    }
                }
            }'
        else
            # Mouse: match barcode exactly
            awk_script='
            BEGIN {
                while ((getline line < barcode_file) > 0) {
                    barcodes[line] = 1
                }
                close(barcode_file)
            }
            /^@/ { print; next }
            {
                for (i = 12; i <= NF; i++) {
                    if ($i ~ /^CB:Z:/) {
                        split($i, arr, ":")
                        cb = arr[3]
                        if (cb in barcodes) {
                            print
                            break
                        }
                    }
                }
            }'
        fi
        
        samtools view -h "$combined_bam" | \
            awk -F'\t' -v barcode_file="$barcode_file" "$awk_script" > "$temp_sam"
        
        # Convert SAM to BAM, sort, and index
        samtools view -bS "$temp_sam" | samtools sort -@ $THREADS -o "$out_bam"
        samtools index "$out_bam"
        rm -f "$temp_sam"
        
        local read_count=$(samtools view -c "$out_bam")
        echo "    ✅ $cell: $read_count reads"
    done
    
    echo "✅ Splitting completed for $histone"
    echo ""
    
    # Step 2: Create input BAMs
    create_input_bams_for_histone "$dataset" "$histone"
}

# Create input BAMs by merging all other cell types
create_input_bams_for_histone() {
    local dataset="$1"
    local histone="$2"
    
    print_header "Create Input BAMs | Dataset: $dataset | Histone: $histone"
    
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${SPLIT_BAM_HUMAN_DIR}"
        OUTPUT_SUBDIR="split_celltype_bams_scCTpeak"
    else
        BASE_DIR="${SPLIT_BAM_MOUSE_DIR}"
        OUTPUT_SUBDIR="split_celltype_bams_scCTpeak"
    fi
    
    local histone_dir="${BASE_DIR}/${histone}"
    local split_dir="${histone_dir}/${OUTPUT_SUBDIR}"
    
    if [[ ! -d "$split_dir" ]]; then
        echo "❌ Split BAM directory not found: $split_dir"
        return 1
    fi
    
    cd "$split_dir" || return 1
    
    # Get all cell type BAMs (exclude input_* and temp files)
    local bam_files=($(ls ${histone}_*.bam 2>/dev/null | grep -v "temp_bam"))
    if [[ ${#bam_files[@]} -eq 0 ]]; then
        echo "❌ No cell type BAM files found in $split_dir"
        return 1
    fi
    
    # Extract cell type names from filenames (remove histone_ prefix and .bam suffix)
    local cell_types=()
    for bam in "${bam_files[@]}"; do
        cell=$(basename "$bam" .bam | sed "s/^${histone}_//")
        cell_types+=("$cell")
    done
    
    echo "Found cell types: ${cell_types[*]}"
    
    # For each cell type, create input BAM by merging all others
    for cell in "${cell_types[@]}"; do
        echo "  Creating input for $cell ..."
        
        local input_bam="input_${cell}.bam"
        local temp_merged="temp_merged_${cell}.bam"
        
        # Collect BAMs of all other cell types
        local other_bams=()
        for other in "${cell_types[@]}"; do
            [[ "$other" == "$cell" ]] && continue
            other_bams+=("${histone}_${other}.bam")
        done
        
        # Check that all files exist and are valid
        local valid_bams=()
        for bam in "${other_bams[@]}"; do
            if [[ -f "$bam" ]] && samtools quickcheck "$bam" 2>/dev/null; then
                valid_bams+=("$bam")
            else
                echo "    ⚠️ Skipping invalid/missing: $bam"
            fi
        done
        
        if [[ ${#valid_bams[@]} -eq 0 ]]; then
            echo "    ❌ No valid BAM files to merge for $cell"
            continue
        fi
        
        # Merge and sort
        samtools merge -f "$temp_merged" "${valid_bams[@]}"
        samtools sort -@ $THREADS -o "$input_bam" "$temp_merged"
        samtools index "$input_bam"
        rm -f "$temp_merged"
        
        local read_count=$(samtools view -c "$input_bam")
        echo "    ✅ Created $input_bam ($read_count reads)"
    done
    
    echo "✅ Input BAM creation completed for $histone"
    echo ""
}

# Batch split for all histones in a dataset
run_split_bam_all() {
    local dataset="$1"
    
    print_header "Batch Split BAM | Dataset: $dataset"
    
    if [[ "$dataset" == "human" ]]; then
        # All human histones (combine broad and narrow lists)
        local histones=(${DATASET_CONFIGS[human_broad_marks]} ${DATASET_CONFIGS[human_narrow_marks]})
        # Remove duplicates
        histones=($(echo "${histones[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    else
        # Mouse histones and TFs
        local histones=(${DATASET_CONFIGS[mouse_broad_marks]} ${DATASET_CONFIGS[mouse_narrow_marks]})
        histones=($(echo "${histones[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    fi
    
    for histone in "${histones[@]}"; do
        echo ""
        echo "═══════════════════════════════════════════════════════════════════"
        echo "  Processing $histone"
        echo "═══════════════════════════════════════════════════════════════════"
        run_split_bam "$dataset" "$histone"
        
        # Small delay between runs
        sleep 2
    done
    
    echo "✅ Batch split BAM completed for $dataset"
    echo ""
}

# ===========================================================================================
# CORE PEAK CALLING FUNCTIONS WITH UNIFIED OUTPUT
# ===========================================================================================

# -------------------------------------------------------------------------------------------
# DROMPAplus parse2wig+ Function
# -------------------------------------------------------------------------------------------
run_drompa_parse2wig() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    
    print_header "DROMPAplus parse2wig+ | Dataset: $dataset | Histone: $histone | Cell: $cell_type"
    
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME_FILE="${DATASET_CONFIGS[human_genome_file]}"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME_FILE="${DATASET_CONFIGS[mouse_genome_file]}"
    fi
    
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    
    if [[ ${#CELL_TYPES[@]} -eq 0 ]] || [[ -z "${CELL_TYPES[0]}" ]]; then
        echo "❌ Error: No valid cell types found for $histone with cell_type=$cell_type"
        return 1
    fi
    
    OUT_DIR="${BASE_DIR}/DROMPAplus_peakbed_test/${histone}"
    mkdir -p "$OUT_DIR"
    cd "$OUT_DIR"
    
    print_step "Base directory: $BASE_DIR"
    print_step "Output directory: $OUT_DIR"
    print_step "Cell types to process: ${CELL_TYPES[*]}"
    
    for cell in "${CELL_TYPES[@]}"; do
        print_step "Processing ${histone} treatment - $cell"
        
        if ! validate_bam_files "$BASE_DIR" "$histone" "$cell"; then
            continue
        fi
        
        "${DROMPA_DIR}/parse2wig+" -i "${BASE_DIR}/${histone}/${histone}_${cell}.bam" \
                   -o ${histone}.${cell} \
                   --pair \
                   --gt "$GENOME_FILE" \
                   -n GR
        
        print_step "Processing input - $cell"
        "${DROMPA_DIR}/parse2wig+" -i "${BASE_DIR}/${histone}/input_${cell}.bam" \
                   -o input.${cell} \
                   --pair \
                   --gt "$GENOME_FILE" \
                   -n GR
    done
    
    echo "✅ Completed parse2wig+ for $histone (Cells: ${CELL_TYPES[*]})"
    echo ""
}

# -------------------------------------------------------------------------------------------
# DROMPAplus Peak Calling with Unified BED Output
# -------------------------------------------------------------------------------------------
run_drompa_peakcalling() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    local mode="$4"       # broad | narrow
    local control="$5"    # with_input | without_input
    
    print_header "DROMPAplus Peak Calling | $dataset | $histone | Cell: $cell_type | $mode | $control"
    
    # ----------------------------
    # Set DROMPA mode & thresholds
    # ----------------------------
    local DROMPA_MODE P_INT P_ENR
    if [[ "$mode" == "broad" ]]; then
        DROMPA_MODE="BROAD"
        P_INT=4
        P_ENR=3
    else
        DROMPA_MODE="SHARP"
        P_INT=5
        P_ENR=4
    fi
    
    # ----------------------------
    # Base directories & genome
    # ----------------------------
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME_FILE="${DATASET_CONFIGS[human_genome_file]}"
        GENE_ANNOT="${DATASET_CONFIGS[human_gene_annot]}"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME_FILE="${DATASET_CONFIGS[mouse_genome_file]}"
        GENE_ANNOT="${DATASET_CONFIGS[mouse_gene_annot]}"
    fi
    
    [[ ! -f "$GENOME_FILE" ]] && { echo "❌ Genome missing: $GENOME_FILE"; return 1; }
    [[ ! -f "$GENE_ANNOT" ]] && { echo "❌ Gene annotation missing: $GENE_ANNOT"; return 1; }
    
    # ----------------------------
    # Get cell types
    # ----------------------------
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    [[ ${#CELL_TYPES[@]} -eq 0 ]] && { echo "❌ No valid cell types for $histone"; return 1; }
    
    PARSE_DIR="${BASE_DIR}/DROMPAplus_peakbed_test/${histone}/parse2wigdir+"
    OUT_DIR="${BASE_DIR}/DROMPAplus_peakbed_test/${histone}/peakbed_${control}_${mode}/"
    UNIFIED_DIR="${BASE_DIR}/unified_peaks/${control}/${mode}/"
    mkdir -p "$OUT_DIR" "$UNIFIED_DIR"
    [[ ! -d "$PARSE_DIR" ]] && { echo "❌ parse2wig+ missing: $PARSE_DIR"; return 1; }
    
    # ----------------------------
    # Prepare DROMPAplus input args
    # ----------------------------
    local -a input_args=()
    local file_count=0
    for cell in "${CELL_TYPES[@]}"; do
        local label
        if [[ "$histone" == "H3K27ac" ]]; then
            [[ "$mode" == "broad" ]] && label="${histone}_b_${cell}" || label="${histone}_s_${cell}"
        else
            label="${histone}_${cell}"
        fi
        
        TREATMENT_BW="${PARSE_DIR}/${histone}.${cell}.100.bw"
        [[ ! -f "$TREATMENT_BW" ]] && { echo "⚠️ Missing treatment: $TREATMENT_BW"; continue; }
        
        if [[ "$control" == "with_input" ]]; then
            CONTROL_BW="${PARSE_DIR}/input.${cell}.100.bw"
            [[ ! -f "$CONTROL_BW" ]] && { echo "⚠️ Missing control: $CONTROL_BW"; continue; }
            input_args+=("-i" "${TREATMENT_BW},${CONTROL_BW},${label},,,100")
        else
            input_args+=("-i" "${TREATMENT_BW},,${label},,,100")
        fi
        ((file_count++))
    done
    [[ $file_count -eq 0 ]] && { echo "❌ No valid input files"; return 1; }
    
    # ----------------------------
    # Run DROMPAplus
    # ----------------------------
    cmd_args=("${DROMPA_DIR}/drompa+" "PC_${DROMPA_MODE}" "${input_args[@]}" \
              "-o" "$OUT_DIR" "--gt" "$GENOME_FILE" "-g" "$GENE_ANNOT" \
              "--lpp" "5" "--showitag" "1" "--callpeak" "--pthre_internal" "$P_INT")
    [[ "$control" == "with_input" ]] && cmd_args+=("--pthre_enrich" "$P_ENR")
    
    LOG_FILE="${OUT_DIR}/drompa_peakcalling.log"
    echo "Command: ${cmd_args[*]}" | tee "$LOG_FILE"
    echo "Started at: $(date)" | tee -a "$LOG_FILE"
    
    if "${cmd_args[@]}" >> "$LOG_FILE" 2>&1; then
        echo "✅ DROMPAplus completed"
        
        # --------------------------------------
        # Sort, filter, and unified BED creation
        # --------------------------------------
        for cell in "${CELL_TYPES[@]}"; do
            # 1️⃣ Try predicted peak filename
            if [[ "$histone" == "H3K27ac" ]]; then
                [[ "$mode" == "broad" ]] && peak_file="${OUT_DIR}/${histone}_b_${cell}.100.bw.peak.bed"
                [[ "$mode" != "broad" ]] && peak_file="${OUT_DIR}/${histone}_s_${cell}.100.bw.peak.bed"
            else
                peak_file="${OUT_DIR}/${histone}_${cell}.100.bw.peak.bed"
            fi
            
            # 2️⃣ Fallback: find in OUT_DIR by cell name
            if [[ ! -f "$peak_file" ]]; then
                peak_file=$(find "$OUT_DIR" -maxdepth 1 -name "*${cell}*.peak.bed" | head -n1)
            fi
            
            if [[ ! -f "$peak_file" ]]; then
                echo "⚠️ Peak missing for $cell. Check DROMPAplus output in $OUT_DIR"
                continue
            fi
            
            # 3️⃣ Define HIST_TAG for naming
            if [[ "$histone" == "H3K27ac" ]]; then
                [[ "$mode" == "broad" ]] && HIST_TAG="H3K27ac-b" || HIST_TAG="H3K27ac-s"
            else
                HIST_TAG="$histone"
            fi
            
            # 4️⃣ Create unified BED
            final_bed="${UNIFIED_DIR}/DROMPAplus_${HIST_TAG}_${cell}.bed"
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' "$peak_file" | \
            sort -k1,1 -k2,2n | grep -E '^chr([0-9]+|X|Y|M)\b' > "$final_bed"
            
            echo "✔ Unified BED: $final_bed"
        done
    else
        echo "❌ DROMPAplus failed. Check $LOG_FILE"
        return 1
    fi
    
    echo ""
}

# -------------------------------------------------------------------------------------------
# MACS2 Functions with Unified Output
# -------------------------------------------------------------------------------------------
run_macs2() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    local mode="$4"
    local control="$5"

    print_header "MACS2 | $dataset | $histone | Cell: $cell_type | $mode | $control"

    local QVAL=0.05
    local ARGS="-q $QVAL --keep-dup all"
    [[ "$mode" == "broad" ]] && ARGS="$ARGS --broad --broad-cutoff $QVAL"

    # Base directories
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME="${DATASET_CONFIGS[human_genome_size]}"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME="${DATASET_CONFIGS[mouse_genome_size]}"
    fi

    # Get cell types
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    [[ ${#CELL_TYPES[@]} -eq 0 ]] && { echo "❌ No valid cell types for $histone"; return 1; }

    # Output directories
    local OUT_DIR="${BASE_DIR}/MACS2_peakbed_test"
    [[ "$control" == "with_input" ]] && OUT_DIR="${OUT_DIR}/With_input" || OUT_DIR="${OUT_DIR}/Without_input"
    local UNIFIED_DIR="${BASE_DIR}/unified_peaks/${control}/${mode}/"
    mkdir -p "$OUT_DIR" "$UNIFIED_DIR"

    for cell in "${CELL_TYPES[@]}"; do
        print_step "Processing: ${histone}_${cell}"

        local TREAT="${BASE_DIR}/${histone}/${histone}_${cell}.bam"
        [[ ! -f "$TREAT" ]] && { echo "⚠️ Missing treatment BAM: $TREAT"; continue; }

        # MACS2 command
        local cmd="macs2 callpeak -t \"$TREAT\" -f BAMPE -g \"$GENOME\" \
                   -n \"${histone}_${cell}_${mode}\" --outdir \"$OUT_DIR\" $ARGS"

        if [[ "$control" == "with_input" ]]; then
            local CTRL="${BASE_DIR}/${histone}/input_${cell}.bam"
            if [[ -f "$CTRL" ]]; then
                cmd="macs2 callpeak -t \"$TREAT\" -c \"$CTRL\" -f BAMPE -g \"$GENOME\" \
                     -n \"${histone}_${cell}_${mode}\" --outdir \"$OUT_DIR\" $ARGS"
            else
                echo "⚠️ Missing control BAM: $CTRL"; continue
            fi
        fi

        print_step "Running MACS2..."
        eval "$cmd"

        if [[ $? -eq 0 ]]; then
            echo "  ✅ MACS2 finished: ${histone}_${cell}"

            # Determine peak file extension
            local peak_ext=".narrowPeak"
            [[ "$mode" == "broad" ]] && peak_ext=".broadPeak"
            local peak_file="${OUT_DIR}/${histone}_${cell}_${mode}_peaks${peak_ext}"

            if [[ -f "$peak_file" ]]; then
                # Naming rule for H3K27ac
                local HISTONE_TAG="$histone"
                if [[ "$histone" == "H3K27ac" ]]; then
                    [[ "$mode" == "broad" ]] && HISTONE_TAG="H3K27ac-b" || HISTONE_TAG="H3K27ac-s"
                fi

                # Unified BED creation
                local unified_bed="${UNIFIED_DIR}/MACS2_${HISTONE_TAG}_${cell}.bed"
                cp "$peak_file" "$unified_bed"

                # Sort BED
                sort -k1,1 -k2,2n "$unified_bed" -o "$unified_bed"
                echo "  ✅ Unified BED sorted: $(basename "$unified_bed")"
            else
                echo "  ⚠️ Peak file not found: $peak_file"
            fi
        else
            echo "  ❌ MACS2 failed: ${histone}_${cell}"
        fi
    done
    echo ""
}


# -------------------------------------------------------------------------------------------
# Genrich Functions with Unified Output
# -------------------------------------------------------------------------------------------
run_genrich() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    local mode="$4"          # broad | narrow
    local control="$5"       # with_input | without_input
    
    print_header "Genrich | $dataset | $histone | Cell: $cell_type | $mode | $control"
    
    local SORT_TAG="qname_sorted_${mode:0:1}"
    
    if [[ "$mode" == "broad" ]]; then
        GENRICH_ARGS="-a 100 -l 500 -g 1000 -p 0.05 -f BAM"
        EXT=".broadPeak"
    else
        GENRICH_ARGS="-a 200 -l 100 -g 100 -p 0.01 -f BAM"
        EXT=".narrowPeak"
    fi
    
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
    fi
    
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    [[ ${#CELL_TYPES[@]} -eq 0 ]] && { echo "❌ No valid cell types found"; return 1; }
    
    local MARK_DIR="${BASE_DIR}/${histone}"
    local SORTED_DIR="${MARK_DIR}/${SORT_TAG}"
    
    local OUT_BASE="${BASE_DIR}/Genrich_peakbed_test"
    local OUT_DIR="${OUT_BASE}/${histone}_${mode}/${control}"
    
    # Unified output directory
    UNIFIED_DIR="${BASE_DIR}/unified_peaks/${control}/${mode}"
    
    mkdir -p "$OUT_DIR" "$SORTED_DIR" "$UNIFIED_DIR"
    
    for cell in "${CELL_TYPES[@]}"; do
        print_step "Processing: $histone ($cell)"
        
        local RAW_TREATMENT="${MARK_DIR}/${histone}_${cell}.bam"
        local TREATMENT="${SORTED_DIR}/${histone}_${cell}_qnamesorted.bam"
        local OUTPUT_FILE="${OUT_DIR}/${cell}_${histone}_${mode}_peaks${EXT}"
        local LOG_FILE="${OUT_DIR}/${cell}_${histone}_${mode}_genrich.log"
        
        [[ ! -f "$RAW_TREATMENT" ]] && { echo "❌ Missing BAM: $RAW_TREATMENT"; continue; }
        
        # Sort BAM if needed
        if [[ ! -f "$TREATMENT" ]]; then
            samtools sort -n -@ $THREADS -o "$TREATMENT" "$RAW_TREATMENT" >> "$LOG_FILE" 2>&1
        fi
        
        local cmd="Genrich -t \"$TREATMENT\" -o \"$OUTPUT_FILE\" -j -r -v $GENRICH_ARGS"
        
        if [[ "$control" == "with_input" ]]; then
            local RAW_CONTROL="${MARK_DIR}/input_${cell}.bam"
            local CONTROL="${SORTED_DIR}/input_${cell}_qnamesorted.bam"
            
            [[ ! -f "$RAW_CONTROL" ]] && { echo "❌ Missing CTRL: $RAW_CONTROL"; continue; }
            
            if [[ ! -f "$CONTROL" ]]; then
                samtools sort -n -@ $THREADS -o "$CONTROL" "$RAW_CONTROL" >> "$LOG_FILE" 2>&1
            fi
            
            cmd="Genrich -t \"$TREATMENT\" -c \"$CONTROL\" -o \"$OUTPUT_FILE\" -j -r -v $GENRICH_ARGS"
        fi
        
        print_step "Running Genrich..."
        eval "$cmd" >> "$LOG_FILE" 2>&1
        
        if [[ $? -eq 0 ]]; then
            echo "  ✅ Genrich completed: $cell"
            
            # ----------------------------
            # Unified BED (Correct naming)
            # ----------------------------
            if [[ -f "$OUTPUT_FILE" ]]; then
                # Naming rule
                local HISTONE_TAG="$histone"
                if [[ "$histone" == "H3K27ac" ]]; then
                    [[ "$mode" == "broad" ]] && HISTONE_TAG="H3K27ac-b" || HISTONE_TAG="H3K27ac-s"
                fi
                
                local FINAL_BED="${UNIFIED_DIR}/Genrich_${HISTONE_TAG}_${cell}.bed"
                
                awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' "$OUTPUT_FILE" | \
                sort -k1,1 -k2,2n | \
                grep -E '^chr([0-9]+|X|Y|M)\b' > "$FINAL_BED"
                
                if [[ -s "$FINAL_BED" ]]; then
                    echo "  ✅ Unified BED created: $(basename "$FINAL_BED")"
                else
                    echo "  ⚠️ Unified BED empty: $FINAL_BED"
                fi
            else
                echo "  ⚠️ Peak file missing: $OUTPUT_FILE"
            fi
        else
            echo "  ❌ Genrich failed: $cell"
        fi
    done
    echo ""
}

# -------------------------------------------------------------------------------------------
# GoPeaks Functions with Unified Output
# -------------------------------------------------------------------------------------------
run_gopeaks() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    local mode="$4"          # broad | narrow
    local control="$5"       # with_input | without_input
    
    print_header "GoPeaks | $dataset | $histone | Cell: $cell_type | $mode | $control"
    
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        CHROMSIZE="${DATASET_CONFIGS[human_chromsizes]}"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        CHROMSIZE="${DATASET_CONFIGS[mouse_chromsizes]}"
    fi
    
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    [[ ${#CELL_TYPES[@]} -eq 0 ]] && { echo "❌ No valid cell types found"; return 1; }
    
    local OUT_BASE="${BASE_DIR}/GoPeaks_peakbed_test"
    local OUT_DIR="${OUT_BASE}/${histone}_${mode}/${control}"
    
    # Unified output directory
    UNIFIED_DIR="${BASE_DIR}/unified_peaks/${control}/${mode}"
    
    mkdir -p "$OUT_DIR" "$UNIFIED_DIR"
    
    for cell in "${CELL_TYPES[@]}"; do
        print_step "Processing: $histone ($cell)"
        
        local BAM="${BASE_DIR}/${histone}/${histone}_${cell}.bam"
        local PREFIX="${OUT_DIR}/${cell}_${histone}_${mode}"
        local PEAK_FILE="${PREFIX}_peaks.bed"
        
        [[ ! -f "$BAM" ]] && { echo "⚠️ Missing BAM: $BAM"; continue; }
        
        # ----------------------------
        # Build GoPeaks command
        # ----------------------------
        local cmd="gopeaks -b \"$BAM\" -s \"$CHROMSIZE\" -o \"$PREFIX\""
        [[ "$mode" == "broad" ]] && cmd="$cmd --broad"
        
        if [[ "$control" == "with_input" ]]; then
            local CTRL="${BASE_DIR}/${histone}/input_${cell}.bam"
            [[ ! -f "$CTRL" ]] && { echo "⚠️ Missing CTRL: $CTRL"; continue; }
            
            cmd="gopeaks -b \"$BAM\" -c \"$CTRL\" -s \"$CHROMSIZE\" -o \"$PREFIX\""
            [[ "$mode" == "broad" ]] && cmd="$cmd --broad"
        fi
        
        print_step "Running GoPeaks..."
        eval "$cmd"
        
        if [[ $? -eq 0 ]]; then
            echo "  ✅ GoPeaks completed: $cell"
            
            # ----------------------------
            # Unified BED (Correct naming)
            # ----------------------------
            if [[ -f "$PEAK_FILE" ]]; then
                # Naming rule (same as Genrich)
                local HISTONE_TAG="$histone"
                if [[ "$histone" == "H3K27ac" ]]; then
                    [[ "$mode" == "broad" ]] && HISTONE_TAG="H3K27ac-b" || HISTONE_TAG="H3K27ac-s"
                fi
                
                local FINAL_BED="${UNIFIED_DIR}/GoPeaks_${HISTONE_TAG}_${cell}.bed"
                
                awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' "$PEAK_FILE" | \
                sort -k1,1 -k2,2n | \
                grep -E '^chr([0-9]+|X|Y|M)\b' > "$FINAL_BED"
                
                if [[ -s "$FINAL_BED" ]]; then
                    echo "  ✅ Unified BED created: $(basename "$FINAL_BED")"
                else
                    echo "  ⚠️ Unified BED empty: $FINAL_BED"
                fi
            else
                echo "  ⚠️ Peak file missing: $PEAK_FILE"
            fi
        else
            echo "  ❌ GoPeaks failed: $cell"
        fi
    done
    echo ""
}

# -------------------------------------------------------------------------------------------
# HOMER Functions with Unified Output
# -------------------------------------------------------------------------------------------
run_homer() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    local mode="$4"        # narrow | broad
    local control="$5"     # with_input | without_input
    
    print_header "HOMER | $dataset | $histone | Cell: $cell_type | $mode | $control"
    
    local STYLE="factor"
    [[ "$mode" == "broad" ]] && STYLE="histone"
    
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
    fi
    
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    [[ ${#CELL_TYPES[@]} -eq 0 ]] && { echo "❌ No valid cell types"; return 1; }
    
    local OUT_BASE="${BASE_DIR}/HOMER_peakbed_test/${control}/${mode}"
    local TAG_BASE="${OUT_BASE}/tag_directories"
    local UNIFIED_DIR="${BASE_DIR}/unified_peaks/${control}/${mode}"
    
    mkdir -p "$OUT_BASE" "$TAG_BASE" "$UNIFIED_DIR"
    
    for cell in "${CELL_TYPES[@]}"; do
        print_step "Processing: ${histone}_${cell}"
        
        local TREAT_BAM="${BASE_DIR}/${histone}/${histone}_${cell}.bam"
        [[ ! -f "$TREAT_BAM" ]] && { echo "⚠️ Missing BAM: $TREAT_BAM"; continue; }
        
        local TREAT_TAG="${TAG_BASE}/${histone}_${cell}_treat"
        local PEAK_TXT="${OUT_BASE}/${histone}_${cell}_${mode}.txt"
        local RAW_BED="${OUT_BASE}/${histone}_${cell}_${mode}.raw.bed"
        
        # -------------------------------
        # Create tag directory
        # -------------------------------
        makeTagDirectory "$TREAT_TAG" "$TREAT_BAM" 2>/dev/null
        
        local cmd="findPeaks \"$TREAT_TAG\" -style \"$STYLE\" -o \"$PEAK_TXT\""
        
        if [[ "$control" == "with_input" ]]; then
            local CTRL_BAM="${BASE_DIR}/${histone}/input_${cell}.bam"
            local CTRL_TAG="${TAG_BASE}/${histone}_${cell}_ctrl"
            [[ ! -f "$CTRL_BAM" ]] && { echo "⚠️ Missing CTRL: $CTRL_BAM"; continue; }
            makeTagDirectory "$CTRL_TAG" "$CTRL_BAM" 2>/dev/null
            cmd="findPeaks \"$TREAT_TAG\" -style \"$STYLE\" -i \"$CTRL_TAG\" -o \"$PEAK_TXT\""
        fi
        
        print_step "Running HOMER..."
        eval "$cmd" 2>/dev/null
        [[ ! -s "$PEAK_TXT" ]] && { echo "⚠️ No peaks detected"; continue; }
        
        # -------------------------------
        # TXT → BED6 (remove header)
        # -------------------------------
        awk 'BEGIN{OFS="\t"} NR>1 {print $2,$3,$4,$1,$8,"."}' \
        "$PEAK_TXT" > "$RAW_BED"
        
        # -------------------------------
        # Unified naming logic
        # -------------------------------
        local FINAL_BED=""
        
        if [[ "$histone" == "H3K27ac" ]]; then
            if [[ "$mode" == "narrow" ]]; then
                FINAL_BED="${UNIFIED_DIR}/HOMER_H3K27ac-s_${cell}.bed"
            else
                FINAL_BED="${UNIFIED_DIR}/HOMER_H3K27ac-b_${cell}.bed"
            fi
        else
            FINAL_BED="${UNIFIED_DIR}/HOMER_${histone}_${cell}.bed"
        fi
        
        # -------------------------------
        # Sort + clean (header-free BED6)
        # -------------------------------
        sort -k1,1V -k2,2n "$RAW_BED" | \
        awk 'BEGIN{OFS="\t"}
        $1 ~ /^chr([0-9]+|X|Y|M)$/ && NF==6 {
            print $1,$2,$3,$4,$5,$6
        }' > "$FINAL_BED"

        rm -f "$RAW_BED"
        
        echo "✅ HOMER BED created (sorted, header-free): $(basename "$FINAL_BED")"
    done
    
    echo ""
}

# -------------------------------------------------------------------------------------------
# SEACR Functions with Unified Output
# -------------------------------------------------------------------------------------------
run_seacr() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    local mode="$4"        # narrow | broad
    local control="$5"     # with_input | without_input
    
    print_header "SEACR | $dataset | $histone | Cell: $cell_type | $mode | $control"
    
    # -------------------------------
    # Threshold logic
    # -------------------------------
    local THRESHOLD="relaxed"
    [[ "$mode" == "broad" ]] && THRESHOLD="stringent"
    
    # -------------------------------
    # Base directories
    # -------------------------------
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
    fi
    
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    [[ ${#CELL_TYPES[@]} -eq 0 ]] && { echo "❌ No valid cell types"; return 1; }
    
    local OUT_BASE="${BASE_DIR}/SEACR_peakbed_test/${control}/${mode}"
    local TMP_DIR="${OUT_BASE}/bedgraphs"
    local UNIFIED_DIR="${BASE_DIR}/unified_peaks/${control}/${mode}"
    
    mkdir -p "$OUT_BASE" "$TMP_DIR" "$UNIFIED_DIR"
    
    for cell in "${CELL_TYPES[@]}"; do
        print_step "Processing: ${histone}_${cell}"
        
        local TREAT_BAM="${BASE_DIR}/${histone}/${histone}_${cell}.bam"
        [[ ! -f "$TREAT_BAM" ]] && { echo "⚠️ Missing BAM: $TREAT_BAM"; continue; }
        
        local TREAT_BG="${TMP_DIR}/${histone}_${cell}_treat.bedgraph"
        local PREFIX="${OUT_BASE}/${histone}_${cell}"
        
        # -------------------------------
        # Generate treatment bedGraph
        # -------------------------------
        print_step "Generating treatment bedGraph..."
        bedtools genomecov -bg -ibam "$TREAT_BAM" 2>/dev/null | \
        LC_ALL=C sort -k1,1V -k2,2n > "$TREAT_BG"
        
        # -------------------------------
        # SEACR command
        # -------------------------------
        local cmd="bash \"$SEACR_SCRIPT\" \"$TREAT_BG\" 0.01 non \"$THRESHOLD\" \"$PREFIX\""
        
        if [[ "$control" == "with_input" ]]; then
            local CTRL_BAM="${BASE_DIR}/${histone}/input_${cell}.bam"
            [[ ! -f "$CTRL_BAM" ]] && { echo "⚠️ Missing control BAM: $CTRL_BAM"; continue; }
            
            local CTRL_BG="${TMP_DIR}/${histone}_${cell}_ctrl.bedgraph"
            print_step "Generating control bedGraph..."
            bedtools genomecov -bg -ibam "$CTRL_BAM" 2>/dev/null | \
            LC_ALL=C sort -k1,1V -k2,2n > "$CTRL_BG"
            
            cmd="bash \"$SEACR_SCRIPT\" \"$TREAT_BG\" \"$CTRL_BG\" norm \"$THRESHOLD\" \"$PREFIX\""
        fi
        
        # -------------------------------
        # Run SEACR
        # -------------------------------
        print_step "Running SEACR..."
        eval "$cmd" 2>/dev/null || { echo "❌ SEACR failed"; continue; }
        
        local RAW_BED="${PREFIX}.${THRESHOLD}.bed"
        [[ ! -f "$RAW_BED" ]] && { echo "⚠️ Missing peak BED"; continue; }
        
        # -------------------------------
        # Unified file naming logic
        # -------------------------------
        local FINAL_BED=""
        
        if [[ "$histone" == "H3K27ac" ]]; then
            if [[ "$mode" == "narrow" ]]; then
                FINAL_BED="${UNIFIED_DIR}/SEACR_H3K27ac-s_${cell}.bed"
            else
                FINAL_BED="${UNIFIED_DIR}/SEACR_H3K27ac-b_${cell}.bed"
            fi
        else
            FINAL_BED="${UNIFIED_DIR}/SEACR_${histone}_${cell}.bed"
        fi
        
        # -------------------------------
        # Clean → BED6 → sort → remove header
        # -------------------------------
        awk -v prefix="SEACR_${histone}_${cell}" '
        BEGIN{OFS="\t"}
        $1 ~ /^chr([0-9]+|X|Y|M)$/ {
            print $1,$2,$3,prefix"-"NR,$4,"."
        }' "$RAW_BED" | \
        sort -k1,1V -k2,2n > "$FINAL_BED"
        
        echo "✅ Unified SEACR BED created: $(basename "$FINAL_BED")"
    done
    
    echo ""
}

# -------------------------------------------------------------------------------------------
# SICER2 Functions with Unified Output
# -------------------------------------------------------------------------------------------
run_sicer2() {
    local dataset="$1"
    local histone="$2"
    local cell_type="$3"
    local mode="$4"        # narrow | broad
    local control="$5"
    
    print_header "SICER2 | $dataset | $histone | Cell: $cell_type | $mode | $control"
    
    local FDR=0.05
    local WINDOW GAP
    if [[ "$mode" == "broad" ]]; then
        WINDOW=100
        GAP=200
    else
        WINDOW=50
        GAP=100
    fi
    
    if [[ "$dataset" == "human" ]]; then
        BASE_DIR="${HUMAN_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME="${DATASET_CONFIGS[human_genome]}"
    else
        BASE_DIR="${MOUSE_DIR}/simulation_with_scCTpeak/celltypeswisebam_50bp"
        GENOME="${DATASET_CONFIGS[mouse_genome]}"
    fi
    
    IFS=' ' read -r -a CELL_TYPES <<< "$(get_cell_types "$dataset" "$histone" "$cell_type")"
    if [[ ${#CELL_TYPES[@]} -eq 0 ]] || [[ -z "${CELL_TYPES[0]}" ]]; then
        echo "❌ Error: No valid cell types found for $histone with cell_type=$cell_type"
        return 1
    fi
    
    local OUT_BASE="${BASE_DIR}/SICER2_peakbed_test"
    local OUT_DIR="${OUT_BASE}/${histone}_${mode}"
    [[ "$control" == "with_input" ]] && OUT_DIR="${OUT_DIR}/With_input" || OUT_DIR="${OUT_DIR}/Without_input"
    
    local UNIFIED_DIR="${BASE_DIR}/unified_peaks/${control}/${mode}"
    mkdir -p "$OUT_DIR" "$UNIFIED_DIR"
    
    for cell in "${CELL_TYPES[@]}"; do
        print_step "Processing: ${histone}_${cell}"
        local CELL_OUT="${OUT_DIR}/${histone}_${cell}"
        mkdir -p "$CELL_OUT"
        
        local TREATMENT_BAM="${BASE_DIR}/${histone}/${histone}_${cell}.bam"
        if [[ ! -f "$TREATMENT_BAM" ]]; then
            echo " ⚠️ Missing treatment BAM: $TREATMENT_BAM"
            continue
        fi
        
        local cmd="sicer --treatment_file \"$TREATMENT_BAM\" --species \"$GENOME\" --fragment_size 150 --window_size $WINDOW --gap_size $GAP --effective_genome_fraction 0.80 --false_discovery_rate $FDR --redundancy_threshold 1 --output_directory \"$CELL_OUT\" --cpu $THREADS --significant_reads"
        
        if [[ "$control" == "with_input" ]]; then
            local CONTROL_BAM="${BASE_DIR}/${histone}/input_${cell}.bam"
            if [[ -f "$CONTROL_BAM" ]]; then
                cmd="sicer --treatment_file \"$TREATMENT_BAM\" --control_file \"$CONTROL_BAM\" --species \"$GENOME\" --fragment_size 150 --window_size $WINDOW --gap_size $GAP --effective_genome_fraction 0.80 --false_discovery_rate $FDR --redundancy_threshold 1 --output_directory \"$CELL_OUT\" --cpu $THREADS --significant_reads"
            else
                echo " ⚠️ Missing control BAM: $CONTROL_BAM"
                continue
            fi
        else
            cmd="$cmd --e_value 1000"
        fi
        
        print_step "Running SICER2..."
        eval "$cmd" > "${CELL_OUT}/sicer.log" 2>&1
        
        # FIXED: Find SICER2 output files correctly
        local SCORE_FILE=""
        
        # First, try to find any relevant peak file
        if [[ "$control" == "with_input" ]]; then
            # For with_input: look for -island.bed files
            SCORE_FILE=$(find "$CELL_OUT" -name "*-island.bed" | head -n 1)
        else
            # For without_input: FIRST try -FDR-island.bed, then .scoreisland
            SCORE_FILE=$(find "$CELL_OUT" -name "*-FDR${FDR}-island.bed" | head -n 1)
            
            # If -FDR-island.bed not found, use .scoreisland
            if [[ ! -f "$SCORE_FILE" ]]; then
                SCORE_FILE=$(find "$CELL_OUT" -name "*-W${WINDOW}-G${GAP}.scoreisland" | head -n 1)
                
                # If specific pattern not found, try any .scoreisland
                if [[ ! -f "$SCORE_FILE" ]]; then
                    SCORE_FILE=$(find "$CELL_OUT" -name "*.scoreisland" | head -n 1)
                fi
            fi
        fi
        
        # Debug: Show what file was found
        if [[ -f "$SCORE_FILE" ]]; then
            echo "  Found SICER2 output: $(basename "$SCORE_FILE")"
        else
            echo "  ⚠️ No peak file found. Looking in: $CELL_OUT"
            echo "  Files present:"
            ls -la "$CELL_OUT/" 2>/dev/null || echo "    Directory not found"
            continue
        fi
        
        if [[ -f "$SCORE_FILE" ]]; then
            OUT_BED="${OUT_DIR}/SICER2_${histone}_${cell}_peaks.sorted.bed"
            
            # Create BED6 format from the scoreisland file
            awk -v prefix="SICER2_${histone}_${cell}" 'BEGIN{OFS="\t"} {
                if ($1 ~ /^chr/) {
                    print $1, $2, $3, prefix"-"NR, $5, "."
                }
            }' "$SCORE_FILE" | sort -k1,1 -k2,2n > "$OUT_BED"
            
            local peak_count=$(wc -l < "$OUT_BED" 2>/dev/null || echo 0)
            echo "✔ Created sorted BED: $OUT_BED ($peak_count peaks)"
            
            # Filter only standard chromosomes (chr1-22, X, Y, M)
            FILTERED_BED="${OUT_DIR}/SICER2_${histone}_${cell}.bed"
            grep -E '^chr([0-9]+|X|Y|M)\b' "$OUT_BED" > "$FILTERED_BED"
            
            local filtered_count=$(wc -l < "$FILTERED_BED" 2>/dev/null || echo 0)
            echo "✔ Filtered BED: $FILTERED_BED ($filtered_count peaks after filtering)"
            
            # naming histone
            local hist_name="$histone"
            if [[ "$histone" == "H3K27ac" ]]; then
                if [[ "$mode" == "broad" ]]; then
                    hist_name="H3K27ac-b"
                else
                    hist_name="H3K27ac-s"
                fi
            fi
            
            # Unified files name
            local unified_name="SICER2_${hist_name}_${cell}.bed"
            local unified_file="${UNIFIED_DIR}/${unified_name}"
            
            # Copy filtered BED to unified location
            cp "$FILTERED_BED" "$unified_file"
            
            # Count peaks
            local unified_count=$(wc -l < "$unified_file" 2>/dev/null || echo 0)
            
            if [[ $unified_count -gt 0 ]]; then
                echo "✅ Created unified BED: $unified_name ($unified_count peaks)"
            else
                echo "⚠️ Created but empty: $unified_name"
            fi
        else
            echo "✘ SICER2 peak file missing: ${histone}_${cell}"
        fi
    done
    echo ""
}

# ===========================================================================================
# MAIN scCTpeak FUNCTION
# ===========================================================================================
scCTpeak() {
    local command="$1"
    shift
    
    # Show welcome message for first run
    if [[ -z "$SCCTPEAK_FIRST_RUN" ]]; then
        echo ""
        echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
        echo "║                        scCTpeak v$SCCTPEAK_VERSION                               ║"
        echo "║          Complete Benchmark Suite for scCUT&Tag Peak-Calling                     ║"
        echo "║        WITH DATA PROCESSING MODULE & UNIFIED OUTPUT FORMAT                       ║"
        echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
        echo ""
        export SCCTPEAK_FIRST_RUN=1
    fi
    
    if [[ -z "$command" ]]; then
        echo "Usage: scCTpeak <command> [options]"
        echo "Use 'scCTpeak help' for detailed usage"
        return 1
    fi
    
    case "$command" in
        # =================================================================
        # DATA PROCESSING COMMANDS
        # =================================================================
        "process")
            local dataset="$1"
            local histone="$2"
            
            if [[ -z "$dataset" ]] || [[ -z "$histone" ]]; then
                echo "Usage: scCTpeak process <dataset> <histone>"
                echo "       scCTpeak process human H3K27ac"
                echo "       scCTpeak process mouse H3K4me3"
                echo "       scCTpeak process mouse Olig2"
                echo ""
                echo "Run 'scCTpeak process_all' to process all histones/TFs"
                return 1
            fi
            
            run_data_processing "$dataset" "$histone"
            ;;
        
        "process_all")
            local dataset="$1"
            
            if [[ -z "$dataset" ]]; then
                echo "Usage: scCTpeak process_all <dataset>"
                echo "       scCTpeak process_all human"
                echo "       scCTpeak process_all mouse"
                return 1
            fi
            
            run_batch_data_processing "$dataset"
            ;;
        
        "verify")
            local dataset="$1"
            local histone="$2"
            local cell="$3"
            
            if [[ -z "$dataset" ]] || [[ -z "$histone" ]] || [[ -z "$cell" ]]; then
                echo "Usage: scCTpeak verify <dataset> <histone> <cell>"
                echo "Example: scCTpeak verify human H3K27ac B"
                return 1
            fi
            
            verify_processed_data "$dataset" "$histone" "$cell"
            ;;
        
        # =================================================================
        # SPLIT BAM COMMANDS
        # =================================================================
        "split_bam")
            local dataset="$1"
            local histone="$2"
            
            if [[ -z "$dataset" ]] || [[ -z "$histone" ]]; then
                echo "Usage: scCTpeak split_bam <dataset> <histone>"
                echo "Example: scCTpeak split_bam human H3K27ac"
                echo "         scCTpeak split_bam mouse Olig2"
                return 1
            fi
            
            run_split_bam "$dataset" "$histone"
            ;;
        
        "split_bam_all")
            local dataset="$1"
            
            if [[ -z "$dataset" ]]; then
                echo "Usage: scCTpeak split_bam_all <dataset>"
                echo "Example: scCTpeak split_bam_all human"
                echo "         scCTpeak split_bam_all mouse"
                return 1
            fi
            
            run_split_bam_all "$dataset"
            ;;
        
        # =================================================================
        # PEAK CALLING COMMANDS
        # =================================================================
        "run")
            local dataset="$1"
            local tool="$2"
            local histone="$3"
            local cell_type="$4"
            local mode="$5"
            local control="$6"
            
            case "$tool" in
                drompa)
                    run_drompa_peakcalling "$dataset" "$histone" "$cell_type" "$mode" "$control"
                    ;;
                macs2)
                    run_macs2 "$dataset" "$histone" "$cell_type" "$mode" "$control"
                    ;;
                genrich)
                    run_genrich "$dataset" "$histone" "$cell_type" "$mode" "$control"
                    ;;
                gopeaks)
                    run_gopeaks "$dataset" "$histone" "$cell_type" "$mode" "$control"
                    ;;
                homer)
                    run_homer "$dataset" "$histone" "$cell_type" "$mode" "$control"
                    ;;
                seacr)
                    run_seacr "$dataset" "$histone" "$cell_type" "$mode" "$control"
                    ;;
                sicer2)
                    run_sicer2 "$dataset" "$histone" "$cell_type" "$mode" "$control"
                    ;;
                *)
                    echo "Error: tool must be one of: drompa, genrich, gopeaks, macs2, homer, seacr, sicer2"
                    return 1
                    ;;
            esac
            ;;
        
        "parse2wig")
            local dataset="$1"
            local histone="$2"
            local cell_type="$3"
            
            if [[ -z "$dataset" ]] || [[ -z "$histone" ]] || [[ -z "$cell_type" ]]; then
                echo "Usage: scCTpeak parse2wig <dataset> <histone> <cell_type>"
                echo "Example: scCTpeak parse2wig human H3K27ac B"
                return 1
            fi
            run_drompa_parse2wig "$dataset" "$histone" "$cell_type"
            ;;
        
        "batch")
            local batch_dataset="$1"
            local batch_tool="$2"
            local batch_cell_type="$3"
            local batch_control="$4"
            
            if [[ -z "$batch_dataset" ]] || [[ -z "$batch_tool" ]] || [[ -z "$batch_cell_type" ]] || [[ -z "$batch_control" ]]; then
                echo "Usage: scCTpeak batch <dataset> <tool> <cell_type> <control>"
                echo "Example: scCTpeak batch human macs2 all with_input"
                return 1
            fi
            
            print_header "BATCH PROCESSING: $batch_dataset | $batch_tool | Cell: $batch_cell_type | $batch_control"
            
            if [[ "$batch_dataset" == "human" ]]; then
                BROAD_MARKS=(${DATASET_CONFIGS[human_broad_marks]})
                NARROW_MARKS=(${DATASET_CONFIGS[human_narrow_marks]})
            else
                BROAD_MARKS=(${DATASET_CONFIGS[mouse_broad_marks]})
                NARROW_MARKS=(${DATASET_CONFIGS[mouse_narrow_marks]})
            fi
            
            # Process broad marks
            for histone in "${BROAD_MARKS[@]}"; do
                if [[ "$histone" == "H3K27ac" ]]; then
                    scCTpeak run "$batch_dataset" "$batch_tool" "$histone" "$batch_cell_type" "broad" "$batch_control"
                    scCTpeak run "$batch_dataset" "$batch_tool" "$histone" "$batch_cell_type" "narrow" "$batch_control"
                else
                    scCTpeak run "$batch_dataset" "$batch_tool" "$histone" "$batch_cell_type" "broad" "$batch_control"
                fi
            done
            
            # Process narrow marks (skip H3K27ac since already processed)
            for histone in "${NARROW_MARKS[@]}"; do
                [[ "$histone" == "H3K27ac" ]] && continue
                scCTpeak run "$batch_dataset" "$batch_tool" "$histone" "$batch_cell_type" "narrow" "$batch_control"
            done
            
            echo "✅ BATCH PROCESSING COMPLETED!"
            echo ""
            ;;
        
        "all_tools")
            local all_dataset="$1"
            local all_cell_type="$2"
            local all_control="$3"
            
            if [[ -z "$all_dataset" ]] || [[ -z "$all_cell_type" ]] || [[ -z "$all_control" ]]; then
                echo "Usage: scCTpeak all_tools <dataset> <cell_type> <control>"
                echo "Example: scCTpeak all_tools human B with_input"
                return 1
            fi
            
            print_header "RUNNING ALL TOOLS: $all_dataset | Cell: $all_cell_type | $all_control"
            
            local TOOLS=("drompa" "genrich" "gopeaks" "macs2" "homer" "seacr" "sicer2")
            
            for tool in "${TOOLS[@]}"; do
                echo ""
                echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
                echo "║                                TOOL: $tool                                        ║"
                echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
                echo ""
                
                if [[ "$tool" == "drompa" ]]; then
                    if [[ "$all_dataset" == "human" ]]; then
                        HISTONES=(${DATASET_CONFIGS[human_broad_marks]} ${DATASET_CONFIGS[human_narrow_marks]})
                    else
                        HISTONES=(${DATASET_CONFIGS[mouse_broad_marks]} ${DATASET_CONFIGS[mouse_narrow_marks]})
                    fi
                    
                    for histone in "${HISTONES[@]}"; do
                        echo "Running parse2wig for $histone..."
                        scCTpeak parse2wig "$all_dataset" "$histone" "$all_cell_type"
                    done
                fi
                
                scCTpeak batch "$all_dataset" "$tool" "$all_cell_type" "$all_control"
            done
            
            echo "✅ ALL TOOLS COMPLETED SUCCESSFULLY!"
            echo ""
            ;;
        
        # =================================================================
        # ORGANIZATION COMMANDS
        # =================================================================
        "organize")
            local dataset="$1"
            local control="$2"
            
            if [[ -z "$dataset" ]] || [[ -z "$control" ]]; then
                echo "Usage: scCTpeak organize <dataset> <control>"
                echo "Example: scCTpeak organize human with_input"
                return 1
            fi
            
            # This function would organize unified outputs
            echo "Organization function would go here"
            ;;
        
        "summary")
            local dataset="$1"
            local control="$2"
            
            if [[ -z "$dataset" ]] || [[ -z "$control" ]]; then
                echo "Usage: scCTpeak summary <dataset> <control>"
                echo "Example: scCTpeak summary human with_input"
                return 1
            fi
            
            # This function would show summary
            echo "Summary function would go here"
            ;;
        
        # =================================================================
        # HELP COMMAND
        # =================================================================
        "help"|"-h"|"--help")
            echo ""
            echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
            echo "║                      scCTpeak v$SCCTPEAK_VERSION - HELP                          ║"
            echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
            echo ""
            echo "COMPLETE WORKFLOW:"
            echo "------------------"
            echo "  1. Process data (R-based):"
            echo "     scCTpeak process human H3K27ac"
            echo "     scCTpeak process_all human"
            echo "     scCTpeak process_all mouse"
            echo ""
            echo "  2. Verify processed data:"
            echo "     scCTpeak verify human H3K27ac B"
            echo ""
            echo "  3. Split real BAM files (if you have combined BAMs):"
            echo "     scCTpeak split_bam human H3K27ac"
            echo "     scCTpeak split_bam_all mouse"
            echo ""
            echo "  4. Run peak calling:"
            echo "     scCTpeak run human macs2 H3K27ac B broad with_input"
            echo "     scCTpeak all_tools human all with_input"
            echo ""
            echo "  5. Organize results:"
            echo "     scCTpeak organize human with_input"
            echo "     scCTpeak summary human with_input"
            echo ""
            echo "DATA PROCESSING COMMANDS:"
            echo "------------------------"
            echo "  process <dataset> <histone>"
            echo "      Process single histone/TF to BAM files"
            echo "      Examples:"
            echo "        scCTpeak process human H3K27ac"
            echo "        scCTpeak process mouse H3K4me3"
            echo ""
            echo "  process_all <dataset>"
            echo "      Process all histones/TFs for dataset"
            echo "      Examples:"
            echo "        scCTpeak process_all human"
            echo "        scCTpeak process_all mouse"
            echo ""
            echo "  verify <dataset> <histone> <cell>"
            echo "      Verify generated BAM files"
            echo "      Example: scCTpeak verify human H3K27ac B"
            echo ""
            echo "SPLIT BAM COMMANDS:"
            echo "------------------"
            echo "  split_bam <dataset> <histone>"
            echo "      Split a combined BAM file into cell‑type BAMs using barcode lists,"
            echo "      then create input (control) BAMs by merging other cell types."
            echo "      Example: scCTpeak split_bam human H3K27ac"
            echo ""
            echo "  split_bam_all <dataset>"
            echo "      Run split_bam for all histones/TFs in the dataset."
            echo "      Example: scCTpeak split_bam_all human"
            echo ""
            echo "PEAK CALLING COMMANDS:"
            echo "---------------------"
            echo "  run <dataset> <tool> <histone> <cell_type> <mode> <control>"
            echo "      Run specific tool with parameters"
            echo "      Example: scCTpeak run human macs2 H3K27ac B broad with_input"
            echo ""
            echo "  parse2wig <dataset> <histone> <cell_type>"
            echo "      Run DROMPAplus parse2wig+ preprocessing"
            echo "      Example: scCTpeak parse2wig human H3K27ac B"
            echo ""
            echo "  batch <dataset> <tool> <cell_type> <control>"
            echo "      Batch process all histone marks for a tool"
            echo "      Example: scCTpeak batch human macs2 all with_input"
            echo ""
            echo "  all_tools <dataset> <cell_type> <control>"
            echo "      Run all tools for a dataset"
            echo "      Example: scCTpeak all_tools human B with_input"
            echo ""
            echo "DATASET OPTIONS:"
            echo "---------------"
            echo "  human: H3K27ac, H3K27me3, H3K4me1, H3K4me2, H3K4me3, H3K9me3"
            echo "  mouse: H3K27ac, H3K27me3, H3K36me3, H3K4me3, Olig2, Rad21"
            echo ""
            echo "CELL TYPE OPTIONS:"
            echo "-----------------"
            echo "  human: B, CD4T, CD8T, DC, Mono, NK, otherT, other, all"
            echo "  mouse: Astrocytes, Microglia, mOL, Neurons1, Neurons2, Neurons3, OEC, OPC, VLMC, Unknown, all"
            echo ""
            ;;
        
        "version"|"-v"|"--version")
            echo "scCTpeak v$SCCTPEAK_VERSION"
            echo "A Complete Benchmark Suite for scCUT&Tag Peak-Calling"
            echo "With Data Processing Module & Unified Output Format"
            echo "Author: $SCCTPEAK_AUTHOR"
            ;;
        
        *)
            echo "Error: Unknown command '$command'"
            echo "Use 'scCTpeak help' for usage information"
            return 1
            ;;
    esac
}

# ===========================================================================================
# AUTO-COMPLETE FUNCTIONALITY
# ===========================================================================================
_scCTpeak_complete() {
    local cur prev words cword
    _init_completion || return
    
    case "${cword}" in
        1)
            COMPREPLY=($(compgen -W "process process_all verify split_bam split_bam_all run parse2wig batch all_tools organize summary help version" -- "$cur"))
            ;;
        2)
            case "${words[1]}" in
                process|process_all|verify|split_bam|split_bam_all|run|batch|all_tools|organize|summary)
                    COMPREPLY=($(compgen -W "human mouse" -- "$cur"))
                    ;;
                parse2wig)
                    COMPREPLY=($(compgen -W "human mouse" -- "$cur"))
                    ;;
            esac
            ;;
        3)
            case "${words[1]}" in
                process)
                    if [[ "${words[2]}" == "human" ]]; then
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K4me1 H3K4me2 H3K4me3 H3K9me3" -- "$cur"))
                    else
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K36me3 H3K4me3 Olig2 Rad21" -- "$cur"))
                    fi
                    ;;
                verify)
                    if [[ "${words[2]}" == "human" ]]; then
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K4me1 H3K4me2 H3K4me3 H3K9me3" -- "$cur"))
                    else
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K36me3 H3K4me3 Olig2 Rad21" -- "$cur"))
                    fi
                    ;;
                split_bam)
                    if [[ "${words[2]}" == "human" ]]; then
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K4me1 H3K4me2 H3K4me3 H3K9me3" -- "$cur"))
                    else
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K36me3 H3K4me3 Olig2 Rad21" -- "$cur"))
                    fi
                    ;;
                run)
                    COMPREPLY=($(compgen -W "drompa genrich gopeaks macs2 homer seacr sicer2" -- "$cur"))
                    ;;
                batch)
                    COMPREPLY=($(compgen -W "drompa genrich gopeaks macs2 homer seacr sicer2" -- "$cur"))
                    ;;
                parse2wig)
                    if [[ "${words[2]}" == "human" ]]; then
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K4me1 H3K4me2 H3K4me3 H3K9me3" -- "$cur"))
                    else
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K36me3 H3K4me3 Olig2 Rad21" -- "$cur"))
                    fi
                    ;;
                organize|summary)
                    COMPREPLY=($(compgen -W "with_input without_input" -- "$cur"))
                    ;;
            esac
            ;;
        4)
            case "${words[1]}" in
                verify)
                    if [[ "${words[2]}" == "human" ]]; then
                        COMPREPLY=($(compgen -W "B CD4T CD8T DC Mono NK otherT other" -- "$cur"))
                    else
                        COMPREPLY=($(compgen -W "Astrocytes Microglia mOL Neurons1 Neurons2 Neurons3 OEC OPC VLMC Unknown" -- "$cur"))
                    fi
                    ;;
                run)
                    if [[ "${words[2]}" == "human" ]]; then
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K4me1 H3K4me2 H3K4me3 H3K9me3" -- "$cur"))
                    else
                        COMPREPLY=($(compgen -W "H3K27ac H3K27me3 H3K36me3 H3K4me3 Olig2 Rad21" -- "$cur"))
                    fi
                    ;;
                batch|all_tools)
                    COMPREPLY=($(compgen -W "all B CD4T CD8T DC Mono NK otherT other Astrocytes Microglia mOL Neurons1 Neurons2 Neurons3 OEC OPC VLMC Unknown" -- "$cur"))
                    ;;
            esac
            ;;
    esac
}

# Register auto-completion if available
if declare -F _init_completion >/dev/null 2>&1; then
    complete -F _scCTpeak_complete scCTpeak
fi

# ===========================================================================================
# EXPORT THE FUNCTION
# ===========================================================================================
export -f scCTpeak

# ===========================================================================================
# MAIN EXECUTION (if script is run directly)
# ===========================================================================================
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    scCTpeak "$@"
fi
