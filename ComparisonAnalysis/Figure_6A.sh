##############################################################################################
## Figure 6 A 
##############################################################################################



###============================================================================================================================
### ChromHMM Part B
### Using Roadmap data for all blood cell
### Without Input Case
### Tools wised output 
### **bash script
###============================================================================================================================

#!/bin/bash
# Define base paths
BASE_PEAK_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/all_without_input_peakbed_corrected/chromhmm_peakbed"
OUTPUT_BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_ChromHMM_roadmapdata/intersect_ChromHMM_without_input"
CHROMHMM_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/ChromHMM_data/all_hg38lift.mnemonics.bedFiles"

# All EIDs for all cell types
EID=("E062" "E034" "E045" "E033" "E044" "E043" "E039" "E040" "E037" "E048" "E038" "E047" "E029" "E032" "E031" "E046" "E030")

tools=("DROMPAplus" "Genrich" "GoPeaks" "HOMER" "MACS2" "SEACR" "SICER2")
# Define all histone modifications
histones=("H3K27ac-b" "H3K27ac-s" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")

# Define all cell types
celltypes=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Log file setup
LOG_FILE="${OUTPUT_BASE_DIR}/intersection_analysis.log"
mkdir -p "$OUTPUT_BASE_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

# Function to log messages with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Check if base directories exist
if [ ! -d "$BASE_PEAK_DIR" ]; then
    log_message "ERROR: Base peak directory does not exist: $BASE_PEAK_DIR"
    exit 1
fi

if [ ! -d "$CHROMHMM_DIR" ]; then
    log_message "ERROR: ChromHMM directory does not exist: $CHROMHMM_DIR"
    exit 1
fi

# Check if bedtools is available
if ! command -v bedtools &> /dev/null; then
    log_message "ERROR: bedtools is not installed or not in PATH"
    exit 1
fi

log_message "Starting ChromHMM intersection analysis - ALL EIDs for ALL cell types"

# Initialize counters
total_processed=0
total_success=0
total_errors=0

# Loop through each histone modification
for hm in "${histones[@]}"; do
    log_message "=========================================="
    log_message "Processing histone modification: $hm"
    log_message "=========================================="
    
    # Set specific directories for this histone modification
    PEAK_DIR="$BASE_PEAK_DIR/${hm}_peakbed"
    OUTPUT_DIR="$OUTPUT_BASE_DIR/$hm"
    
    # Check if the peak directory exists
    if [ ! -d "$PEAK_DIR" ]; then
        log_message "WARNING: Peak directory does not exist: $PEAK_DIR, skipping histone $hm"
        continue
    fi
    
    # Create output directory for this histone modification
    mkdir -p "$OUTPUT_DIR"
    
    # Loop through each cell type
    for cell_type in "${celltypes[@]}"; do
        log_message "Processing cell type: $cell_type"
        
        # Loop through each EID for this cell type
        for eid in "${EID[@]}"; do
            log_message "Processing EID: $eid"
            
            # Define ChromHMM file for this EID
            CHROMHMM_FILE="${CHROMHMM_DIR}/${eid}_15_coreMarks_hg38lift_mnemonics.bed"
            
            # Check if ChromHMM file exists
            if [ ! -f "$CHROMHMM_FILE" ]; then
                log_message "WARNING: ChromHMM file not found: $CHROMHMM_FILE, skipping..."
                ((total_errors++))
                continue
            fi
            
            # Loop through each tool
            for tool in "${tools[@]}"; do
                log_message "Processing tool: $tool"
                
                # Define input peak file path using the histone-specific directory
                PEAK_FILE="${PEAK_DIR}/${tool}_${hm}_${cell_type}.bed"
                
                # Check if peak file exists
                if [ ! -f "$PEAK_FILE" ]; then
                    log_message "WARNING: Peak file not found: $PEAK_FILE, skipping..."
                    ((total_errors++))
                    continue
                fi
                
                # Create output directory with organized structure
                OUTPUT_SUB_DIR="${OUTPUT_DIR}/${tool}"
                mkdir -p "$OUTPUT_SUB_DIR"
                
                # Define output file
                OUTPUT_FILE="${OUTPUT_SUB_DIR}/${tool}_${cell_type}_${hm}_${eid}_intersect.bed"
                
                # Perform intersection using bedtools with 20% overlap
                log_message "Intersecting: $(basename "$PEAK_FILE") with $(basename "$CHROMHMM_FILE")"
                
                if bedtools intersect -a "$PEAK_FILE" -b "$CHROMHMM_FILE" -wo -f 0.20 > "$OUTPUT_FILE"; then
                    intersection_count=$(wc -l < "$OUTPUT_FILE" 2>/dev/null || echo "0")
                    log_message "SUCCESS: $intersection_count intersections -> $(basename "$OUTPUT_FILE")"
                    ((total_success++))
                else
                    log_message "ERROR: Intersection failed for $PEAK_FILE with $CHROMHMM_FILE"
                    # Remove empty output file if creation failed
                    [[ -f "$OUTPUT_FILE" && ! -s "$OUTPUT_FILE" ]] && rm -f "$OUTPUT_FILE"
                    ((total_errors++))
                fi
                
                ((total_processed++))
            done
        done
    done
    
    log_message "Completed processing $hm"
    log_message ""
done

# Final summary
log_message "=========================================="
log_message "ANALYSIS COMPLETED"
log_message "Total combinations attempted: $total_processed"
log_message "Successful intersections: $total_success"
log_message "Errors/Skipped: $total_errors"
log_message "Success rate: $((total_success * 100 / total_processed))%"
log_message "=========================================="

