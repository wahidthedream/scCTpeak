###########################################################################################################################
#### calculation peak count, peak width, and frip 
###########################################################################################################################


#!/bin/bash
#=========================================================================================================================#
# Calculation of Number of Peaks and Peak Widths for HumanPBMC BED files (with input)
# Works for multiple histones and transcription factors
#=========================================================================================================================#

# Define Histones/TFs and Methods
Histones=("H3K27ac-b" "H3K27ac-s" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
Methods=("DROMPAplus" "Genrich" "GoPeaks" "HOMER" "MACS2" "SEACR" "SICER2")

# Input directory containing BED files
PEAK_DIR="~/HumanPBMC_peakbed/all_with_input_peakbed"

# Output CSV file
output_file="~/result_frip/HumanPBMC_peak_summary_with_input.csv"

# Write header
echo "Histone/TF,Method,Sample,Number_of_Peaks,Mean_Peak_Width" > "$output_file"

# Loop through each Histone/TF
for histone in "${Histones[@]}"; do
    # Loop through each Method
    for method in "${Methods[@]}"; do
        # Find all BED files matching the method and histone
        for bed in "$PEAK_DIR"/${method}_${histone}*.bed; do
            # Skip if file does not exist
            [[ -f "$bed" ]] || continue
            
            # Extract sample name
            sample=$(basename "$bed" .bed)
            
            # Calculate number of peaks
            num_peaks=$(wc -l < "$bed")
            
            # Calculate mean peak width
            mean_width=$(awk '{sum += $3 - $2} END {if (NR>0) print int(sum / NR); else print 0}' "$bed")
            
            # Append to output CSV
            echo "${histone},${method},${sample},${num_peaks},${mean_width}" >> "$output_file"
        done
    done
done

echo "Peak summary calculation completed. Output saved to $output_file"

#!/bin/bash
#=========================================================================================================================#
# Calculation of Number of Peaks and Peak Widths for Mouse Brain BED files (with input)
# Works for multiple histones and transcription factors
#=========================================================================================================================#

# Define Histones/TFs and Methods
Histones=("H3K27ac-b" "H3K27ac-s" "H3K27me3" "H3K36me3" "H3K4me3" "Olig2" "Rad21")
Methods=("DROMPAplus" "Genrich" "GoPeaks" "HOMER" "MACS2" "SEACR" "SICER2")

# Input directory containing BED files
PEAK_DIR="~/MouseBrain_peakbed/all_with_input_peakbed_corrected"

# Output CSV file
output_file="~/result_frip/MouseBrain_peak_summary_with_input.csv"

# Write header
echo "Histone/TF,Method,Sample,Number_of_Peaks,Mean_Peak_Width" > "$output_file"

# Loop through each Histone/TF
for histone in "${Histones[@]}"; do
    # Loop through each Method
    for method in "${Methods[@]}"; do
        # Find all BED files matching the method and histone
        for bed in "$PEAK_DIR"/${method}_${histone}*.bed; do
            # Skip if file does not exist
            [[ -f "$bed" ]] || continue
            
            # Extract sample name
            sample=$(basename "$bed" .bed)
            
            # Calculate number of peaks
            num_peaks=$(wc -l < "$bed")
            
            # Calculate mean peak width
            mean_width=$(awk '{sum += $3 - $2} END {if (NR>0) print int(sum / NR); else print 0}' "$bed")
            
            # Append to output CSV
            echo "${histone},${method},${sample},${num_peaks},${mean_width}" >> "$output_file"
        done
    done
done

echo "Peak summary calculation completed. Output saved to $output_file"


#!/bin/bash
#=========================================================================================================================#
##     Calculation FRiPs for Real Bam data for HumanPBMC_with input
#=========================================================================================================================#

# Define arrays (must be of equal length)
HM=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
Histone=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
Method=("DROMPAplus" "Genrich" "GoPeaks" "HOMER" "MACS2" "SEACR" "SICER2")

BAM_BASE="~/project_scCTpeak/BAM"
PEAK_BASE="~/project_scCTpeak/HumanPBMC_peakbed"

# Output file
OUTPUT_FILE="~/result_frip/HumanPBMC_frip_with_input.csv"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Write CSV header
echo "HM,Histone,Method,Sample,Total_Reads,Reads_in_Peaks,FRiP" > "$OUTPUT_FILE"

# Iterate and append results
for i in "${!HM[@]}"; do
    hm=${HM[$i]}
    histone=${Histone[$i]}
    
    BAM_DIR="${BAM_BASE}/${hm}/split_celltype_bams"
    
    # Check if BAM directory exists
    if [[ ! -d "$BAM_DIR" ]]; then
        echo "Warning: BAM directory $BAM_DIR does not exist"
        continue
    fi
    
    for method in "${Method[@]}"; do
        PEAK_DIR="${PEAK_BASE}/all_with_input_peakbed"
        
        # Find all BAM files for this histone type
        for bam in "$BAM_DIR"/${histone}_*.bam; do
            # Skip if no files found
            [[ ! -f "$bam" ]] && continue
            
            base=$(basename "$bam" .bam)
            
            # Extract celltype from filename (remove the histone prefix and underscore)
            celltype="${base#${histone}_}"
            
            # For H3K27ac, process both -b and -s BED files
            if [[ "$histone" == "H3K27ac" ]]; then
                # Process H3K27ac-b BED file
                peak_b="${PEAK_DIR}/${method}_${histone}-b_${celltype}.bed"
                if [[ -s "$peak_b" ]]; then
                    total_reads=$(samtools view -c -F 4 "$bam" 2>/dev/null)
                    reads_in_peaks=$(bedtools intersect -u -a "$bam" -b "$peak_b" 2>/dev/null | samtools view -c 2>/dev/null)
                    
                    if [[ $total_reads -gt 0 ]]; then
                        frip=$(echo "scale=6; $reads_in_peaks / $total_reads" | bc 2>/dev/null)
                        [[ -z "$frip" ]] && frip=0
                    else
                        frip=0
                    fi
                    echo "${hm},${histone},${method},${histone}-b_${celltype},${total_reads},${reads_in_peaks},${frip}" >> "$OUTPUT_FILE"
                else
                    echo "${hm},${histone},${method},${histone}-b_${celltype},NA,NA,0" >> "$OUTPUT_FILE"
                fi
                
                # Process H3K27ac-s BED file
                peak_s="${PEAK_DIR}/${method}_${histone}-s_${celltype}.bed"
                if [[ -s "$peak_s" ]]; then
                    total_reads=$(samtools view -c -F 4 "$bam" 2>/dev/null)
                    reads_in_peaks=$(bedtools intersect -u -a "$bam" -b "$peak_s" 2>/dev/null | samtools view -c 2>/dev/null)
                    
                    if [[ $total_reads -gt 0 ]]; then
                        frip=$(echo "scale=6; $reads_in_peaks / $total_reads" | bc 2>/dev/null)
                        [[ -z "$frip" ]] && frip=0
                    else
                        frip=0
                    fi
                    echo "${hm},${histone},${method},${histone}-s_${celltype},${total_reads},${reads_in_peaks},${frip}" >> "$OUTPUT_FILE"
                else
                    echo "${hm},${histone},${method},${histone}-s_${celltype},NA,NA,0" >> "$OUTPUT_FILE"
                fi
                
            else
                # For other histones (H3K27me3, H3K36me3, H3K4me3, Olig2, Rad21)
                peak="${PEAK_DIR}/${method}_${histone}_${celltype}.bed"
                
                if [[ ! -s "$peak" ]]; then
                    echo "${hm},${histone},${method},${base},NA,NA,0" >> "$OUTPUT_FILE"
                    continue
                fi
                
                total_reads=$(samtools view -c -F 4 "$bam" 2>/dev/null)
                reads_in_peaks=$(bedtools intersect -u -a "$bam" -b "$peak" 2>/dev/null | samtools view -c 2>/dev/null)
                
                if [[ $total_reads -gt 0 ]]; then
                    frip=$(echo "scale=6; $reads_in_peaks / $total_reads" | bc 2>/dev/null)
                    [[ -z "$frip" ]] && frip=0
                else
                    frip=0
                fi
                
                echo "${hm},${histone},${method},${base},${total_reads},${reads_in_peaks},${frip}" >> "$OUTPUT_FILE"
            fi
        done
    done
done

echo "FRiP calculation completed. Results saved to: $OUTPUT_FILE"


#!/bin/bash
#=========================================================================================================================#
##     Calculation FRiPs for Real Bam data for MouseBrain_with input
#=========================================================================================================================#

# Define arrays (must be of equal length)
HM=("H3K27ac" "H3K27me3" "H3K36me3" "H3K4me3" "Olig2" "Rad21")
Histone=("H3K27ac" "H3K27me3" "H3K36me3" "H3K4me3" "Olig2" "Rad21")
Method=("DROMPAplus" "Genrich" "GoPeaks" "HOMER" "MACS2" "SEACR" "SICER2")

BAM_BASE="~/project_scCTpeak/BAM"
PEAK_BASE="~/project_scCTpeak/MouseBrain_peakbed"

# Output file
OUTPUT_FILE="~/result_frip/MouseBrain_frip_with_input.csv"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Write CSV header
echo "HM,Histone,Method,Sample,Total_Reads,Reads_in_Peaks,FRiP" > "$OUTPUT_FILE"

# Iterate and append results
for i in "${!HM[@]}"; do
    hm=${HM[$i]}
    histone=${Histone[$i]}
    
    BAM_DIR="${BAM_BASE}/${hm}/split_celltype_bams"
    
    # Check if BAM directory exists
    if [[ ! -d "$BAM_DIR" ]]; then
        echo "Warning: BAM directory $BAM_DIR does not exist"
        continue
    fi
    
    for method in "${Method[@]}"; do
        PEAK_DIR="${PEAK_BASE}/all_with_input_peakbed_corrected"
        
        # Find all BAM files for this histone type
        for bam in "$BAM_DIR"/${histone}_*.bam; do
            # Skip if no files found
            [[ ! -f "$bam" ]] && continue
            
            base=$(basename "$bam" .bam)
            
            # Extract celltype from filename (remove the histone prefix and underscore)
            celltype="${base#${histone}_}"
            
            # For H3K27ac, process both -b and -s BED files
            if [[ "$histone" == "H3K27ac" ]]; then
                # Process H3K27ac-b BED file
                peak_b="${PEAK_DIR}/${method}_${histone}-b_${celltype}.bed"
                if [[ -s "$peak_b" ]]; then
                    total_reads=$(samtools view -c -F 4 "$bam" 2>/dev/null)
                    reads_in_peaks=$(bedtools intersect -u -a "$bam" -b "$peak_b" 2>/dev/null | samtools view -c 2>/dev/null)
                    
                    if [[ $total_reads -gt 0 ]]; then
                        frip=$(echo "scale=6; $reads_in_peaks / $total_reads" | bc 2>/dev/null)
                        [[ -z "$frip" ]] && frip=0
                    else
                        frip=0
                    fi
                    echo "${hm},${histone},${method},${histone}-b_${celltype},${total_reads},${reads_in_peaks},${frip}" >> "$OUTPUT_FILE"
                else
                    echo "${hm},${histone},${method},${histone}-b_${celltype},NA,NA,0" >> "$OUTPUT_FILE"
                fi
                
                # Process H3K27ac-s BED file
                peak_s="${PEAK_DIR}/${method}_${histone}-s_${celltype}.bed"
                if [[ -s "$peak_s" ]]; then
                    total_reads=$(samtools view -c -F 4 "$bam" 2>/dev/null)
                    reads_in_peaks=$(bedtools intersect -u -a "$bam" -b "$peak_s" 2>/dev/null | samtools view -c 2>/dev/null)
                    
                    if [[ $total_reads -gt 0 ]]; then
                        frip=$(echo "scale=6; $reads_in_peaks / $total_reads" | bc 2>/dev/null)
                        [[ -z "$frip" ]] && frip=0
                    else
                        frip=0
                    fi
                    echo "${hm},${histone},${method},${histone}-s_${celltype},${total_reads},${reads_in_peaks},${frip}" >> "$OUTPUT_FILE"
                else
                    echo "${hm},${histone},${method},${histone}-s_${celltype},NA,NA,0" >> "$OUTPUT_FILE"
                fi
                
            else
                # For other histones (H3K27me3, H3K36me3, H3K4me3, Olig2, Rad21)
                peak="${PEAK_DIR}/${method}_${histone}_${celltype}.bed"
                
                if [[ ! -s "$peak" ]]; then
                    echo "${hm},${histone},${method},${base},NA,NA,0" >> "$OUTPUT_FILE"
                    continue
                fi
                
                total_reads=$(samtools view -c -F 4 "$bam" 2>/dev/null)
                reads_in_peaks=$(bedtools intersect -u -a "$bam" -b "$peak" 2>/dev/null | samtools view -c 2>/dev/null)
                
                if [[ $total_reads -gt 0 ]]; then
                    frip=$(echo "scale=6; $reads_in_peaks / $total_reads" | bc 2>/dev/null)
                    [[ -z "$frip" ]] && frip=0
                else
                    frip=0
                fi
                
                echo "${hm},${histone},${method},${base},${total_reads},${reads_in_peaks},${frip}" >> "$OUTPUT_FILE"
            fi
        done
    done
done

echo "FRiP calculation completed. Results saved to: $OUTPUT_FILE"
