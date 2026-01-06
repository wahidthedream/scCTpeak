####################################################################################################################
### Application DROMPAplus on HumanPBMC data
### Flexible function for multiple histones
####################################################################################################################

run_drompa_for_histone() {
    local histone="$1"
    
    echo "####################################################################################################################"
    echo "### Processing $histone"
    echo "####################################################################################################################"
    
    # Create and change to output directory
    mkdir -p ~/project_scHMTF/GSE195725_processed_data/simulation_analysis_50bp/DROMPAplus_peakbed/${histone}
    cd ~/project_scHMTF/GSE195725_processed_data/simulation_analysis_50bp/DROMPAplus_peakbed/${histone}/
    
    # Cell types array
    local cell_types=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")
    # Process treatment BAMs
    for cell in "${cell_types[@]}"; do
        echo "Processing ${histone} treatment - $cell"
        parse2wig+ -i ~/project_scHMTF/GSE195725_processed_data/BAM/celltypeswisebam_50bp/${histone}/${histone}_${cell}.bam \
                   -o ${histone}.${cell} \
                   --pair \
                   --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
                   -n GR
    done
    
    # Process input BAMs  
    for cell in "${cell_types[@]}"; do
        echo "Processing input - $cell"
        parse2wig+ -i ~/project_scHMTF/GSE195725_processed_data/BAM/celltypeswisebam_50bp/${histone}/input_${cell}.bam \
                   -o input.${cell} \
                   --pair \
                   --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
                   -n GR
    done
    
    echo "✓ Completed $histone processing"
    echo ""
}

# Main function to run for all histones
run_drompa_all_histones() {
    local histones=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
    
    for histone in "${histones[@]}"; do
        run_drompa_for_histone "$histone"
    done
    
    echo "####################################################################################################################"
    echo "### ALL HISTONES PROCESSING COMPLETED!"
    echo "####################################################################################################################"
}

# AUTO-RUN WHEN SCRIPT IS EXECUTED
run_drompa_all_histones




#!/bin/bash
####################################################################################################################
### Application DROMPAplus on REAL BAMMouse Brain data
### Flexible function for multiple histones 
####################################################################################################################

run_drompa_for_histone() {
    local histone="$1"
    
    echo "####################################################################################################################"
    echo "### Processing $histone"
    echo "####################################################################################################################"

    local outdir=/home/wahid/project_scHMTF/GSE157637_processed_data/simulation_analysis_50bp/DROMPAplus_peakbed/${histone}
    mkdir -p "$outdir"
    cd "$outdir" || { echo "❌ Failed to enter $outdir"; exit 1; }

    # Define cell types depending on histone
    local cell_types=()
    case "$histone" in
        H3K27ac) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC" "VLMC") ;;
        H3K27me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons3" "OEC" "OPC" "VLMC") ;;
        H3K36me3) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC") ;;
        H3K4me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons2" "Neurons3" "OEC" "OPC" "VLMC") ;;
        Olig2) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
        Rad21) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
        *)
            echo "⚠️  Unknown histone type: $histone"
            return
            ;;
    esac

    # Reference genome file
    local genome_ref=/home/wahid/project_scHMTF/GSE157637_processed_data/ref/mm10.genome_file.txt

    # Process treatment BAMs
    for cell in "${cell_types[@]}"; do
        echo "Processing ${histone} treatment - $cell"
        parse2wig+ -i /home/wahid/project_scHMTF/GSE157637_processed_data/result/celltypeswisebam_50bp/${histone}/${histone}_${cell}.bam \
                   -o ${histone}.${cell} \
                   --pair \
                   --gt "$genome_ref" \
                   -n GR
    done

    # Process input BAMs
    for cell in "${cell_types[@]}"; do
        echo "Processing input - $cell"
        parse2wig+ -i /home/wahid/project_scHMTF/GSE157637_processed_data/result/celltypeswisebam_50bp/${histone}/input_${cell}.bam \
                   -o input.${cell} \
                   --pair \
                   --gt "$genome_ref" \
                   -n GR
    done

    echo "✓ Completed $histone processing"
    echo ""
}

# Run for all histones
run_drompa_all_histones() {
    local histones=("H3K27ac" "H3K27me3" "H3K36me3" "H3K4me3" "Olig2" "Rad21")
    
    for histone in "${histones[@]}"; do
        run_drompa_for_histone "$histone"
    done
    
    echo "####################################################################################################################"
    echo "### ALL HISTONES PROCESSING COMPLETED!"
    echo "####################################################################################################################"
}

# AUTO-RUN WHEN SCRIPT IS EXECUTED
run_drompa_all_histones



#!/bin/bash
###############################################################################################################
### DROMPAplus Peak Calling WITH Control (Input bigWig files)
### Human PBMC real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
SHARP_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base directories
#------------------------------------------
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/simulation_analysis_50bp"
REF_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/ref"
GENOME_FILE="${REF_DIR}/genome_file.txt"
GENE_ANNOT="${REF_DIR}/refFlat.dupremoved.txt"

THREADS=32

###############################################################################################################
# Function to run DROMPAplus
###############################################################################################################
run_drompa() {
    local histone="$1"
    local mode="$2"
    local P_INT P_ENR

    if [[ "$mode" == "BROAD" ]]; then
        P_INT=4
        P_ENR=3
    else
        P_INT=5
        P_ENR=4
    fi

    echo "===================================================================="
    echo "Processing $histone ($mode)"
    echo "===================================================================="

    local HISTONE_DIR="${BASE_DIR}/DROMPAplus_peakbed/${histone}"
    local PARSE_DIR="${HISTONE_DIR}/parse2wigdir+"
    local OUT_DIR="${HISTONE_DIR}/peakbed_with_input_${mode,,}/DROMPAplus"
    mkdir -p "$OUT_DIR"

    local INPUTS=""
    for cell in "${CELLTYPES[@]}"; do
        local label=""

        # H3K27ac special naming
        if [[ "$histone" == "H3K27ac" ]]; then
            if [[ "$mode" == "BROAD" ]]; then
                label="${histone}-b_${cell}"
            else
                label="${histone}-s_${cell}"
            fi
        else
            label="${histone}_${cell}"
        fi

        INPUTS+=" -i ${PARSE_DIR}/${histone}.${cell}.100.bw,${PARSE_DIR}/input.${cell}.100.bw,${label},,,100"
    done

    drompa+ PC_${mode} \
        ${INPUTS} \
        -o "${OUT_DIR}" \
        --gt "${GENOME_FILE}" \
        -g "${GENE_ANNOT}" \
        --lpp 5 --showitag 1 --callpeak \
        --pthre_internal "${P_INT}" \
        --pthre_enrich "${P_ENR}"

    echo "✅ Done: $histone ($mode)"
    echo
}

###############################################################################################################
# Main loop over histones
###############################################################################################################
for histone in "${BROAD_MARKS[@]}" "${SHARP_MARKS[@]}"; do
    if [[ "$histone" == "H3K27ac" ]]; then
        # Run both broad and sharp for H3K27ac
        run_drompa "$histone" "BROAD"
        run_drompa "$histone" "SHARP"
    elif [[ " ${BROAD_MARKS[@]} " =~ " ${histone} " ]]; then
        run_drompa "$histone" "BROAD"
    else
        run_drompa "$histone" "SHARP"
    fi
done


#!/bin/bash
###############################################################################################################
### DROMPAplus Peak Calling WITHOUT Control (No Input bigWig files)
### Human PBMC real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
SHARP_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base directories
#------------------------------------------
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/simulation_analysis_50bp"
REF_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/ref"
GENOME_FILE="${REF_DIR}/genome_file.txt"
GENE_ANNOT="${REF_DIR}/refFlat.dupremoved.txt"

THREADS=32

###############################################################################################################
# Function to run DROMPAplus
###############################################################################################################
run_drompa() {
    local histone="$1"
    local mode="$2"
    local P_INT P_ENR

    if [[ "$mode" == "BROAD" ]]; then
        P_INT=4

    else
        P_INT=5

    fi

    echo "===================================================================="
    echo "Processing $histone ($mode)"
    echo "===================================================================="

    local HISTONE_DIR="${BASE_DIR}/DROMPAplus_peakbed/${histone}"
    local PARSE_DIR="${HISTONE_DIR}/parse2wigdir+"
    local OUT_DIR="${HISTONE_DIR}/peakbed_without_input_${mode,,}/DROMPAplus"
    mkdir -p "$OUT_DIR"

    local INPUTS=""
    for cell in "${CELLTYPES[@]}"; do
        local label=""

        # H3K27ac special naming
        if [[ "$histone" == "H3K27ac" ]]; then
            if [[ "$mode" == "BROAD" ]]; then
                label="${histone}-b_${cell}"
            else
                label="${histone}-s_${cell}"
            fi
        else
            label="${histone}_${cell}"
        fi

        INPUTS+=" -i ${PARSE_DIR}/${histone}.${cell}.100.bw,,${label},,,100"
    done

    drompa+ PC_${mode} \
        ${INPUTS} \
        -o "${OUT_DIR}" \
        --gt "${GENOME_FILE}" \
        -g "${GENE_ANNOT}" \
        --lpp 5 --showitag 1 --callpeak \
        --pthre_internal "${P_INT}" 


    echo "✅ Done: $histone ($mode)"
    echo
}

###############################################################################################################
# Main loop over histones
###############################################################################################################
for histone in "${BROAD_MARKS[@]}" "${SHARP_MARKS[@]}"; do
    if [[ "$histone" == "H3K27ac" ]]; then
        # Run both broad and sharp for H3K27ac
        run_drompa "$histone" "BROAD"
        run_drompa "$histone" "SHARP"
    elif [[ " ${BROAD_MARKS[@]} " =~ " ${histone} " ]]; then
        run_drompa "$histone" "BROAD"
    else
        run_drompa "$histone" "SHARP"
    fi
done


#!/bin/bash
###############################################################################################################
### DROMPAplus Peak Calling WITH Control (Input bigWig files)
### MouseBrain real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
SHARP_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

#------------------------------------------
# Function: define cell types per histone mark
#------------------------------------------
get_celltypes_for_mark() {
  local MARK="$1"
  case "$MARK" in
    H3K27ac) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC" "VLMC") ;;
    H3K27me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons3" "OEC" "OPC" "VLMC") ;;
    H3K36me3) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC") ;;
    H3K4me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons2" "Neurons3" "OEC" "OPC" "VLMC") ;;
    Olig2) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
    Rad21) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
    *) echo "Unknown mark: $MARK"; CELLTYPES=() ;;
  esac
}

#------------------------------------------
# Base directories
#------------------------------------------
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/simulation_analysis_50bp/"
REF_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/ref"
GENOME_FILE="${REF_DIR}/mm10.genome_file.txt"
GENE_ANNOT="${REF_DIR}/refFlat.dupremoved.txt"

THREADS=32

###############################################################################################################
# Function to run DROMPAplus
###############################################################################################################
run_drompa() {
    local histone="$1"
    local mode="$2"
    local P_INT P_ENR

    if [[ "$mode" == "BROAD" ]]; then
        P_INT=4
        P_ENR=3
    else
        P_INT=5
        P_ENR=4
    fi

    echo "===================================================================="
    echo "Processing $histone ($mode)"
    echo "===================================================================="

    local HISTONE_DIR="${BASE_DIR}/DROMPAplus_peakbed/${histone}"
    local PARSE_DIR="${HISTONE_DIR}/parse2wigdir+"
    local OUT_DIR="${HISTONE_DIR}/peakbed_with_input_${mode,,}/DROMPAplus"
    mkdir -p "$OUT_DIR"

    local INPUTS=""
    for cell in "${CELLTYPES[@]}"; do
        local label=""

        # H3K27ac special naming
        if [[ "$histone" == "H3K27ac" ]]; then
            if [[ "$mode" == "BROAD" ]]; then
                label="${histone}-b_${cell}"
            else
                label="${histone}-s_${cell}"
            fi
        else
            label="${histone}_${cell}"
        fi

        INPUTS+=" -i ${PARSE_DIR}/${histone}.${cell}.100.bw,${PARSE_DIR}/input.${cell}.100.bw,${label},,,100"
    done

    drompa+ PC_${mode} \
        ${INPUTS} \
        -o "${OUT_DIR}" \
        --gt "${GENOME_FILE}" \
        -g "${GENE_ANNOT}" \
        --lpp 5 --showitag 1 --callpeak \
        --pthre_internal "${P_INT}" \
        --pthre_enrich "${P_ENR}"

    echo "✅ Done: $histone ($mode)"
    echo
}

###############################################################################################################
# Main loop over histones
###############################################################################################################
for histone in "${BROAD_MARKS[@]}" "${SHARP_MARKS[@]}"; do
    get_celltypes_for_mark "$histone"

    if [[ "$histone" == "H3K27ac" ]]; then
        run_drompa "$histone" "BROAD"
        run_drompa "$histone" "SHARP"
    elif [[ " ${BROAD_MARKS[@]} " =~ " ${histone} " ]]; then
        run_drompa "$histone" "BROAD"
    else
        run_drompa "$histone" "SHARP"
    fi
done

#!/bin/bash
###############################################################################################################
### DROMPAplus Peak Calling WITHOUT Control (No Input bigWig files)
### MouseBrain real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
SHARP_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

#------------------------------------------
# Function: define cell types per histone mark
#------------------------------------------
get_celltypes_for_mark() {
  local MARK="$1"
  case "$MARK" in
    H3K27ac) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC" "VLMC") ;;
    H3K27me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons3" "OEC" "OPC" "VLMC") ;;
    H3K36me3) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC") ;;
    H3K4me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons2" "Neurons3" "OEC" "OPC" "VLMC") ;;
    Olig2|Rad21) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
    *) echo "Unknown mark: $MARK"; CELLTYPES=() ;;
  esac
}

#------------------------------------------
# Base directories
#------------------------------------------
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/simulation_analysis_50bp"
REF_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/ref"
GENOME_FILE="${REF_DIR}/mm10.genome_file.txt"
GENE_ANNOT="${REF_DIR}/refFlat.dupremoved.txt"

THREADS=32

###############################################################################################################
# Function to run DROMPAplus
###############################################################################################################
run_drompa_noctrl() {
    local histone="$1"
    local mode="$2"
    local P_INT

    if [[ "$mode" == "BROAD" ]]; then
        P_INT=4
    else
        P_INT=5
    fi

    echo "===================================================================="
    echo "Processing $histone ($mode) WITHOUT control"
    echo "===================================================================="

    local HISTONE_DIR="${BASE_DIR}/DROMPAplus_peakbed/${histone}"
    local PARSE_DIR="${HISTONE_DIR}/parse2wigdir+"
    local OUT_DIR="${HISTONE_DIR}/peakbed_without_input_${mode,,}/DROMPAplus"
    mkdir -p "$OUT_DIR"

    # Get cell types for this histone
    get_celltypes_for_mark "$histone"

    local INPUTS=""
    for cell in "${CELLTYPES[@]}"; do
        local label=""

        # H3K27ac special naming
        if [[ "$histone" == "H3K27ac" ]]; then
            if [[ "$mode" == "BROAD" ]]; then
                label="${histone}_b_${cell}"
            else
                label="${histone}_s_${cell}"
            fi
        else
            label="${histone}_${cell}"
        fi

        INPUTS+=" -i ${PARSE_DIR}/${histone}.${cell}.100.bw,,${label},,,100"
    done

    drompa+ PC_${mode} \
        ${INPUTS} \
        -o "${OUT_DIR}" \
        --gt "${GENOME_FILE}" \
        -g "${GENE_ANNOT}" \
        --lpp 5 --showitag 1 --callpeak \
        --pthre_internal "${P_INT}"

    echo "✅ Done: $histone ($mode) WITHOUT control"
    echo
}

###############################################################################################################
# Main loop over histones
###############################################################################################################
for histone in "${BROAD_MARKS[@]}" "${SHARP_MARKS[@]}"; do
    if [[ "$histone" == "H3K27ac" ]]; then
        # Run both broad and sharp for H3K27ac
        run_drompa_noctrl "$histone" "BROAD"
        run_drompa_noctrl "$histone" "SHARP"
    elif [[ " ${BROAD_MARKS[@]} " =~ " ${histone} " ]]; then
        run_drompa_noctrl "$histone" "BROAD"
    else
        run_drompa_noctrl "$histone" "SHARP"
    fi
done

