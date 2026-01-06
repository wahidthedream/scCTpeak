#!/bin/bash
###############################################################################################################
### SEACR Peak Calling WITH Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

#------------------------------------------
# Histone mark definitions
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Path configurations
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/SEACR_peakbed_corrected"
SEACR_SCRIPT="/home/wahid/tools/SEACR/SEACR_1.3.sh"

mkdir -p "${OUT_BASE}"

###############################################################################################################
### Function: Run SEACR with Control
###############################################################################################################
run_seacr_with_control() {
    local MARK="$1"
    local MODE="$2"   # "broad" or "narrow"
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    local TMP_DIR="${OUT_DIR}/bedgraphs"

    mkdir -p "$OUT_DIR" "$TMP_DIR"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CONTROL_BAM="${MARK_DIR}/input_${CELL}.bam"

        if [[ ! -f "$TREATMENT_BAM" || ! -f "$CONTROL_BAM" ]]; then
            echo "⚠️ Missing BAM for ${MARK} (${CELL}). Skipping..."
            continue
        fi

        local TREATMENT_BG="${TMP_DIR}/${MARK}_${CELL}_treat.bedgraph"
        local CONTROL_BG="${TMP_DIR}/${MARK}_${CELL}_ctrl.bedgraph"

        echo "========== Converting BAM to BEDGRAPH for: ${MARK} (${CELL}) =========="
        bedtools genomecov -bg -ibam "$TREATMENT_BAM" | sort -k1,1 -k2,2n > "$TREATMENT_BG"
        bedtools genomecov -bg -ibam "$CONTROL_BAM"   | sort -k1,1 -k2,2n > "$CONTROL_BG"

        echo "========== Running SEACR for: ${MARK} (${CELL}) =========="
        local OUTPUT_PREFIX="${OUT_DIR}/${CELL}_${MARK}"

        # Broad marks → stringent threshold, Narrow → relaxed
        if [[ "$MODE" == "broad" ]]; then
            THRESHOLD="stringent"
        else
            THRESHOLD="relaxed"
        fi

        # ✅ Correct SEACR command
        bash "$SEACR_SCRIPT" "$TREATMENT_BG" "$CONTROL_BG" norm "$THRESHOLD" "$OUTPUT_PREFIX"

        echo "✅ Finished: ${MARK}_${CELL} (SEACR ${THRESHOLD})"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SEACR Peak Summary for ${MARK} (${MODE}) ========="
    for CELL in "${CELLTYPES[@]}"; do
        local BED_FILE=$(find "$OUT_DIR" -type f -name "${CELL}_${MARK}*.bed" | head -n 1)
        if [[ -f "$BED_FILE" ]]; then
            local COUNT=$(wc -l < "$BED_FILE")
            echo "${MARK}_${CELL}: ${COUNT} peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

###############################################################################################################
### Run SEACR for Broad and Narrow Marks
###############################################################################################################
for MARK in "${BROAD_MARKS[@]}"; do
    run_seacr_with_control "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_seacr_with_control "$MARK" "narrow"
done




#!/bin/bash
###############################################################################################################
### SEACR Peak Calling WITHOUT Control ( No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

#------------------------------------------
# Histone mark definitions
#------------------------------------------
BROAD_MARKS=("H3K27ac" )
NARROW_MARKS=("H3K27ac" )
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Path configurations
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/SEACR_peakbed_corrected_evaluation"
SEACR_SCRIPT="/home/wahid/tools/SEACR/SEACR_1.3.sh"

mkdir -p "${OUT_BASE}"

###############################################################################################################
### Function: Run SEACR WITHOUT Control
###############################################################################################################
run_seacr_without_control() {
    local MARK="$1"
    local MODE="$2"   # "broad" or "narrow"
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected_evaluation"
    local TMP_DIR="${OUT_DIR}/bedgraphs"

    mkdir -p "$OUT_DIR" "$TMP_DIR"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        if [[ ! -f "$TREATMENT_BAM" ]]; then
            echo "⚠️ Missing BAM for ${MARK} (${CELL}). Skipping..."
            continue
        fi

        local TREATMENT_BG="${TMP_DIR}/${MARK}_${CELL}_treat.bedgraph"

        echo "========== Converting BAM to BEDGRAPH for: ${MARK} (${CELL}) =========="
        bedtools genomecov -bg -ibam "$TREATMENT_BAM" | sort -k1,1 -k2,2n > "$TREATMENT_BG"

        echo "========== Running SEACR for: ${MARK} (${CELL}) =========="
        local OUTPUT_PREFIX="${OUT_DIR}/${CELL}_${MARK}"

        # Broad marks → stringent threshold, Narrow → relaxed
        if [[ "$MODE" == "broad" ]]; then
            THRESHOLD="stringent"
        else
            THRESHOLD="relaxed"
        fi

        # ✅ Correct SEACR command for no control (use "non")
        bash "$SEACR_SCRIPT" "$TREATMENT_BG" "0.01" "non" "$THRESHOLD" "$OUTPUT_PREFIX"

        echo "✅ Finished: ${MARK}_${CELL} (SEACR ${THRESHOLD})"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SEACR Peak Summary for ${MARK} (${MODE}) ========="
    for CELL in "${CELLTYPES[@]}"; do
        local BED_FILE=$(find "$OUT_DIR" -type f -name "${CELL}_${MARK}*.bed" | head -n 1)
        if [[ -f "$BED_FILE" ]]; then
            local COUNT=$(wc -l < "$BED_FILE")
            echo "${MARK}_${CELL}: ${COUNT} peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

###############################################################################################################
### Run SEACR for Broad and Narrow Marks
###############################################################################################################
for MARK in "${BROAD_MARKS[@]}"; do
    run_seacr_without_control "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_seacr_without_control "$MARK" "narrow"
done



#!/bin/bash
###############################################################################################################
### SEACR Peak Calling WITH Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Mouse Brain real BAMs
###############################################################################################################

#------------------------------------------
# Histone mark definitions
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

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
    *)
      echo "⚠️ Unknown mark: $MARK"
      CELLTYPES=()
      ;;
  esac
}

#------------------------------------------
# Path configurations
#------------------------------------------
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/SEACR_peakbed"
SEACR_SCRIPT="/home/wahid/tools/SEACR/SEACR_1.3.sh"

mkdir -p "${OUT_BASE}"

###############################################################################################################
### Function: Run SEACR WITH Control
###############################################################################################################
run_seacr_with_control() {
    local MARK="$1"
    local MODE="$2"   # "broad" or "narrow"
    get_celltypes_for_mark "$MARK"

    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    local TMP_DIR="${OUT_DIR}/bedgraphs"

    mkdir -p "$OUT_DIR"
    rm -rf "$TMP_DIR" && mkdir -p "$TMP_DIR"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CONTROL_BAM="${MARK_DIR}/input_${CELL}.bam"

        if [[ ! -f "$TREATMENT_BAM" || ! -f "$CONTROL_BAM" ]]; then
            echo "⚠️ Missing BAM for ${MARK} (${CELL}). Skipping..."
            continue
        fi

        local TREATMENT_BG="${TMP_DIR}/${MARK}_${CELL}_treat.bedgraph"
        local CONTROL_BG="${TMP_DIR}/${MARK}_${CELL}_ctrl.bedgraph"

        echo "========== Converting BAM → BEDGRAPH for: ${MARK} (${CELL}) =========="
        bedtools genomecov -bg -ibam "$TREATMENT_BAM" | LC_ALL=C sort -k1,1 -k2,2n > "$TREATMENT_BG"
        bedtools genomecov -bg -ibam "$CONTROL_BAM"   | LC_ALL=C sort -k1,1 -k2,2n > "$CONTROL_BG"

        echo "========== Running SEACR for: ${MARK} (${CELL}) =========="
        local OUTPUT_PREFIX="${OUT_DIR}/${CELL}_${MARK}"

        # Broad marks → stringent threshold, Narrow → relaxed
        if [[ "$MODE" == "broad" ]]; then
            THRESHOLD="stringent"
        else
            THRESHOLD="relaxed"
        fi

        # ✅ Correct SEACR syntax for WITH control
        bash "$SEACR_SCRIPT" "$TREATMENT_BG" "$CONTROL_BG" norm "$THRESHOLD" "$OUTPUT_PREFIX"

        echo "✅ Finished: ${MARK}_${CELL} (SEACR ${THRESHOLD})"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SEACR Peak Summary for ${MARK} (${MODE}) ========="
    for CELL in "${CELLTYPES[@]}"; do
        local BED_FILE
        BED_FILE=$(find "$OUT_DIR" -type f -name "${CELL}_${MARK}*.bed" | head -n 1)
        if [[ -f "$BED_FILE" ]]; then
            local COUNT
            COUNT=$(wc -l < "$BED_FILE")
            echo "${MARK}_${CELL}: ${COUNT} peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

###############################################################################################################
### Run SEACR for Broad and Narrow Marks
###############################################################################################################
for MARK in "${BROAD_MARKS[@]}"; do
    run_seacr_with_control "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_seacr_with_control "$MARK" "narrow"
done




#!/bin/bash
###############################################################################################################
### SEACR Peak Calling WITHOUT Control (No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Mouse Brain real BAMs
###############################################################################################################

#------------------------------------------
# Histone mark definitions
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

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
    *)
      echo "⚠️ Unknown mark: $MARK"
      CELLTYPES=()
      ;;
  esac
}

#------------------------------------------
# Path configurations
#------------------------------------------
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/SEACR_peakbed"
SEACR_SCRIPT="/home/wahid/tools/SEACR/SEACR_1.3.sh"

mkdir -p "${OUT_BASE}"

###############################################################################################################
### Function: Run SEACR WITHOUT Control
###############################################################################################################
run_seacr_without_control() {
    local MARK="$1"
    local MODE="$2"   # "broad" or "narrow"
    get_celltypes_for_mark "$MARK"

    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
    local TMP_DIR="${OUT_DIR}/bedgraphs"

    mkdir -p "$OUT_DIR"
    rm -rf "$TMP_DIR" && mkdir -p "$TMP_DIR"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"

        if [[ ! -f "$TREATMENT_BAM" ]]; then
            echo "⚠️ Missing BAM for ${MARK} (${CELL}). Skipping..."
            continue
        fi

        local TREATMENT_BG="${TMP_DIR}/${MARK}_${CELL}_treat.bedgraph"

        echo "========== Converting BAM → BEDGRAPH for: ${MARK} (${CELL}) =========="
        bedtools genomecov -bg -ibam "$TREATMENT_BAM" | LC_ALL=C sort -k1,1 -k2,2n > "$TREATMENT_BG"

        echo "========== Running SEACR for: ${MARK} (${CELL}) =========="
        local OUTPUT_PREFIX="${OUT_DIR}/${CELL}_${MARK}"

        # Broad marks → stringent threshold, Narrow → relaxed
        if [[ "$MODE" == "broad" ]]; then
            THRESHOLD="stringent"
        else
            THRESHOLD="relaxed"
        fi

        # ✅ Correct SEACR syntax for WITHOUT control
        bash "$SEACR_SCRIPT" "$TREATMENT_BG" "0.01" "non" "$THRESHOLD" "$OUTPUT_PREFIX"

        echo "✅ Finished: ${MARK}_${CELL} (SEACR ${THRESHOLD})"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SEACR Peak Summary for ${MARK} (${MODE}) ========="
    for CELL in "${CELLTYPES[@]}"; do
        local BED_FILE
        BED_FILE=$(find "$OUT_DIR" -type f -name "${CELL}_${MARK}*.bed" | head -n 1)
        if [[ -f "$BED_FILE" ]]; then
            local COUNT
            COUNT=$(wc -l < "$BED_FILE")
            echo "${MARK}_${CELL}: ${COUNT} peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

###############################################################################################################
### Run SEACR for Broad and Narrow Marks
###############################################################################################################
for MARK in "${BROAD_MARKS[@]}"; do
    run_seacr_without_control "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_seacr_without_control "$MARK" "narrow"
done

