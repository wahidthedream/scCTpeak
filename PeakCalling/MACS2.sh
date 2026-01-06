


#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling with Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

# Define cell types
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Base directories
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/MACS2_peakbed"
GENOME=2.7e9  # Human genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark
# -------------------------------------------------------------------------------------------
run_macs2_all_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "${OUT_DIR}"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CTRL="${MARK_DIR}/input_${CELL}.bam"

        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" || ! -f "$CTRL" ]]; then
            echo "⚠️ Missing treatment or control BAM for ${CELL}. Skipping..."
            continue
        fi

        # Set parameters
        local QVAL=0.05
        local ARGS=" -q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -c "$CTRL" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "✅ Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "narrow"
done


#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling (Without Control)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

# Define cell types
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/MACS2_peakbed"
GENOME=2.7e9  # Human genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark (no control)
# -------------------------------------------------------------------------------------------
run_macs2_no_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
    mkdir -p "${OUT_DIR}"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"

        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" ]]; then
            echo "⚠️ Missing treatment BAM for ${CELL}. Skipping..."
            continue
        fi

        # Parameters
        local QVAL=0.05
        local ARGS="-q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "✅ Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "narrow"
done


#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling with Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

# Function: define cell types per histone mark
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

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/MACS2_peakbed_corrected"
GENOME=1.87e9  # Mouse genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark
# -------------------------------------------------------------------------------------------
run_macs2_all_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/With_input_corrected"
    mkdir -p "${OUT_DIR}"

    get_celltypes_for_mark "$MARK"
    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CTRL="${MARK_DIR}/input_${CELL}.bam"

        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" || ! -f "$CTRL" ]]; then
            echo "⚠️ Missing treatment or control BAM for ${CELL}. Skipping..."
            continue
        fi

        # Set parameters
        local QVAL=0.05
        local ARGS="-q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -c "$CTRL" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "✅ Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "broad"
done 

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "narrow"
done



#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling WITHOUT Control (NO Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

# Function: define cell types per histone mark
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
      echo " Unknown mark: $MARK"
      CELLTYPES=()
      ;;
  esac
}

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/MACS2_peakbed_corrected"
GENOME=1.87e9  # Mouse genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark (without control)
# -------------------------------------------------------------------------------------------
run_macs2_no_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/Without_input_corrected"
    mkdir -p "${OUT_DIR}"

    get_celltypes_for_mark "$MARK"
    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"
        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" ]]; then
            echo "Missing treatment BAM for ${CELL}. Skipping..."
            continue
        fi

        # Set parameters
        local QVAL=0.05
        local ARGS="-q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "✅ Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "narrow"
done



