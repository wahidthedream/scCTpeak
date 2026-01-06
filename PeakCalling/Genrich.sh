
###############################################################################################################
### Genrich peak calling for multiple histone marks (broad and narrow)
### Supports dual processing for H3K27ac (both broad and narrow)
### Real Bam: HumanPMC data
###############################################################################################################

#!/bin/bash
#####=========================================================================================================####
###
### for with input
#####=========================================================================================================####
# Define histone mark groups
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "other" "otherT")
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/Genrich_peakbed"

###############################################################################################################
run_genrich_all_ctrl() {
  local MARK="$1"
  local MARK_TYPE="$2" # "broad" or "narrow"
  local SORT_TAG="qname_sorted_${MARK_TYPE:0:1}"  # b for broad, s for narrow

  # If same mark has two versions (e.g., H3K27ac_broad and H3K27ac_narrow)
  local MARK_LABEL="${MARK}_${MARK_TYPE}"

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
  local SORTED_DIR="${MARK_DIR}/${SORT_TAG}"
  local OUT_DIR="${OUT_BASE}/${MARK_LABEL}/With_input"

  mkdir -p "$OUT_DIR" "$SORTED_DIR" || {
    echo "❌ Failed to create output directories for $MARK_LABEL"
    return 1
  }

  # Adjust Genrich parameters based on mark type
  if [[ "$MARK_TYPE" == "broad" ]]; then
    GENRICH_ARGS="-a 100 -l 500 -g 1000 -p 0.05 -f BAM"
    EXT=".broadPeak"
  else
    GENRICH_ARGS="-a 200 -l 100 -g 100 -p 0.01 -f BAM"
    EXT=".narrowPeak"
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo -e "\n========== Processing $MARK_LABEL ($CELL) =========="

    local RAW_TREATMENT="${MARK_DIR}/${MARK}_${CELL}.bam"
    local RAW_CONTROL="${MARK_DIR}/input_${CELL}.bam"
    local TREATMENT="${SORTED_DIR}/${MARK}_${CELL}_qnamesorted.bam"
    local CONTROL="${SORTED_DIR}/input_${CELL}_qnamesorted.bam"
    local OUTPUT_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_peaks${EXT}"
    local LOG_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_genrich.log"

    if [[ ! -f "$RAW_TREATMENT" ]]; then
      echo "❌ Missing treatment BAM: $RAW_TREATMENT"
      continue
    fi
    if [[ ! -f "$RAW_CONTROL" ]]; then
      echo "❌ Missing control BAM: $RAW_CONTROL"
      continue
    fi

    if [[ ! -f "$TREATMENT" ]]; then
      echo "🔹 Sorting treatment BAM..."
      samtools sort -n -@ $THREADS -o "$TREATMENT" "$RAW_TREATMENT" >> "$LOG_FILE" 2>&1
    fi
    if [[ ! -f "$CONTROL" ]]; then
      echo "🔹 Sorting control BAM..."
      samtools sort -n -@ $THREADS -o "$CONTROL" "$RAW_CONTROL" >> "$LOG_FILE" 2>&1
    fi

    echo "🔹 Running Genrich peak calling..."
    Genrich -t "$TREATMENT" -c "$CONTROL" -o "$OUTPUT_FILE" -j -r -v $GENRICH_ARGS >> "$LOG_FILE" 2>&1

    echo "✅ Completed: $CELL ($MARK_LABEL)"
  done

  echo -e "\n========= Peak Summary for $MARK_LABEL ========="
  for peakfile in "${OUT_DIR}"/*${EXT}; do
    [[ -f "$peakfile" ]] && echo "$(basename "$peakfile"): $(wc -l < "$peakfile") peaks"
  done
}

#####=========================================================================================================####
# Run for all marks
for MARK in "${BROAD_MARKS[@]}"; do
  run_genrich_all_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  run_genrich_all_ctrl "$MARK" "narrow"
done
###############################################################################################################


#!/bin/bash
#####=========================================================================================================####
### Without input
#####=========================================================================================================####

# Define histone mark groups
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "other" "otherT")
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/Genrich_peakbed"

#####=========================================================================================================####
run_genrich_all_noctrl() {
  local MARK="$1"
  local MARK_TYPE="$2"  # "broad" or "narrow"
  local SORT_TAG="qname_sorted_${MARK_TYPE:0:1}"  # b for broad, s for narrow
  local MARK_LABEL="${MARK}_${MARK_TYPE}"

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
  local SORTED_DIR="${MARK_DIR}/${SORT_TAG}"
  local OUT_DIR="${OUT_BASE}/${MARK_LABEL}/Without_input"

  mkdir -p "$OUT_DIR" "$SORTED_DIR" || {
    echo "❌ Failed to create output directories for $MARK_LABEL"
    return 1
  }

  # Adjust Genrich parameters based on peak type
  if [[ "$MARK_TYPE" == "broad" ]]; then
    GENRICH_ARGS="-a 100 -l 500 -g 1000 -p 0.05 -f BAM"
    EXT=".broadPeak"
  else
    GENRICH_ARGS="-a 200 -l 100 -g 100 -p 0.01 -f BAM"
    EXT=".narrowPeak"
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo -e "\n========== Processing $MARK_LABEL ($CELL) =========="

    local RAW_TREATMENT="${MARK_DIR}/${MARK}_${CELL}.bam"
    local TREATMENT="${SORTED_DIR}/${MARK}_${CELL}_qnamesorted.bam"
    local OUTPUT_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_peaks${EXT}"
    local LOG_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_genrich.log"

    if [[ ! -f "$RAW_TREATMENT" ]]; then
      echo "❌ Missing treatment BAM: $RAW_TREATMENT"
      continue
    fi

    # Sort BAM by queryname if not already done
    if [[ ! -f "$TREATMENT" ]]; then
      echo "🔹 Sorting treatment BAM..."
      samtools sort -n -@ $THREADS -o "$TREATMENT" "$RAW_TREATMENT" >> "$LOG_FILE" 2>&1
    fi

    # Run Genrich without control
    echo "🔹 Running Genrich peak calling (no control)..."
    Genrich -t "$TREATMENT" -o "$OUTPUT_FILE" -j -r -v $GENRICH_ARGS >> "$LOG_FILE" 2>&1

    if [[ $? -eq 0 ]]; then
      echo "✅ Completed: $CELL ($MARK_LABEL)"
    else
      echo "❌ Genrich failed: $CELL ($MARK_LABEL)"
    fi
  done

  # Peak summary
  echo -e "\n========= Peak Summary for $MARK_LABEL ========="
  for peakfile in "${OUT_DIR}"/*${EXT}; do
    [[ -f "$peakfile" ]] && echo "$(basename "$peakfile"): $(wc -l < "$peakfile") peaks"
  done
}

#####=========================================================================================================####
# Run for all broad and narrow histone marks
for MARK in "${BROAD_MARKS[@]}"; do
  run_genrich_all_noctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  run_genrich_all_noctrl "$MARK" "narrow"
done
###############################################################################################################




#!/bin/bash
###############################################################################################################
### Genrich peak calling for multiple histone marks (broad and narrow)
### Supports dual processing for H3K27ac (both broad and narrow)
### Real Bam MouseBrain Data
###############################################################################################################

#####=========================================================================================================####
### for with input
#####=========================================================================================================####

# Define histone mark groups
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/MouseBrain_peakbed/Genrich_peakbed"

#####=========================================================================================================####

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

#####=========================================================================================================####

run_genrich_all_ctrl() {
  local MARK="$1"
  local MARK_TYPE="$2" # "broad" or "narrow"
  local SORT_TAG="qname_sorted_${MARK_TYPE:0:1}"  # b for broad, s for narrow
  local MARK_LABEL="${MARK}_${MARK_TYPE}"

  get_celltypes_for_mark "$MARK"
  if [[ ${#CELLTYPES[@]} -eq 0 ]]; then
    echo "❌ No cell types defined for $MARK"
    return
  fi

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
  local SORTED_DIR="${MARK_DIR}/${SORT_TAG}"
  local OUT_DIR="${OUT_BASE}/${MARK_LABEL}/With_input"

  mkdir -p "$OUT_DIR" "$SORTED_DIR" || {
    echo "❌ Failed to create output directories for $MARK_LABEL"
    return 1
  }

  # Adjust Genrich parameters based on mark type
  if [[ "$MARK_TYPE" == "broad" ]]; then
    GENRICH_ARGS="-a 100 -l 500 -g 1000 -p 0.05 -f BAM"
    EXT=".broadPeak"
  else
    GENRICH_ARGS="-a 200 -l 100 -g 100 -p 0.01 -f BAM"
    EXT=".narrowPeak"
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo -e "\n========== Processing $MARK_LABEL ($CELL) =========="

    local RAW_TREATMENT="${MARK_DIR}/${MARK}_${CELL}.bam"
    local RAW_CONTROL="${MARK_DIR}/input_${CELL}.bam"
    local TREATMENT="${SORTED_DIR}/${MARK}_${CELL}_qnamesorted.bam"
    local CONTROL="${SORTED_DIR}/input_${CELL}_qnamesorted.bam"
    local OUTPUT_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_peaks${EXT}"
    local LOG_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_genrich.log"

    if [[ ! -f "$RAW_TREATMENT" ]]; then
      echo "❌ Missing treatment BAM: $RAW_TREATMENT"
      continue
    fi
    if [[ ! -f "$RAW_CONTROL" ]]; then
      echo "❌ Missing control BAM: $RAW_CONTROL"
      continue
    fi

    if [[ ! -f "$TREATMENT" ]]; then
      echo "🔹 Sorting treatment BAM..."
      samtools sort -n -@ $THREADS -o "$TREATMENT" "$RAW_TREATMENT" >> "$LOG_FILE" 2>&1
    fi
    if [[ ! -f "$CONTROL" ]]; then
      echo "🔹 Sorting control BAM..."
      samtools sort -n -@ $THREADS -o "$CONTROL" "$RAW_CONTROL" >> "$LOG_FILE" 2>&1
    fi

    echo "🔹 Running Genrich peak calling..."
    Genrich -t "$TREATMENT" -c "$CONTROL" -o "$OUTPUT_FILE" -j -r -v $GENRICH_ARGS >> "$LOG_FILE" 2>&1

    echo "✅ Completed: $CELL ($MARK_LABEL)"
  done

  echo -e "\n========= Peak Summary for $MARK_LABEL ========="
  for peakfile in "${OUT_DIR}"/*${EXT}; do
    [[ -f "$peakfile" ]] && echo "$(basename "$peakfile"): $(wc -l < "$peakfile") peaks"
  done
}

###############################################################################################################
# Run for all marks
for MARK in "${BROAD_MARKS[@]}"; do
  run_genrich_all_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  run_genrich_all_ctrl "$MARK" "narrow"
done
###############################################################################################################


#!/bin/bash
###############################################################################################################
### Genrich peak calling for multiple histone marks (WITHOUT input/control)
### Works for both broad and narrow marks
### Optimized and robust version
###############################################################################################################

#--------------------------------------------
# Define histone mark groups
#--------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/MouseBrain_peakbed/Genrich_peakbed"

#--------------------------------------------
# Function: define cell types per histone mark
#--------------------------------------------
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

#--------------------------------------------
# Function: Run Genrich (no control)
#--------------------------------------------
run_genrich_all_noctrl() {
  local MARK="$1"
  local MARK_TYPE="$2"  # "broad" or "narrow"
  local SORT_TAG="qname_sorted_${MARK_TYPE:0:1}"  # b or s
  local MARK_LABEL="${MARK}_${MARK_TYPE}"

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
  local SORTED_DIR="${MARK_DIR}/${SORT_TAG}"
  local OUT_DIR="${OUT_BASE}/${MARK_LABEL}/Without_input"

  mkdir -p "$OUT_DIR" "$SORTED_DIR" || {
    echo "❌ Failed to create output directories for $MARK_LABEL"
    return 1
  }

  # Load cell types
  get_celltypes_for_mark "$MARK"

  # Adjust Genrich parameters
  if [[ "$MARK_TYPE" == "broad" ]]; then
    GENRICH_ARGS="-a 100 -l 500 -g 1000 -p 0.05 -f BAM"
    EXT=".broadPeak"
  else
    GENRICH_ARGS="-a 200 -l 100 -g 100 -p 0.01 -f BAM"
    EXT=".narrowPeak"
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo -e "\n========== Processing $MARK_LABEL ($CELL) =========="

    local RAW_TREATMENT="${MARK_DIR}/${MARK}_${CELL}.bam"
    local TREATMENT="${SORTED_DIR}/${MARK}_${CELL}_qnamesorted.bam"
    local OUTPUT_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_peaks${EXT}"
    local LOG_FILE="${OUT_DIR}/${CELL}_${MARK_LABEL}_genrich.log"

    if [[ ! -f "$RAW_TREATMENT" ]]; then
      echo "❌ Missing treatment BAM: $RAW_TREATMENT"
      continue
    fi

    # Sort BAM by queryname if needed
    if [[ ! -f "$TREATMENT" ]]; then
      echo "🔹 Sorting treatment BAM..."
      samtools sort -n -@ $THREADS -o "$TREATMENT" "$RAW_TREATMENT" >> "$LOG_FILE" 2>&1
    fi

    # Backup old results (if exist)
    [[ -f "$OUTPUT_FILE" ]] && mv "$OUTPUT_FILE" "${OUTPUT_FILE}.bak_$(date +%s)"

    # Run Genrich
    echo "🔹 Running Genrich peak calling (no control)... $(date)" >> "$LOG_FILE"
    Genrich -t "$TREATMENT" -o "$OUTPUT_FILE" -j -r -v $GENRICH_ARGS >> "$LOG_FILE" 2>&1

    if [[ $? -eq 0 ]]; then
      echo "✅ Completed: $CELL ($MARK_LABEL)"
    else
      echo "❌ Genrich failed: $CELL ($MARK_LABEL)"
    fi
  done

  # Summary
  echo -e "\n========= Peak Summary for $MARK_LABEL ========="
  total=0
  for peakfile in "${OUT_DIR}"/*${EXT}; do
    if [[ -f "$peakfile" ]]; then
      count=$(wc -l < "$peakfile")
      echo "$(basename "$peakfile"): $count peaks"
      total=$((total + count))
    fi
  done
  echo "📊 Total peaks across all cell types for $MARK_LABEL: $total"
}

#--------------------------------------------
# Run for all marks
#--------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
  run_genrich_all_noctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  run_genrich_all_noctrl "$MARK" "narrow"
done

###############################################################################################################
