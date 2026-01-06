#!/bin/bash
###############################################################################################################
### SICER2 Peak Calling WITH Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/SICER2_peakbed_corrected"
mkdir -p "${OUT_BASE}"

# Human genome parameters
GENOME="hg38"
FDR=0.05

###############################################################################################################
# Function: Run SICER2 with Control
###############################################################################################################
run_sicer2_all_input() {
    local MARK="$1"
    local MODE="$2"  # broad or narrow

    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected/"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "$OUT_DIR"

    # Adjust window and gap for broad vs narrow
    if [[ "$MODE" == "broad" ]]; then
        WINDOW=100
        GAP=200
    else
        WINDOW=50
        GAP=100
    fi

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running SICER2 for: $MARK (${CELL}) =========="
        CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
        mkdir -p "$CELL_OUT"

        # Check BAM files
        if [[ ! -f "${MARK_DIR}/${MARK}_${CELL}.bam" ]]; then
            echo "Treatment BAM missing: ${MARK}_${CELL}.bam. Skipping..."
            continue
        fi
        if [[ ! -f "${MARK_DIR}/input_${CELL}.bam" ]]; then
            echo "Control BAM missing: input_${CELL}.bam. Skipping..."
            continue
        fi

        # Run SICER2
        sicer \
          --treatment_file "${MARK_DIR}/${MARK}_${CELL}.bam" \
          --control_file "${MARK_DIR}/input_${CELL}.bam" \
          --species "$GENOME" \
          --fragment_size 150 \
          --window_size $WINDOW \
          --gap_size $GAP \
          --effective_genome_fraction 0.80 \
          --false_discovery_rate $FDR \
          --redundancy_threshold 1 \
          --output_directory "$CELL_OUT" \
          --cpu $THREADS \
          --significant_reads \
          > "${CELL_OUT}/sicer.log" 2>&1

        echo "✅ Finished: ${MARK}_${CELL} (SICER2)"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SICER2 Peak Summary for $MARK ($MODE) ========="
    for CELL in "${CELLTYPES[@]}"; do
        DIR="${OUT_DIR}/${MARK}_${CELL}"
        PEAK_FILE=$(find "$DIR" -type f -name "*-FDR0.05-island.bed" | head -n 1)
        if [[ -f "$PEAK_FILE" ]]; then
            COUNT=$(wc -l < "$PEAK_FILE")
            echo "${MARK}_${CELL}: $COUNT peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

###############################################################################################################
# Run for Broad and Narrow Marks
###############################################################################################################
for MARK in "${BROAD_MARKS[@]}"; do
    run_sicer2_all_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_sicer2_all_input "$MARK" "narrow"
done

#!/bin/bash
###############################################################################################################
### SICER2 Peak Calling WITHOUT Control (No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/SICER2_peakbed_corrected"
mkdir -p "${OUT_BASE}"

GENOME="hg38"
FDR=0.05

###############################################################################################################
# Function: Run SICER2 without Control
###############################################################################################################
run_sicer2_no_input() {
  local MARK="$1"
  local MODE="$2"  # broad or narrow

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
  local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
  mkdir -p "$OUT_DIR"

  # Adjust window and gap for broad vs narrow
  if [[ "$MODE" == "broad" ]]; then
    WINDOW=100; GAP=200
  else
    WINDOW=50; GAP=100
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo "========== Running SICER2 for: $MARK (${CELL}) =========="

    CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
    mkdir -p "$CELL_OUT"

    BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
    if [[ ! -f "$BAM" ]]; then
      echo "Missing BAM file: $BAM"
      continue
    fi

    sicer \
      --treatment_file "$BAM" \
      --species "$GENOME" \
      --fragment_size 150 \
      --window_size "$WINDOW" \
      --gap_size "$GAP" \
      --effective_genome_fraction 0.80 \
      --false_discovery_rate "$FDR" \
      --redundancy_threshold 1 \
      --output_directory "$CELL_OUT" \
      --cpu "$THREADS" \
      --e_value 1000 \
      --significant_reads \
      > "${CELL_OUT}/sicer.log" 2>&1

    # Convert .scoreisland to .bed
    SCORE_FILE=$(find "$CELL_OUT" -name "*scoreisland" | head -n 1)
    if [[ -f "$SCORE_FILE" ]]; then
      BED_FILE="${SCORE_FILE%.scoreisland}-FDR${FDR}-island.bed"
      awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $5, "."}' "$SCORE_FILE" > "$BED_FILE"
    fi

    echo "✅ Finished: ${MARK}_${CELL}"
    echo
  done

  echo "========= SICER2 Peak Summary for $MARK ========="
  for CELL in "${CELLTYPES[@]}"; do
    DIR="${OUT_DIR}/${MARK}_${CELL}"
    PEAK_FILE=$(find "$DIR" -type f -name "*-FDR${FDR}-island.bed" | head -n 1)
    if [[ -f "$PEAK_FILE" ]]; then
      count=$(wc -l < "$PEAK_FILE")
      echo "${MARK}_${CELL}: $count peaks"
    else
      echo "${MARK}_${CELL}: 0 peaks"
    fi
  done
  echo
}

#------------------------------------------
# Run SICER2 for Broad and Narrow Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
  run_sicer2_no_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  run_sicer2_no_input "$MARK" "narrow"
done



#!/bin/bash
###############################################################################################################
### SICER2 Peak Calling WITH Control (with Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
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
    Olig2|Rad21) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
    *) echo "Unknown mark: $MARK"; CELLTYPES=() ;;
  esac
}

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/SICER2_peakbed"
mkdir -p "${OUT_BASE}"

# Mouse genome
GENOME="mm10"
FDR=0.05

#------------------------------------------
# Function: Run SICER2 with Control
#------------------------------------------
run_sicer2_all_input() {
    local MARK="$1"
    local MODE="$2"  # broad or narrow

    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "$OUT_DIR"

    # Adjust window and gap for broad vs narrow
    if [[ "$MODE" == "broad" ]]; then
        WINDOW=100
        GAP=200
    else
        WINDOW=50
        GAP=100
    fi

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running SICER2 for: $MARK (${CELL}) =========="
        CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
        mkdir -p "$CELL_OUT"

        # Check BAM files
        if [[ ! -f "${MARK_DIR}/${MARK}_${CELL}.bam" ]]; then
            echo "Treatment BAM missing: ${MARK}_${CELL}.bam. Skipping..."
            continue
        fi
        if [[ ! -f "${MARK_DIR}/input_${CELL}.bam" ]]; then
            echo "Control BAM missing: input_${CELL}.bam. Skipping..."
            continue
        fi

        # Run SICER2
        sicer \
          --treatment_file "${MARK_DIR}/${MARK}_${CELL}.bam" \
          --control_file "${MARK_DIR}/input_${CELL}.bam" \
          --species "$GENOME" \
          --fragment_size 150 \
          --window_size $WINDOW \
          --gap_size $GAP \
          --effective_genome_fraction 0.80 \
          --false_discovery_rate $FDR \
          --redundancy_threshold 1 \
          --output_directory "$CELL_OUT" \
          --cpu $THREADS \
          --significant_reads \
          > "${CELL_OUT}/sicer.log" 2>&1

        echo "✅ Finished: ${MARK}_${CELL} (SICER2)"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SICER2 Peak Summary for $MARK ($MODE) ========="
    for CELL in "${CELLTYPES[@]}"; do
        DIR="${OUT_DIR}/${MARK}_${CELL}"
        PEAK_FILE=$(find "$DIR" -type f -name "*-FDR0.05-island.bed" | head -n 1)
        if [[ -f "$PEAK_FILE" ]]; then
            COUNT=$(wc -l < "$PEAK_FILE")
            echo "${MARK}_${CELL}: $COUNT peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

#------------------------------------------
# Run SICER2 for Broad Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_sicer2_all_input "$MARK" "broad"
done

#------------------------------------------
# Run SICER2 for Narrow Marks
#------------------------------------------
for MARK in "${NARROW_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_sicer2_all_input "$MARK" "narrow"
done



#!/bin/bash
###############################################################################################################
### SICER2 Peak Calling WITHOUT Control (No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

#------------------------------------------
# Histone marks
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
    Olig2|Rad21) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
    *) echo "Unknown mark: $MARK"; CELLTYPES=() ;;
  esac
}

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/SICER2_peakbed"
GENOME="mm10"
FDR=0.05
mkdir -p "$OUT_BASE"

#------------------------------------------
# Function: Run SICER2 without Control
#------------------------------------------
run_sicer2_no_input() {
  local MARK="$1"
  local MODE="$2"  # broad or narrow

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
  local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
  mkdir -p "$OUT_DIR"

  # Adjust window and gap
  if [[ "$MODE" == "broad" ]]; then
    WINDOW=100; GAP=200
  else
    WINDOW=50; GAP=100
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo "========== Running SICER2 for: $MARK (${CELL}) =========="

    CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
    mkdir -p "$CELL_OUT"

    BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
    if [[ ! -f "$BAM" ]]; then
      echo "⚠️ Missing BAM file: $BAM"
      continue
    fi

    sicer \
      --treatment_file "$BAM" \
      --species "$GENOME" \
      --fragment_size 150 \
      --window_size "$WINDOW" \
      --gap_size "$GAP" \
      --effective_genome_fraction 0.80 \
      --false_discovery_rate "$FDR" \
      --redundancy_threshold 1 \
      --output_directory "$CELL_OUT" \
      --cpu "$THREADS" \
      --e_value 1000 \
      --significant_reads \
      > "${CELL_OUT}/sicer.log" 2>&1

    # Convert .scoreisland to .bed
    SCORE_FILE=$(find "$CELL_OUT" -name "*scoreisland" | head -n 1)
    if [[ -f "$SCORE_FILE" ]]; then
      BED_FILE="${SCORE_FILE%.scoreisland}-FDR${FDR}-island.bed"
      awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $5, "."}' "$SCORE_FILE" > "$BED_FILE"
    fi

    echo "✅ Finished: ${MARK}_${CELL}"
    echo
  done

  echo "========= SICER2 Peak Summary for $MARK ========="
  for CELL in "${CELLTYPES[@]}"; do
    DIR="${OUT_DIR}/${MARK}_${CELL}"
    PEAK_FILE=$(find "$DIR" -type f -name "*-FDR${FDR}-island.bed" | head -n 1)
    if [[ -f "$PEAK_FILE" ]]; then
      count=$(wc -l < "$PEAK_FILE")
      echo "${MARK}_${CELL}: $count peaks"
    else
      echo "${MARK}_${CELL}: 0 peaks"
    fi
  done
  echo
}

#------------------------------------------
# Run SICER2 for Broad and Narrow Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
  get_celltypes_for_mark "$MARK"
  run_sicer2_no_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  get_celltypes_for_mark "$MARK"
  run_sicer2_no_input "$MARK" "narrow"
done

