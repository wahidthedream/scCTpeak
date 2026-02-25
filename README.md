<img width="100" height="100" alt="image" src="https://github.com/user-attachments/assets/af2565e3-bcc2-4adf-8f5d-3fddda8f51ca" />





# Comprehensive multiomics benchmarking of peak-calling methods for scCUT&Tag with scCTpeak

**scCTpeak** is a fully reproducible benchmarking framework designed to systematically evaluate seven widely used peak-calling tools across two distinct biological systems (Human PBMC and Mouse Brain) and multiple histone modifications of single-cell CUT&Tag (scCUT&Tag) datasets. It integrates an end‑to‑end data processing module with two complementary input modes: (i) spliting a combined BAM files into cell‑type‑specific treatment and input (control) BAMs using cell barcodes, and (ii) generating cell‑type‑specific BAMs directly from Seurat objects (RDS) and fragment files. The framework assesses performance with and without control samples, providing insights into tool robustness for single-cell epigenomics.
The pipeline benchmarks the following tools:


**Framework:**

<img width="1200" height="500" alt="benchmark_framework" src="https://github.com/user-attachments/assets/9652a027-f050-44ee-acf5-e591a57b6a92" />


**Command-Line Interface:** 

<img width="1200" height="1100" alt="scCTpeak_commandline" src="https://github.com/user-attachments/assets/79dd21f5-b76a-4003-a9fc-c9b5134c8049" />







## Using Instruction of *scCTpeak* 

### Step 1: ***Environment Setup***

First, create a working directory and install required tools in your conda environment:
To use scCTpeak, ensure you have Python 3.9 or higher. Install the required dependencies using:
````bash
pip install -r requirements.txt

# Create and activate conda environment
conda create -n scCTpeak python=3.9
conda activate scCTpeak
````


### Step 2: ***Run scCTpeak***

Ensure the scCTpeak function available in your shell session:

```bash

# Source the scCTpeak script or add it to your bash profile

source scCTpeak.sh  # or add to ~/.bashrc

````

### Step 3: ***Explore scCTpeak Commands***

Use the help command to see all available options:

```bash

scCTpeak help
```

### Command Reference for Data Processing

#### Split a combined BAM file into cell‑type BAMs using barcode lists for real data

````bash
  split_bam <dataset> <histone>

````
Split a combined BAM file into cell‑type BAMs using barcode lists, then create input (control) BAMs by merging other cell types.
Example: 

````bash
scCTpeak split_bam human H3K27ac
````

Run split_bam for all histones/TFs in the dataset.
Example: 

````bash
scCTpeak split_bam_all human
````

#### Data simulation data processing from RDS:

Process single histone/TF to BAM files
Examples:

````bash
scCTpeak process human H3K27ac
scCTpeak process mouse H3K4me3
````
Process all histones/TFs for dataset
Examples:

````bash
scCTpeak process_all human
scCTpeak process_all mouse
````


### Command Reference for Peak Calling

***run - Run specific tool***

Execute a single peak-calling tool with specific parameters.

Syntax:

```bash

scCTpeak run <dataset> <tool> <histone> <cell_type> <mode> <control>
```

Example:

````bash

# Run MACS2 for H3K27ac in B cells with broad peak calling using input control

scCTpeak run human macs2 H3K27ac B broad with_input

`````

***parse2wig - DROMPAplus preprocessing***

Required preprocessing step for DROMPAplus before peak calling.

Syntax:

````bash
scCTpeak parse2wig <dataset> <histone> <cell_type>
````

Example:

````bash
# Preprocess H3K27ac data for B cells
scCTpeak parse2wig human H3K27ac B
````

***batch - Batch process all histone marks***

Run all histone marks for a specific tool.

Syntax:

````bash

scCTpeak batch <dataset> <tool> <cell_type> <control>
````

Example:

````bash

# Process all histone marks for MACS2 in all cell types with input control
scCTpeak batch human macs2 all with_input

````

***all_tools - Run complete benchmark***
Execute all supported tools for comprehensive analysis.

Syntax:

````bash
scCTpeak all_tools <dataset> <cell_type> <control>
````

Example:

````bash
# Run all tools on human PBMC data for B cells with input control
scCTpeak all_tools human B with_input
````


***organize - Organize unified outputs***
Collect all unified BED files into a single directory.

Syntax:

````bash

scCTpeak organize <dataset> <control>
````
Example:

````bash
# Organize all with_input control files
scCTpeak organize human with_input

````


***summary - Generate summary report***
Create a comprehensive summary of all unified outputs.

Syntax:

```bash
scCTpeak summary <dataset> <control>
```

Example:

```bash
# Generate summary for with_input control files
scCTpeak summary human with_input
````

# Key Features
## Multi-Tool Evaluation (Seven Tools)

**DROMPAplus** - Broad/Sharp peak detection with specialized parameters

**Genrich** - ATAC-seq optimized, works with BAM processing

**GoPeaks** - Simple, direct BAM processing

**HOMER** - Comprehensive motif analysis pipeline

**MACS2** - Industry standard for ChIP-seq

**SEACR** - Sparse enrichment analysis for CUT&RUN/Tag

**SICER2** - Spatial clustering approach for broad domains


## Experimental Conditions

With Input Control: For each cell type, we constructed input by pooling sequencing fragments from all other cell types.

Without Input Control: Peak calling on treatment-only data

Broad vs Sharp: Parameter optimization per mark type
### Conserved and divergent peak-type patterns across HMs and TFs
| Category | Human PBMC | Mouse Brain | Conservation Status |
|----------|------------|-------------|---------------------|
| **Dual-Mode (Sharp & Broad)** | • H3K27ac | • H3K27ac | Conserved |
| **Sharp Peaks** | • H3K4me1<br>• H3K4me2<br>• H3K4me3 | • H3K4me3<br>• Olig2 (TF)<br>• Rad21 (cohesin) | Context-dependent |
| **Broad Peaks** | • H3K27me3<br>• H3K9me3 | • H3K27me3<br>• H3K36me3 | Conserved |


## Biological System & Data Specifications

| System | Accession | Epigenomic Targets | Cell Types (n) | Genome Assembly |
|--------|-----------|-------------------|----------------|-----------------|
| **Human PBMC** | GSE195725 | **Histone modifications**<br>• H3K27ac<br>• H3K27me3<br>• H3K4me1<br>• H3K4me2<br>• H3K4me3<br>• H3K9me3 | **8 Cell types**<br>B cells, CD4 T, CD8 T, DC, Mono, NK, other T, other | GRCh38 (hg38) |
| **Mouse Brain** | GSE157637 | **Histone modifications**<br>• H3K27ac<br>• H3K27me3<br>• H3K36me3<br>• H3K4me3<br>**Transcription factors**<br>• Olig2<br>• Rad21 | **Cell-types**<br>• H3K27ac: Astrocytes, mOL, OEC, OPC, VLMC (n=5)<br>• H3K27me3: Astrocytes, Microglia, mOL, Neurons1, Neurons3, OEC, OPC, VLMC (n=8)<br>• H3K36me3: Astrocytes, mOL, OEC, OPC (n=4)<br>• H3K4me3: Astrocytes, Microglia, mOL, Neurons1–3, OEC, OPC, VLMC (n=9)<br>• Olig2: Astrocytes, mOL, OEC, Unknown (n=4)<br>• Rad21: Astrocytes, mOL, OEC, Unknown (n=4) | GRCm38 (mm10) |

## Tool-Specific Configurations

| Tool | Version | Broad Parameters | Sharp Parameters | Key Features |
|------|---------|------------------|------------------|--------------|
| **DROMPAplus** | v1.20.1 | `p_int=4, p_enr=3` | `p_int=5, p_enr=4` | parse2wig+ preprocessing, paired-end optimized |
| **Genrich** | v0.6.1 | `-a 100 -l 500 -g 1000 -p 0.05` | `-a 200 -l 100 -g 100 -p 0.01` | Requires BAMs, ATAC-seq optimized |
| **GoPeaks** | v1.0.0 | `--broad` | Default sharp | ATAC-seq optimized |
| **HOMER** | v5.1 | `histone` | `factor` | motif analysis |
| **MACS2** | v2.2.9.1 | `--broad --broad-cutoff 0.05` | Default sharp settings | BAMPE mode, FDR=0.05, --keep dup all |
| **SEACR** | v1.3 | Threshold `stringent` | Threshold `relaxed` | bedGraph input from coverage, CUT&RUN/Tag specialized |
| **SICER2** | v1.0.3 | `window=100, gap=200` | `window=50, gap=100` | FDR=0.05, spatial clustering for broad domains |


## Unified Input Strategy
For each cell type, we constructed a pseudo-input by pooling fragments from all remaining cell types.
In this design, the input control for a given cell type is derived from the complementary cell populations.
This cross-cell-type input strategies:
- Provides balanced background signal
- Reduces cell-type-specific bias
- Improves peak-calling robustness in scCUT&Tag data

## Quick Start

```bash
# Clone repository
git clone https://github.com/wahidthedream/scCTpeak
cd scCTpeak

