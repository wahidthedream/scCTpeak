<img width="833" height="811" alt="image" src="https://github.com/user-attachments/assets/af2565e3-bcc2-4adf-8f5d-3fddda8f51ca" />





# scCTPEAK: a benchmark suite for scCUT&Tag peak-calling across diverse epigenomic contexts

**scCTPEAK** is a fully reproducible benchmarking framework designed to systematically evaluate seven widely used peak-calling algorithms across two distinct biological systems (Human PBMC and Mouse Brain) and multiple histone modifications of single-cell CUT&Tag (scCUT&Tag) datasets. The framework assesses performance with and without control samples, providing insights into tool robustness for single-cell epigenomics.
The pipeline benchmarks the following tools:

# Key Features
## 1. Multi-Tool Evaluation (7 Algorithms)

**DROMPAplus** - Broad/Sharp peak detection with specialized parameters

**Genrich** - ATAC-seq optimized, works with BAM processing

**GoPeaks** - Simple, direct BAM processing

**HOMER** - Comprehensive motif analysis pipeline

**MACS2** - Industry standard for ChIP-seq

**SEACR** - Sparse enrichment analysis for CUT&RUN/Tag

**SICER2** - Spatial clustering approach for broad domains

## 2. Dual Biological Systems

Human PBMC (GSE195725): 8 cell types, 6 histone marks

Mouse Brain (GSE157637): 4-9 cell types, 6 histone marks + transcription factors

## 3. Experimental Conditions

With Control: Using matched input BAMs

Without Control: Peak calling on treatment-only data

Broad vs Sharp: Parameter optimization per mark type

Dual-mode for H3K27ac: Processed as both broad and sharp

# Dataset Specifications

## Human PBMC Dataset

Histone Marks: H3K27ac, H3K27me3, H3K9me3, H3K4me1, H3K4me2, H3K4me3

Cell Types: B, CD4T, CD8T, DC, Mono, NK, otherT, other (8 types)

Genome: hg38

Total Conditions: 6 marks × 8 cells × 2 (with /without input control) = 96 per tool

## Mouse Brain Dataset

Targets: H3K27ac, H3K27me3, H3K36me3, H3K4me3, Olig2, Rad21

Cell Types: Astrocytes, mOL, OEC, OPC, VLMC, Microglia, Neurons(1-3), Unknown

Genome: mm10

Cell-type specific: Different cell populations per histone mark

# Tool-Specific Configurations

### DROMPAplus
```bash
Broad marks: p_int=4, p_enr=3
Sharp marks: p_int=5, p_enr=4
Preprocessing: parse2wig+ for bigWig generation
````
### Genrich
```bash
Broad: -a 100 -l 500 -g 1000 -p 0.05
Sharp: -a 200 -l 100 -g 100 -p 0.01
Requirement: Queryname-sorted BAMs
````
### MACS2
```bash
Genome size: Human=2.7e9, Mouse=1.87e9
Mode: BAMPE for paired-end
FDR: 0.05 for both broad and narrow (sharp)
--keep-dup all
```
### SEACR
```bash
Threshold: Stringent (broad) vs Relaxed (narrow)
Input: bedGraph files from BAM coverage
```
### SICER2
```bash
Broad: window=100, gap=200
Narrow: window=50, gap=100
FDR: 0.05 across all
```
## Key Features
1. Unified peak-calling interface
For each cell type, we constructed a pseudo-input by pooling fragments from all remaining cell types.
In this design, the input control for a given cell type is derived from the complementary cell populations. This cross–cell-type input strategy provides a balanced background signal, reduces cell-type–specific bias, and improves the robustness of peak calling in scCUT&Tag data.
