<img width="833" height="811" alt="image" src="https://github.com/user-attachments/assets/af2565e3-bcc2-4adf-8f5d-3fddda8f51ca" />





# scCTPEAK: a benchmark suite for scCUT&Tag peak-calling across diverse epigenomic contexts

**scCTPEAK** is a fully reproducible benchmarking framework designed to systematically evaluate seven widely used peak-calling algorithms across two distinct biological systems (Human PBMC and Mouse Brain) and multiple histone modifications of single-cell CUT&Tag (scCUT&Tag) datasets. The framework assesses performance with and without control samples, providing insights into tool robustness for single-cell epigenomics.
The pipeline benchmarks the following tools:

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
git clone https://github.com/wahidthedream/scCTPEAK
cd scCTPEAK

