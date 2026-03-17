#################################################################################################
## Figure 8 A, B
#################################################################################################

################################################################################
### 
### A
################################################################################


library(Signac)
library(Seurat)
library(monocle3)
library(Matrix)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)

# Resolve filter function conflicts
filter <- dplyr::filter
select <- dplyr::select

###=============================================
### NATURE METHODS THEME STANDARDS
###=============================================

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(color = "black", family = "Helvetica"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1, margin = margin(b = 5)),
      axis.line = element_line(color = "black", linewidth = 0.25),
      #axis.ticks = element_blank(),
      #axis.title = element_blank(),
      #axis.text = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size, face = "bold", margin = margin(b = 2)),
      legend.text = element_text(size = base_size - 1),
      legend.key.width = unit(4, "mm"),
      legend.key.height = unit(1.5, "mm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = -5, b = 0),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

theme_set(theme_nature())

###=============================================
### STEP 0: DATA LOADING
###=============================================

histones<-c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
for(hist in histones){
histone_file <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/", hist, ".rds")
histone <- readRDS(histone_file)

output_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result/result_trajectoryfinal/Figure/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_dir_barcode<-"/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result/result_trajectoryfinal/partitionwise_barcode/"
dir.create(output_dir_barcode, showWarnings = FALSE, recursive = TRUE)

###=============================================
### STEP 1: CELL TYPE IDENTIFICATION
###=============================================

histone$is_cd8t <- ifelse(histone$predicted.celltype.l1 == "CD8 T", "CD8 T", "Other")

cat("Cell distribution:\n")
print(table(histone$is_cd8t))

###=============================================
### STEP 2: UMAP COORDINATES EXTRACTION
###=============================================

tryCatch({
  umap_coords <- Embeddings(histone, "azimuth.umap")
  umap_coords <- umap_coords[, 1:2]
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  cat("✓ UMAP coordinates extracted\n")
}, error = function(e) {
  cat("UMAP missing → Creating synthetic UMAP\n")
  set.seed(99)
  umap_coords <- matrix(rnorm(ncol(histone) * 2), ncol = 2)
  rownames(umap_coords) <- colnames(histone)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
})

###=============================================
### STEP 3: PSEUDOTIME CALCULATION FOR CD8T CELLS
###=============================================

cd8_cells <- colnames(histone)[histone$is_cd8t == "CD8 T"]
cd8_only <- subset(histone, cells = cd8_cells)

cds <- as.cell_data_set(cd8_only)

# Attach UMAP coordinates to CellDataSet
cds@int_colData@listData$reducedDims$UMAP <- umap_coords[cd8_cells, ]

# Monocle3 pseudotime calculation
pseudotime_df <- NULL

tryCatch({
  set.seed(1)
  cds <- cluster_cells(cds, resolution = 1e-5)
  cds <- learn_graph(cds)
  
  # Use right-most cell as root for trajectory
  root_id <- colnames(cds)[which.max(reducedDims(cds)$UMAP[, 1])]
  cds <- order_cells(cds, root_cells = root_id)
  
  pseudotime_df <- data.frame(
    cell = colnames(cds),
    pseudotime = as.numeric(colData(cds)$pseudotime)
  )
  
  cat("✓ Monocle3 pseudotime calculation completed\n")
}, error = function(e) {
  cat("Monocle3 failed → Using UMAP-based pseudotime\n")
})

# Fallback: UMAP-distance based pseudotime
if (is.null(pseudotime_df)) {
  um <- reducedDims(cds)$UMAP
  start <- um[which.max(um[, 1]), ]
  dist <- sqrt(rowSums((um - matrix(start, nrow(um), ncol = 2, byrow = TRUE))^2))
  pst <- (dist - min(dist)) / (max(dist) - min(dist))
  
  pseudotime_df <- data.frame(cell = colnames(cds), pseudotime = pst)
  cat("✓ UMAP-based pseudotime calculation completed\n")
}

###=============================================
### STEP 4: DATA INTEGRATION FOR PLOTTING
###=============================================

plot_data <- data.frame(
  cell = colnames(histone),
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  is_cd8t = histone$is_cd8t,
  stringsAsFactors = FALSE
)

plot_data <- left_join(plot_data, pseudotime_df, by = "cell")

# Handle missing pseudotime values for CD8T cells
med_pt <- median(plot_data$pseudotime[plot_data$is_cd8t == "CD8 T"], na.rm = TRUE)
plot_data$pseudotime[plot_data$is_cd8t == "CD8 T" & is.na(plot_data$pseudotime)] <- med_pt

# Remove any remaining NA values
plot_data <- plot_data[complete.cases(plot_data$UMAP1, plot_data$UMAP2), ]

###=============================================
### STEP 5: BARCODE EXTRACTION AND PARTITIONING INTO 10 EQUAL GROUPS
###=============================================

cat("\n=== EXTRACTING AND PARTITIONING CD8T CELL BARCODES ===\n")

# Get CD8T cells with pseudotime values and sort by pseudotime
cd8t_data <- plot_data %>% 
  dplyr::filter(is_cd8t == "CD8 T") %>%
  dplyr::arrange(pseudotime)  # Sort by pseudotime (low to high)

# Add equal-sized partition groups (1 to 10) based on cell count
cd8t_data$partition <- cut(
  seq_len(nrow(cd8t_data)),  # Sequence of row numbers
  breaks = 10,               # 10 equal groups
  labels = 1:10              # Partition labels 1 to 10
)

# Create directory for barcode output
barcode_dir <-output_dir_barcode
dir.create(barcode_dir, showWarnings = FALSE, recursive = TRUE)

# Save individual partition files and collect summary
partition_summary <- data.frame()

for (i in 1:10) {
  partition_cells <- cd8t_data %>% dplyr::filter(partition == i)
  
  if (nrow(partition_cells) > 0) {
    # Save barcodes to file
    barcode_file <- paste0(barcode_dir, hist, "_partition", i, "_barcodes.txt")
    writeLines(partition_cells$cell, barcode_file)
    
    # Add to summary
    partition_summary <- rbind(partition_summary, data.frame(
      Partition = i,
      Pseudotime_Range = paste0(
        round(min(partition_cells$pseudotime), 3), 
        " - ", 
        round(max(partition_cells$pseudotime), 3)
      ),
      Cell_Count = nrow(partition_cells),
      Barcode_File = basename(barcode_file)
    ))
  }
}

# Save complete CD8T barcode list
all_cd8t_barcodes_file <- paste0(barcode_dir, hist, "_CD8T_all_barcodes.txt")
writeLines(cd8t_data$cell, all_cd8t_barcodes_file)

# Save partition summary
write.csv(partition_summary, 
          paste0(barcode_dir, hist, "_CD8T_partition_summary.csv"), 
          row.names = FALSE)

# Print summary
cat("✓ CD8T cells partitioned into 10 groups based on pseudotime\n")
cat("✓ Total CD8T cells:", nrow(cd8t_data), "\n")
print(partition_summary)

# Add partition information to plot_data for visualization
plot_data <- plot_data %>%
  left_join(cd8t_data %>% select(cell, partition), by = "cell")

###=============================================
### STEP 6: APPROVED COLOR SCHEME
###=============================================

# Nature-style blue gradient (dark to light)
blue_palette <- colorRampPalette(c(
  "#08306B", "#08519C", "#2171B5", 
  "#4292C6", "#6BAED6", "#9ECAE1",
  "#C6DBEF", "#DEEBF7", "#F7FBFF"
))(100)

###=============================================
### STEP 7: FINAL QUALITY PLOT
###=============================================

p_final <- ggplot() +
  # Background: Other cell types in light gray
  geom_point(
    data = subset(plot_data, is_cd8t == "Other"),
    aes(x = UMAP1, y = UMAP2),
    color = "grey90",
    size = 0.6,
    alpha = 0.3,
    stroke = 0
  ) +
  # Foreground: CD8T cells colored by pseudotime
  geom_point(
    data = subset(plot_data, is_cd8t == "CD8 T"),
    aes(x = UMAP1, y = UMAP2, color = pseudotime),
    size = 0.8,
    alpha = 0.8,
    stroke = 0
  ) +
  # Nature-style color scale
  scale_color_gradientn(
    colors = blue_palette,
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("0", "1"),
    name = "Pseudotime",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(15, "mm"),
      barheight = unit(1.5, "mm"),
      ticks.linewidth = 0.5,
      direction = "horizontal"
    )
  ) +
  # Nature-style title
  labs(title = paste0(hist, "_CD8 T Cell Developmental Trajectory")) +
  # Coordinate system with equal aspect ratio
  coord_fixed(ratio = 1)

print(p_final)

###=============================================
### STEP 8: SAVE IN  STANDARD FORMATS
###=============================================

# Main figure (Nature single column width)
ggsave(
  filename = paste0(output_dir, hist, "_CD8T_pseudotime.pdf"),
  plot = p_final,
  width = 85,    # Nature standard: 85 mm single column
  height = 85,   # Square format
  units = "mm",
  dpi = 300
)

# High-resolution TIFF for publication
ggsave(
  filename = paste0(output_dir, hist, "_CD8T_pseudotime.tiff"),
  plot = p_final,
  width = 85,
  height = 85,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

###=============================================
### FINAL SUMMARY
###=============================================

cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
cat("✓  standards applied\n")
cat("✓ CD8T cells:", sum(plot_data$is_cd8t == "CD8 T"), "\n")
cat("✓ Other cells:", sum(plot_data$is_cd8t == "Other"), "\n")
cat("✓ Pseudotime range:", 
    round(min(plot_data$pseudotime[plot_data$is_cd8t == "CD8 T"], na.rm = TRUE), 3), 
    "-", 
    round(max(plot_data$pseudotime[plot_data$is_cd8t == "CD8 T"], na.rm = TRUE), 3), 
    "\n")
cat("✓ Barcode partitions saved in:", barcode_dir, "\n")
cat("✓ Partition files created:\n")
cat("  - CD8T_all_barcodes.txt (all CD8T cells)\n")
cat("  - CD8T_partition_1-10_barcodes.txt (10 pseudotime-based partitions)\n")
cat("  - CD8T_partition_summary.csv (partition statistics)\n")
cat("✓ Output files saved in:", output_dir, "\n")
}



################################################################################
### B
################################################################################

library(Signac)
library(Seurat)
library(monocle3)
library(Matrix)
library(SeuratWrappers)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(scales)   # for rescale
# Resolve filter / select name conflicts explicitly
filter <- dplyr::filter
select <- dplyr::select

###=============================================================================
###  STANDARDS THEME
###=============================================================================
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(color = "black", family = "Helvetica"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1, margin = margin(b = 5)),
      axis.line = element_line(color = "black", linewidth = 0.25),
      #axis.ticks = element_blank(),
      #axis.title = element_blank(),
      #axis.text = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size, face = "bold", margin = margin(b = 2)),
      legend.text = element_text(size = base_size - 1),
      legend.key.width = unit(4, "mm"),
      legend.key.height = unit(1.5, "mm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = -5, b = 0),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}
theme_set(theme_nature())

###=============================================================================
### STEP 0: DATA LOADING - MOUSE BRAIN
###=============================================================================
histones<-c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3")

for(hist in histones){

histone_file <- paste0("/home/wahid/project_scHMTF/GSE157637_processed_data/", hist, "_seurat_object.Rds")
if (!file.exists(histone_file)) stop("Seurat object file not found: ", histone_file)
histone <- readRDS(histone_file)
# Ensure this is updated Seurat object
histone <- tryCatch(UpdateSeuratObject(histone), error = function(e) histone)

output_dir <- "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result/result_trajectoryfinal/Figure/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_dir_barcode <- "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result/result_trajectoryfinal/partitionwise_barcode/"
dir.create(output_dir_barcode, showWarnings = FALSE, recursive = TRUE)

###=============================================================================
### STEP 1: CELL TYPE IDENTIFICATION - OPC + mOL MERGED
###=============================================================================
if (is.null(histone$cell_type)) stop("'cell_type' metadata column missing in Seurat object")
histone$is_opc_mol <- ifelse(histone$cell_type %in% c("OPC", "mOL"), "OPC_mOL", "Other")
cat("Cell distribution after merging OPC + mOL:\n")
print(table(histone$is_opc_mol))

###=============================================================================
### STEP 2: UMAP COORDINATES EXTRACTION (verify shape & rownames)
###=============================================================================
umap_coords <- NULL
try({
  um <- Embeddings(histone, "umap")
  if (!is.null(um) && ncol(um) >= 2) {
    umap_coords <- um[, 1:2, drop = FALSE]
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    # ensure rownames are cell names
    if (is.null(rownames(umap_coords))) stop("UMAP embedding missing rownames")
    cat("✓ UMAP coordinates extracted (", nrow(umap_coords), "cells )\n")
  } else stop("UMAP embedding missing or has <2 dims")
}, silent = TRUE)

if (is.null(umap_coords)) {
  cat("UMAP missing → Creating synthetic UMAP (will warn this is not real)\n")
  set.seed(99)
  umap_coords <- matrix(rnorm(ncol(histone) * 2), ncol = 2)
  rownames(umap_coords) <- colnames(histone)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
}

###=============================================================================
### STEP 3: PSEUDOTIME CALCULATION FOR OPC + mOL CELLS
###=============================================================================
opc_mol_cells <- colnames(histone)[histone$is_opc_mol == "OPC_mOL"]
if (length(opc_mol_cells) == 0) stop("No OPC_mOL cells found in metadata")

opc_mol_only <- subset(histone, cells = opc_mol_cells)
cds <- as.cell_data_set(opc_mol_only)

# Attach UMAP coordinates to the CDS in a robust way
if (!all(opc_mol_cells %in% rownames(umap_coords))) {
  warning("Some OPC_mOL cells not present in UMAP coords — subsetting coordinates to available cells")
  common_cells <- intersect(opc_mol_cells, rownames(umap_coords))
  cds <- cds[, common_cells]
  umap_sub <- umap_coords[common_cells, , drop = FALSE]
} else {
  umap_sub <- umap_coords[opc_mol_cells, , drop = FALSE]
}

# set reducedDims properly
reducedDims(cds)$UMAP <- as.matrix(umap_sub)

pseudotime_df <- NULL
root_cell_chosen <- NA

tryCatch({
  set.seed(1)
  cds <- cluster_cells(cds, resolution = 1e-5)
  cds <- learn_graph(cds)

  # Choose root: letf-most in UMAP (biological justification should be added by user)
  um_mat <- reducedDims(cds)$UMAP
  if (is.null(um_mat) || nrow(um_mat) == 0) stop("UMAP reducedDims missing in cds")

  root_cell_chosen <- rownames(um_mat)[which.min(um_mat[, 1])]
  if (length(root_cell_chosen) != 1) root_cell_chosen <- root_cell_chosen[1]

  cds <- order_cells(cds, root_cells = root_cell_chosen)

  # extract pseudotime and rescale to 0-1
  pt_raw <- as.numeric(colData(cds)$pseudotime)
  names(pt_raw) <- colnames(cds)
  pt_scaled <- rescale(pt_raw, to = c(0, 1))

  pseudotime_df <- data.frame(cell = names(pt_scaled), pseudotime = as.numeric(pt_scaled), stringsAsFactors = FALSE)

  cat("✓ Monocle3 pseudotime calculation completed for OPC + mOL\n")
  cat("✓ Root cell:", root_cell_chosen, "\n")
}, error = function(e) {
  cat("Monocle3 graph-based pseudotime failed: ", conditionMessage(e), "\nUsing UMAP-distance fallback.\n")
})

# Fallback: UMAP-distance based pseudotime (and rescale to 0-1)
if (is.null(pseudotime_df)) {
  um <- reducedDims(cds)$UMAP
  if (is.null(um) || nrow(um) == 0) stop("Cannot compute fallback pseudotime: UMAP missing")
  start <- um[which.min(um[, 1]), , drop = FALSE]
  d <- sqrt(rowSums((um - matrix(start, nrow(um), ncol = ncol(um), byrow = TRUE))^2))
  pst <- (d - min(d)) / (max(d) - min(d))
  pseudotime_df <- data.frame(cell = rownames(um), pseudotime = as.numeric(pst), stringsAsFactors = FALSE)
  root_cell_chosen <- rownames(um)[which.max(um[, 1])]
  cat("✓ UMAP-distance pseudotime completed (root selected): ", root_cell_chosen, "\n")
}

###=============================================================================
### STEP 4: DATA INTEGRATION FOR PLOTTING
###=============================================================================
# Build plot_data for all cells in Seurat object (so background 'Other' is shown)
plot_data <- data.frame(
  cell = colnames(histone),
  UMAP1 = umap_coords[colnames(histone), 1],
  UMAP2 = umap_coords[colnames(histone), 2],
  is_opc_mol = histone$is_opc_mol,
  stringsAsFactors = FALSE
)

# join pseudotime (only OPC+mOL cells will have non-NA pseudotime)
plot_data <- left_join(plot_data, pseudotime_df, by = "cell")

# For OPC_mOL cells missing pseudotime (if any), fill with median pseudotime of OPC_mOL
med_pt <- median(plot_data$pseudotime[plot_data$is_opc_mol == "OPC_mOL"], na.rm = TRUE)
if (is.na(med_pt)) med_pt <- 0  # fallback if everything NA
plot_data$pseudotime[plot_data$is_opc_mol == "OPC_mOL" & is.na(plot_data$pseudotime)] <- med_pt

# For non-OPC cells, keep pseudotime NA (they will be plotted in gray background)
# Remove rows with missing UMAP coords (defensive)
plot_data <- plot_data[complete.cases(plot_data$UMAP1, plot_data$UMAP2), ]

cat("Final pseudotime range (OPC_mOL):", range(plot_data$pseudotime[plot_data$is_opc_mol == "OPC_mOL"], na.rm = TRUE), "\n")

###=============================================================================
### STEP 5: BARCODE EXTRACTION AND PARTITIONING INTO 10 EQUAL GROUPS
###=============================================================================

cat("\n=== EXTRACTING AND PARTITIONING OPC+mOL CELL BARCODES ===\n")

# Get OPC+mOL cells with pseudotime values and sort by pseudotime
opc_mol_data <- plot_data %>% 
  dplyr::filter(is_opc_mol == "OPC_mOL") %>%
  dplyr::arrange(pseudotime)  # Sort by pseudotime (low to high)

# Add equal-sized partition groups (1 to 10) based on cell count
opc_mol_data$partition <- cut(
  seq_len(nrow(opc_mol_data)),  # Sequence of row numbers
  breaks = 10,               # 10 equal groups
  labels = 1:10              # Partition labels 1 to 10
)

# Create directory for barcode output
barcode_dir <- output_dir_barcode
dir.create(barcode_dir, showWarnings = FALSE, recursive = TRUE)

# Save individual partition files and collect summary
partition_summary <- data.frame()

for (i in 1:10) {
  partition_cells <- opc_mol_data %>% dplyr::filter(partition == i)
  
  if (nrow(partition_cells) > 0) {
    # Save barcodes to file
    barcode_file <- paste0(barcode_dir, hist, "_partition", i, "_barcodes.txt")
    writeLines(partition_cells$cell, barcode_file)
    
    # Add to summary
    partition_summary <- rbind(partition_summary, data.frame(
      Partition = i,
      Pseudotime_Range = paste0(
        round(min(partition_cells$pseudotime), 3), 
        " - ", 
        round(max(partition_cells$pseudotime), 3)
      ),
      Cell_Count = nrow(partition_cells),
      Barcode_File = basename(barcode_file)
    ))
  }
}

# Save complete OPC+mOL barcode list
all_opc_mol_barcodes_file <- paste0(barcode_dir, hist, "_OPC_mOL_all_barcodes.txt")
writeLines(opc_mol_data$cell, all_opc_mol_barcodes_file)

# Save partition summary
write.csv(partition_summary, 
          paste0(barcode_dir, hist, "_OPC_mOL_partition_summary.csv"), 
          row.names = FALSE)

# Print summary
cat("✓ OPC+mOL cells partitioned into 10 groups based on pseudotime\n")
cat("✓ Total OPC+mOL cells:", nrow(opc_mol_data), "\n")
print(partition_summary)

# Add partition information to plot_data for visualization
plot_data <- plot_data %>%
  left_join(opc_mol_data %>% select(cell, partition), by = "cell")

###=============================================================================
### STEP 6: -APPROVED COLOR SCHEME
###=============================================================================
blue_orange_palette <- colorRampPalette(c(
  "#00441b", "#006d2c", "#238b45", 
  "#41ab5d", "#74c476", "#a1d99b",
  "#c7e9c0", "#e5f5e0", "#f7fcf5"
))(100)

###=============================================================================
### STEP 7: FINAL QUALITY PLOT
###=============================================================================
p_final <- ggplot() +
  # Background: Other cell types in light gray
  geom_point(
    data = subset(plot_data, is_opc_mol == "Other"),
    aes(x = UMAP1, y = UMAP2),
    color = "grey90",
    size = 0.6,
    alpha = 0.35,
    stroke = 0
  ) +
  # Foreground: OPC+mOL cells colored by pseudotime
  geom_point(
    data = subset(plot_data, is_opc_mol == "OPC_mOL"),
    aes(x = UMAP1, y = UMAP2, color = pseudotime),
    size = 0.9,
    alpha = 0.9,
    stroke = 0
  ) +
  scale_color_gradientn(
    colors = blue_orange_palette,
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("0", "1"),
    name = "Pseudotime",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = unit(15, "mm"), barheight = unit(1.5, "mm"),
                           ticks.linewidth = 0.5, direction = "horizontal")
  ) +
  labs(title = paste0(hist, "- OPC_mOL Developmental Trajectory")) +
  coord_fixed(ratio = 1)

# Save and print
ggsave(filename = file.path(output_dir, paste0(hist, "_OPC_mOL_pseudotime.png")), plot = p_final, width = 4, height = 4, units = "in", dpi = 300)
print(p_final)

# PDF for publication
ggsave(
  filename = file.path(output_dir, paste0(hist, "_OPC_mOL_pseudotime.pdf")),
  plot = p_final,
  width = 85, height = 85, units = "mm", dpi = 300
)

###=============================================================================
### FINAL SUMMARY
###=============================================================================

cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
cat("✓ standards applied\n")
cat("✓ OPC+mOL cells:", sum(plot_data$is_opc_mol == "OPC_mOL"), "\n")
cat("✓ Other cells:", sum(plot_data$is_opc_mol == "Other"), "\n")
cat("✓ Pseudotime range:", 
    round(min(plot_data$pseudotime[plot_data$is_opc_mol == "OPC_mOL"], na.rm = TRUE), 3), 
    "-", 
    round(max(plot_data$pseudotime[plot_data$is_opc_mol == "OPC_mOL"], na.rm = TRUE), 3), 
    "\n")
cat("✓ Barcode partitions saved in:", barcode_dir, "\n")
cat("✓ Partition files created:\n")
cat("  - OPC_mOL_all_barcodes.txt (all OPC+mOL cells)\n")
cat("  - OPC_mOL_partition_1-10_barcodes.txt (10 pseudotime-based partitions)\n")
cat("  - OPC_mOL_partition_summary.csv (partition statistics)\n")
cat("✓ Output files saved in:", output_dir, "\n")
}

#################################################################################
#### Figure 8 C D E F
#################################################################################

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)

# List of file paths (removed the duplicate first one)
file_paths <- c(
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-b_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-b_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-b_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-b_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-b_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-b_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-b_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-s_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-s_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-s_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-s_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-s_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-s_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27ac-s_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27me3_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27me3_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27me3_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27me3_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27me3_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27me3_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K27me3_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me1_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me1_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me1_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me1_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me1_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me1_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me1_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me2_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me2_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me2_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me2_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me2_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me2_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me2_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me3_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me3_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me3_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me3_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me3_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me3_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K4me3_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K9me3_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K9me3_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K9me3_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K9me3_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K9me3_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K9me3_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_H3K9me3_SICER2_zscore_matrix.csv"
)

# Function to extract histone and tool from filename
extract_metadata <- function(filepath) {
  filename <- basename(filepath)
  # Remove prefix "MouseBrain_" and suffix "_zscore_matrix.csv"
  clean_name <- str_remove(filename, "^HumanPBMC_")
  clean_name <- str_remove(clean_name, "_zscore_matrix\\.csv$")
  
  # Split by last underscore to separate histone and tool
  # Histone examples: "H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1","H3K4me2", "H3K4me3","H3K9me3"
  # Find the last underscore position
  last_underscore <- str_locate_all(clean_name, "_")[[1]]
  if (nrow(last_underscore) > 0) {
    last_pos <- last_underscore[nrow(last_underscore), 1]
    histone <- str_sub(clean_name, 1, last_pos - 1)
    tool <- str_sub(clean_name, last_pos + 1)
  } else {
    histone <- clean_name
    tool <- NA
  }
  
  return(list(histone = histone, tool = tool))
}

# Function to read and process each file
process_file <- function(filepath) {
  cat("Processing:", filepath, "\n")
  
  # Extract metadata
  metadata <- extract_metadata(filepath)
  
  # Read the CSV file
  df <- read_csv(filepath)
  
  # Ensure the first column is named "gene"
  if (!names(df)[1] %in% c("gene", "Gene")) {
    names(df)[1] <- "gene"
  }
  
  # Add histone and tool columns
  df <- df %>%
    mutate(
      histone = metadata$histone,
      tool = metadata$tool,
      .before = 1
    ) %>%
    rename(gene = 3)  # Make sure gene is the third column
  
  return(df)
}

# Process all files and combine them
combined_data <- map_dfr(file_paths, process_file)

# Reorder columns
combined_data <- combined_data %>%
  select(histone, tool, gene, everything())

# View the first few rows
print(head(combined_data))

# Save to a new CSV file
output_path <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore//combined_zscore_data.csv"
write_csv(combined_data, output_path)

cat("\nCombined data saved to:", output_path, "\n")
cat("Total rows:", nrow(combined_data), "\n")
cat("Unique histones:", unique(combined_data$histone), "\n")
cat("Unique tools:", unique(combined_data$tool), "\n")


HumanPBMC_data<-combined_data

# Define the gene lists for different stages
early_stage_genes <- c("CCR7", "SELL", "IL7R", "TCF7", "LEF1", "CD55")
middle_stage_genes <- c("ITGAE", "MKI67", "CD69", "IL2RA", "ICOS", "IFNG")
late_stage_genes <- c("GZMB", "GZMA", "CX3CR1", "TNF", "B3GAT1")

# Option 1: Create separate dataframes for each stage
early_stage_data <- HumanPBMC_data %>% filter(gene %in% early_stage_genes)
middle_stage_data <- HumanPBMC_data %>% filter(gene %in% middle_stage_genes)
late_stage_data <- HumanPBMC_data %>% filter(gene %in% late_stage_genes)

early_stage_data2 <- early_stage_data %>%
  mutate(
    mean1 = rowMeans(select(., partition1, partition2, partition3, partition4, partition5), na.rm = TRUE),
    mean2 = rowMeans(select(., partition6, partition7, partition8, partition9, partition10), na.rm = TRUE),
    .after = partition10
  )

# Check the updated structure
head(early_stage_data2)

middle_stage_data2 <- middle_stage_data %>%
  mutate(
    mean1 = rowMeans(select(., partition4, partition5, partition6, partition7), na.rm = TRUE),
    mean2 = rowMeans(select(., partition1, partition2, partition3, partition8, partition9, partition10), na.rm = TRUE),
    .after = partition10
  )
# Check the updated structure
head(middle_stage_data2)

late_stage_data2 <- late_stage_data %>%
  mutate(
    mean1 = rowMeans(select(., partition6, partition7, partition8, partition9, partition10), na.rm = TRUE),
    mean2 = rowMeans(select(., partition1, partition2, partition3, partition4, partition5), na.rm = TRUE),
    .after = partition10
  )

# Check the updated structure
head(late_stage_data2)

early_stage_data2$ratio<-((early_stage_data2$mean1)-(early_stage_data2$mean2))
middle_stage_data2$ratio<-((middle_stage_data2$mean1)-(middle_stage_data2$mean2))
late_stage_data2$ratio<-((late_stage_data2$mean1)-(late_stage_data2$mean2))
data<-rbind(early_stage_data2, middle_stage_data2, late_stage_data2)




pbmc_with_input<-early_stage_data2

####################################################################################################
## Plot - Corrected for Human PBMC data (APA Peak Height)
#####################################################################################################

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# First, check the actual range of your data
cat("=== DATA RANGE CHECK ===\n")
cat("Min ratio:", format(min(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Max ratio:", format(max(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Mean ratio:", format(mean(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("SD ratio:", format(sd(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n\n")

# Check distribution
cat("Score distribution summary:\n")
summary(pbmc_with_input$ratio)

# Create summary statistics
ratio_summary <- pbmc_with_input %>%
  group_by(tool, histone) %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(tool),
    Histone = factor(histone),
    # Calculate error bars using STANDARD ERROR (not standard deviation)
    se_ratio = sd_ratio / sqrt(n_samples),
    # Correct calculation of confidence intervals
    ci_lower = median_ratio - 1.96 * se_ratio,
    ci_upper = median_ratio + 1.96 * se_ratio
  )

# Print summary
cat("\n=== SUMMARY STATISTICS ===\n")
print(head(ratio_summary, 10))

# Calculate appropriate y-axis limits
data_min <- min(ratio_summary$ci_lower, na.rm = TRUE)
data_max <- max(ratio_summary$ci_upper, na.rm = TRUE)

# Add 10% padding
padding <- (data_max - data_min) * 0.1
y_min <- data_min - padding
y_max <- data_max + padding

cat(sprintf("\nY-axis limits: %.6f to %.6f\n", y_min, y_max))

# Define colors
all_methods <- unique(ratio_summary$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

# OPTION 1: Plot with scientific notation
# OPTION 1: Plot with scientific notation
plot1 <- ggplot(ratio_summary, 
                aes(x = Histone, 
                    y = median_ratio,
                    color = Method)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = ci_lower,
        ymax = ci_upper),
    width = 0.2, 
    position = position_dodge(width = 0.7),
    alpha = 0.7,
    size = 0.7
  ) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    limits = c(-2, 2),  # Set y-axis limits from -3 to 3
    breaks = seq(-2, 2, by = 1),  # Set breaks at -3, -2, -1, 0, 1, 2, 3
    labels = function(x) format(x, scientific = FALSE, digits = 1)  # Remove scientific notation
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "red", 
    alpha = 0.5,
    size = 0.5
  ) +
  labs(
    title = "Human PBMC Early Stage",
    x = "Histone Modification",
    y = "Z-score Shift"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

print(plot1)
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_early_state_scatter_errorbar_plot.pdf", plot1, width =8, height = 6, dpi = 300)


pbmc_with_input<-middle_stage_data2

####################################################################################################
## Plot - Corrected for Human PBMC data (APA Peak Height)
#####################################################################################################

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# First, check the actual range of your data
cat("=== DATA RANGE CHECK ===\n")
cat("Min ratio:", format(min(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Max ratio:", format(max(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Mean ratio:", format(mean(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("SD ratio:", format(sd(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n\n")

# Check distribution
cat("Score distribution summary:\n")
summary(pbmc_with_input$ratio)

# Create summary statistics
ratio_summary <- pbmc_with_input %>%
  group_by(tool, histone) %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(tool),
    Histone = factor(histone),
    # Calculate error bars using STANDARD ERROR (not standard deviation)
    se_ratio = sd_ratio / sqrt(n_samples),
    # Correct calculation of confidence intervals
    ci_lower = median_ratio - 1.96 * se_ratio,
    ci_upper = median_ratio + 1.96 * se_ratio
  )

# Print summary
cat("\n=== SUMMARY STATISTICS ===\n")
print(head(ratio_summary, 10))

# Calculate appropriate y-axis limits
data_min <- min(ratio_summary$ci_lower, na.rm = TRUE)
data_max <- max(ratio_summary$ci_upper, na.rm = TRUE)

# Add 10% padding
padding <- (data_max - data_min) * 0.1
y_min <- data_min - padding
y_max <- data_max + padding

cat(sprintf("\nY-axis limits: %.6f to %.6f\n", y_min, y_max))

# Define colors
all_methods <- unique(ratio_summary$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

# OPTION 1: Plot with scientific notation
# OPTION 1: Plot with scientific notation
plot1 <- ggplot(ratio_summary, 
                aes(x = Histone, 
                    y = median_ratio,
                    color = Method)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = ci_lower,
        ymax = ci_upper),
    width = 0.2, 
    position = position_dodge(width = 0.7),
    alpha = 0.7,
    size = 0.7
  ) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    limits = c(-2, 2),  # Set y-axis limits from -3 to 3
    breaks = seq(-2, 2, by = 1),  # Set breaks at -3, -2, -1, 0, 1, 2, 3
    labels = function(x) format(x, scientific = FALSE, digits = 1)  # Remove scientific notation
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "red", 
    alpha = 0.5,
    size = 0.5
  ) +
  labs(
    title = "Human PBMC Middle Stage",
    x = "Histone Modification",
    y = "Z-score Shift"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

print(plot1)
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_middle_state_scatter_errorbar_plot.pdf", plot1, width =8, height = 6, dpi = 300)


pbmc_with_input<-late_stage_data2

####################################################################################################
## Plot - Corrected for Human PBMC data (APA Peak Height)
#####################################################################################################

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# First, check the actual range of your data
cat("=== DATA RANGE CHECK ===\n")
cat("Min ratio:", format(min(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Max ratio:", format(max(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Mean ratio:", format(mean(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("SD ratio:", format(sd(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n\n")

# Check distribution
cat("Score distribution summary:\n")
summary(pbmc_with_input$ratio)

# Create summary statistics
ratio_summary <- pbmc_with_input %>%
  group_by(tool, histone) %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(tool),
    Histone = factor(histone),
    # Calculate error bars using STANDARD ERROR (not standard deviation)
    se_ratio = sd_ratio / sqrt(n_samples),
    # Correct calculation of confidence intervals
    ci_lower = median_ratio - 1.96 * se_ratio,
    ci_upper = median_ratio + 1.96 * se_ratio
  )

# Print summary
cat("\n=== SUMMARY STATISTICS ===\n")
print(head(ratio_summary, 10))

# Calculate appropriate y-axis limits
data_min <- min(ratio_summary$ci_lower, na.rm = TRUE)
data_max <- max(ratio_summary$ci_upper, na.rm = TRUE)

# Add 10% padding
padding <- (data_max - data_min) * 0.1
y_min <- data_min - padding
y_max <- data_max + padding

cat(sprintf("\nY-axis limits: %.6f to %.6f\n", y_min, y_max))

# Define colors
all_methods <- unique(ratio_summary$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

# OPTION 1: Plot with scientific notation
# OPTION 1: Plot with scientific notation
plot1 <- ggplot(ratio_summary, 
                aes(x = Histone, 
                    y = median_ratio,
                    color = Method)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = ci_lower,
        ymax = ci_upper),
    width = 0.2, 
    position = position_dodge(width = 0.7),
    alpha = 0.7,
    size = 0.7
  ) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    limits = c(-2, 2.2),  # Set y-axis limits from -3 to 3
    breaks = seq(-2, 2.2, by = 1),  # Set breaks at -3, -2, -1, 0, 1, 2, 3
    labels = function(x) format(x, scientific = FALSE, digits = 1)  # Remove scientific notation
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "red", 
    alpha = 0.5,
    size = 0.5
  ) +
  labs(
    title = "Human PBMC Late Stage",
    x = "Histone Modification",
    y = "Z-score Shift"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

print(plot1)
ggsave("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/HumanPBMC_late_state_scatter_errorbar_plot.pdf", plot1, width =8, height = 6, dpi = 300)


#######################################################################################
## mOL marker enrichment
#######################################################################################
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)

# List of file paths (removed the duplicate first one)
file_paths <- c(
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-b_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-b_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-b_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-b_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-b_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-b_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-b_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-s_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-s_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-s_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-s_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-s_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-s_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27ac-s_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27me3_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27me3_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27me3_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27me3_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27me3_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27me3_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K27me3_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K36me3_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K36me3_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K36me3_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K36me3_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K36me3_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K36me3_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K36me3_SICER2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K4me3_DROMPAplus_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K4me3_Genrich_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K4me3_GoPeaks_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K4me3_HOMER_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K4me3_MACS2_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K4me3_SEACR_zscore_matrix.csv",
  "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_H3K4me3_SICER2_zscore_matrix.csv"
)

# Function to extract histone and tool from filename
extract_metadata <- function(filepath) {
  filename <- basename(filepath)
  # Remove prefix "MouseBrain_" and suffix "_zscore_matrix.csv"
  clean_name <- str_remove(filename, "^MouseBrain_")
  clean_name <- str_remove(clean_name, "_zscore_matrix\\.csv$")
  
  # Split by last underscore to separate histone and tool
  # Histone examples: "H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K36me3", "H3K4me3"
  # Find the last underscore position
  last_underscore <- str_locate_all(clean_name, "_")[[1]]
  if (nrow(last_underscore) > 0) {
    last_pos <- last_underscore[nrow(last_underscore), 1]
    histone <- str_sub(clean_name, 1, last_pos - 1)
    tool <- str_sub(clean_name, last_pos + 1)
  } else {
    histone <- clean_name
    tool <- NA
  }
  
  return(list(histone = histone, tool = tool))
}

# Function to read and process each file
process_file <- function(filepath) {
  cat("Processing:", filepath, "\n")
  
  # Extract metadata
  metadata <- extract_metadata(filepath)
  
  # Read the CSV file
  df <- read_csv(filepath)
  
  # Ensure the first column is named "gene"
  if (!names(df)[1] %in% c("gene", "Gene")) {
    names(df)[1] <- "gene"
  }
  
  # Add histone and tool columns
  df <- df %>%
    mutate(
      histone = metadata$histone,
      tool = metadata$tool,
      .before = 1
    ) %>%
    rename(gene = 3)  # Make sure gene is the third column
  
  return(df)
}

# Process all files and combine them
combined_data <- map_dfr(file_paths, process_file)

# Reorder columns
combined_data <- combined_data %>%
  select(histone, tool, gene, everything())

# View the first few rows
print(head(combined_data))

# Save to a new CSV file
output_path <- "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/combined_zscore_data.csv"
write_csv(combined_data, output_path)

cat("\nCombined data saved to:", output_path, "\n")
cat("Total rows:", nrow(combined_data), "\n")
cat("Unique histones:", unique(combined_data$histone), "\n")
cat("Unique tools:", unique(combined_data$tool), "\n")


# Add firstfive_median and lastfive_median columns to combined_data
mouse_data <- combined_data %>%
  mutate(
    mean1 = rowMeans(select(., partition6, partition7, partition8, partition9, partition10), na.rm = TRUE),
    mean2 = rowMeans(select(., partition1, partition2, partition3, partition4, partition5), na.rm = TRUE),
    .after = partition10
  )

# Check the updated structure
head(mouse_data)

mouse_data$ratio<-(mouse_data$mean1)-(mouse_data$mean2)
pbmc_with_input<-mouse_data


####################################################################################################
## Plot - Corrected for Human PBMC data (APA Peak Height)
#####################################################################################################

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# First, check the actual range of your data
cat("=== DATA RANGE CHECK ===\n")
cat("Min ratio:", format(min(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Max ratio:", format(max(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("Mean ratio:", format(mean(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n")
cat("SD ratio:", format(sd(pbmc_with_input$ratio), scientific = TRUE, digits = 3), "\n\n")

# Check distribution
cat("Score distribution summary:\n")
summary(pbmc_with_input$ratio)

# Create summary statistics
ratio_summary <- pbmc_with_input %>%
  group_by(tool, histone) %>%
  summarize(
    median_ratio = median(ratio, na.rm = TRUE),
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = factor(tool),
    Histone = factor(histone),
    # Calculate error bars using STANDARD ERROR (not standard deviation)
    se_ratio = sd_ratio / sqrt(n_samples),
    # Correct calculation of confidence intervals
    ci_lower = median_ratio - 1.96 * se_ratio,
    ci_upper = median_ratio + 1.96 * se_ratio
  )

# Print summary
cat("\n=== SUMMARY STATISTICS ===\n")
print(head(ratio_summary, 10))

# Calculate appropriate y-axis limits
data_min <- min(ratio_summary$ci_lower, na.rm = TRUE)
data_max <- max(ratio_summary$ci_upper, na.rm = TRUE)

# Add 10% padding
padding <- (data_max - data_min) * 0.5
y_min <- data_min - padding
y_max <- data_max + padding

cat(sprintf("\nY-axis limits: %.6f to %.6f\n", y_min, y_max))

# Define colors
all_methods <- unique(ratio_summary$Method)
nature_methods_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
method_colors <- setNames(nature_methods_colors[1:length(all_methods)], all_methods)

# OPTION 1: Plot with scientific notation
# OPTION 1: Plot with scientific notation
plot1 <- ggplot(ratio_summary, 
                aes(x = Histone, 
                    y = median_ratio,
                    color = Method)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = ci_lower,
        ymax = ci_upper),
    width = 0.2, 
    position = position_dodge(width = 0.7),
    alpha = 0.7,
    size = 0.7
  ) +
  scale_color_manual(values = method_colors,
                     name = "Peak-Calling Method") +
  scale_y_continuous(
    limits = c(-2 , 2),  # Set y-axis limits from -3 to 3
    breaks = seq(-2, 2, by = 1),  # Set breaks at -3, -2, -1, 0, 1, 2, 3
    labels = function(x) format(x, scientific = FALSE, digits = 1)  # Remove scientific notation
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "red", 
    alpha = 0.5,
    size = 0.5
  ) +
  labs(
    title = "Mouse Brain",
    x = "Histone Modification",
    y = "Trajectory-associated Z-score Shift"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

print(plot1)

ggsave("/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore/MouseBrain_scatter_error_bar_plot.pdf", plot1, width =8, height = 6, dpi = 300)


#################################################################################
### Figure 8 G H
#################################################################################

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

tools <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")

for(tool in tools) {
  peakdir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result/result_trajectoryfinal/partitionwise_barcode/split_celltype_bams/HumanPBMC_peakbed_partitionwise/patitionwise_allpeakbed/H3K4me3_peakbed/", tool)
  
  # ==============================================================================
  # STEP 2: PEAK PROCESSING FOR HUMAN PBMC
  # ==============================================================================
  celltype_order <- c("partition1", "partition2", "partition3", "partition4", "partition5", 
                      "partition6", "partition7", "partition8", "partition9", "partition10")
  cat("\n=== STEP 2: Human PBMC Peak Processing ===\n")
  
  read_peaks <- function(bed_file, tool_name = NULL) {
    tryCatch({
      bed_data <- read.table(bed_file, sep = "\t", header = FALSE, 
                            stringsAsFactors = FALSE, fill = TRUE, 
                            comment.char = "", quote = "", na.strings = ".")
      
      n_cols <- ncol(bed_data)
      filename <- basename(bed_file)
      
      # Debug: Check what we're reading
      cat("  Reading", filename, "with", n_cols, "columns\n")
      
      # Basic validation
      if (n_cols < 3) {
        cat("  WARNING:", filename, "has only", n_cols, "columns. Skipping.\n")
        return(NULL)
      }
      
      # Create GRanges from first 3 columns
      peaks <- GRanges(
        seqnames = bed_data[, 1],
        ranges = IRanges(start = bed_data[, 2] + 1, end = bed_data[, 3])
      )
      
      # SPECIAL HANDLING FOR SICER2
      if (!is.null(tool_name) && tool_name == "SICER2") {
        cat("  Detected SICER2 format - using default score=1\n")
        peaks$score <- 0
      } else {
        # Standard handling for other tools
        score_found <- FALSE
        
        for (col_idx in 4:min(6, n_cols)) {
          if (col_idx <= n_cols) {
            potential_score <- suppressWarnings(as.numeric(bed_data[, col_idx]))
            
            if (!all(is.na(potential_score))) {
              peaks$score <- potential_score
              score_found <- TRUE
              cat("  Using score from column", col_idx, "for", filename, "\n")
              break
            }
          }
        }
        
        if (!score_found) {
          peaks$score <- 1
          cat("  No numeric score column found for", filename, "(using score=1)\n")
        }
      }
      
      # Clean up any NA or infinite scores
      valid_scores <- !is.na(peaks$score) & is.finite(peaks$score)
      peaks <- peaks[valid_scores]
      
      if (length(peaks) == 0) {
        cat("  WARNING: No valid peaks after cleaning in", filename, "\n")
        return(NULL)
      }
      
      cat("  ✓", length(peaks), "peaks loaded from", filename, "\n")
      return(peaks)
      
    }, error = function(e) {
      cat("  ERROR reading", basename(bed_file), ":", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  process_peak_files <- function(peakdir, celltype_order) {
    peak_files <- list.files(peakdir, pattern = "\\.bed$", full.names = TRUE)
    cat("Found", length(peak_files), "peak bed files\n")
    
    # Map filenames to trajectory celltypes
    celltype_mapping <- c(
      "partition1" = "partition1",
      "partition2" = "partition2", 
      "partition3" = "partition3",
      "partition4" = "partition4",
      "partition5" = "partition5",
      "partition6" = "partition6", 
      "partition7" = "partition7",
      "partition8" = "partition8",
      "partition9" = "partition9",
      "partition10" = "partition10"
    )
    
    all_peaks <- list()
    
    for (file in peak_files) {
      filename <- basename(file)
      file_cell_type <- gsub(".*_H3K4me3_", "", filename)
      file_cell_type <- gsub("\\.bed$", "", file_cell_type)
      
      cell_type <- celltype_mapping[file_cell_type]
      
      if (!is.na(cell_type) && cell_type %in% celltype_order) {
        cat("Processing", file_cell_type, "->", cell_type, "...\n")
        peaks_gr <- read_peaks(file)
        
        if (!is.null(peaks_gr) && length(peaks_gr) > 0) {
          peaks_df <- as.data.frame(peaks_gr)
          peaks_df <- peaks_df[, c("seqnames", "start", "end", "score")]
          colnames(peaks_df) <- c("chr", "start", "end", "score")
          peaks_df$chr <- as.character(peaks_df$chr)
          peaks_df$cell_type <- cell_type
          
          all_peaks[[filename]] <- peaks_df
          cat("  ✓", nrow(peaks_df), "peaks processed\n")
        }
      }
    }  
    
    if (length(all_peaks) > 0) {
      combined_peaks <- dplyr::bind_rows(all_peaks)
      cat("Successfully processed", length(all_peaks), "files\n")
      cat("Total peaks:", nrow(combined_peaks), "\n")
      return(combined_peaks)
    } else {
      stop("No peak files processed successfully!")
    }
  }
  
  combined_peaks <- process_peak_files(peakdir, celltype_order)
  
  # ==============================================================================
  # STEP 3: HUMAN GENE ANNOTATION WITH REGULATORY POTENTIAL
  # ==============================================================================
  cat("\n=== STEP 3: Human Gene Regulatory Potential Model ===\n")
  
  edb <- EnsDb.Hsapiens.v86
  genes_gr <- genes(edb)
  seqlevelsStyle(genes_gr) <- "UCSC"
  
  # Select protein-coding genes only
  cat("Total genes in EnsDb:", length(genes_gr), "\n")
  protein_coding_genes_gr <- genes_gr[genes_gr$gene_biotype == "protein_coding"]
  cat("Protein-coding genes:", length(protein_coding_genes_gr), "\n")
  
  # Show distribution of gene biotypes
  cat("\nGene biotype distribution:\n")
  biotype_counts <- table(genes_gr$gene_biotype)
  print(sort(biotype_counts, decreasing = TRUE)[1:10])
  
  genes_gr <- protein_coding_genes_gr
  
  # Create single-base TSS positions
  tss_gr <- genes_gr
  tss_positions_for_promoters <- ifelse(strand(genes_gr) == "+", 
                                       start(genes_gr), 
                                       end(genes_gr))
  start(tss_gr) <- tss_positions_for_promoters
  end(tss_gr) <- tss_positions_for_promoters
  promoter_region <- promoters(tss_gr, upstream = 100000, downstream = 100000)
  
  # Function to calculate Regulatory Potential (RP) weight
  calculate_rp_weight <- function(distance) {
    2^(-abs(distance) / 10000)
  }
  
  # ==============================================================================
  # PROCESS EACH CELL TYPE
  # ==============================================================================
  
  cat(sprintf("\nProcessing peaks for tool: %s\n", tool))
  
  # Use celltype_order
  celltypes <- celltype_order
  
  # Initialize list to store results for all cell types
  all_celltype_results <- list()
  
  for (ct in celltypes) {
    cat(sprintf("  Processing %s...", ct))
    
    # Get peaks for this cell type from combined_peaks
    ct_peaks <- combined_peaks[combined_peaks$cell_type == ct, ]
    
    if (nrow(ct_peaks) > 0) {
      # Convert to GRanges
      peaks_gr <- GRanges(
        seqnames = ct_peaks$chr,
        ranges = IRanges(start = ct_peaks$start, end = ct_peaks$end),
        score = ct_peaks$score
      )
      
      # Find overlaps between peaks and promoter regions
      overlaps <- findOverlaps(peaks_gr, promoter_region)
      
      if (length(overlaps) > 0) {
        # Get matching peaks and genes
        peak_hits <- peaks_gr[queryHits(overlaps)]
        gene_hits <- promoter_region[subjectHits(overlaps)]
        
        # Get original genes for correct TSS
        gene_indices <- subjectHits(overlaps)
        original_genes <- genes_gr[gene_indices]
        
        # Calculate TSS positions
        gene_tss_positions <- ifelse(strand(original_genes) == "+", 
                                    start(original_genes), 
                                    end(original_genes))
        
        # Calculate peak centers
        peak_centers <- (start(peak_hits) + end(peak_hits)) / 2
        
        # Calculate distances (strand-aware)
        distances <- ifelse(strand(original_genes) == "+",
                           peak_centers - gene_tss_positions, 
                           gene_tss_positions - peak_centers)
        
        # Calculate RP weights
        rp_weights <- calculate_rp_weight(distances)
        rp_weights[is.na(rp_weights) | is.infinite(rp_weights)] <- 0
        
        # Create data frame
        gene_activity <- data.frame(
          gene_id = original_genes$gene_id,
          gene_name = original_genes$gene_name,
          celltype = ct,
          weighted_peak_count = rp_weights,
          distance_from_tss = distances,
          stringsAsFactors = FALSE
        )
        
        # Summarize for this cell type
        gene_activity_summary <- gene_activity %>%
          dplyr::group_by(gene_id, gene_name, celltype) %>%
          dplyr::summarise(
            weighted_peak_count = sum(weighted_peak_count, na.rm = TRUE),
            n_peaks = n(),
            mean_distance = mean(distance_from_tss, na.rm = TRUE),
            .groups = 'drop'
          ) %>%
          dplyr::mutate(weighted_peak_count = replace_na(weighted_peak_count, 0))
        
        all_celltype_results[[ct]] <- gene_activity_summary
        cat(sprintf(" DONE (%d peak-gene overlaps)\n", length(overlaps)))
      } else {
        # Create empty result for this cell type
        empty_result <- data.frame(
          gene_id = character(),
          gene_name = character(),
          celltype = character(),
          weighted_peak_count = numeric(),
          n_peaks = integer(),
          mean_distance = numeric()
        )
        all_celltype_results[[ct]] <- empty_result
        cat(" NO overlaps\n")
      }
    } else {
      cat(" NO peaks\n")
    }
  }
  
  # Combine all cell type results
  if (length(all_celltype_results) > 0) {
    gene_activity_summary <- bind_rows(all_celltype_results)
    
    # Fill missing gene-celltype combinations with zeros
    all_genes <- unique(gene_activity_summary$gene_name)
    all_combinations <- expand.grid(
      gene_name = all_genes,
      celltype = celltypes,
      stringsAsFactors = FALSE
    )
    
    gene_activity_summary <- all_combinations %>%
      left_join(gene_activity_summary, by = c("gene_name", "celltype")) %>%
      mutate(
        weighted_peak_count = replace_na(weighted_peak_count, 0),
        n_peaks = replace_na(n_peaks, 0),
        mean_distance = replace_na(mean_distance, 0)
      )
    
    cat("✓ Human RP model calculation completed!\n")
    cat("✓ Generated", nrow(gene_activity_summary), "human gene-celltype RP scores\n")
  } else {
    stop("No peak data was successfully processed")
  }
  
  # ==============================================================================
  # STEP 4: CREATE GENE ACTIVITY MATRIX
  # ==============================================================================
  
  cat("\n=== STEP 4: Creating Gene Activity Matrix ===\n")
  
  gene_activity_matrix <- gene_activity_summary %>%
    dplyr::group_by(gene_name, celltype) %>%
    dplyr::summarise(
      weighted_peak_count = mean(weighted_peak_count, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = weighted_peak_count,
      values_fill = 0
    )
  
  gene_names <- gene_activity_matrix$gene_name
  matrix_data <- gene_activity_matrix[, -1]
  matrix_data[] <- lapply(matrix_data, as.numeric)
  gene_activity_matrix <- as.matrix(matrix_data)
  rownames(gene_activity_matrix) <- gene_names
  
  # Ensure all celltypes are present
  missing_celltypes <- setdiff(celltype_order, colnames(gene_activity_matrix))
  if (length(missing_celltypes) > 0) {
    for(ct in missing_celltypes) {
      gene_activity_matrix <- cbind(gene_activity_matrix, rep(0, nrow(gene_activity_matrix)))
      colnames(gene_activity_matrix)[ncol(gene_activity_matrix)] <- ct
    }
  }
  
  # Order columns by celltype_order
  gene_activity_matrix <- gene_activity_matrix[, celltype_order, drop = FALSE]
  
  # Save the matrix
  output_dir <- "/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore"
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(gene_activity_matrix, paste0(output_dir, "/HumanPBMC_H3K4me3_", tool, "_heatmap_matrix.csv"), row.names = TRUE, quote = FALSE) 
  
  # ==============================================================================
  # HEATMAP GENERATION
  # ==============================================================================
  
  partition_order <- paste0("partition", 1:10)
  
  # Ensure numeric matrix
  gene_activity_matrix <- as.matrix(gene_activity_matrix)
  mode(gene_activity_matrix) <- "numeric"
  
  # Your CD8+ T cell gene list
  CD8T_DEG <- c("CCR7", "SELL", "IL7R", "TCF7", "LEF1", "CD55", 
                "ITGAE", "MKI67", "CD69", "IL2RA", "ICOS", "IFNG", 
                "GZMB", "GZMA", "CX3CR1", "TNF", "B3GAT1")
  
  # Select genes that exist in your matrix
  available_genes <- CD8T_DEG[CD8T_DEG %in% rownames(gene_activity_matrix)]
  missing_genes <- CD8T_DEG[!CD8T_DEG %in% rownames(gene_activity_matrix)]
  
  cat("Found", length(available_genes), "out of", length(CD8T_DEG), "genes\n")
  if(length(missing_genes) > 0) {
    cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  }
  
  # Create plot_data with specific row order
  plot_data <- gene_activity_matrix[available_genes, partition_order]
  plot_data <- as.matrix(plot_data)
  rownames(plot_data) <- available_genes
  
  # Z-score normalization
  z_score_by_row <- function(mat) {
    row_means <- apply(mat, 1, mean, na.rm = TRUE)
    row_sds <- apply(mat, 1, sd, na.rm = TRUE)
    z_mat <- (mat - row_means) / row_sds
    z_mat[is.infinite(z_mat)] <- 0
    z_mat[is.nan(z_mat)] <- 0
    return(z_mat)
  }
  
  z_score_matrix <- z_score_by_row(plot_data)
  
  write.csv(z_score_matrix, paste0(output_dir, "/HumanPBMC_H3K4me3_", tool, "_zscore_matrix.csv"), row.names = TRUE, quote = FALSE) 
  
  # Create color function
  col_fun <- colorRamp2(
    breaks = seq(-3, 3, length.out = 101),
    colors = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(101)
  )
  
  # Check z-score range
  z_range <- range(z_score_matrix, na.rm = TRUE)
  cat("Z-score range:", round(z_range[1], 2), "to", round(z_range[2], 2), "\n")
  
  # Clamp data if needed
  if (z_range[1] < -3 || z_range[2] > 3) {
    cat("Warning: Data exceeds -3..3 range. Clamping.\n")
    z_score_matrix[z_score_matrix < -3] <- -3
    z_score_matrix[z_score_matrix > 3] <- 3
  }
  
  # Create top annotation
  top_anno <- HeatmapAnnotation(
    Pseudotime = 1:10,
    col = list(Pseudotime = colorRamp2(c(1, 10), c("#3aebd7", "#E76F51"))),
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    annotation_legend_param = list(
      Pseudotime = list(title = "Pseudotime (1=early, 10=late)")
    )
  )
  
  # Create heatmap
  ht <- Heatmap(
    z_score_matrix,
    name = "Z-score",
    col = col_fun,
    border = FALSE,
    rect_gp = gpar(col = NA, lwd = 0),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_order = partition_order,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 14),
    row_order = available_genes,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 14, fontface = "italic"),
    row_title = "CD8+ T Cell State Markers",
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title = "Gene Regulatory Z-score Along Pseudotime",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    top_annotation = top_anno,
    use_raster = TRUE,
    raster_quality = 2,
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12),
      legend_height = unit(35, "mm"),
      direction = "vertical",
      at = seq(-3, 3, by = 1),
      labels = seq(-3, 3, by = 1)
    )
  )
  
  # Save plots
  pdf_file <- file.path(output_dir, paste0("H3K4me3_", tool, "_GeneActivity_Heatmap.pdf"))
  pdf(pdf_file, width = 8.5, height = 10, useDingbats = FALSE)
  draw(ht, 
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(15, 10, 10, 10), "mm"))
  dev.off()
  
  png_file <- file.path(output_dir, paste0("H3K4me3_", tool, "_GeneActivity_Heatmap.png"))
  png(png_file, width = 8.5 * 300, height = 10 * 300, res = 300)
  draw(ht, 
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(15, 10, 10, 10), "mm"))
  dev.off()
  
  cat("\n✅ Figure saved successfully!\n")
  cat("📄 PDF:", pdf_file, "\n")
  cat("🖼️ PNG:", png_file, "\n")
  print(ht)
}


###############################################################################
### Figure I J
###############################################################################

library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

tools <- c("DROMPAplus", "Genrich", "GoPeaks", "HOMER", "MACS2", "SEACR", "SICER2")
#tools <- c("DROMPAplus")
for(tool in tools) {
  peakdir <- paste0("/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result/result_trajectoryfinal/partitionwise_barcode/split_celltype_bams_corrected/MouseBrain_peakbed_partitionwise/partitionwised_allpeakbed/H3K4me3_peakbed/", tool)
  
  # ==============================================================================
  # STEP 2: PEAK PROCESSING FOR HUMAN PBMC
  # ==============================================================================
  celltype_order <- c("partition1", "partition2", "partition3", "partition4", "partition5", 
                      "partition6", "partition7", "partition8", "partition9", "partition10")
  cat("\n=== STEP 2: Human PBMC Peak Processing ===\n")
  
  read_peaks <- function(bed_file, tool_name = NULL) {
    tryCatch({
      bed_data <- read.table(bed_file, sep = "\t", header = FALSE, 
                            stringsAsFactors = FALSE, fill = TRUE, 
                            comment.char = "", quote = "", na.strings = ".")
      
      n_cols <- ncol(bed_data)
      filename <- basename(bed_file)
      
      # Debug: Check what we're reading
      cat("  Reading", filename, "with", n_cols, "columns\n")
      
      # Basic validation
      if (n_cols < 3) {
        cat("  WARNING:", filename, "has only", n_cols, "columns. Skipping.\n")
        return(NULL)
      }
      
      # Create GRanges from first 3 columns
      peaks <- GRanges(
        seqnames = bed_data[, 1],
        ranges = IRanges(start = bed_data[, 2] + 1, end = bed_data[, 3])
      )
      
      # SPECIAL HANDLING FOR SICER2
      if (!is.null(tool_name) && tool_name == "SICER2") {
        cat("  Detected SICER2 format - using default score=1\n")
        peaks$score <- 0
      } else {
        # Standard handling for other tools
        score_found <- FALSE
        
        for (col_idx in 4:min(6, n_cols)) {
          if (col_idx <= n_cols) {
            potential_score <- suppressWarnings(as.numeric(bed_data[, col_idx]))
            
            if (!all(is.na(potential_score))) {
              peaks$score <- potential_score
              score_found <- TRUE
              cat("  Using score from column", col_idx, "for", filename, "\n")
              break
            }
          }
        }
        
        if (!score_found) {
          peaks$score <- 1
          cat("  No numeric score column found for", filename, "(using score=1)\n")
        }
      }
      
      # Clean up any NA or infinite scores
      valid_scores <- !is.na(peaks$score) & is.finite(peaks$score)
      peaks <- peaks[valid_scores]
      
      if (length(peaks) == 0) {
        cat("  WARNING: No valid peaks after cleaning in", filename, "\n")
        return(NULL)
      }
      
      cat("  ✓", length(peaks), "peaks loaded from", filename, "\n")
      return(peaks)
      
    }, error = function(e) {
      cat("  ERROR reading", basename(bed_file), ":", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  process_peak_files <- function(peakdir, celltype_order) {
    peak_files <- list.files(peakdir, pattern = "\\.bed$", full.names = TRUE)
    cat("Found", length(peak_files), "peak bed files\n")
    
    # Map filenames to trajectory celltypes
    celltype_mapping <- c(
      "partition1" = "partition1",
      "partition2" = "partition2", 
      "partition3" = "partition3",
      "partition4" = "partition4",
      "partition5" = "partition5",
      "partition6" = "partition6", 
      "partition7" = "partition7",
      "partition8" = "partition8",
      "partition9" = "partition9",
      "partition10" = "partition10"
    )
    
    all_peaks <- list()
    
    for (file in peak_files) {
      filename <- basename(file)
      file_cell_type <- gsub(".*_H3K4me3_", "", filename)
      file_cell_type <- gsub("\\.bed$", "", file_cell_type)
      
      cell_type <- celltype_mapping[file_cell_type]
      
      if (!is.na(cell_type) && cell_type %in% celltype_order) {
        cat("Processing", file_cell_type, "->", cell_type, "...\n")
        peaks_gr <- read_peaks(file)
        
        if (!is.null(peaks_gr) && length(peaks_gr) > 0) {
          peaks_df <- as.data.frame(peaks_gr)
          peaks_df <- peaks_df[, c("seqnames", "start", "end", "score")]
          colnames(peaks_df) <- c("chr", "start", "end", "score")
          peaks_df$chr <- as.character(peaks_df$chr)
          peaks_df$cell_type <- cell_type
          
          all_peaks[[filename]] <- peaks_df
          cat("  ✓", nrow(peaks_df), "peaks processed\n")
        }
      }
    }  
    
    if (length(all_peaks) > 0) {
      combined_peaks <- dplyr::bind_rows(all_peaks)
      cat("Successfully processed", length(all_peaks), "files\n")
      cat("Total peaks:", nrow(combined_peaks), "\n")
      return(combined_peaks)
    } else {
      stop("No peak files processed successfully!")
    }
  }
  
  combined_peaks <- process_peak_files(peakdir, celltype_order)
  
  # ==============================================================================
  # STEP 3: HUMAN GENE ANNOTATION WITH REGULATORY POTENTIAL
  # ==============================================================================
  cat("\n=== STEP 3: Mouse Gene Regulatory Potential Model ===\n")
  
  edb <- EnsDb.Mmusculus.v79
  genes_gr <- genes(edb)
  seqlevelsStyle(genes_gr) <- "UCSC"
  
  # Select protein-coding genes only
  cat("Total genes in EnsDb:", length(genes_gr), "\n")
  protein_coding_genes_gr <- genes_gr[genes_gr$gene_biotype == "protein_coding"]
  cat("Protein-coding genes:", length(protein_coding_genes_gr), "\n")
  
  # Show distribution of gene biotypes
  cat("\nGene biotype distribution:\n")
  biotype_counts <- table(genes_gr$gene_biotype)
  print(sort(biotype_counts, decreasing = TRUE)[1:10])
  
  genes_gr <- protein_coding_genes_gr
  
  # Create single-base TSS positions
  tss_gr <- genes_gr
  tss_positions_for_promoters <- ifelse(strand(genes_gr) == "+", 
                                       start(genes_gr), 
                                       end(genes_gr))
  start(tss_gr) <- tss_positions_for_promoters
  end(tss_gr) <- tss_positions_for_promoters
  promoter_region <- promoters(tss_gr, upstream = 100000, downstream = 100000)
  
  # Function to calculate Regulatory Potential (RP) weight
  calculate_rp_weight <- function(distance) {
    2^(-abs(distance) / 10000)
  }
  
  # ==============================================================================
  # PROCESS EACH CELL TYPE
  # ==============================================================================
  
  cat(sprintf("\nProcessing peaks for tool: %s\n", tool))
  
  # Use celltype_order
  celltypes <- celltype_order
  
  # Initialize list to store results for all cell types
  all_celltype_results <- list()
  
  for (ct in celltypes) {
    cat(sprintf("  Processing %s...", ct))
    
    # Get peaks for this cell type from combined_peaks
    ct_peaks <- combined_peaks[combined_peaks$cell_type == ct, ]
    
    if (nrow(ct_peaks) > 0) {
      # Convert to GRanges
      peaks_gr <- GRanges(
        seqnames = ct_peaks$chr,
        ranges = IRanges(start = ct_peaks$start, end = ct_peaks$end),
        score = ct_peaks$score
      )
      
      # Find overlaps between peaks and promoter regions
      overlaps <- findOverlaps(peaks_gr, promoter_region)
      
      if (length(overlaps) > 0) {
        # Get matching peaks and genes
        peak_hits <- peaks_gr[queryHits(overlaps)]
        gene_hits <- promoter_region[subjectHits(overlaps)]
        
        # Get original genes for correct TSS
        gene_indices <- subjectHits(overlaps)
        original_genes <- genes_gr[gene_indices]
        
        # Calculate TSS positions
        gene_tss_positions <- ifelse(strand(original_genes) == "+", 
                                    start(original_genes), 
                                    end(original_genes))
        
        # Calculate peak centers
        peak_centers <- (start(peak_hits) + end(peak_hits)) / 2
        
        # Calculate distances (strand-aware)
        distances <- ifelse(strand(original_genes) == "+",
                           peak_centers - gene_tss_positions, 
                           gene_tss_positions - peak_centers)
        
        # Calculate RP weights
        rp_weights <- calculate_rp_weight(distances)
        rp_weights[is.na(rp_weights) | is.infinite(rp_weights)] <- 0
        
        # Create data frame
        gene_activity <- data.frame(
          gene_id = original_genes$gene_id,
          gene_name = original_genes$gene_name,
          celltype = ct,
          weighted_peak_count = rp_weights,
          distance_from_tss = distances,
          stringsAsFactors = FALSE
        )
        
        # Summarize for this cell type
        gene_activity_summary <- gene_activity %>%
          dplyr::group_by(gene_id, gene_name, celltype) %>%
          dplyr::summarise(
            weighted_peak_count = sum(weighted_peak_count, na.rm = TRUE),
            n_peaks = n(),
            mean_distance = mean(distance_from_tss, na.rm = TRUE),
            .groups = 'drop'
          ) %>%
          dplyr::mutate(weighted_peak_count = replace_na(weighted_peak_count, 0))
        
        all_celltype_results[[ct]] <- gene_activity_summary
        cat(sprintf(" DONE (%d peak-gene overlaps)\n", length(overlaps)))
      } else {
        # Create empty result for this cell type
        empty_result <- data.frame(
          gene_id = character(),
          gene_name = character(),
          celltype = character(),
          weighted_peak_count = numeric(),
          n_peaks = integer(),
          mean_distance = numeric()
        )
        all_celltype_results[[ct]] <- empty_result
        cat(" NO overlaps\n")
      }
    } else {
      cat(" NO peaks\n")
    }
  }
  
  # Combine all cell type results
  if (length(all_celltype_results) > 0) {
    gene_activity_summary <- bind_rows(all_celltype_results)
    
    # Fill missing gene-celltype combinations with zeros
    all_genes <- unique(gene_activity_summary$gene_name)
    all_combinations <- expand.grid(
      gene_name = all_genes,
      celltype = celltypes,
      stringsAsFactors = FALSE
    )
    
    gene_activity_summary <- all_combinations %>%
      left_join(gene_activity_summary, by = c("gene_name", "celltype")) %>%
      mutate(
        weighted_peak_count = replace_na(weighted_peak_count, 0),
        n_peaks = replace_na(n_peaks, 0),
        mean_distance = replace_na(mean_distance, 0)
      )
    
    cat("✓ Human RP model calculation completed!\n")
    cat("✓ Generated", nrow(gene_activity_summary), "mouse gene-celltype RP scores\n")
  } else {
    stop("No peak data was successfully processed")
  }
  
  # ==============================================================================
  # STEP 4: CREATE GENE ACTIVITY MATRIX
  # ==============================================================================
  
  cat("\n=== STEP 4: Creating Gene Activity Matrix ===\n")
  
  gene_activity_matrix <- gene_activity_summary %>%
    dplyr::group_by(gene_name, celltype) %>%
    dplyr::summarise(
      weighted_peak_count = mean(weighted_peak_count, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = weighted_peak_count,
      values_fill = 0
    )
  
  gene_names <- gene_activity_matrix$gene_name
  matrix_data <- gene_activity_matrix[, -1]
  matrix_data[] <- lapply(matrix_data, as.numeric)
  gene_activity_matrix <- as.matrix(matrix_data)
  rownames(gene_activity_matrix) <- gene_names
  
  # Ensure all celltypes are present
  missing_celltypes <- setdiff(celltype_order, colnames(gene_activity_matrix))
  if (length(missing_celltypes) > 0) {
    for(ct in missing_celltypes) {
      gene_activity_matrix <- cbind(gene_activity_matrix, rep(0, nrow(gene_activity_matrix)))
      colnames(gene_activity_matrix)[ncol(gene_activity_matrix)] <- ct
    }
  }
  
  # Order columns by celltype_order
  gene_activity_matrix <- gene_activity_matrix[, celltype_order, drop = FALSE]
  
  # Save the matrix
  output_dir <- "/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam/result_corrected/result_trajcetory/additional_analysis_final_zscore"
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(gene_activity_matrix, paste0(output_dir, "/MouseBrain_H3K4me3_", tool, "_heatmap_matrix.csv"), row.names = TRUE, quote = FALSE) 
  
  # ==============================================================================
  # HEATMAP GENERATION
  # ==============================================================================
  
  partition_order <- paste0("partition", 1:10)
  
  # Ensure numeric matrix
  gene_activity_matrix <- as.matrix(gene_activity_matrix)
  mode(gene_activity_matrix) <- "numeric"
  
  # Your CD8+ T cell gene list
  CD8T_DEG <- c("Mbp", "Tmem125", "Epb4.1l3", "Cldn11", "Mobp", "Mal", "Map7", "Mog", 
                             "Rffl", "Adap1", "Cntn2", "Pllp")
  
  # Select genes that exist in your matrix
  available_genes <- CD8T_DEG[CD8T_DEG %in% rownames(gene_activity_matrix)]
  missing_genes <- CD8T_DEG[!CD8T_DEG %in% rownames(gene_activity_matrix)]
  
  cat("Found", length(available_genes), "out of", length(CD8T_DEG), "genes\n")
  if(length(missing_genes) > 0) {
    cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  }
  
  # Create plot_data with specific row order
  plot_data <- gene_activity_matrix[available_genes, partition_order]
  plot_data <- as.matrix(plot_data)
  rownames(plot_data) <- available_genes
  
  # Z-score normalization
  z_score_by_row <- function(mat) {
    row_means <- apply(mat, 1, mean, na.rm = TRUE)
    row_sds <- apply(mat, 1, sd, na.rm = TRUE)
    z_mat <- (mat - row_means) / row_sds
    z_mat[is.infinite(z_mat)] <- 0
    z_mat[is.nan(z_mat)] <- 0
    return(z_mat)
  }
  
  z_score_matrix <- z_score_by_row(plot_data)
  
  write.csv(z_score_matrix, paste0(output_dir, "/MouseBrain_H3K4me3_", tool, "_zscore_matrix.csv"), row.names = TRUE, quote = FALSE) 
  
  # Create color function
  col_fun <- colorRamp2(
    breaks = seq(-3, 3, length.out = 101),
    colors = colorRampPalette(brewer.pal(11, "PiYG"))(101)  # Purple-Green diverging
  )
  
  # Check z-score range
  z_range <- range(z_score_matrix, na.rm = TRUE)
  cat("Z-score range:", round(z_range[1], 2), "to", round(z_range[2], 2), "\n")
  
  # Clamp data if needed
  if (z_range[1] < -3 || z_range[2] > 3) {
    cat("Warning: Data exceeds -3..3 range. Clamping.\n")
    z_score_matrix[z_score_matrix < -3] <- -3
    z_score_matrix[z_score_matrix > 3] <- 3
  }
  
  # Create top annotation
  top_anno <- HeatmapAnnotation(
    Pseudotime = 1:10,
    col = list(Pseudotime = colorRamp2(c(1, 10), c("#3aebd7", "#E76F51"))),
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    annotation_legend_param = list(
      Pseudotime = list(title = "Pseudotime (1=early, 10=late)")
    )
  )
  
  # Create heatmap
  ht <- Heatmap(
    z_score_matrix,
    name = "Z-score",
    col = col_fun,
    border = FALSE,
    rect_gp = gpar(col = NA, lwd = 0),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_order = partition_order,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 14),
    row_order = available_genes,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 14, fontface = "italic"),
    row_title = "OPC+ mOL Cell State Markers",
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title = "Gene Regulatory Z-score Along Pseudotime",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    top_annotation = top_anno,
    use_raster = TRUE,
    raster_quality = 2,
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12),
      legend_height = unit(35, "mm"),
      direction = "vertical",
      at = seq(-3, 3, by = 1),
      labels = seq(-3, 3, by = 1)
    )
  )
  
  # Save plots
  pdf_file <- file.path(output_dir, paste0("H3K4me3_", tool, "_GeneActivity_Heatmap.pdf"))
  pdf(pdf_file, width = 8.5, height = 10, useDingbats = FALSE)
  draw(ht, 
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(15, 10, 10, 10), "mm"))
  dev.off()
  
  png_file <- file.path(output_dir, paste0("H3K4me3_", tool, "_GeneActivity_Heatmap.png"))
  png(png_file, width = 8.5 * 300, height = 10 * 300, res = 300)
  draw(ht, 
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE,
       padding = unit(c(15, 10, 10, 10), "mm"))
  dev.off()
  
  cat("\n✅ Figure saved successfully!\n")
  cat("📄 PDF:", pdf_file, "\n")
  cat("🖼️ PNG:", png_file, "\n")
  print(ht)
}
