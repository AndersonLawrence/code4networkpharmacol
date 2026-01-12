# ==============================================================================
# 0. Environment Setup & Dependencies
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Install Seurat v5 and core bioinformatics packages
pacman::p_load(
  Seurat, SeuratWrappers, Matrix, dplyr, patchwork, qs, 
  harmony, ggplot2, AUCell, CellChat, remotes, devtools, tibble
)

# Install Monocle3 and dependencies if not present
if (!requireNamespace("monocle3", quietly = TRUE)) {
  if (!requireNamespace("leidenbase", quietly = TRUE)) {
    devtools::install_github("cole-trapnell-lab/leidenbase")
  }
  devtools::install_github("cole-trapnell-lab/monocle3")
}

set.seed(123)

# ==============================================================================
# 1. Configuration & I/O
# ==============================================================================

# Update working directory path as needed
# setwd("C:/Users/whisk/Desktop/AS/GSE159677")

# Create output directories
dir.create("Results/QC", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Clustering", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Expression", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/CellChat", recursive = TRUE, showWarnings = FALSE)

sample_folders <- c("GSM4837523", "GSM4837524", "GSM4837525", 
                    "GSM4837526", "GSM4837527", "GSM4837528")

sample_map <- c(
  "GSM4837523" = "AC", "GSM4837524" = "PA", "GSM4837525" = "AC",
  "GSM4837526" = "PA", "GSM4837527" = "AC", "GSM4837528" = "PA"
)

# ==============================================================================
# 2. Data Loading & Merging (Seurat v5)
# ==============================================================================

seurat_list <- lapply(sample_folders, function(sf) {
  message("Loading Sample: ", sf)
  counts <- Read10X(data.dir = sf)
  obj <- CreateSeuratObject(
    counts = counts, 
    project = sf, 
    min.cells = 3, 
    min.features = 200
  )
  # Add sample prefix to barcodes to prevent duplication
  obj <- RenameCells(obj, add.cell.id = sf)
  return(obj)
})

# Merge samples into a single object
sce.all <- Reduce(function(a, b) merge(a, y = b), seurat_list)

# Assign metadata based on sample ID
cell_prefix <- sub("_.*$", "", colnames(sce.all))
sce.all$sample         <- cell_prefix
sce.all$orig.ident     <- cell_prefix
sce.all$classification <- unname(sample_map[cell_prefix])

# Calculate mitochondrial percentage
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")

# Save raw merged object
qsave(sce.all, "Results/sce_raw.qs")

# ==============================================================================
# 3. Quality Control (QC)
# ==============================================================================

# Visualize QC metrics before filtering
pdf("Results/QC/QC_Pre_Filtering.pdf", width = 12, height = 6)
VlnPlot(sce.all, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        group.by = "orig.ident", pt.size = 0, ncol = 3) + 
  plot_annotation(title = "Pre-QC Metrics")
dev.off()

# Apply filtering criteria
sce.filt <- subset(sce.all, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Visualize QC metrics after filtering
pdf("Results/QC/QC_Post_Filtering.pdf", width = 12, height = 6)
VlnPlot(sce.filt, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        group.by = "orig.ident", pt.size = 0, ncol = 3) + 
  plot_annotation(title = "Post-QC Metrics")
dev.off()

# ==============================================================================
# 4. Normalization, Integration, and Dimensionality Reduction
# ==============================================================================

DefaultAssay(sce.filt) <- "RNA"
sce.filt <- NormalizeData(sce.filt)
sce.filt <- FindVariableFeatures(sce.filt, selection.method = "vst", nfeatures = 2000)
sce.filt <- ScaleData(sce.filt)
sce.filt <- RunPCA(sce.filt, npcs = 30, verbose = FALSE)

# Integration using Harmony to correct batch effects
sce.filt <- RunHarmony(sce.filt, group.by.vars = "orig.ident", dims.use = 1:30)

# Clustering and Visualization
sce.filt <- FindNeighbors(sce.filt, reduction = "harmony", dims = 1:20)
sce.filt <- FindClusters(sce.filt, resolution = 0.5)
sce.filt <- RunUMAP(sce.filt, reduction = "harmony", dims = 1:20)
sce.filt <- RunTSNE(sce.filt, reduction = "harmony", dims = 1:20)

# Save visualization plots
pdf("Results/Clustering/UMAP_tSNE_Overview.pdf", width = 12, height = 10)
(DimPlot(sce.filt, reduction = "umap", group.by = "orig.ident") + ggtitle("UMAP by Sample")) +
  (DimPlot(sce.filt, reduction = "umap", group.by = "classification") + ggtitle("UMAP by Class")) /
  (DimPlot(sce.filt, reduction = "tsne", group.by = "orig.ident") + ggtitle("tSNE by Sample")) +
  (DimPlot(sce.filt, reduction = "tsne", group.by = "classification") + ggtitle("tSNE by Class"))
dev.off()

# ==============================================================================
# 5. Cell Type Annotation (AUCell)
# ==============================================================================

marker_list <- list(
  Macrophages   = c("AIF1", "CD14"),
  Endothelial   = c("PECAM1", "ECSCR", "CALD1"),
  VSMCs         = c("MYL9", "TAGLN"),
  NK_T          = c("NKG7", "XCL1", "CTSW"),
  T_Lymphocytes = c("CD2", "TRAC", "CD69"),
  B_Lymphocytes = c("CD79A", "MS4A1", "IGKC")
)

# Run AUCell scoring
cells_rankings <- AUCell_buildRankings(GetAssayData(sce.filt, slot = "data"))
cells_AUC      <- AUCell_calcAUC(marker_list, cells_rankings)
sce.filt$celltype <- apply(getAUC(cells_AUC), 2, function(x) names(which.max(x)))
Idents(sce.filt)  <- "celltype"

# Plot annotated clusters
pdf("Results/Clustering/CellType_Annotation.pdf", width = 10, height = 8)
DimPlot(sce.filt, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("Cell Type Annotation (AUCell)")
dev.off()

# ==============================================================================
# 6. Gene Expression Visualization
# ==============================================================================

# FeaturePlot for specific genes
pdf("Results/Expression/FeaturePlot_Markers.pdf", width = 12, height = 4)
FeaturePlot(sce.filt, features = c("TNF", "RXRA", "MMP9"), reduction = "umap", 
            order = TRUE, min.cutoff = "q10", max.cutoff = "q90", ncol = 3)
dev.off()

# DotPlot for broad markers
pdf("Results/Expression/DotPlot_Markers.pdf", width = 10, height = 6)
DotPlot(sce.filt, features = unique(unlist(marker_list)), group.by = "celltype") + 
  RotatedAxis() + scale_color_viridis_c() + ggtitle("Marker Gene Expression")
dev.off()

# ==============================================================================
# 7. Cell-Cell Communication (CellChat)
# ==============================================================================

# Downsample for computational efficiency if necessary
sce.comm <- if (ncol(sce.filt) > 15000) {
  subset(sce.filt, cells = sample(Cells(sce.filt), 15000))
} else sce.filt

# Initialize CellChat
cellchat <- createCellChat(object = GetAssayData(sce.comm, slot = "data"), 
                           meta = sce.comm@meta.data, group.by = "celltype")
cellchat@DB <- CellChatDB.human # Use Human Database
cellchat <- subsetData(cellchat) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(raw.use = FALSE) %>% 
  filterCommunication(min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

# Circle Plots: Interaction Count and Strength
mat.count  <- as.matrix(cellchat@net$count)
mat.weight <- as.matrix(cellchat@net$weight)
groupSize  <- as.numeric(table(cellchat@idents))

pdf("Results/CellChat/Interaction_Network_Circle.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
netVisual_circle(mat.count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Number of Interactions")
netVisual_circle(mat.weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Interaction Strength")
dev.off()

# Signaling Heatmap
pdf("Results/CellChat/Interaction_Heatmap.pdf", width = 10, height = 8)
netVisual_heatmap(cellchat, signaling = NULL, color.heatmap = c("white", "#EA3F35"),
                  cluster.rows = TRUE, cluster.cols = TRUE)
dev.off()

# Final Save
qsave(sce.filt, "Results/sce_final_annotated.qs")