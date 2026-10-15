#!/usr/bin/env Rscript

# --- General Set up ---

# Load Libraries
library(argparse)
library(Seurat)
library(ggplot2)

# Source helper functions
source("workflow/scripts/seurat/function_helpers.R")

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("--slx_id", required = TRUE)
parser$add_argument("--output_dir", required = TRUE)
xargs <- parser$parse_args()

# Extract SLX-ID  & output path
SLX_ID <- xargs$slx_id
output_dir <- xargs$output_dir

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set seed
set.seed(241202)

# Find QC filtered and SCTransformed Seurat Objects for the requested SLX-ID
assigned_seurats <- list.files(
  path = file.path("results", SLX_ID),
  pattern = "_all_assigned_seurat\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(assigned_seurats) == 0) {
  stop(paste0("No assigned Seurat RDS files found for ", SLX_ID, ". Please run clone_assignment.smk"))
}

message("Found ", length(assigned_seurats), " samples for ", SLX_ID)

# --- QC various inherent mathematical properties of the SCTransformed data ---

# Derive cell and gene attributes
attributes_list <- extract_cell_and_gene_attributes(assigned_seurats, SLX_ID)

# Extract cell attributes
cell_attributes <- attributes_list$cell_attributes

# 1) Visualize relationship between library size (nCount_RNA) and number of genes detected (nFeature_RNA) before SCTransform normalization
before_normalization <- visualize_library_size_detected_genes_relationship(cell_attributes, SLX_ID, normalized = FALSE)
ggsave(filename = file.path(output_dir, paste0(SLX_ID, "_combined_lib_size_detectedgenes_relationship_before_norm.pdf")), plot = before_normalization, width = 16, height = 10)

# 2) Visualize relationship between library size (nCount_RNA) and number of genes detected (nFeature_RNA) after SCTransform normalization
after_normalization  <- visualize_library_size_detected_genes_relationship(cell_attributes, SLX_ID, normalized = TRUE)
ggsave(filename = file.path(output_dir, paste0(SLX_ID, "_combined_lib_size_detectedgenes_relationship_after_norm.pdf")), plot = after_normalization, width = 16, height = 10)

# Extract gene attributes
gene_attributes <- attributes_list$gene_attributes

# 3) Visualize mean–variance relationship
mean_variance_plot <- visualize_mean_variance_relationship(gene_attributes,SLX_ID)
ggsave(filename = file.path(output_dir, paste0(SLX_ID, "_combined_mean_variance_relationship_before_norm.pdf")), plot = mean_variance_plot, width = 16, height = 10)

# 4) Visualize mean-detection-rate relationship 
mean_detection_plot <- visualize_mean_detection_rate_relationship(gene_attributes, SLX_ID)
ggsave(filename = file.path(output_dir, paste0(SLX_ID, "_combined_mean_detection_vs_poisson_before_norm.pdf")), plot = mean_detection_plot, width = 16,height = 10)

# Extract residuals
residual_attributes <- attributes_list$residuals

# 5) Visualize residual variance vs geometric mean (homoskedasticity)
residual_plot <- visualize_residual_variance_vs_mean(residual_attributes, SLX_ID)
ggsave(filename = file.path(output_dir, paste0(SLX_ID, "_combined_residual_variance_vs_geomean_after_norm.pdf")), plot = residual_plot, width = 16, height = 10)


message("QC plots are saved for ", SLX_ID)



# # --- Identification of highly variable features (feature selection) ---

# seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# # Identify the top 10 most highly variable genes
# head(VariableFeatures(seurat), 10)

# # Plot Variable Features with and without labels
# LabelPoints(plot = VariableFeaturePlot(seurat, cols = c("#D5D5D3", "#B40F20")) + labs(y = "Variance", x = "Mean Expression"),
#   points = head(VariableFeatures(seurat), 10), repel = TRUE)

# # --- Dimensionality Reduction ---

# # Perform Principal Component Analysis linear dimensionality reduction
# seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), verbose = FALSE)

# # Examine PCA results in numerous ways

# # Visualize PCs using Screeplot
# visualize_principal_components(seurat, pca_method = "Seurat", y_variable_to_plot = "cumvar", num_PCs = 30)

# # Examine what features (genes) is associated with which PCs (loadings)
# print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

# # Visualize the features (loadings) associated with each PCs
# VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

# # Visualize the PCs profiles for each cell on Scatterplot
# DimPlot(seurat, reduction = "pca") # not sure why this is coloring by cell cycle?

# # Visualize feature expression for each cell using a Heatmap
# DimHeatmap(seurat, dims = 1:10, cells = 500, balanced = TRUE)

# # Apply a graph-based clustering approach (K-nearest neighbour graph)

# # Construct KNN graph based on euclidean distacne in the PCA space
# seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)

# # Apply Louvain optimization method to cluster
# seurat <- FindClusters(seurat, resolution = 0.5, verbose = FALSE)

# # Run dimensionality reduction (suggest using the same PCs as input to clustering analysis)
# seurat <- RunUMAP(seurat, dims = 1:30, n.neighbors = 15, min.dist = 0.01, seed.use = 241202, verbose = FALSE)  # Sharp well-separated clusters, for tight groupings

# # Visualize single cell clustering on UMAP embeddings
# DimPlot(seurat, reduction = "umap", label = TRUE)

# # Define color palette 
# color_palette <- c(
#   "#CCCCCC","#001219", "#003049", "#005F73", "#0A9396", "#1B998B",
#   "#94D2BD", "#E9D8A6", "#FFB703", "#EE9B00", "#CA6702", 
#   "#BB3E03", "#ae2f12", "#AE2012", "#9B2226", "#780000", 
#   "#000000")

# # Compute proportions per barcode
# barcode_proportions <- sort(round(prop.table(table(seurat$clone)) * 100, 2), decreasing = FALSE)
# barcode_order <- names(barcode_proportions)

# # Visualize single cell clustering on UMAP embeddings with Wildseq barcodes "populations"
# pdf("results/SLX-26635/SITTC1/seurat_analyses/wildseq_barcodes_populations_umap_overlay.pdf", w = 10, h = 10)
# DimPlot(seurat, reduction = "umap", group.by = "clone", cols = color_palette, order = barcode_order)
# dev.off()

# # Plot only the top 10 barcodes
# top_barcodes <- names(sort(round(prop.table(table(seurat$clone)) * 100, 2), decreasing = TRUE))[1:10]
# top_barcodes_order <- rev(top_barcodes)
# color_ten_palette <- c(
#   "#808080","#77030b","#5d4b20","#469374","#9341b3",
#   "#e3427d","#e68653","#ebe0b0","#edfbba","#303638")
# DimPlot(subset(seurat, subset = clone %in% top_barcodes_order), reduction = "umap", group.by = "clone", order = top_barcodes_order, cols = color_ten_palette)













