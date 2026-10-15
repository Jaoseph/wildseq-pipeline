#!/usr/bin/env Rscript
library(argparse)
library(Seurat)
library(tidyverse)
library(glmGamPoi)
library(future)
plan("multisession", workers = 4) 
options(future.globals.maxSize = 20000 * 1024^2)

# parse data
parser <- argparse::ArgumentParser(description = 'clone_assignment')
parser$add_argument('--slx_id', help = "SLX_ID")
parser$add_argument('--sample', help = 'sample_index')
parser$add_argument('--filtered_feature_bc_matrix', help = 'filtered_feature_bc_matrix path')
parser$add_argument('--barcode_assigned_table', help = 'barcode_assigned_table.txt')
parser$add_argument('--filter_nFeatureRNA_lower', type = "integer", help = 'Lower threshold for nFeature_RNA')
parser$add_argument('--filter_nFeatureRNA_upper', type = "integer", help = 'Upper threshold for nFeature_RNA')
parser$add_argument('--filter_percent_mt', type = "double", help = 'Mitochondrial percentage cutoff')
parser$add_argument('--full_seurat', help = 'Output path for full Seurat object')
parser$add_argument('--all_assigned_seurat', help = 'Output path for assigned Seurat object')
parser$add_argument('--cell_numbers_table', help = 'Output path for cell numbers table')
xargs <- parser$parse_args()

# Set seed
set.seed(8290)

# Write function to print title
gp <- function(...) print(glue::glue(...))

# # If directory to store plots does not exist, create the directory
# if (!dir.exists(xargs$plot_dir)) {
#   dir.create(xargs$plot_dir, recursive = TRUE)
# }

# --- Mitochondria and nFEATURERNA quality control pre-processing

gp("{xargs$slx_id} {xargs$sample}: importing data")
# Import cellranger output matrix
seurat <- Seurat::Read10X(data.dir = xargs$filtered_feature_bc_matrix)
# Create Seurat object
seurat <- Seurat::CreateSeuratObject(counts = seurat, min.cells = 5)
# Add sample index for Seurate Object
seurat$sample <- xargs$sample
# Add SLX index for Seurate Object
seurat$SLX_ID <- xargs$slx_id

# Compute the percentage of transcripts that map to mitochondrial genes
seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^MT-")

# Save the full (unfiltered) Seurat object
saveRDS(seurat, xargs$full_seurat)

# apply mitochondria and nfeature (UMI) thresholds
gp("{xargs$slx_id} {xargs$sample}: filtering")
seurat <- subset(seurat, subset = nFeature_RNA > xargs$filter_nFeatureRNA_lower & 
                                  nFeature_RNA < xargs$filter_nFeatureRNA_upper & 
                                  percent.mt < xargs$filter_percent_mt)

# --- Prevent WILDseq barcodes from influencing clustering by storing in Seurat object metadata

# Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor
seurat <- Seurat::NormalizeData(seurat)
# Extract noramlized WildSeq "expression" given that WildSeq barcode is coded as another gene in the ref genome
barcode <- Seurat::FetchData(object = seurat, vars = "WILDseq")
# Remove the WildSeq barcode from feature matrix to prevent the barcode from contributing the clustering algorithms (i.e. PC, UMAP, etc)
seurat <- subset(seurat, features = row.names(seurat)[-which(row.names(seurat) == "WILDseq")])
# Store extracted WildSeq barcode in metadata instead of feature matrix
seurat <- Seurat::AddMetaData(object = seurat, metadata = barcode, col.name=c("barcode_transcript"))

# --- Calculate cell cycle phase scores based on canonical markers to mitigate effects of cell cycle heterogeneity in scRNA-seq by regress these out of the data during pre-processing

# Cell cycle scoring
gp("{xargs$slx_id} {xargs$sample}: cell cycle scoring")
# Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor
seurat <- Seurat::NormalizeData(seurat)
# Extract genes active during S phase
s.genes <- Seurat::cc.genes$s.genes
# Extract genes active during G2/m phase
g2m.genes <- Seurat::cc.genes$g2m.genes
# Calculate cell cycle phase scores
seurat <- Seurat::CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


# --- Fit a Negative Binomial Model to perform normalization and variance stabilization of scRNA-seq datasets 

gp("{xargs$slx_id} {xargs$sample}: SCT transforming")
# Normalize raw counts using a Negative Binominal Model and then correct for Cell cycle, precent.mt, nFeature_RNA
# NB: if we statistically remove the contribution of percent.mt and nFeature_RNA, then automated quantile at 0.05 or 0.95 will work to simply remove outliers.
seurat <- Seurat::SCTransform(seurat, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt","nFeature_RNA"), 
  verbose = TRUE, conserve.memory = TRUE, ncells = ncol(seurat))

# --- Collate Wildseq Barcode Information with Seurat results

gp("{xargs$slx_id} {xargs$sample}: assigning barcodes")
# Read Wildseq Barcodes extracted from PCR enrichemnt?
barcode.table <- readr::read_tsv(xargs$barcode_assigned_table)
# Extract Cell Barcode, WildSeq Barcode (Clone), and where the barcode was extracted (10x or PCR-enrichment) 
barcode.table <- as.data.frame(barcode.table[,c("CBC","clone", "origin")])
# Set Rownames as Cell Barcode (each cell should have an unique cell barcode)
barcode.table <- barcode.table |>
  tibble::column_to_rownames(var = "CBC")
# Add Wildseq barcode information to metadata in Seurat Object
seurat <- Seurat::AddMetaData(object = seurat, metadata = barcode.table, col.name=c("clone","identification.method"))
# If the cell does not belong in any clone (aka. NA), set as clone as unknown
if (any(is.na(seurat@meta.data$clone))) {
  seurat@meta.data$clone[is.na(seurat@meta.data$clone)] <- "Unknown"
}
# If the cell does not belong in any clone (aka. NA), set as origin as none
if (any(is.na(seurat@meta.data$identification.method))) {
  seurat@meta.data$identification.method[is.na(seurat@meta.data$identification.method)] <- "None"
}

# Save Seurat object with barcode information
saveRDS(seurat, xargs$all_assigned_seurat)

# --- Examine clonal diversity

gp("{xargs$slx_id} {xargs$sample}: determining clones")
# Compute the frequency of each unique clone
cell_numbers <- as.data.frame(table(seurat$clone))
# Rename columns
colnames(cell_numbers) <- c("clone", "count")
# Compute the percent of each unique column 
cell_numbers = cell_numbers |>
  dplyr::filter(!clone %in% c("Unknown", "non-match")) |>
  dplyr::mutate(pct = count*100/sum(count))

# Save clonal diversity table and percentage
readr::write_tsv(cell_numbers,xargs$cell_numbers_table)

