# ------------------------------
# --- Functions to help pre-processing of Seurat Objects
# ------------------------------


# --- Function to apply nFeature (UMI) and Mitochondrial percentage thresholds
subset_seurat <- function(seurat_obj, nFeature_lower, nFeature_upper, percent_mt) {
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > nFeature_lower & 
                                            nFeature_RNA < nFeature_upper & 
                                            percent.mt < percent_mt)
  
  return(seurat_obj)
}

# --- Function to reformat WILDseq barcodes to metadata
barcodes_to_metadata <- function(seurat_obj) {
  # Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  # Extract noramlized WildSeq "expression" given that WildSeq barcode is coded as another gene in the ref genome
  barcode <- Seurat::FetchData(object = seurat_obj, vars = "WILDseq")
  # Remove the WildSeq barcode from feature matrix to prevent the barcode from contributing the clustering algorithms (i.e. PC, UMAP, etc)
  seurat_obj <- subset(seurat_obj, features = row.names(seurat_obj)[-which(row.names(seurat_obj) == "WILDseq")])
  # Store extracted WildSeq barcode in metadata instead of feature matrix
  seurat_obj <- Seurat::AddMetaData(object = seurat_obj, metadata = barcode, col.name = c("barcode_transcript"))
  
  return(seurat_obj)
}

# --- Function to score cell cycle
score_cell_cycle <- function(seurat_obj) {
  # Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  # Extract genes active during S phase
  s.genes <- Seurat::cc.genes$s.genes
  # Extract genes active during G2/m phase
  g2m.genes <- Seurat::cc.genes$g2m.genes
  # Calculate cell cycle phase scores
  seurat_obj <- Seurat::CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  return(seurat_obj)
}

add_wildseq_metadata <- function(seurat_obj, barcode_table_path) {
  # Read Wildseq barcode assignments
  barcode.table <- readr::read_tsv(barcode_table_path)
  barcode.table <- as.data.frame(barcode.table[, c("CBC", "clone", "origin")])
  barcode.table <- tibble::column_to_rownames(barcode.table, var = "CBC")

  # Add metadata to Seurat object
  seurat_obj <- Seurat::AddMetaData(
    object = seurat_obj,
    metadata = barcode.table,
    col.name = c("clone", "identification.method")
  )

  # Fill in NA values
  if (any(is.na(seurat_obj@meta.data$clone))) {
    seurat_obj@meta.data$clone[is.na(seurat_obj@meta.data$clone)] <- "Unknown"
  }
  if (any(is.na(seurat_obj@meta.data$identification.method))) {
    seurat_obj@meta.data$identification.method[is.na(
      seurat_obj@meta.data$identification.method
    )] <- "None"
  }

  return(seurat_obj)
}