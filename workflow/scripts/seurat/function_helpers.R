# -------- Function to visualize principal components

visualize_principal_components <- function(
  pca_object,
  pca_method = "PCA",
  y_variable_to_plot = "variance",
  num_PCs = 10
) {
  # If PCA method is normal PCA from prcomp
  if (pca_method == "PCA") {
    # Extract eigenvalues from the PCA object
    eigenvalues <- pca_object$sdev^2
  } else if (pca_method == "ordPCA") {
    # Extract eigenvalues from the PCA object
    eigenvalues <- pca_object$pca$sdev^2
  } else if (pca_method == "Seurat"){
    eigenvalues <- pca_object@reductions$pca@stdev^2
  } else {
    stop("Unknown PCA method: ", pca_method)
  }

  # Compute variance table
  pca_plot <- data.frame(
    PC = paste0("PC", 1:length(eigenvalues)),
    eigenvalue = eigenvalues,
    variance = eigenvalues / sum(eigenvalues)
  ) |>
    dplyr::mutate(cumvar = cumsum(variance) * 100)

  # Dynamic y-axis label
  y_label <- dplyr::case_when(
    y_variable_to_plot == "cumvar" ~ "Cumulative Variance Explained (%)",
    y_variable_to_plot == "variance" ~ "Variance Explained",
    y_variable_to_plot == "eigenvalue" ~ "Eigenvalues",
    TRUE ~ "Value")

  # Generate Scree plot
  ggplot2::ggplot(
    pca_plot,
    ggplot2::aes(
      x = as.numeric(gsub("PC", "", PC)),
      y = .data[[y_variable_to_plot]]
    )
  ) +
    ggplot2::geom_line(group = 1, color = "black") +
    ggplot2::geom_point(color = "black") +
    ggplot2::labs(
      x = "Principal Components",
      y = y_label
    ) +
    ggplot2::theme_classic()

  # Print Variance table
  # print(head(pca_plot, n = num_PCs))
}


# --- Function to extract and collate cell attributes for SLX ID

# extract_cell_attributes <- function(assigned_seurats_file_path_list, SLX_ID) {
#   # Intialize cell attributes storage list
#   all_cell_attributes <- list()

#   # Loop over all the assigned Seurat files
#   for (file in assigned_seurats_file_path_list) {
#     # Extract sample index from the filename
#     sample_index <- sub(
#       paste0(".*", SLX_ID, "_(.*)_all_assigned_seurat\\.rds"),
#       "\\1",
#       file
#     )

#     # Read Seurat object
#     seurat <- readRDS(file)

#     # Extract cell-level metadata & add sample_index
#     cell_df <- seurat@meta.data |>
#       as.data.frame() |>
#       dplyr::mutate(sample_index = sample_index)

#     # Store it in list
#     all_cell_attributes[[sample_index]] <- cell_df
#   }

#   # Collate cell attributes across all samples
#   all_cell_attributes <- dplyr::bind_rows(all_cell_attributes)

#   # Return cell attributes as collated dataframe
#   return(all_cell_attributes)
# }


# --- Function to extract and collate cell and gene attributes for SLX ID

extract_cell_and_gene_attributes <- function(assigned_seurats_file_path_list, SLX_ID) {
  
  # Initialize storage lists
  all_cell_attributes <- list()
  all_gene_attributes <- list()
  all_residual_attributes <- list()
  
  # Loop over all Seurat files
  for (file in assigned_seurats_file_path_list) {
    # Extract sample index from the filename
    sample_index <- sub(
      paste0(".*", SLX_ID, "_(.*)_all_assigned_seurat\\.rds"),
      "\\1",
      file
    )

    # Load Seurat object
    seurat <- readRDS(file)

    ## -----------------
    ## CELL-LEVEL METADATA
    ## -----------------

    cell_df <- seurat@meta.data |>
      as.data.frame() |>
      dplyr::mutate(sample_index = sample_index)

    all_cell_attributes[[sample_index]] <- cell_df

    ## -----------------
    ## GENE-LEVEL METADATA
    ## -----------------

    counts_matrix <- seurat@assays$RNA@layers$counts

    gene_df <- data.frame(
      mean = rowMeans(counts_matrix),
      detection_rate = rowMeans(counts_matrix > 0),
      var = matrixStats::rowVars(as.matrix(counts_matrix))
    ) |>
      dplyr::mutate(
        log_mean = log10(mean),
        log_var = log10(var),
        sample_index = sample_index
      )

    all_gene_attributes[[sample_index]] <- gene_df

    ## -----------------
    ## RESIDUALS
    ## -----------------

    residual_df <- seurat@assays$SCT@SCTModel.list$counts@feature.attributes |>
    as.data.frame() |>
    dplyr::mutate(sample_index = sample_index)

    all_residual_attributes[[sample_index]] <- residual_df
  }
  
  # Combine all samples into single dataframes
  cell_attributes_all <- dplyr::bind_rows(all_cell_attributes)
  gene_attributes_all <- dplyr::bind_rows(all_gene_attributes)
  residual_attributes_all <- dplyr::bind_rows(all_residual_attributes)
  
  # Return as a list
  return(list(
    cell_attributes = cell_attributes_all,
    gene_attributes = gene_attributes_all,
    residuals = residual_attributes_all
  ))
}

# --- Function to explore and visualize the relationship between library size and number of genes detected

visualize_library_size_detected_genes_relationship <- function(collated_cell_attributes_dataframe, SLX_ID, normalized = TRUE) {

  if (normalized == TRUE) {
    p <- ggplot(collated_cell_attributes_dataframe, aes(nCount_SCT, nFeature_SCT)) +
          geom_point(alpha = 0.3, shape = 16) +
          geom_density_2d(size = 0.3) +
          facet_wrap(~ sample_index, scales = "free") +
          theme_classic() +
          labs(title = paste("Library size vs gene detection after SCTransform (", SLX_ID, ")", sep = ""),
          x = "Library Size", y = "Detected genes")
  } else if (normalized == FALSE) {
    p <- ggplot(collated_cell_attributes_dataframe, aes(nCount_RNA, nFeature_RNA)) +
          geom_point(alpha = 0.3, shape = 16) +
          geom_density_2d(size = 0.3) +
          facet_wrap(~ sample_index, scales = "free") +
          theme_classic() +
          labs(title = paste("Library size vs gene detection prior to SCTransform (", SLX_ID, ")", sep = ""),
              x = "Library Size", y = "Number of genes detected")
  } else {
    stop("Please select a normalized argument, either TRUE or FALSE.")
  }

  return(p)

}

# --- Function to explore and visualize the relationship between mean and variance

visualize_mean_variance_relationship <- function(gene_attributes_df, SLX_ID) {
  ggplot(gene_attributes_df, aes(log_mean, log_var)) +
    geom_point(alpha = 0.3, shape = 16) +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    facet_wrap(~ sample_index, scales = "free") +
    labs(
      title = paste("Mean–variance relationship (", SLX_ID, ")", sep = ""),
      x = "log10(mean expression)",
      y = "log10(variance)"
    ) +
    theme_classic()
}

# --- Function to explore and visualize mean-detection-rate relationship 

visualize_mean_detection_rate_relationship <- function(gene_attributes_df, SLX_ID) {
  
  # Poisson model curve
  x_vals <- seq(from = -3, to = 2, length.out = 1000)
  poisson_model <- data.frame(
    log_mean = x_vals,
    detection_rate = 1 - dpois(0, lambda = 10^x_vals)
  )
  
  ggplot(gene_attributes_df, aes(log_mean, detection_rate)) +
    geom_point(alpha = 0.3, shape = 16) +
    geom_line(data = poisson_model, aes(log_mean, detection_rate), color = 'red') +
    facet_wrap(~ sample_index, scales = "free") +
    theme_classic() +
    labs(
      title = paste("Mean–detection rate vs Poisson model (", SLX_ID, ")", sep = ""),
      x = "log10(mean expression)",
      y = "Detection rate"
    )
}

# --- Function to visualize residual variance vs geometric mean for Homoskedasticity

visualize_residual_variance_vs_mean <- function(residual_df, SLX_ID) {
  ggplot(residual_df, aes(x = log(gmean), y = residual_variance)) +
    geom_point(alpha = 0.3, shape = 16) +
    facet_wrap(~ sample_index, scales = "free") +
    xlim(0, 10) +
    theme_classic() +
    labs(
      title = paste("Residual variance vs geometric mean (", SLX_ID, ")", sep = ""),
      x = "log(geometric mean)",
      y = "Residual variance"
    )
}


# --- Function to Visualize WILDseq barcodes on UMAP embeddings for a specific Sample Index (SITT*)

plot_top_barcodes_umap <- function(SLX_ID, sample_index, top_n = 6, seed = 241202) {
  
  # Build path to Seurat object
  seurat_path <- file.path("results", SLX_ID, sample_index, "clone_assignment",
                           paste0(SLX_ID, "_", sample_index, "_all_assigned_seurat.rds"))
  
  if (!file.exists(seurat_path)) {
    stop("Seurat object not found: ", seurat_path)
  }
  message("Loading Seurat object: ", seurat_path)
  
  # Load assigned Seurat object
  seurat <- readRDS(seurat_path)
  
  # Identify the top 2000 most variable features
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  # Perform Principal Component Analysis linear dimensionality reduction
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), verbose = FALSE)
  
  # Run dimensionality reduction
  seurat <- RunUMAP(seurat, dims = 1:30, n.neighbors = 15, min.dist = 0.01,
                    seed.use = seed, verbose = FALSE)
  
  # Extract the Top 6 Barcodes Indexes (Under the assumption that the First highest proportion is Unknown)
  top_barcodes <- names(sort(round(prop.table(table(seurat$clone)) * 100, 2), decreasing = TRUE))[1:top_n]
  
  # Print top_n barcodes + proportions
  message("Top barcodes for ", sample_index, ":")
  print(sort(round(prop.table(table(seurat$clone)) * 100, 2), decreasing = TRUE))

  # Define color palette 
  color_palette <- c("#CCCCCC", "#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")[1:top_n]
  
  # Visualize WILDsesq barcodes on UMAP embeddings
  seurat_top <- subset(seurat, subset = clone %in% top_barcodes)

  # Visualize WILDsesq barcodes on UMAP embeddings
  p <- DimPlot(subset(seurat, subset = clone %in% top_barcodes), reduction = "umap", group.by = "clone", order = rev(top_barcodes), cols = color_palette) +
  labs(title = paste0(sample_index, " (", SLX_ID, ")"), x = "UMAP1", y = "UMAP2")

  return(p)
}

# --- Function to save WILDseq barcodes on UMAP embeddings for a specified number of sample indices

save_multiple_top_barcodes_umap <- function(SLX_ID, sample_indices, top_n = 6) {
  
  # Define output directory for saving
  outdir <- file.path("results", SLX_ID, "seurat_analyses", "Barcode_UMAPs")
  
  # If output directory does not exist, create the output directory
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # Loop across all defined sample indices
  for (sample in sample_indices){
  message("Processing sample: ", sample)
  
  # Generate Plot for Sample 
  p <- plot_top_barcodes_umap(SLX_ID = SLX_ID, sample_index = sample, top_n = top_n)

  # Save plot
  ggsave(file.path(outdir, paste0(SLX_ID, "_", sample, "_top_clones_umap.pdf")), plot = p, width = 8, height = 6) 
}
  message("Plots are saved under ", outdir)
}

