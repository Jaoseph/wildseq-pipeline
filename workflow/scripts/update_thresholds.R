library(Seurat)

samples <- c("SITTE11", "SITTF11", "SITTG11", "SITTH11")
base_path <- "/scratche/slow/jblab/chong02/wildseq/results/SLX-28281"
metadata_path <- "/scratche/slow/jblab/chong02/wildseq/data/metadata/sample_metadata_SLX-28281.csv"

seurat_list <- lapply(samples, function(s) {
  data <- Seurat::Read10X(data.dir = file.path(base_path, s, "cellranger_count/outs/filtered_feature_bc_matrix"))
  obj <- Seurat::CreateSeuratObject(counts = data)
  obj$sample <- s
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})
names(seurat_list) <- samples

for (i in seq_along(samples)) {
  obj <- seurat_list[[i]]
  cat("\n===", samples[i], "===\n")
  cat("nFeature_RNA 5%:", quantile(obj$nFeature_RNA, 0.05), "\n")
  cat("nFeature_RNA 95%:", quantile(obj$nFeature_RNA, 0.95), "\n")
  cat("percent.mt 95%:", quantile(obj$percent.mt, 0.95), "\n")
}

metadata <- readr::read_csv(metadata_path)
thresholds <- data.frame(
  sample_index = samples,
  lower = sapply(seq_along(samples), function(i) round(quantile(seurat_list[[i]]$nFeature_RNA, 0.05))),
  upper = sapply(seq_along(samples), function(i) round(quantile(seurat_list[[i]]$nFeature_RNA, 0.95))),
  mt    = sapply(seq_along(samples), function(i) round(quantile(seurat_list[[i]]$percent.mt,  0.95), 2))
)
metadata <- metadata |>
  dplyr::select(sample, sample_index, enrich_index) |>
  dplyr::left_join(thresholds, by = "sample_index")

readr::write_csv(metadata, metadata_path)
cat("\nUpdated metadata:\n")
print(metadata)
