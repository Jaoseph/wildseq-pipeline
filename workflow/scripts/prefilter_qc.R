#!/usr/bin/env Rscript
# library(tidyverse)
library(argparse)


#parse data
parser <- ArgumentParser(description= 'prefilter qc plot')
parser$add_argument('--slx_id', help = "SLX_ID")
parser$add_argument('--sample_index', help = 'sample_index')
parser$add_argument("--filtered_feature_bc_matrix")
xargs<- parser$parse_args()

gp <- function(...) print(glue::glue(...))

outdir = gp("results/{xargs$slx_id}/{xargs$sample_index}/prefilter_qc_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


#import matrix and make Seurat object
gp("{xargs$slx_id} {xargs$sample}: importing data")
seurat <- Seurat::Read10X(data.dir = xargs$filtered_feature_bc_matrix)
seurat <- Seurat::CreateSeuratObject(counts = seurat, min.cells = 5)
seurat$sample = xargs$sample
seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^MT-")


plot_title = gp("{xargs$sample_index} prefilter")

gp("{xargs$slx_id} {xargs$sample}: plotting percent.mt")
filename = gp("{outdir}/{xargs$sample_index}_featureScatter_percent.mt_prefilter.png")
png(filename, height=6,width=6,units="in",res=300)
print(Seurat::FeatureScatter(seurat, 
        feature1 = "nCount_RNA", feature2 = "percent.mt", 
        pt.size=0.1, group.by = "sample", cols = scales::alpha("black",0.5)) + 
    ggplot2::theme(legend.position="none") +
    ggplot2::labs(title=plot_title))
dev.off()
        
gp("{xargs$slx_id} {xargs$sample}: plotting nFeature_RNA")
filename = gp("{outdir}/{xargs$sample_index}_featureScatter_nFeature_RNA_prefilter.png")
png(filename, height=6,width=6,units="in",res=300)
print(Seurat::FeatureScatter(seurat, 
        feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
        pt.size=0.1, group.by = "sample", cols = scales::alpha("black",0.5)) + 
    ggplot2::theme(legend.position="none") + 
    ggplot2::labs(title=plot_title)) 
dev.off()









# # ------------ SINGLE RUN ------------ #


# #!/usr/bin/env Rscript
# #library(tidyverse)
# library(argparse)

# gp <- function(...) print(glue::glue(...))

# seurat <- Seurat::Read10X(data.dir = "results/SLX-27045/SITTE11/cellranger_count/outs/filtered_feature_bc_matrix")
# seurat <- Seurat::CreateSeuratObject(counts = seurat)
# seurat$sample = "SITTE11"
# seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^MT-")

# quantile(seurat$nFeature_RNA, 0.05)
#   5% 
# 1337 
# quantile(seurat$nFeature_RNA, 0.95)
#    95% 
# 9380.2 
# quantile(seurat$percent.mt, 0.95)
#      95% 
# 8.493693 


# plot_title = gp("SITTC1 prefilter")


# pdf("demo_results_wiki/lib_size_vs_nfeature_SITTC1_scatterplot.pdf", w=15,h=10)
# FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     ggplot2::theme(legend.position="none") +
#     ggplot2::labs(title="SITTC1 Library size (nCount_RNA) vs Number of genes detected (nFeature_RNA)")
# dev.off()

# pdf("demo_results_wiki/lib_size_vs_nfeature_SITTC1_scatterplot_with_thresholds.pdf", w=15,h=10)
# FeatureScatter(seurat,
#         feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
#         pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     theme(legend.position="none") +
#     labs(title="SITTC1 Library size (nCount_RNA) vs Number of genes detected (nFeature_RNA)") +
#     geom_hline(yintercept = quantile(seurat$nFeature_RNA, 0.05), linetype = "dashed", color = "red") +
#     geom_hline(yintercept = quantile(seurat$nFeature_RNA, 0.95), linetype = "dashed", color = "red")
# dev.off()

# pdf("demo_results_wiki/lib_size_vs_percent_mt_SITTC1_scatterplot.pdf", w=15,h=10)
# FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     ggplot2::theme(legend.position="none") +
#     ggplot2::labs(title="SITTC1 Library size (nCount_RNA) vs percentage of Mitochondrial transcripts (percent.mt)")
# dev.off()

# pdf("demo_results_wiki/lib_size_vs_percent_mt_SITTC1_scatterplot_with_thresholds.pdf", w=15,h=10)
# FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     theme(legend.position="none") +
#     labs(title="SITTC1 Library size (nCount_RNA) vs percentage of Mitochondrial transcripts (percent.mt)") +
#     geom_hline(yintercept = quantile(seurat$percent.mt, 0.95), linetype = "dashed", color = "red")
# dev.off()

# # ----

# png("demo_results_wiki/lib_size_vs_nfeature_SITTC1_scatterplot.png", width = 1500, height = 1000, res = 150, units = "px")
# FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     ggplot2::theme(legend.position="none") +
#     ggplot2::labs(title="SITTC1 Library size (nCount_RNA) vs Number of genes detected (nFeature_RNA)")
# dev.off()

# png("demo_results_wiki/lib_size_vs_nfeature_SITTC1_scatterplot_with_thresholds.png", width = 1500, height = 1000, res = 150, units = "px")
# print(Seurat::FeatureScatter(seurat,
#         feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
#         pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     ggplot2::theme(legend.position="none") +
#     ggplot2::labs(title="SITTC1 Library size (nCount_RNA) vs Number of genes detected (nFeature_RNA)") +
#     ggplot2::geom_hline(yintercept = quantile(seurat$nFeature_RNA, 0.05), linetype = "dashed", color = "red") +
#     ggplot2::geom_hline(yintercept = quantile(seurat$nFeature_RNA, 0.95), linetype = "dashed", color = "red"))
# dev.off()

# png("demo_results_wiki/lib_size_vs_percent_mt_SITTC1_scatterplot.png", width = 1500, height = 1000, res = 150, units = "px")
# FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     ggplot2::theme(legend.position="none") +
#     ggplot2::labs(title="SITTC1 Library size (nCount_RNA) vs percentage of Mitochondrial transcripts (percent.mt)")
# dev.off()

# png("demo_results_wiki/lib_size_vs_percent_mt_SITTC1_scatterplot_with_thresholds.png", width = 1500, height = 1000, res = 150, units = "px")
# FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.1, cols = scales::alpha("black",0.5)) + 
#     ggplot2::theme(legend.position="none") +
#     ggplot2::labs(title="SITTC1 Library size (nCount_RNA) vs percentage of Mitochondrial transcripts (percent.mt)") +
#     ggplot2::geom_hline(yintercept = quantile(seurat$percent.mt, 0.95), linetype = "dashed", color = "red")
# dev.off()
