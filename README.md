# WILDSeq Singlecell-RNAseq Barcoding Pipeline

## Description

### *Summary*

Perform lineage tracing and clonal analyses using the WILDseq barcoding approach ["Clonal transcriptomics identifies mechanisms of chemoresistance and empowers rational design of combination therapies"](https://doi.org/10.7554/eLife.80981)

### *Detailed description* ###

This pipeline processes raw fastq files from CRUK-CI (with the possibility of starting from outside of the CRUK-CI sequencing file nomenclature using standard illumina file naming by omitting the file conversion rule `crukci_to_illumina`). Standard 10X CellRanger is used to generate gene expression matrices. WILDSeq barcodes are extracted and assigned through a comprehensive dual-validation approach combining transriptiomic analysis and PCR enrichment data.


## Pipeline set up

### Step 1 Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository

### Step 2 Ensure you have the correct environment

It is important to note that all the required packages and functions needed for both pre-processing and downstream analyses can be obtained through a conda environment unless stated otherwise below. 

Firstly, make sure that either mamba or conda are installed on your system. The recommendation is installing miniforge which is a optimized lighter version of conda. Please see[miniforge](https://github.com/conda-forge/miniforge) for installation guide.

A conda.yaml is available at the root of the directory and can be installed by running 

```
conda env create -f conda.yaml
```
### Step 3 Installing softwares

Certain softwares will need to be manually installed such as [Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in) and Bowtie into the software/ directorty. Currently, cellranger v9.0.1 and bowtie v1.3.0 are used. Important that it is Bowtie not Bowtie2.

### Step 4 Obtaining Theoretical barcodes...

### Step 5 Preparing metadata

#### Metadata sheet

The pipeline requires generate a metadata table in the data/metadata/ folder named sample_metadata_SLX-XXXXX.csv. This metadata should contain the following columns: `sample` (sample name), `sample_index` (scRNA run index), `enrich_index` (PCR run index), `lower` (lower nFeature threshold), `upper` (upper nFeature threshold), and `mt` (mitocondrial percentage threshold).

When you run **Step 1** of the pipeline, you should leave upper, lower and mt blank as currently, the pipeline is ran under manual user review of the prefilter_QC metrics.

The table below is an example of the csv table

|sample|sample_index|enrich_index|lower|upper          |mt
|----------|---------|--------|------|--------------|----------|
|CIOV123456_BASELINE      |SITTE15     |SITTA1    |   BLANK    | BLANK         |BLANK          |
|CIOV123456_16H_RUCAP     |SITTE16     |SITTB1    |   BLANK    | BLANK         |BLANK          |


##### main.yaml 

The main.yaml (`config/main.yaml`) contains the necessary information for the pipeline which you wish to run. This includes the SLX-ID and the path to the metadata generated (see above) which includes the sample and enrich indices with QC thresholds.


### Workflow management

The workflow profile configuration file (`workflow/profiles/config.v8+.yaml`) has the necessary information to configure job submissions developed for SLURM work load manager under snakemake. 

Currently, this configuration file is adapted for JBLAB members for the CRUK CI HPC cluster. Under `set-resources`, users can specify the cpus and memory required for their specific jobs. Please note that rules `cellranger_count`, `enrich_pcr`, and `clone_assignment` should be already optimized and should not be touched unless explicitly shown an error regarding too little `mem_mb`. If so, increase as necessary, by manually editting the YAML file.

### Step 6

Quite simply, have your raw fastq file data be located in the `data/raw_data/SLX-XXXX` directory. Each SLX-runs should be under a different directory in `raw_data/`.


## Running the Pipeline
### Step 7, Stage 1

Once you have the metadata set up, the conda environment set up as well as the necessary softwares, the first stage is ready to run. At the root of the directory (with your conda environment activated), run 

```
sbatch submit_slurm_step1.sh 
```

to SLURM workload manager (CRUK CI) or if you are manually running soley using snakemake and not SLURM, run 

```
snakemake --workflow-profile workflow/profiles/ \
    --use-conda -p -k --snakefile Snakefile \
    --until barcode_assignment prefilter_qc 
```

***IMPORTANT:*** If you are not downloading the sequencing data via CRUK CI clarity tools and/or already have te fastq files with the correct cellranger naming convention (if the name does not match the required name for cellranger, it will not run), then you can go to `Snakefile` and comment out `include: "workflow/rules/crukci_to_illumina.smk"`



### Step 7, Stage 2

At the end of Stage 1, you will have two plots under `results/{SLX-XXXXX}/{Sample_Index}/prefilter_qc_plots/`: **{SAMPLE_INDEX}_featureScatter_nFeature_RNA_prefilter.png** and **{SAMPLE_INDEX}_featureScatter_percent.mt_prefilter.png**.

Using these plots you can determine empirically the thresholds for both nFeature and percent mitochondria. 
If you wish, you can load `results/{SLX-XXXXX}/{Sample_Index}/cellranger_count/outs/filtered_feature_bc_matrix` and run 

```
seurat <- Seurat::Read10X(data.dir = "results/{SLX-XXXXX}/{Sample_Index}/cellranger_count/outs/filtered_feature_bc_matrix")
seurat <- Seurat::CreateSeuratObject(counts = seurat)
seurat$sample = "{Sample_Index}"
seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^MT-")

quantile(seurat$nFeature_RNA, 0.05)
#   5% 
# 1337 
quantile(seurat$nFeature_RNA, 0.95)
#    95% 
# 9380.2 
quantile(seurat$percent.mt, 0.95)
#      95% 
# 8.493693 

```
to determine the thresholds by quantile. An example script is available and located at the bottom of the prefilter_qc.R script.

Following this, you will need to manually update the metadata sheet in the columns `upper`, `lower`, `mt` corresponding to your chosen thresholds. 

Finally, once you have selected and updated the thresholds, at the root of the directory (with your conda environment activated), run 

```
sbatch submit_slurm_step2.sh 
```

to SLURM workload manager (CRUK CI) or if you are manually running soley using snakemake and not SLURM, run 

```
snakemake --workflow-profile workflow/profiles/ \
    --use-conda -p -k --snakefile Snakefile \
    --allowed-rules clone_assignment
```

## Output files

The final output files generated are located in the `results/` directory under its own respective SLX-ID and indices. 

- `results/{SLX_ID}/{sample_index}/clone_assignment/{SLX_ID}_{sample_index}_cellNumbers.tsv`: tab-separated file containing  your barcodes counts and proportion in the scRNA-seq run
- `results/{SLX_ID}/{sample_index}/clone_assignment/{SLX_ID}_{sample_index}_all_assigned_seurat.rds`: Seurat Object containing both your raw and SCTransformed expression matrices, along with metadata and most importantly your WILDSeq barcodes which can be found in seurat@meta.data$clones.

## Authors

* Jao Lok Chong 