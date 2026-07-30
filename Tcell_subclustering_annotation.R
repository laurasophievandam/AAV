# ==============================================================================
# T-cell subclustering and annotation
# ==============================================================================
# Starting from the merged, clustered GEX+CITE object, after splitting T- and 
# B-lineage clusters,(GEX_CITE_processing.R)
# re-processes and re-clusters the T-cell subset, removes residual
# low-quality/contaminating clusters, assigns final cell-type labels and merges 
# with TCR.
# ==============================================================================

INPUT_RDS    <- "<path/to/merged_clustered.rds>"
SAMPLE_SHEET <- "<path/to/sample_sheet.csv>"  # same sample sheet as GEX_CITE_processing.R,
                                               # extended with a `tcr_dir` column (path to
                                               # each sample's filtered_contig_annotations.csv,
                                               # relative to DATA_ROOT)
DATA_ROOT    <- "<path/to/raw/data>"
OUTPUT_DIR   <- "<path/to/output/directory>"

library(Seurat)
library(dplyr)
library(scRepertoire)

merged_T <- readRDS(INPUT_RDS)

# ==============================================================================
# 1. Re-process and cluster the T-cell subset
# ==============================================================================
merged_T <- NormalizeData(merged_T, normalization.method = "LogNormalize", scale.factor = 10000)
merged_T <- FindVariableFeatures(merged_T, selection.method = "vst", nfeatures = 3000)
merged_T <- ScaleData(merged_T)
merged_T <- RunPCA(merged_T, features = VariableFeatures(merged_T), npcs = 50)

merged_T <- RunUMAP(merged_T, reduction = "harmony", dims = 1:20)
merged_T <- FindNeighbors(merged_T, reduction = "harmony", dims = 1:20)
merged_T <- FindClusters(merged_T, reduction = "harmony", resolution = 1)

saveRDS(merged_T, file = file.path(OUTPUT_DIR, "Tcells_clustered.rds"))

# ==============================================================================
# 3. QC pass: remove residual non-T contaminating cells
# ==============================================================================
# Removes cells with detectable CD19 expression (B-lineage contamination not
# captured by the cluster-level QC pass above).
merged_T <- subset(merged_T, CD19 > 0, invert = TRUE)

# ==============================================================================
# 4. Cell-type annotation
# ==============================================================================
DefaultAssay(merged_T) <- "RNA"
Idents(merged_T) <- "seurat_clusters"

new_cluster_ids <- c(
  "CD4 Naive", "CD4 CM", "CD4 Naive", "CD4/CD8 EM CCR6+",
  "CD8 EMRA GZMB+/NKT like", "CD4 EM", "CD4/8 CM", "CD8 EMRA GZMB/GZMK+",
  "CD8 Naive", "CD4 CM CD69+", "Activated T cells CD44+", "CD4 Treg",
  "CD8 EM GZMK+", "CD4 TFH", "NKT CD1D+", "CD8 EM ITGA4+", "NKT",
  "CD4 CTLA4+", "CD4 Naive", "CD4 CTLA4+"
)
names(new_cluster_ids) <- levels(merged_T)
merged_T <- RenameIdents(merged_T, new_cluster_ids)

cluster_order <- c(
  "CD4 Naive", "CD8 Naive", "CD4 CM CD69+", "CD4 CTLA4+", "CD4 CM",
  "CD4/8 CM", "CD4 EM", "CD4/CD8 EM CCR6+", "CD4 Treg", "CD4 TFH", "NKT",
  "CD8 EM GZMK+", "CD8 EMRA GZMB/GZMK+", "CD8 EMRA GZMB+/NKT like",
  "NKT CD1D+", "CD8 EM ITGA4+", "Activated T cells CD44+"
)
Idents(merged_T) <- factor(Idents(merged_T), levels = cluster_order)
merged_T$annotated_clusters <- Idents(merged_T)

saveRDS(merged_T, file = file.path(OUTPUT_DIR, "Tcells_annotated.rds"))

# ==============================================================================
# 5. TCR integration
# ==============================================================================
# Loads each sample's filtered_contig_annotations.csv (standard CellRanger VDJ
# output), merges them into one clonotype-level object, and attaches
# clonotype calls to the annotated T-cell object as new metadata columns.
sample_sheet <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE)
stopifnot("tcr_dir" %in% colnames(sample_sheet))
sample_sheet <- sample_sheet[!is.na(sample_sheet$tcr_dir) & nzchar(sample_sheet$tcr_dir), ]

contig_list <- lapply(sample_sheet$tcr_dir, function(d) {
  df <- read.csv(file.path(DATA_ROOT, d))
  # some contig files already carry a "sample" column from an earlier
  # processing step, which conflicts with combineTCR()'s own sample labeling
  if ("sample" %in% colnames(df)) df <- select(df, -sample)
  df
})
names(contig_list) <- sample_sheet$sample_id

combined_tcr <- combineTCR(contig_list,
                           samples = sample_sheet$sample_id,
                           removeNA = TRUE,
                           removeMulti = FALSE,
                           filterMulti = TRUE)

combined_tcr <- addVariable(combined_tcr, name = "patient_id",     variables = sample_sheet$patient_id)
combined_tcr <- addVariable(combined_tcr, name = "group",          variables = sample_sheet$group)
combined_tcr <- addVariable(combined_tcr, name = "disease_status", variables = sample_sheet$disease_status)

saveRDS(combined_tcr, file = file.path(OUTPUT_DIR, "Tcells_TCR_combined.rds"))

merged_T <- combineExpression(combined_tcr, merged_T,
                              cloneCall = "gene",
                              group.by = "sample_id",
                              proportion = FALSE,
                              cloneTypes = c(Single = 1, Small = 5, Medium = 20,
                                             Large = 100, Hyperexpanded = 500,
                                             Ultraexpanded = 2000))

saveRDS(merged_T, file = file.path(OUTPUT_DIR, "Tcells_annotated_TCR.rds"))
