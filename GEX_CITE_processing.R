# ==============================================================================
# Single-cell RNA-seq + CITE-seq processing pipeline: AAV samples
# ==============================================================================
# Per-sample loading, QC filtering, optional CITE-seq (ADT) attachment,
# merging across samples, batch correction (Harmony), and clustering.
#
# Expects one CellRanger GEX output (filtered_feature_bc_matrix) per sample,
# and optionally a matching CITE-seq antibody-derived-tag (ADT) UMI count
# directory (e.g. from CITE-seq-Count). Sample paths and metadata are supplied
# via an external sample sheet (see below).
# ==============================================================================

# ---- Configuration: edit these three paths for your own data ---------------
DATA_ROOT    <- "<path/to/raw/data>"          # root directory containing all sample subfolders
SAMPLE_SHEET <- "<path/to/sample_sheet.csv>"  # see "Sample sheet" section below
OUTPUT_DIR   <- "<path/to/output/directory>"  # where .rds checkpoints are written

library(Seurat)
library(Matrix)
library(harmony)

# ---- Sample sheet ------------------------------------------------------------
# CSV with one row per sample and the following columns:
#   sample_id       - unique sample identifier
#   gex_dir         - path to the CellRanger filtered_feature_bc_matrix directory,
#                      relative to DATA_ROOT
#   adt_dir         - path to the CITE-seq ADT UMI count directory, relative to
#                      DATA_ROOT; leave blank if the sample has no ADT data
#   patient_id      - patient/subject identifier
#   group           - study group (AAV/HC)
#   timepoint       - sample timepoint label
#   ANCA            - HC/MPO/PR3 label
sample_sheet <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE)
required_cols <- c("sample_id", "gex_dir", "adt_dir", "patient_id", "group",
                    "timepoint", "disease_status", "treatment", "relapse_status")
stopifnot(all(required_cols %in% colnames(sample_sheet)))

# ==============================================================================
# Helper functions
# ==============================================================================

# Reads a CITE-seq ADT count directory manually rather than via Read10X().
# Read10X() parses features.tsv with data.table::fread(), which auto-detects a
# field delimiter; antibody panels mixing space-free names with names
# containing spaces/parentheses can break that delimiter inference and
# silently truncate the parsed features vector. readLines() has no delimiter
# to infer, so it parses any antibody-naming convention correctly.
read_adt_manual <- function(data_dir) {
  mat      <- Matrix::readMM(gzfile(file.path(data_dir, "matrix.mtx.gz")))
  barcodes <- readLines(gzfile(file.path(data_dir, "barcodes.tsv.gz")))
  features <- readLines(gzfile(file.path(data_dir, "features.tsv.gz")))
  stopifnot(length(features) == nrow(mat), length(barcodes) == ncol(mat))
  rownames(mat) <- features
  colnames(mat) <- barcodes
  as(mat, "CsparseMatrix")
}

# Attaches a CITE-seq ADT assay to a GEX Seurat object, restricted to the
# barcodes present in both, and drops any "unmapped" control feature.
attach_adt <- function(obj, adt_counts) {
  colnames(adt_counts) <- paste0(colnames(adt_counts), "-1")
  rownames(adt_counts) <- sub("-[^-]*$", "", rownames(adt_counts))
  adt_obj <- CreateSeuratObject(counts = adt_counts)
  joint_bcs <- intersect(colnames(obj), colnames(adt_obj))
  obj     <- obj[, joint_bcs]
  adt_obj <- adt_obj[, joint_bcs]
  adt_data <- GetAssayData(adt_obj)
  adt_data <- adt_data[!grepl("unmapped", rownames(adt_data), ignore.case = TRUE), colnames(obj), drop = FALSE]
  obj[["ADT"]] <- CreateAssayObject(counts = adt_data)
  obj
}

# QC filtering: fixed nFeature/percent.mt thresholds.
qc_filter <- function(obj) {
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6240 & percent.mt < 20)
}

# Loads one sample: GEX counts, QC filtering, optional ADT attachment, and
# sample-level metadata from the sample sheet. NormalizeData()/
# FindVariableFeatures() are deliberately NOT run here - both are only ever
# computed on the full merged object downstream (see step 3), so running them
# per-sample here would be redundant, discarded computation.
load_sample <- function(row) {
  gex_counts <- Read10X(data.dir = file.path(DATA_ROOT, row$gex_dir))
  obj <- CreateSeuratObject(counts = gex_counts, project = paste0(row$sample_id, "_GEX"), min.features = 200)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- qc_filter(obj)

  if (!is.na(row$adt_dir) && nzchar(row$adt_dir)) {
    adt_counts <- read_adt_manual(file.path(DATA_ROOT, row$adt_dir))
    obj <- attach_adt(obj, adt_counts)
  }

  obj$sample_id      <- row$sample_id
  obj$patient_id     <- row$patient_id
  obj$group          <- row$group
  obj$timepoint      <- row$timepoint
  obj$ANCA           <- row$ANCA
  obj
}

# Merges a list of Seurat objects in chunks rather than one merge() call
# across all samples at once, to keep peak memory manageable on large
# sample sets. Reduce chunk_size if memory is still limiting.
chunked_merge <- function(obj_list, ids, chunk_size = 8, project_prefix = "chunk") {
  n <- length(obj_list)
  starts <- seq(1, n, by = chunk_size)
  intermediates <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    idx <- starts[k]:min(starts[k] + chunk_size - 1, n)
    if (length(idx) == 1) {
      intermediates[[k]] <- obj_list[[idx]]
    } else {
      intermediates[[k]] <- merge(obj_list[[idx[1]]], y = obj_list[idx[-1]],
                                   add.cell.ids = ids[idx],
                                   project = paste0(project_prefix, k), merge.data = TRUE)
    }
    message("merged chunk ", k, " of ", length(starts), " (", length(idx), " samples)")
    gc()
  }
  intermediates
}

# ==============================================================================
# 1. Load all samples
# ==============================================================================
sample_objects <- setNames(
  lapply(seq_len(nrow(sample_sheet)), function(i) load_sample(sample_sheet[i, ])),
  sample_sheet$sample_id
)

# ==============================================================================
# 2. Merge into one object
# ==============================================================================
intermediates <- chunked_merge(sample_objects, names(sample_objects), chunk_size = 8, project_prefix = "chunk")
rm(sample_objects); gc()

merged <- if (length(intermediates) == 1) {
  intermediates[[1]]
} else {
  merge(intermediates[[1]], y = intermediates[-1], project = "merged_ALL", merge.data = TRUE)
}
rm(intermediates); gc()

saveRDS(merged, file = file.path(OUTPUT_DIR, "merged_raw.rds"))

# ==============================================================================
# 3. Normalize, PCA, batch correction (Harmony), clustering
# ==============================================================================
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = 50)
merged <- RunHarmony(merged, group.by.vars = "orig.ident", plot_convergence = TRUE)

saveRDS(merged, file = file.path(OUTPUT_DIR, "merged_harmony.rds"))

# Inspect an elbow plot (ElbowPlot(merged)) before committing to the number of
# PCs used below - 1:15 reflects the variance structure of our own dataset and
# should be re-evaluated for a different dataset.
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:20)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:20)
merged <- FindClusters(merged, reduction = "harmony", resolution = 1.0)

saveRDS(merged, file = file.path(OUTPUT_DIR, "merged_clustered.rds"))

# ==============================================================================
# 4. Split T- and B-lineage clusters
# ==============================================================================
# Cluster-to-lineage assignment based on marker-gene inspection (not
# re-derived here; see Methods for the marker panel used).
T_CLUSTERS <- c("1", "2", "4", "5", "7", "9", "10", "13", "18", "25", "26", "28", "29")
B_CLUSTERS <- c("0", "3", "6", "8", "11", "12", "15", "17", "18", "19", "20", "21", "22", "24", "27")

merged_T <- subset(merged, idents = T_CLUSTERS)
merged_B <- subset(merged, idents = B_CLUSTERS)

saveRDS(merged_T, file = file.path(OUTPUT_DIR, "Tcells_raw.rds"))
saveRDS(merged_B, file = file.path(OUTPUT_DIR, "Bcells_raw.rds"))
