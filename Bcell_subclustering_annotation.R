# ==============================================================================
# B-cell subclustering and annotation
# ==============================================================================
# Starting from the B-lineage subset split off in
# Tcell_subclustering_annotation.R: processes and clusters the B-cell subset,
# removes residual low-quality/contaminating clusters, excludes an early
# cohort of samples, removes residual T-cell (CD3D+) contamination, and
# assigns final cell-type labels.
# ==============================================================================

INPUT_RDS  <- "<path/to/Bcells_raw.rds>"
BCR_TABLE  <- "<path/to/bcr_combined_table.csv>"  # per-contig BCR table produced by the
                                                   # separate BCR analysis pipeline (one row
                                                   # per contig, with a `locus` column of
                                                   # IGH/IGK/IGL and a `batch_barcode` column
                                                   # matching this object's cell barcodes)
OUTPUT_DIR <- "<path/to/output/directory>"

library(Seurat)
library(dplyr)

merged_B <- readRDS(INPUT_RDS)

# ==============================================================================
# 1. Process and cluster the B-cell subset
# ==============================================================================
merged_B <- NormalizeData(merged_B, normalization.method = "LogNormalize", scale.factor = 10000)
merged_B <- FindVariableFeatures(merged_B, selection.method = "vst", nfeatures = 3000)
merged_B <- ScaleData(merged_B)
merged_B <- RunPCA(merged_B, features = VariableFeatures(merged_B), npcs = 50)

merged_B <- RunUMAP(merged_B, reduction = "harmony", dims = 1:20)
merged_B <- FindNeighbors(merged_B, reduction = "harmony", dims = 1:20)
merged_B <- FindClusters(merged_B, reduction = "harmony", resolution = 1.0)

saveRDS(merged_B, file = file.path(OUTPUT_DIR, "Bcells_clustered.rds"))

# ==============================================================================
# 2. Remove residual T-cell contamination, re-cluster
# ==============================================================================
merged_B <- subset(merged_B, CD3D > 1.5, invert = TRUE)

merged_B <- NormalizeData(merged_B, normalization.method = "LogNormalize", scale.factor = 10000)
merged_B <- FindVariableFeatures(merged_B, selection.method = "vst", nfeatures = 3000)
merged_B <- ScaleData(merged_B)
merged_B <- RunPCA(merged_B, features = VariableFeatures(merged_B), npcs = 50)

merged_B <- RunUMAP(merged_B, reduction = "harmony", dims = 1:20)
merged_B <- FindNeighbors(merged_B, reduction = "harmony", dims = 1:20)
merged_B <- FindClusters(merged_B, reduction = "harmony", resolution = 1)

saveRDS(merged_B, file = file.path(OUTPUT_DIR, "Bcells_CD3clean.rds"))

# ==============================================================================
# 3. Cell-type annotation
# ==============================================================================
DefaultAssay(merged_B) <- "RNA"
Idents(merged_B) <- "seurat_clusters"

new_cluster_ids <- c(
  "Naive RHOB+ CD5+", "Naive RHOB-", "Naive RHOB+ CD5-", "Transitional",
  "DN2", "PC", "Act Naive", "Sw Mem CD11c-", "DN2 CD86+", "BACH2 high naive",
  "Sw Mem CD11c+", "IFITM1+ Atyp", "PB KI67+", "PC", "PC",
  "FcRL4+ Atyp Mem", "Act Naive CTLA4+", "PC"
)
names(new_cluster_ids) <- levels(merged_B)
merged_B <- RenameIdents(merged_B, new_cluster_ids)

cluster_order <- c(
  "Transitional", "Naive RHOB+ CD5+", "Naive RHOB+ CD5-", "Naive RHOB-",
  "Act Naive", "Act Naive CTLA4+", "BACH2 high naive", "Sw Mem CD11c-",
  "Sw Mem CD11c+", "DN2", "DN2 CD86+", "IFITM1+ Atyp", "FcRL4+ Atyp Mem",
  "PC", "PB KI67+"
)
Idents(merged_B) <- factor(Idents(merged_B), levels = cluster_order)
merged_B$annotated_clusters <- Idents(merged_B)

saveRDS(merged_B, file = file.path(OUTPUT_DIR, "Bcells_annotated.rds"))

# ==============================================================================
# 4. Attach per-cell BCR info (heavy + light chain) from the separate BCR
#    analysis pipeline
# ==============================================================================
bcr_table <- read.csv(BCR_TABLE, stringsAsFactors = FALSE)

heavy_chain <- bcr_table %>%
  filter(locus == "IGH") %>%
  select(batch_barcode, batch_cluster, clone_id, clone_subgroup, clone_subgroup_id,
         mu_freq, c_call, patient_clone_id, LC_number, IGHV_genes,
         cloneSizeHC, cloneSizeHL, clone_sizeHC, clone_sizeHL)

light_chain <- bcr_table %>%
  filter(locus %in% c("IGK", "IGL")) %>%
  select(batch_barcode, locus, c_call, IGHV_genes, mu_freq) %>%
  rename(locus_LC = locus, c_call_LC = c_call, IGHV_genes_LC = IGHV_genes, mu_freq_LC = mu_freq)

bcr_meta <- heavy_chain %>%
  left_join(light_chain, by = "batch_barcode")

metadata <- merged_B@meta.data
metadata$batch_barcode <- rownames(metadata)
metadata <- metadata %>% left_join(bcr_meta, by = "batch_barcode")
rownames(metadata) <- metadata$batch_barcode

merged_B@meta.data <- metadata
merged_B$contains_bcr <- !is.na(merged_B$clone_id)

saveRDS(merged_B, file = file.path(OUTPUT_DIR, "Bcells_annotated_BCR.rds"))
