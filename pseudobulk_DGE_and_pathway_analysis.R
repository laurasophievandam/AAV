# ==============================================================================
# Pseudobulk differential gene expression
# ==============================================================================
# Works on either the annotated T-cell or B-cell object (see
# Tcell_subclustering_annotation.R / Bcell_subclustering_annotation.R) - both
# carry the same `annotated_clusters` and `expansion_status` metadata columns.
# Subsets to one cell-type cluster (and, by default, clonally expanded cells),
# aggregates to pseudobulk per patient, and runs an edgeR-based differential
# expression test between two groups.
#
# Dependency note: convertBulkFilter() and edgeRtest() come from the internal
# lab packages `bioscrna`/`scgtest` (thin wrappers around standard
# single-cell-to-pseudobulk aggregation and edgeR, respectively) - a reader
# without access to those packages cannot run this script as-is; the
# aggregation is equivalent to summing raw counts per patient per gene, and
# edgeRtest() to a standard edgeR::glmLRT() workflow with TMM normalization.
#
# Usage: Rscript pseudobulk_DGE.R <cell_type> <design_var> <contrast>
#   e.g. Rscript pseudobulk_DGE.R "PB KI67+" group AAV-HC
# ==============================================================================

SEURAT_RDS   <- "<path/to/annotated_object.rds>"  # final annotated T- or B-cell Seurat object
SCE_RDS      <- "<path/to/output/directory>/sce_raw.rds"  # converted SCE, written by step 0,
                                                           # read by step 1
PATIENT_META <- "<path/to/output/directory>/metadata.csv"  # full per-cell metadata, written by
                                                             # step 0, deduped to one row per
                                                             # patient in step 2
OUTPUT_DIR   <- "<path/to/output/directory>"

cmd_args   <- commandArgs(trailingOnly = TRUE)
CELL_TYPE  <- if (length(cmd_args) >= 1) cmd_args[1] else "PB KI67+"
DESIGN_VAR <- if (length(cmd_args) >= 2) cmd_args[2] else "group"
CONTRAST   <- if (length(cmd_args) >= 3) cmd_args[3] else "AAV-HC"

library(SingleCellExperiment)
library(Seurat)
library(edgeR)
library(bioscrna)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(forcats)

# ==============================================================================
# 0. Convert the annotated Seurat object to a SingleCellExperiment
# ==============================================================================
# Only needs to be run once per Seurat object - re-run whenever the upstream
# annotated object is regenerated.
seurat <- readRDS(SEURAT_RDS)
seurat <- JoinLayers(seurat)

counts <- seurat@assays[["RNA"]]@layers[["counts"]]
colnames(counts) <- colnames(seurat)
rownames(counts) <- rownames(seurat)

logcounts <- seurat@assays[["RNA"]]@layers[["data"]]
colnames(logcounts) <- colnames(seurat)
rownames(logcounts) <- rownames(seurat)

seurat_metadata <- seurat@meta.data
stopifnot(identical(rownames(seurat_metadata), colnames(counts)))
write.csv(seurat_metadata, file = PATIENT_META, row.names = TRUE)

sce <- SingleCellExperiment(assays = list(counts = counts, logcounts = logcounts))
colData(sce) <- cbind(colData(sce), seurat_metadata)
sce$batch <- sce$sample_id

reducedDim(sce, "PCA")     <- seurat@reductions[["pca"]]@cell.embeddings
reducedDim(sce, "HARMONY") <- seurat@reductions[["harmony"]]@cell.embeddings
reducedDim(sce, "UMAP")    <- seurat@reductions[["umap"]]@cell.embeddings

saveRDS(sce, file = SCE_RDS)

# ==============================================================================
# 1. Subset to one cell type, aggregate to pseudobulk per patient
# ==============================================================================
sce <- readRDS(SCE_RDS)

sce <- sce[, sce$expansion_status == "Expanded"]
sce <- sce[, sce$annotated_clusters == CELL_TYPE]

bulk <- convertBulkFilter(sce, group = sce$patient_id, filter = FALSE)
counts <- bulk$counts
metadata <- bulk$metadata

# ==============================================================================
# 2. Attach patient-level metadata
# ==============================================================================
patient_meta <- read.csv(PATIENT_META, row.names = 1)
patient_meta <- patient_meta[!duplicated(patient_meta$patient_id),
                             c("patient_id", "group", "timepoint", "disease_status", "treatment")]
patient_meta <- patient_meta[match(colnames(counts), patient_meta$patient_id), ]
stopifnot(identical(patient_meta$patient_id, colnames(counts)))

metadata <- cbind(metadata, patient_meta)

# ==============================================================================
# 3. Differential expression (edgeR, TMM normalization, LRT)
# ==============================================================================
design <- model.matrix(as.formula(paste0("~ 0 + ", DESIGN_VAR)), data = metadata)
colnames(design) <- gsub(DESIGN_VAR, "", colnames(design))

markers <- edgeRtest(counts, metadata, group = DESIGN_VAR, design, contrasts = CONTRAST,
                     method = "lrt", normalize = "TMM", filter = TRUE)

top <- markers$top

# Remove uninformative/noisy gene classes: antisense/lncRNA/pseudogene
# families (AC*, AL*, LINC*), ribosomal and mitochondrial genes, sex-linked
# XIST/TSIX, and immune-receptor loci (which reflect clonotype identity
# rather than a regulated expression program).
noise_patterns <- c("^AC", "^RP[SL]", "^MT-", "XIST", "TSIX", "^LINC", "^AL", "^IGH", "^IGK", "^IGL", "^TRBV")
for (pattern in noise_patterns) {
  top <- top[grep(pattern, rownames(top), invert = TRUE), ]
}

write.csv(top, file = file.path(OUTPUT_DIR, paste0("DE_", gsub("[^A-Za-z0-9]", "_", CELL_TYPE), "_", CONTRAST, ".csv")))

# ==============================================================================
# 4. Volcano plot
# ==============================================================================
top_sig <- top[top$FDR < 0.05, ]
up   <- head(top_sig[order(-top_sig$logFC), ], 20)
down <- head(top_sig[order(top_sig$logFC), ], 20)
top_labels <- c(rownames(up), rownames(down))

gra <- EnhancedVolcano(top, lab = rownames(top), x = "logFC", y = "FDR",
                       selectLab = top_labels,
                       xlab = bquote(~Log[2] ~ "fold change"),
                       pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.0, labSize = 3,
                       labCol = "black", labFace = "bold", boxedLabels = FALSE,
                       colAlpha = 4 / 5, legendPosition = "right", legendLabSize = 14,
                       legendIconSize = 4.0, drawConnectors = TRUE, widthConnectors = 0.5,
                       colConnectors = "black", max.overlaps = 60, subtitle = NULL) +
  theme_classic() +
  ggtitle(paste0(CELL_TYPE, ": ", CONTRAST)) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13),
        legend.text = element_text(size = 13), legend.title = element_text(size = 14),
        panel.grid = element_blank(), title = element_text(face = "bold"),
        legend.position = "none")

ggsave(filename = file.path(OUTPUT_DIR, paste0("volcano_", gsub("[^A-Za-z0-9]", "_", CELL_TYPE), "_", CONTRAST, ".pdf")),
       plot = gra, width = 9, height = 5)

# ==============================================================================
# 5. Hallmark pathway enrichment (fgsea)
# ==============================================================================
set.seed(1)
NPERM <- 10000

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

gene_rank <- top$logFC
names(gene_rank) <- rownames(top)
gene_rank <- sort(gene_rank, decreasing = TRUE)

gsea_res <- fgsea(hallmark_list, stats = gene_rank, minSize = 10, eps = 0,
                  nPermSimple = NPERM, scoreType = "std")
gsea_res$leadingEdge <- vapply(gsea_res$leadingEdge, paste, collapse = ", ", FUN.VALUE = character(1))
gsea_res <- gsea_res %>% arrange(padj) %>% as.data.frame()
gsea_res$pathway <- gsub("^HALLMARK_", "", gsea_res$pathway)

write.csv(gsea_res,
         file = file.path(OUTPUT_DIR, paste0("hallmark_GSEA_", gsub("[^A-Za-z0-9]", "_", CELL_TYPE), "_", CONTRAST, ".csv")),
         row.names = FALSE)

# top 20 by significance, or all significant (padj < 0.05) pathways if more than 20
top_pathways <- gsea_res[gsea_res$padj < 0.05, ]
if (nrow(top_pathways) < 20) top_pathways <- head(gsea_res[order(gsea_res$padj), ], 20)
top_pathways <- top_pathways[order(top_pathways$padj, top_pathways$NES), ]
top_pathways$pathway <- fct_rev(fct_reorder(top_pathways$pathway, top_pathways$NES, .desc = TRUE))
top_pathways$direction <- factor(ifelse(top_pathways$NES > 0, "up", "down"), levels = c("up", "down"))

contrast_groups <- strsplit(CONTRAST, "-")[[1]]

gsea_plot <- ggplot(top_pathways, aes(NES, pathway, label = ifelse(padj < 0.05, "*", " "),
                                      alpha = -log10(pval), fill = direction)) +
  geom_col() +
  geom_text() +
  labs(x = "Normalized enrichment score", fill = CONTRAST,
       title = paste0("Hallmark gene sets - ", CELL_TYPE)) +
  scale_fill_manual(values = c(up = "darkred", down = "darkblue"),
                    labels = c(paste0("Up in ", contrast_groups[1]), paste0("Down in ", contrast_groups[1]))) +
  scale_alpha(guide = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"))

ggsave(filename = file.path(OUTPUT_DIR, paste0("hallmark_GSEA_", gsub("[^A-Za-z0-9]", "_", CELL_TYPE), "_", CONTRAST, ".pdf")),
       plot = gsea_plot, width = 10, height = 6, dpi = 600)
