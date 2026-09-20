# Cytomegalovirus drives autoreactive B cells in ANCA-associated vasculitis

Analysis code accompanying the paper *Cytomegalovirus drives autoreactive B cells in ANCA-associated vasculitis*.

This repository contains the R scripts used for single-cell GEX/CITE-seq processing, T- and B-cell annotation, BCR dual- and multi-light-chain QC, pseudobulk differential expression, and GLIPH2 TCR specificity-group analysis.

Scripts are written as analysis notebooks rather than a single entry-point package. Edit the path variables at the top of each file, then run that script in R (or via `Rscript` where noted). Cluster IDs, cell-type labels, Harmony dimensions, and several QC thresholds are specific to this cohort and should be re-evaluated for a new dataset.

## Companion R packages

These packages were developed with the paper and are used by parts of this pipeline:

- [screceptor](https://github.com/zhangzmr/screceptor): **BCR analysis.** Single-cell BCR/TCR repertoire processing (clonal clustering, germline reconstruction, SHM, merge with GEX). Used to produce the combined BCR tables consumed by the B-cell and light-chain scripts.
- [scgtest](https://github.com/zhangzmr/scgtest): **Pseudobulk DE testing.** Differential expression for single-cell and pseudobulk RNA-seq (edgeR, RUV, pathway enrichment). Provides `edgeRtest()` used in `pseudobulk_DGE_and_pathway_analysis.R`.
- [bioscrna](https://github.com/zhangzmr/bioscrna): Bioconductor-based single-cell RNA-seq pipeline (loading, QC, normalization, clustering, markers). Provides `convertBulkFilter()` used in the pseudobulk script.

Install those packages from their GitHub repositories before running the DE or BCR-merge steps. The remaining scripts use standard CRAN/Bioconductor packages listed under [R dependencies](#r-dependencies).

## Repository layout

| Script | Role |
|---|---|
| `GEX_CITE_processing.R` | Per-sample GEX (+ optional CITE-seq ADT) load, QC, merge, Harmony, clustering; split T- vs B-lineage |
| `Tcell_subclustering_annotation.R` | T-cell re-clustering, annotation, TCR attach via scRepertoire |
| `Bcell_subclustering_annotation.R` | B-cell re-clustering, CD3D cleanup, annotation, BCR metadata attach |
| `dualLC_identification_QC.R` | Cells with exactly 1 heavy + 2 light chains (cell-level dual LC) |
| `multiLC_identification_QC.R` | Clones whose member cells use more than one light chain (clone-level multi LC) |
| `pseudobulk_DGE_and_pathway_analysis.R` | Patient-level pseudobulk DE (edgeR) + Hallmark GSEA |
| `gliph2_merged.R` | Helper functions to prepare GLIPH2 input, annotate output, and plot specificity networks |

License: MIT (see `LICENSE`).

## Analysis order

```
CellRanger GEX (+ optional ADT) ──► GEX_CITE_processing.R
                                      │
                                      ├── Tcells_raw.rds ──► Tcell_subclustering_annotation.R
                                      │                         └── CellRanger TCR contigs
                                      │
                                      └── Bcells_raw.rds ──► Bcell_subclustering_annotation.R
                                                                └── BCR table from screceptor

screceptor BCR pipeline ──► BCR_raw.rds ──────────► dualLC_identification_QC.R
                        └── bcr_combined_table.csv ──► multiLC_identification_QC.R
                                                   └── B-cell annotation (BCR attach)

annotated T or B object ──► pseudobulk_DGE_and_pathway_analysis.R

scRepertoire TCR table ──► gliph2_merged.R (steps 1–2)
                              │
                              ▼
                         GLIPH2 (web tool or binary; not run in R)
                              │
                              ▼
                         gliph2_merged.R (steps 3–5)
```

1. Process GEX/CITE and split lineages (`GEX_CITE_processing.R`).
2. Annotate T cells and attach TCR (`Tcell_subclustering_annotation.R`).
3. Annotate B cells and attach BCR (`Bcell_subclustering_annotation.R`).
4. In parallel, run dual-LC and multi-LC QC on the BCR tables from [screceptor](https://github.com/zhangzmr/screceptor).
5. Run pseudobulk DE / GSEA on an annotated object.
6. Prepare GLIPH2 input, run GLIPH2 outside R, then filter, annotate, and plot.

## Expected inputs

Placeholders such as `<path/to/raw/data>` at the top of each script must be replaced before running.

### Sample sheet (GEX / T-cell scripts)

CSV, one row per sample. `GEX_CITE_processing.R` requires these columns:

- `sample_id` — unique sample identifier
- `gex_dir` — CellRanger `filtered_feature_bc_matrix` directory, relative to `DATA_ROOT`
- `adt_dir` — CITE-seq ADT UMI-count directory, relative to `DATA_ROOT`; leave blank if that sample has no ADT
- `patient_id` — patient / subject identifier
- `group` — study group (e.g. AAV / HC)
- `timepoint` — sample timepoint label
- `disease_status`, `treatment`, `relapse_status` — clinical metadata (required by the column check)

The GEX loader also writes `ANCA` onto each cell from a column named `ANCA` (HC / MPO / PR3). Include that column as well.

For TCR attach, extend the same sheet with:

- `tcr_dir` — path to each sample’s CellRanger `filtered_contig_annotations.csv`, relative to `DATA_ROOT`

### GEX and CITE-seq

- One CellRanger GEX `filtered_feature_bc_matrix` per sample (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`).
- Optional matching ADT count directory (e.g. from CITE-seq-Count), same 10X-style files. Antibody names are read with `readLines()` rather than `Read10X()` so names that contain spaces or parentheses are not truncated.

### BCR (from screceptor, not generated in this repo)

- `BCR_raw.rds` — named list of per-sample AIRR contig data frames (`sample_id` → table). Used by dual-LC QC.
- Combined per-contig BCR table (CSV) with at least: `batch_barcode` or `cell_id`, `locus` (IGH / IGK / IGL), `productive`, `patient_clone_id`, `LC_number`, `sample_id`, `ANCA`, `v_call`, `junction_aa`, `umi_count`, plus the clone / SHM columns listed in `Bcell_subclustering_annotation.R`.

### TCR

- Per-sample CellRanger VDJ `filtered_contig_annotations.csv` for `Tcell_subclustering_annotation.R`.
- For GLIPH2: a per-cell TCR table (e.g. from scRepertoire) with `sample`, `CTaa`, `cdr3_aa1`, `cdr3_aa2`, `TCR2`.
- Optional VDJdb export (with HLA columns) if you want a viral reference / background file.

## Script notes

### 1. `GEX_CITE_processing.R`

Per-sample load → QC → optional ADT attach → chunked merge → log-normalize, HVG, PCA, Harmony (`orig.ident`) → UMAP / neighbors / clusters on Harmony dims 1–20 at resolution 1.0 → split T- and B-lineage clusters.

Edit at the top:

```r
DATA_ROOT    <- "<path/to/raw/data>"
SAMPLE_SHEET <- "<path/to/sample_sheet.csv>"
OUTPUT_DIR   <- "<path/to/output/directory>"
```

QC (fixed in this cohort): `nFeature_RNA > 200`, `nFeature_RNA < 6240`, `percent.mt < 20`.

Outputs written to `OUTPUT_DIR`:

| File | Contents |
|---|---|
| `merged_raw.rds` | Merged object after QC, before normalization |
| `merged_harmony.rds` | After Harmony |
| `merged_clustered.rds` | After UMAP and clustering |
| `Tcells_raw.rds` | T-lineage subset |
| `Bcells_raw.rds` | B-lineage subset |

Lineage split uses cluster IDs from this dataset (inspect markers before changing them):

- T clusters: `1, 2, 4, 5, 7, 9, 10, 13, 18, 25, 26, 28, 29`
- B clusters: `0, 3, 6, 8, 11, 12, 15, 17, 18, 19, 20, 21, 22, 24, 27`

Inspect `ElbowPlot(merged)` before committing to the number of Harmony dimensions. Merge is done in chunks of 8 samples to limit peak memory; lower `chunk_size` if needed.

### 2. `Tcell_subclustering_annotation.R`

Re-normalizes and clusters the T-cell subset, drops residual `CD19+` cells, assigns the paper’s T-cell labels, then attaches clonotypes with scRepertoire.

Set `INPUT_RDS` to the T-lineage object from step 1 (`Tcells_raw.rds`). The placeholder in the script currently says `merged_clustered.rds`; use the T-cell subset so cluster labels match T-cell biology.

TCR attach (`combineTCR` / `combineExpression`):

- `removeNA = TRUE`, `removeMulti = FALSE`, `filterMulti = TRUE`
- Clone call: `gene`
- Size bins: Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500, Ultraexpanded = 2000

Outputs: `Tcells_clustered.rds`, `Tcells_annotated.rds`, `Tcells_TCR_combined.rds`, `Tcells_annotated_TCR.rds`.

Annotation labels (dataset-specific, mapped onto `levels(merged_T)` after clustering): CD4 Naive, CD4 CM, CD4/CD8 EM CCR6+, CD8 EMRA GZMB+/NKT like, CD4 EM, CD4/8 CM, CD8 EMRA GZMB/GZMK+, CD8 Naive, CD4 CM CD69+, Activated T cells CD44+, CD4 Treg, CD8 EM GZMK+, CD4 TFH, NKT CD1D+, CD8 EM ITGA4+, NKT, CD4 CTLA4+.

### 3. `Bcell_subclustering_annotation.R`

Re-clusters `Bcells_raw.rds`, removes residual T-cell contamination (`CD3D > 1.5`), re-clusters, assigns B-cell labels, then left-joins heavy- and light-chain columns from the screceptor BCR table on `batch_barcode`.

Outputs: `Bcells_clustered.rds`, `Bcells_CD3clean.rds`, `Bcells_annotated.rds`, `Bcells_annotated_BCR.rds`.

`contains_bcr` is `TRUE` when `clone_id` is present.

Annotation labels (dataset-specific): Transitional, Naive RHOB+ CD5+, Naive RHOB+ CD5-, Naive RHOB-, Act Naive, Act Naive CTLA4+, BACH2 high naive, Sw Mem CD11c-, Sw Mem CD11c+, DN2, DN2 CD86+, IFITM1+ Atyp, FcRL4+ Atyp Mem, PC, PB KI67+.

UMAP / neighbors / clusters reuse the Harmony reduction already stored on the parent object (Harmony is not re-run in this script).

### 4. `dualLC_identification_QC.R`

Identifies **cells** with exactly one productive IGH contig and two productive light-chain contigs (IGK / IGL) from `BCR_raw.rds`.

QC: `min_umi >= 5` on the weaker light chain **and** `HC_umi >= 5`. There is no UMI-ratio filter.

`SAMPLE_SHEET` here is patient metadata (`patient_id` matching BCR sample names, plus `disease_status`, renamed to `ANCA`).

Outputs:

- `dualLC_cells_QCpass.csv` — QC-passing dual-LC cells (HC UMI, LC min/max UMI, chain pairing, ratio)
- `total_cells_per_sample.csv` — per-sample cell totals for use as a denominator

This is **not** the same as multi-LC (below). Dual-LC is two light chains in one cell; multi-LC is more than one light chain across cells that share a heavy-chain clone.

### 5. `multiLC_identification_QC.R`

Identifies **clonal families** (`patient_clone_id`) whose cells collectively use more than one distinct light-chain V gene, consistent with receptor revision / secondary LC rearrangement within a lineage.

Default thresholds:

- Start from clones with `LC_number >= 2`
- Clone-level median IGH UMI (`HC_umi_med`) ≥ 10
- A light-chain lineage is “confirmed” if `max UMI >= 10` **or** it is seen in ≥ 2 cells
- Keep clones with ≥ 2 confirmed LC V genes

Outputs: `multiLC_clones_QCpass.csv`, `multiLC_cohort_stats.csv` (per-sample % QC-passing multi-LC clones; Kruskal–Wallis across HC / MPO / PR3).

`cohort_comparison()` calls `apply_sample_exclusions()` and `assign_timepoint()`, which are **not defined in this repository**. Provide those helpers in the session (or comment those two lines out) before running the cohort comparison.

### 6. `pseudobulk_DGE_and_pathway_analysis.R`

Works on a final annotated T- or B-cell Seurat object. Converts it once to a `SingleCellExperiment`, subsets to **expanded** cells (`expansion_status == "Expanded"`) in one `annotated_clusters` label, aggregates raw counts per patient (`convertBulkFilter()` from bioscrna), and runs edgeR LRT with TMM normalization (`edgeRtest()` from scgtest).

Command-line usage:

```bash
Rscript pseudobulk_DGE_and_pathway_analysis.R "<cell_type>" <design_var> <contrast>
# example:
Rscript pseudobulk_DGE_and_pathway_analysis.R "PB KI67+" group AAV-HC
```

Defaults if arguments are omitted: cell type `PB KI67+`, design `group`, contrast `AAV-HC`.

After DE, genes matching `AC*`, `RPS`/`RPL`, `MT-`, `XIST`, `TSIX`, `LINC*`, `AL*`, `IGH`/`IGK`/`IGL`, and `TRBV` are dropped as uninformative.

Then:

1. Volcano plot (`EnhancedVolcano`; FDR 0.05, |log2FC| 0.5); labels the top 20 up and 20 down FDR-significant genes.
2. Hallmark GSEA (`msigdbr` category `H`, `fgsea`, 10,000 permutations, `minSize = 10`).

Outputs (names include the cell type and contrast):

- `sce_raw.rds`, `metadata.csv` (from the one-time Seurat → SCE conversion)
- `DE_<cell_type>_<contrast>.csv`
- `volcano_<cell_type>_<contrast>.pdf`
- `hallmark_GSEA_<cell_type>_<contrast>.csv` and `.pdf`

`expansion_status` must already exist on the object. If you only have scRepertoire clone-size bins, add that column first (e.g. treat clone count > 1 as Expanded).

### 7. `gliph2_merged.R`

Function library for the GLIPH2 workflow. GLIPH2 clustering itself is **not** run in R. Use the [GLIPH2 web tool](http://50.255.35.37:8080/) or the standalone binary between steps 2 and 3.

The file concatenates config plus five steps that were originally separate scripts. The later `source("00_config.R")` lines refer to a file that is not in this repo — the shared settings already live at the top of `gliph2_merged.R`. Source or run the relevant section, then call the example usage blocks.

| Step | Function(s) | What it does |
|---|---|---|
| 1 | `prepare_gliph_input()` | Collapse a per-cell TCR table to one row per clonotype per sample in GLIPH2 format (CDR3b, TRBV, TRBJ, CDR3a, `subject:condition`, count) |
| 2 | `prepare_vdjdb_reference()` | Encode VDJdb antigen / HLA metadata into `subject:condition` for a viral background file |
| 3 | `filter_significant_groups()`, `annotate_viral_specificity()`, `split_sample_metadata()`, `summarize_dominant_value()`, `count_antigen_hits()` | Keep significant groups, map CDR3 motifs to curated CMV / EBV / AdV epitopes, parse encoded sample fields |
| 4 | `flag_grouped_clones()`, `compute_expansion()`, `summarize_proportion()`, `plot_proportion()`, `merge_final_annotations()` | Join groups back onto the full TCR table; expansion and proportion plots |
| 5 | `build_specificity_network()`, `plot_specificity_network()` | igraph network of groups that share CDR3b sequences |

Default group filters: Fisher score ≤ 0.05 and ≥ 3 subjects. AAV vs viral reference samples are distinguished with the regex `MPO|PR3|:HC|neg` — change this if your sample IDs differ.

Default directories created if missing: `data/gliph2/` and `results/`.

## R dependencies

```r
# CRAN
install.packages(c(
  "Seurat", "Matrix", "harmony", "dplyr", "purrr", "tidyr",
  "stringr", "tidyverse", "ggplot2", "ggpubr", "forcats", "igraph"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "SingleCellExperiment", "edgeR", "EnhancedVolcano", "fgsea"
))

# other
install.packages("msigdbr")
# scRepertoire: https://github.com/ncborcherding/scRepertoire

# paper-specific (see above)
# remotes::install_github("zhangzmr/screceptor")
# remotes::install_github("zhangzmr/scgtest")
# remotes::install_github("zhangzmr/bioscrna")
```

Typical session: R ≥ 4.x with Seurat v5 (the DE script uses `JoinLayers()` and assay layers).

## Re-using this code on another dataset

These choices are hardcoded for this AAV / CMV cohort and will not transfer automatically:

- QC cutoffs (`nFeature_RNA`, `percent.mt`, CD3D / CD19 filters, UMI floors)
- Harmony group (`orig.ident`), dimensions `1:20`, clustering resolution `1.0`
- T- and B-lineage cluster IDs after the first clustering
- Cell-type label vectors (length must equal the number of clusters you obtain)
- Dual-LC and multi-LC UMI / lineage rules
- GLIPH2 viral epitope motif map (CMV, EBV, AdV)
- Pseudobulk subset to `expansion_status == "Expanded"`

Re-inspect markers, elbow plots, and cluster counts before keeping those values.
