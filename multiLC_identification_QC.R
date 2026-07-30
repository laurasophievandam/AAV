# ==============================================================================
# Multi-LC clonal light-chain diversity: QC pipeline
# ==============================================================================
# Distinct from dual-LC (two light chains within a single CELL - see
# dualLC_identification_QC.R / dualLC_figures.R): this identifies clonal
# families (cells sharing a common heavy-chain lineage, already assigned a
# patient_clone_id upstream) whose member cells collectively use more than one
# distinct light chain - i.e. light-chain diversity across a clone, consistent
# with receptor revision/secondary light-chain rearrangement having occurred at
# different points within the same B-cell lineage.
#
# Built directly from the per-contig combined BCR table - no separate
# pre-aggregated "master" table dependency, since the two clone-level values
# this pipeline actually needs (LC_number, HC_umi_med) are either already a
# column in that table (LC_number) or trivially computed from it (HC_umi_med).

# ==============================================================================

BCR_COMBINED <- "<path/to/bcr_combined_annotated_table.csv>"  # per-contig BCR table with
                                                               # patient_clone_id, LC_number,
                                                               # sample_id, ANCA, locus,
                                                               # productive, v_call,
                                                               # junction_aa, umi_count
OUTPUT_DIR   <- "<path/to/output/directory>"

HC_UMI_CUT    <- 10
LIN_UMI_CUT   <- 10
LIN_MIN_CELLS <- 2

library(dplyr)
library(stringr)
library(tidyr)

build_clone_HC_summary <- function(raw) {
  raw %>%
    filter(productive == TRUE, locus == "IGH") %>%
    group_by(patient_clone_id) %>%
    summarise(HC_umi_med = median(umi_count), .groups = "drop")
}

build_lc_lineages <- function(raw) {
  raw %>%
    filter(productive == TRUE, locus %in% c("IGK", "IGL")) %>%
    group_by(patient_clone_id, junction_aa) %>%
    summarise(lin_umi_max = max(umi_count),
              lin_n_cells = n_distinct(cell_id),
              v_call = first(v_call),
              .groups = "drop") %>%
    mutate(lin_pass = (lin_umi_max >= LIN_UMI_CUT) | (lin_n_cells >= LIN_MIN_CELLS),
           v_gene = str_split(v_call, "\\*", simplify = TRUE)[, 1],
           v_gene = str_split(v_gene, ",", simplify = TRUE)[, 1])
}

confirmed_lc_vgenes_per_clone <- function(lin) {
  lin %>%
    filter(lin_pass) %>%
    group_by(patient_clone_id) %>%
    summarise(n_LC_vgenes_confirmed = n_distinct(v_gene), .groups = "drop")
}

run_pipeline <- function(raw) {
  clones <- raw %>%
    select(patient_clone_id, sample_id, ANCA, LC_number) %>%
    distinct(patient_clone_id, .keep_all = TRUE)

  hc_summary <- build_clone_HC_summary(raw)
  lin <- build_lc_lineages(raw)
  vgenes <- confirmed_lc_vgenes_per_clone(lin)

  df <- clones %>%
    left_join(hc_summary, by = "patient_clone_id") %>%
    left_join(vgenes, by = "patient_clone_id") %>%
    mutate(n_LC_vgenes_confirmed = replace_na(n_LC_vgenes_confirmed, 0L))

  start <- df %>% filter(LC_number >= 2)
  cat(sprintf("multi-LC candidate clones (LC_number>=2):        %5d\n", nrow(start)))

  after_hc <- start %>% filter(HC_umi_med >= HC_UMI_CUT)
  cat(sprintf("+ HC UMI(med) >= %d:                             %5d\n", HC_UMI_CUT, nrow(after_hc)))

  qc <- after_hc %>% filter(n_LC_vgenes_confirmed >= 2)
  cat(sprintf("+ >=2 confirmed LC V-genes (lineage QC):         %5d\n", nrow(qc)))

  qc
}

cohort_comparison <- function(qc, raw) {
  totals <- raw %>%
    group_by(sample_id) %>%
    summarise(n_total_clones = n_distinct(patient_clone_id), .groups = "drop")

  # pull ANCA from the raw table (covers every sample, even ones with zero
  # multi-LC candidate clones) rather than from qc, which only covers
  # samples with >=1 candidate - otherwise patients with zero surviving
  # clones lose their ANCA label and silently drop out of the comparison.
  pat_anca <- raw %>% distinct(sample_id, .keep_all = TRUE) %>% select(sample_id, ANCA)

  n_qc <- qc %>% count(sample_id, name = "n_qc")

  per_pat <- totals %>%
    left_join(pat_anca, by = "sample_id") %>%
    left_join(n_qc, by = "sample_id") %>%
    mutate(n_qc = replace_na(n_qc, 0)) %>%
    filter(!is.na(ANCA)) %>%
    apply_sample_exclusions() %>%
    assign_timepoint() %>%
    mutate(pct_qc = 100 * n_qc / n_total_clones)

  order <- c("HC", "MPO", "PR3")
  groups <- lapply(order, function(c) per_pat$pct_qc[per_pat$ANCA == c])
  p <- kruskal.test(groups)$p.value

  cat("\nn patients per cohort:\n")
  print(table(per_pat$ANCA))
  cat("median % QC-passing multi-LC clones per cohort:\n")
  print(per_pat %>% group_by(ANCA) %>% summarise(median_pct_qc = median(pct_qc), .groups = "drop"))
  cat(sprintf("\nKruskal-Wallis p = %.4g\n", p))

  list(per_pat = per_pat, p = p)
}

# ==============================================================================
# Run
# ==============================================================================
cols <- c("patient_clone_id", "sample_id", "ANCA", "LC_number", "cell_id",
          "locus", "productive", "v_call", "junction_aa", "umi_count")
raw <- read.csv(BCR_COMBINED, stringsAsFactors = FALSE)[, cols]

qc <- run_pipeline(raw)
write.csv(qc, file.path(OUTPUT_DIR, "multiLC_clones_QCpass.csv"), row.names = FALSE)

result <- cohort_comparison(qc, raw)
write.csv(result$per_pat, file.path(OUTPUT_DIR, "multiLC_cohort_stats.csv"), row.names = FALSE)

cat("\nwrote multiLC_clones_QCpass.csv and multiLC_cohort_stats.csv\n")
