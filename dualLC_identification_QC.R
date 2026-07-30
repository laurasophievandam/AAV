# ==============================================================================
# Dual light-chain (LC) cell identification and QC
# ==============================================================================
# Identifies cells with exactly 1 heavy chain (HC) and 2 light chain (LC)
# productive contigs directly from the raw per-sample AIRR contig list
# (BCR_raw.rds) 
#
# QC: LC min_umi >= 5 AND HC_umi >= 5 (symmetric 5-UMI floor per chain). 
# ==============================================================================

BCR_RAW_RDS    <- "<path/to/BCR_raw.rds>"           # named list of per-sample AIRR contig
                                                     # data.frames (one element per sample_id)
SAMPLE_SHEET   <- "<path/to/patient_metadata.csv>"  # one row per patient, with columns
                                                     # patient_id (matching BCR_raw's sample
                                                     # names) and disease_status (ANCA)
OUTPUT_DIR     <- "<path/to/output/directory>"

library(dplyr)
library(purrr)
library(tidyr)

BCR_raw <- readRDS(BCR_RAW_RDS)

# ==============================================================================
# 1. Per-cell HC/LC contig summary across all samples
# ==============================================================================
cell_summary <- imap_dfr(BCR_raw, function(df, sample_id) {
  prod <- df %>% filter(productive == TRUE)

  contigs <- prod %>%
    group_by(cell_id, locus, junction_aa) %>%
    summarise(umi_count = max(umi_count), .groups = "drop")

  hc <- contigs %>% filter(locus == "IGH")
  lc <- contigs %>% filter(locus %in% c("IGK", "IGL"))

  n_h  <- hc %>% count(cell_id, name = "n_H")
  n_lc <- lc %>% count(cell_id, name = "n_LC")

  full_join(n_h, n_lc, by = "cell_id") %>%
    mutate(n_H = replace_na(n_H, 0L),
           n_LC = replace_na(n_LC, 0L),
           sample_id = sample_id)
})

candidates <- cell_summary %>% filter(n_H == 1, n_LC == 2)
cat("cells with exactly 1 HC + 2 LC contigs (from raw BCR data):", nrow(candidates), "\n")

# ==============================================================================
# 2. Pull HC UMI and LC UMI/pairing info for candidates, apply UMI floors
# ==============================================================================
qc_rows <- imap_dfr(BCR_raw, function(df, sample_id) {
  cand_cells <- candidates %>% filter(sample_id == !!sample_id) %>% pull(cell_id)
  if (length(cand_cells) == 0) return(NULL)

  prod <- df %>% filter(productive == TRUE, cell_id %in% cand_cells)
  contigs <- prod %>%
    group_by(cell_id, locus, junction_aa) %>%
    summarise(umi_count = max(umi_count), .groups = "drop")

  hc <- contigs %>% filter(locus == "IGH") %>%
    group_by(cell_id) %>% summarise(HC_umi = max(umi_count), .groups = "drop")

  lc <- contigs %>% filter(locus %in% c("IGK", "IGL")) %>%
    group_by(cell_id) %>%
    summarise(min_umi = min(umi_count), max_umi = max(umi_count),
              chains = paste(sort(locus), collapse = "+"), .groups = "drop") %>%
    mutate(ratio = min_umi / max_umi)

  hc %>% inner_join(lc, by = "cell_id") %>% mutate(sample_id = sample_id)
})

qc_rows <- qc_rows %>%
  filter(min_umi >= 5, HC_umi >= 5)  # symmetric 5-UMI floor per chain, no ratio filter

cat("+ LC min_umi>=5, HC_umi>=5 (no UMI-ratio filter):", nrow(qc_rows), "\n")

write.csv(qc_rows, file.path(OUTPUT_DIR, "dualLC_cells_QCpass.csv"), row.names = FALSE)

# ==============================================================================
# 3. Total cell counts per sample (denominator), with disease status attached
# ==============================================================================
totals <- imap_dfr(BCR_raw, function(df, sample_id) {
  data.frame(sample_id = sample_id, n_total_cells = length(unique(df$cell_id)))
})

sample_sheet <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE)
totals <- totals %>%
  left_join(sample_sheet %>% select(patient_id, disease_status),
            by = c("sample_id" = "patient_id")) %>%
  rename(ANCA = disease_status)

write.csv(totals, file.path(OUTPUT_DIR, "total_cells_per_sample.csv"), row.names = FALSE)
cat("wrote dualLC_cells_QCpass.csv and total_cells_per_sample.csv\n")
