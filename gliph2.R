# ============================================================================
# GLIPH2 TCR specificity-group pipeline — shared configuration
#
# Source this file first from every other script in this pipeline:
#   source("00_config.R")
#
# Pipeline order:
#   01_prepare_gliph_input.R           single-cell TCR table -> GLIPH2 input
#   02_prepare_viral_reference_database.R   VDJdb -> GLIPH2 reference/background
#   03_filter_and_annotate_gliph_output.R   raw GLIPH2 output -> significant,
#                                            annotated specificity groups
#   04_merge_and_plot.R                 join groups back onto the full TCR
#                                        table; clonal-expansion summary plots
#   05_network_analysis.R               specificity-group sharing network
#
# GLIPH2 itself (the clustering algorithm) is run outside R, between scripts
# 01/02 and 03 — either via the GLIPH2 web tool (http://50.255.35.37:8080/)
# or the standalone binary. This pipeline only prepares its input and
# processes its output.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# ---- Paths (edit these for your own project) ------------------------------
# data_dir:    single-cell TCR tables and the VDJdb reference download
# gliph_dir:   GLIPH2 input/output files
# results_dir: final merged tables and figures
data_dir    <- "data"
gliph_dir   <- file.path(data_dir, "gliph2")
results_dir <- "results"

dir.create(gliph_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Default significance thresholds for GLIPH2 specificity groups --------
# GLIPH2 reports two alternative significance scores per group:
#   Fisher_score - enrichment of the CDR3 motif itself
#   vb_score     - enrichment of V-beta gene usage within the group
# Pick whichever your cohort design supports (Fisher_score is the more
# commonly reported default) via the `score_col` argument of
# filter_significant_groups().
default_thresholds <- list(
  score_max      = 0.05,
  min_subjects   = 3,
)

# ---- Curated CDR3-motif -> viral epitope lookup ----------------------------
# Built from manual annotation of GLIPH2 "pattern" motifs against known
# CMV/EBV/AdV epitope-specific CDR3b motifs. Applied in order with
# last-match-wins semantics (a later row overrides an earlier one if a
# pattern matches the same group) — this replicates the original
# chained-mutate annotation exactly; see annotate_viral_specificity()
# in 03_filter_and_annotate_gliph_output.R.
viral_epitope_map <- tribble(
  ~pattern,      ~epitope,
  "S%GQG",       "CMV-pp50_VTE",
  "VTGGT",       "CMV-pp65_YSE",
  "LAPG",        "CMV-pp65_NLV",
  "GQGGV",       "CMV-IE1_VLE",
  "SLIG%S",      "CMV-pp65_TPR",
  "P%RNT",       "CMV-pp65_RPH",
  "V%R%%ANT",    "CMV-1E1_QIK",
  "TDLN",        "CMV-1E1_ELR",
  "S%QGG",       "EBV_LMP2_FLY",
  "DGMN",        "EBV_LMP2_CLG",
  "EGQ%%S",      "EBV_LMP2_CLG",
  "IAL",         "EBV_EBNA3C_LLD",
  "RDRVG",       "EBV_BMLF1_GLC",
  "RDRTG",       "EBV_BMLF1_GLC",
  "G%GGTN",      "EBV_BMLF1_GLC",
  "AIGGS",       "EBV_BRLF1_YVL",
  "RLTG",        "EBV_EBNA3A_RPP",
  "SSSLN",       "EBV_BZLF1_RAK",
  "RDR%%EN",     "EBV_BZLF1_RAK",
  "S%GQA",       "EBV_EBNA3A_FLR",
  "LETA",        "EBV_EBNA3A_QAK",
  "LETG",        "EBV_EBNA3A_QAK",
  "APGQ",        "AdV_HEXON_TDL",
  "VPGQ",        "AdV_HEXON_TDL",
  "AR%GLA",      "AdV_E1A_LLD",
  "INP",         "AdV_HEXON_KPY"
)

# ============================================================================
# Step 1 — build a GLIPH2 input file from a single-cell TCR table
#
# Expects a per-cell TCR table (e.g. combined across samples with
# scRepertoire) containing at minimum:
#   sample    - subject/sample identifier
#   CTaa      - combined alpha_beta clonotype string (used to collapse cells
#               into clonotypes; not written to the GLIPH2 file itself)
#   cdr3_aa1  - TCR-alpha CDR3 amino acid sequence
#   cdr3_aa2  - TCR-beta CDR3 amino acid sequence
#   TCR2      - beta-chain gene usage string, e.g. "TRBV20-1.TRBD1.TRBJ2-1.TRBC2"
# ============================================================================

source("00_config.R")

#' Collapse a per-cell TCR table to one row per clonotype per sample and
#' write it in GLIPH2 input format (CDR3b, TRBV, TRBJ, CDR3a,
#' subject:condition, count).
#'
#' @param tcr_table data frame with columns sample, CTaa, cdr3_aa1, cdr3_aa2, TCR2
#' @param output_path path to write the GLIPH2 input CSV
prepare_gliph_input <- function(tcr_table, output_path) {
  aggregated <- tcr_table %>%
    select(sample, CTaa, cdr3_aa1, cdr3_aa2, TCR2) %>%
    group_by(sample, CTaa) %>%
    summarise(
      count    = n(),
      cdr3_aa1 = dplyr::first(cdr3_aa1),
      cdr3_aa2 = dplyr::first(cdr3_aa2),
      TCR2     = dplyr::first(TCR2),
      .groups  = "drop"
    ) %>%
    separate(TCR2, into = c("TRBV", "col2", "TRBJ", "TRBC"), sep = "\\.") %>%
    select(-col2, -TRBC) %>%
    rename(
      `CDR3b`             = cdr3_aa2,
      `CDR3a`              = cdr3_aa1,
      `subject:condition` = sample
    ) %>%
    relocate(`CDR3b`, TRBV, TRBJ, `CDR3a`, `subject:condition`, count)
  
  write_csv(aggregated, output_path)
  aggregated
}

# ---- Example usage ----------------------------------------------------
# tcr_table <- read_csv(file.path(data_dir, "combined_TCR_table.csv"))
# prepare_gliph_input(tcr_table, file.path(gliph_dir, "gliph2_input.csv"))

# ============================================================================
# Step 2 — build a viral-specificity reference/background file for GLIPH2
# from a VDJdb export
#
# GLIPH2 only carries five fields per TCR (CDR3b, TRBV, TRBJ, CDR3a,
# subject:condition); antigen/HLA metadata has to be smuggled into the
# `subject:condition` field so it survives the clustering step and can be
# recovered afterwards. This step builds that encoded field from a VDJdb
# export that has been annotated with HLA (mhc.a / mhc.b / mhc.class).
# ============================================================================

source("00_config.R")

#' Reformat a VDJdb export into GLIPH2 reference format, encoding
#' antigen + HLA metadata into the `subject:condition` field as
#' "<original subject:condition>:<epitope>_<gene>_<species>_<mhc.a>_<mhc.b>_<mhc.class>".
#'
#' @param vdjdb_table data frame with columns `subject:condition`,
#'   antigen.epitope, antigen.gene, antigen.species, mhc.a, mhc.b, mhc.class,
#'   and `#CDR3b`
#' @param output_path path to write the reformatted reference CSV
prepare_vdjdb_reference <- function(vdjdb_table, output_path) {
  reference <- vdjdb_table %>%
    mutate(
      `subject:condition` = paste(`subject:condition`, antigen.epitope, sep = ":"),
      `subject:condition` = paste(`subject:condition`, antigen.gene, sep = "_"),
      `subject:condition` = paste(`subject:condition`, antigen.species, sep = "_"),
      `subject:condition` = paste(`subject:condition`, mhc.a, sep = "_"),
      `subject:condition` = paste(`subject:condition`, mhc.b, sep = "_"),
      `subject:condition` = paste(`subject:condition`, mhc.class, sep = "_")
    ) %>%
    select(-antigen.epitope, -antigen.gene, -antigen.species, -mhc.a, -mhc.b, -mhc.class) %>%
    filter(!is.na(`#CDR3b`))
  
  write_csv(reference, output_path)
  reference
}

# ---- Example usage ----------------------------------------------------
# vdjdb_table <- read_csv(file.path(data_dir, "vdjdb_HLA_clean.csv"))
# prepare_vdjdb_reference(vdjdb_table, file.path(gliph_dir, "viral_reference.csv"))

# ============================================================================
# Step 3 — filter raw GLIPH2 output to significant specificity groups and
# annotate each group with its dominant viral epitope / HLA / gene / peptide
#
# Expects a raw GLIPH2 output table (one row per TCR per specificity group)
# with at least: index, pattern, TcRb, TcRa, V, J, Sample, Freq,
# Fisher_score, vb_score, number_subject, number_unique_cdr3, final_score.
# `Sample` here is the encoded `subject:condition` field produced upstream
# (see 01_prepare_gliph_input.R / 02_prepare_viral_reference_database.R).
# ============================================================================

source("00_config.R")

#' Keep only statistically significant GLIPH2 specificity groups.
#'
#' @param df raw GLIPH2 output
#' @param score_col which GLIPH2 significance score to filter on —
#'   "Fisher_score" (CDR3 motif enrichment) or "vb_score" (V-beta usage
#'   enrichment). GLIPH2 reports both; pick the one appropriate to your
#'   cohort design.
#' @param score_max maximum value of score_col to retain (default 0.05)
#' @param min_subjects minimum distinct subjects contributing to a group
#' @param min_unique_cdr3 minimum distinct CDR3b sequences in a group
filter_significant_groups <- function(df,
                                      score_col       = "Fisher_score",
                                      score_max       = default_thresholds$score_max,
                                      min_subjects    = default_thresholds$min_subjects,
                                      min_unique_cdr3 = default_thresholds$min_unique_cdr3) {
  df %>%
    filter(.data[[score_col]] <= score_max) %>%
    filter(number_subject >= min_subjects) %>%
    filter(number_unique_cdr3 >= min_unique_cdr3)
}

#' Annotate each row's GLIPH2 `pattern` motif with a viral epitope, using
#' viral_epitope_map (see 00_config.R). Rows are checked against every
#' pattern in the map in order, and a later match overrides an earlier one
#' for the same row (last-match-wins) — this reproduces the original
#' sequential-mutate annotation exactly.
annotate_viral_specificity <- function(df, epitope_map = viral_epitope_map, pattern_col = "pattern") {
  df$annotation <- NA_character_
  for (i in seq_len(nrow(epitope_map))) {
    hits <- grepl(epitope_map$pattern[i], df[[pattern_col]])
    df$annotation[hits] <- epitope_map$epitope[i]
  }
  df
}

#' Split the encoded `Sample` (subject:condition) field back into its
#' component parts, and correct the parsed virus label using anything found
#' in the trailing "extra" component (handles CMV/EBV strings that spilled
#' into extra fields during encoding).
#'
#' @param df table with a `Sample` column
#' @param aav_pattern regex identifying AAV-cohort samples (vs. viral
#'   reference samples) in `Sample` — adjust for your own cohort labels
split_sample_metadata <- function(df, sample_col = "Sample", aav_pattern = "MPO|PR3|:HC|neg") {
  df %>%
    mutate(Sample_orig = .data[[sample_col]],
           TCR_type     = if_else(str_detect(.data[[sample_col]], aav_pattern), "AAV", "Viral")) %>%
    separate(!!sample_col, into = c("ID", "gene", "Virus", "HLA", "B2M", "MHC", "extra"),
             sep = "_", fill = "right", extra = "merge") %>%
    mutate(Virus = case_when(
      is.na(extra)                 ~ Virus,
      str_detect(extra, "CMV")     ~ "CMV",
      str_detect(extra, "EBV")     ~ "EBV",
      TRUE                          ~ Virus
    )) %>%
    mutate(Peptide = if_else(TCR_type == "Viral",
                             str_split_fixed(ID, ":", 2)[, 2],
                             NA_character_))
}

#' For a categorical column (e.g. gene, HLA, B2M, MHC), find the most
#' frequent non-NA value within each specificity group and join it back as
#' `<value_col>_max`. Generalizes the repeated gene_max/HLA_max/B2M_max/
#' MHC_max blocks in the original script into one reusable step.
summarize_dominant_value <- function(df, value_col, group_col = "index") {
  dominant <- df %>%
    filter(!is.na(.data[[value_col]])) %>%
    count(.data[[group_col]], .data[[value_col]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-n) %>%
    rename(!!paste0(value_col, "_max") := .data[[value_col]])
  
  df %>% left_join(dominant, by = group_col)
}

#' Tabulate antigen/virus hits per specificity group and report the
#' dominant antigen for each group. Replaces the original's fixed list of
#' ~28 hardcoded antigen names with a tabulation over whatever values are
#' actually present in `antigen_col` — no antigen list to maintain, but
#' verify your `antigen_col` values are already clean, discrete labels
#' (that's what split_sample_metadata()'s `Virus` column produces).
count_antigen_hits <- function(df, antigen_col = "Virus", group_col = "index") {
  dominant <- df %>%
    filter(!is.na(.data[[antigen_col]])) %>%
    count(.data[[group_col]], .data[[antigen_col]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(dominant_antigen = .data[[antigen_col]], dominant_antigen_n = n)
  
  df %>% left_join(dominant, by = group_col)
}

# ---- Example usage ------------------------------------------------------
# raw <- read_csv(file.path(gliph_dir, "gliph2_raw_output.csv"))
#
# annotated <- raw %>%
#   filter_significant_groups(score_col = "Fisher_score") %>%
#   annotate_viral_specificity() %>%
#   split_sample_metadata() %>%
#   summarize_dominant_value("gene") %>%
#   summarize_dominant_value("HLA") %>%
#   summarize_dominant_value("B2M") %>%
#   summarize_dominant_value("MHC") %>%
#   count_antigen_hits()
#
# aav_groups   <- annotated %>% filter(TCR_type == "AAV")
# viral_groups <- annotated %>% filter(TCR_type == "Viral")
#
# write_csv(annotated, file.path(results_dir, "gliph2_annotated_groups.csv"))

# ============================================================================
# Step 4 — join significant GLIPH2 specificity groups back onto the full
# single-cell TCR table, flag clonal expansion, and summarize/plot
#
# `sig_groups` here is the output of 03_filter_and_annotate_gliph_output.R.
# `tcr_all` is the same per-cell TCR table used to build the GLIPH2 input in
# step 1, plus a `Group` column (e.g. disease vs. control) you've already
# assigned.
# ============================================================================

source("00_config.R")
library(ggpubr)

#' Within each (TcRb, TcRa, sample-ID) combination, keep only the single
#' best-scoring row (lowest final_score). GLIPH2 can report the same TCR in
#' more than one overlapping specificity group; this collapses to one
#' representative group assignment per TCR.
get_representative_hits <- function(df, id_col = "ID") {
  df %>%
    group_by(TcRb, TcRa, .data[[id_col]]) %>%
    slice_min(final_score, n = 1, with_ties = FALSE) %>%
    ungroup()
}

#' Flag cells in the full TCR table whose beta-chain CDR3 belongs to a
#' significant specificity group.
#'
#' @param tcr_all full per-cell TCR table
#' @param grouped_cdr3b character vector of CDR3b sequences that were part
#'   of a significant, non-"unknown" specificity group
#' @param cdr3b_col column in tcr_all holding the beta-chain CDR3 sequence
flag_grouped_clones <- function(tcr_all, grouped_cdr3b, cdr3b_col = "cdr3_aa2") {
  tcr_all %>%
    mutate(grouped = if_else(.data[[cdr3b_col]] %in% grouped_cdr3b, "yes", "no"))
}

#' Add a per-sample clone count and an Expanded/Non-expanded flag.
compute_expansion <- function(tcr_all, sample_col = "sample", clonotype_col = "CTaa") {
  tcr_all %>%
    add_count(.data[[sample_col]], .data[[clonotype_col]], name = "clone_count") %>%
    mutate(Expansion = if_else(clone_count > 1, "Expanded", "Non-expanded"))
}

#' Summarize the proportion of grouped (or any binary flag) cells per
#' sample, optionally stratified by an extra grouping column (e.g.
#' Expansion). Replaces the three near-duplicate proportion_df blocks in
#' the original script with one reusable summary.
summarize_proportion <- function(tcr_all, group_col = "Group", flag_col = "grouped",
                                 extra_group_col = NULL, positive_value = "yes") {
  grouping_cols <- c("sample", group_col, extra_group_col)
  tcr_all %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      total       = n(),
      n_positive  = sum(.data[[flag_col]] == positive_value),
      proportion  = n_positive / total,
      .groups = "drop"
    )
}

#' Bar + jitter plot of a proportion column across a categorical group,
#' with a Wilcoxon comparison, optionally faceted. Replaces the three
#' near-identical ggplot blocks in the original script.
plot_proportion <- function(df, x_col, y_col = "proportion", facet_col = NULL,
                            comparisons = list(c("AAV", "HC")),
                            ylim = NULL, title = "", ylab = "Proportion") {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    stat_summary(fun = mean, geom = "bar", fill = "white", color = "black", width = 0.6) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
    labs(title = title, x = x_col, y = ylab) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif")
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  if (!is.null(facet_col)) p <- p + facet_wrap(vars(.data[[facet_col]]))
  p
}

#' Join annotated significant-group results back onto the full TCR table.
#'
#' The two sides need a matching join key built the same way on both tables
#' (e.g. paste(alpha CDR3, beta CDR3, subject:condition)). Build
#' `tcr_key_fn` / `sig_key_fn` to match your own `subject:condition`
#' encoding — the two key formats MUST agree exactly, so check
#' `mean(tcr_all$TCR_code %in% sig_groups$TCR_code)` is high before trusting
#' the merge.
#'
#' @param tcr_all full per-cell TCR table
#' @param sig_groups annotated significant-group table (one row per TCR)
#' @param tcr_key_fn function(df) -> character vector, join key for tcr_all
#' @param sig_key_fn function(df) -> character vector, join key for sig_groups
merge_final_annotations <- function(tcr_all, sig_groups, tcr_key_fn, sig_key_fn) {
  tcr_all$TCR_code    <- tcr_key_fn(tcr_all)
  sig_groups$TCR_code <- sig_key_fn(sig_groups)
  
  match_rate <- mean(tcr_all$TCR_code %in% sig_groups$TCR_code)
  message(sprintf("TCR_code match rate: %.1f%% of tcr_all rows found a join key in sig_groups",
                  100 * match_rate))
  
  tcr_all %>% left_join(sig_groups, by = "TCR_code")
}

# ---- Example usage ------------------------------------------------------
# tcr_all <- read_csv(file.path(data_dir, "combined_TCR_table.csv")) %>%
#   compute_expansion()
#
# grouped_cdr3b <- sig_groups %>%
#   filter(!is.na(dominant_antigen)) %>%
#   pull(TcRb) %>% unique()
#
# tcr_all_grouped <- flag_grouped_clones(tcr_all, grouped_cdr3b)
#
# prop_grouped <- summarize_proportion(tcr_all_grouped, flag_col = "grouped")
# plot_proportion(prop_grouped, x_col = "Group", ylim = c(0, 0.15),
#                 title = "Virus-specific TCRs", ylab = "Proportion of grouped TCRs")
#
# prop_by_expansion <- summarize_proportion(tcr_all_grouped, flag_col = "grouped",
#                                            extra_group_col = "Expansion")
# plot_proportion(prop_by_expansion, x_col = "Group", facet_col = "Expansion",
#                 ylim = c(0, 0.3), title = "Grouped TCRs by expansion status")
#
# final_table <- merge_final_annotations(
#   tcr_all_grouped, sig_groups,
#   tcr_key_fn = function(df) paste(df$CTaa, df$sample, sep = "_"),
#   sig_key_fn = function(df) paste(df$TcRa, df$TcRb, df$Sample_orig, sep = "_")
# )
# write_csv(final_table, file.path(results_dir, "final_annotated_TCR_table.csv"))

# ============================================================================
# Step 5 — network of GLIPH2 specificity groups that share CDR3b sequences
#
# Each node is a specificity group (GLIPH2 `index`); an edge connects two
# groups that share at least one CDR3b sequence, weighted by the number of
# shared sequences. Node size reflects total clone frequency in the group;
# node color reflects dominant viral antigen.
# ============================================================================

source("00_config.R")
library(igraph)

default_virus_colors <- c(
  InfluenzaA_total = "darkred", InfluenzaB_total = "orange",
  CMV_total = "darkblue", EBV_total = "darkgreen",
  Sars_CoV_2_total = "purple", HCV_total = "pink",
  unknown = "lightgreen", other = "lightgrey"
)

#' Bucket total clone frequency per specificity group into an ordinal size
#' class for plotting (breaks match the original figure's binning).
bucket_total_freq <- function(total_freq) {
  case_when(
    total_freq == 2                        ~ 1,
    total_freq == 3                        ~ 2,
    total_freq == 4                        ~ 3,
    total_freq == 5                        ~ 4,
    total_freq >= 6  & total_freq <= 10     ~ 5,
    total_freq >= 11 & total_freq <= 20     ~ 6,
    total_freq >= 21 & total_freq <= 30     ~ 7,
    total_freq >= 31 & total_freq <= 40     ~ 8,
    total_freq >= 41 & total_freq <= 50     ~ 9,
    total_freq >= 51 & total_freq <= 60     ~ 10,
    total_freq >= 61 & total_freq <= 100    ~ 11,
    total_freq >= 101 & total_freq <= 200   ~ 12,
    total_freq > 200                        ~ 13,
    TRUE                                     ~ NA_real_
  )
}

#' Build an igraph object of specificity groups connected by shared CDR3b
#' sequences.
#'
#' @param data GLIPH2 group table (one row per TCR), must have `index`,
#'   `Freq`, `TcRb`, plus one column per requested node attribute
#' @param node_attr_cols named vector mapping node attribute name -> source
#'   column, e.g. c(type = "dominant_antigen", gene = "gene_max",
#'   peptide = "Peptide_max")
#' @param color_by which node_attr_cols entry to color nodes by (must be a
#'   name in node_attr_cols); pass NULL to skip coloring
#' @param color_map named vector mapping category -> color for the
#'   color_by attribute (falls back to "lightgrey" for unmapped values)
#' @param subsample_n if set, randomly keep this many nodes (for
#'   readability on dense networks) — omit for the full network
#' @param seed random seed used only when subsample_n is set
build_specificity_network <- function(data,
                                      node_attr_cols = c(type = "dominant_antigen"),
                                      color_by = names(node_attr_cols)[1],
                                      color_map = default_virus_colors,
                                      subsample_n = NULL,
                                      seed = NULL) {
  node_freq <- data %>%
    group_by(index) %>%
    summarise(total_freq = sum(Freq), .groups = "drop")
  
  node_attrs <- purrr::imap_dfc(node_attr_cols, function(col, name) {
    data %>%
      filter(!is.na(.data[[col]])) %>%
      distinct(index, .data[[col]]) %>%
      group_by(index) %>%
      slice(1) %>%
      ungroup() %>%
      select(index, !!name := .data[[col]])
  }) %>%
    select(index, everything()) %>%
    distinct(index, .keep_all = TRUE)
  
  vert <- node_freq %>%
    left_join(node_attrs, by = "index") %>%
    mutate(freq_group = bucket_total_freq(total_freq),
           index = as.character(index))
  
  if (!is.null(subsample_n)) {
    if (!is.null(seed)) set.seed(seed)
    vert <- vert[sample(nrow(vert), min(subsample_n, nrow(vert))), ]
  }
  
  if (!is.null(color_by)) {
    vert$color <- color_map[vert[[color_by]]]
    vert$color[is.na(vert$color)] <- "lightgrey"
  }
  
  # group x CDR3b incidence matrix -> shared-CDR3b counts via crossprod,
  # instead of an O(n^2) loop of pairwise intersect() calls
  membership <- data %>%
    filter(index %in% vert$index) %>%
    distinct(index, TcRb) %>%
    mutate(index = as.character(index))
  
  incidence <- table(membership$index, membership$TcRb)
  shared <- incidence %*% t(incidence)
  diag(shared) <- 0
  
  shared_long <- as.data.frame(as.table(shared)) %>%
    filter(Var1 < Var2, Freq > 0) %>%
    rename(from = Var1, to = Var2, number = Freq)
  
  vert <- vert[vert$index %in% unique(c(shared_long$from, shared_long$to)), ]
  
  graph_from_data_frame(shared_long, directed = FALSE, vertices = vert)
}

#' Plot a specificity-group network to PDF with size and color legends.
#'
#' @param graph output of build_specificity_network()
#' @param output_pdf path to write the PDF
#' @param size_by "freq_group" to scale node size by clone frequency, or
#'   "fixed" for uniform small nodes
#' @param color_legend_title legend title for the color categories
plot_specificity_network <- function(graph, output_pdf,
                                     size_by = c("freq_group", "fixed"),
                                     color_legend_title = "Category") {
  size_by <- match.arg(size_by)
  vert <- igraph::as_data_frame(graph, what = "vertices")
  vertex_size <- if (size_by == "freq_group") vert$freq_group else 2
  
  pdf(output_pdf, width = 10, height = 8)
  plot(graph,
       vertex.size = vertex_size, vertex.label = NA,
       vertex.color = vert$color, vertex.frame.color = NA,
       edge.width = igraph::E(graph)$number, edge.color = "black")
  
  if (size_by == "freq_group") {
    sizes <- sort(unique(vert$freq_group))
    legend(x = 1.2, y = 0.8, legend = paste("Size", sizes),
           pt.cex = sizes * 0.25, pch = 21, col = "black", pt.bg = "grey",
           bty = "n", cex = 1.2, text.col = "black")
  }
  
  legend(x = -1.5, y = 1, legend = unique(vert[[color_legend_title]] %||% vert$type),
         col = unique(vert$color), pch = 19, pt.cex = 1, bty = "n",
         cex = 1, text.col = "black")
  dev.off()
}

# ---- Example usage ------------------------------------------------------
# aav_expanded <- aav_groups %>% filter(Freq > 1, !is.na(dominant_antigen), dominant_antigen != "unknown")
#
# net <- build_specificity_network(
#   aav_expanded,
#   node_attr_cols = c(type = "dominant_antigen", gene = "gene_max", peptide = "Peptide_max"),
#   color_by = "type", color_map = default_virus_colors
# )
# plot_specificity_network(net, file.path(results_dir, "specificity_network.pdf"))

