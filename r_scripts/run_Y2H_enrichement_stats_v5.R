#!/usr/bin/env Rscript

# v5 of run_Y2H_enrichement_stats.R: switches DESeq2's size-factor
# normalization from simple total-read-count scaling to DESeq2's native
# "poscounts" estimator (via estimateSizeFactors(dds, type="poscounts")).
#
# poscounts computes each gene's normalization reference using only that
# gene's nonzero values, so a gene is not excluded from the calculation
# merely because it happens to be zero in one or two samples - unlike
# DESeq2's classic median-of-ratios ("ratio") estimator, which requires a
# nonzero count in every sample for a gene to count at all. This matters a
# lot for Y2H selected-condition data, where competitive selection drives
# most genes toward zero.
#
# Verified across a 73-dataset scan: total-read-count normalization's hit
# count (p<0.01) correlates strongly with a bait's own "crush strength"
# (casualty rate) at r=-0.73 - i.e. weakly-crushing ("low-crush"/papa-bear)
# datasets get a systematically inflated hit count purely from failing to
# fully deplete background prey, not from real interactions. Switching to
# poscounts drops that correlation to r=0.18 and compresses the hit-count
# range roughly 10-fold (51x range down to ~5x across the same 73
# datasets). The trade-off: on true-null vector-vs-vector comparisons (no
# real biology, any hit is by definition a false positive), poscounts runs
# a modestly higher empirical false-positive rate than total-count
# normalization (~1.7% vs ~1.1% at nominal p<0.01) - judged an acceptable
# cost given the much larger, systematic crush-strength bias it removes,
# especially since real hits are still expected to be filtered downstream
# for genuine absolute PPM enrichment (Bait-Selected > Bait-Non-Selected),
# which independently screens out the "less de-enriched than vector, but
# still not really enriched" pattern that drives low-crush false calls.
#
# All other logic (raw-count construction, DESeq2 contrasts, specificity
# score, everything_combined.csv contract) is unchanged from v4.
#
# Reads a single collated_input.csv that contains:
#   1) a manifest block (top)
#   2) CONFIG_START block (analysis parameters + sample groupings)
#   3) DATA_START block (wide counts table; first column "GENE")

library(optparse)
suppressPackageStartupMessages({
  library(DESeq2)
library(dplyr)
  library(tidyr)
})
library(readr)
library(stringr)


option_list <- list(
  make_option(c("--collated_file"), type="character", default=NULL,
              help="Path to collated_input.csv (contains manifest + CONFIG_START + DATA_START blocks)", metavar="character"),
  make_option(c("--verbose"), action="store_true", default=TRUE,
              help="Print step-by-step progress messages and overdispersion summaries")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$collated_file)) stop("Must supply --collated_file", call.=FALSE)

collated_file <- opt$collated_file
verbose <- isTRUE(opt$verbose)

if (verbose) message("👍 Doing calculations")

# ---------- helpers ----------
read_block_csv <- function(lines, header_line_idx, data_start_idx, data_end_idx, ...) {
  if (data_end_idx < data_start_idx) {
    return(data.frame())
  }
  block_lines <- c(lines[header_line_idx], lines[data_start_idx:data_end_idx])
  con <- textConnection(block_lines)
  on.exit(close(con), add = TRUE)

  df <- read.csv(con, stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE, ...)

  nm <- names(df)
  bad <- is.na(nm) | nm == ""
  if (any(bad)) {
    nm[bad] <- paste0("X", seq_len(sum(bad)))
  }
  nm <- make.unique(nm, sep = ".")
  names(df) <- nm
  df
}


infer_bait_name <- function(sample_names_selected) {
  s <- strsplit(sample_names_selected, ";", fixed=TRUE)[[1]]
  s <- s[nchar(s) > 0]
  if (length(s) == 0) return(NA_character_)
  # Trailing-digit stripping is only for reconciling multiple replicate
  # names of the *same* bait (e.g. "Vector1S;Vector2S" -> shared name
  # "Vector") down to one common prefix. For a single name (Bait1S,
  # Bait2S - one file each, no replicates) the digit is part of the bait's
  # own identity, not a replicate index - stripping it would collapse
  # Bait1 and Bait2 into the same name "Bait" and silently merge them.
  if (length(s) == 1) {
    return(str_replace(s[[1]], "(S|N)$", ""))
  }
  base <- s %>%
    str_replace("(S|N)$", "") %>%
    str_replace("\\d+$", "")
  u <- unique(base)
  if (length(u) == 1) return(u[[1]])
  lcp <- Reduce(function(a,b) {
    i <- 1
    max_i <- min(nchar(a), nchar(b))
    while (i <= max_i && substr(a,i,i) == substr(b,i,i)) i <- i + 1
    substr(a,1,i-1)
  }, u)
  if (nchar(lcp) > 0) return(lcp)
  paste(u, collapse="_")
}

everything_output <- data.frame(gene=character(), bait=character(), pvalue=numeric(), log2FoldChange=numeric(), stringsAsFactors = FALSE)

# ---------- parse collated file ----------
lines <- lines <- readLines(collated_file, warn=FALSE)

idx_config <- which(str_detect(trimws(lines), "^CONFIG_START"))
idx_data   <- which(str_detect(trimws(lines), "^DATA_START"))

if (length(idx_config) == 0 || length(idx_data) == 0) {
  stop("collated_input.csv must contain CONFIG_START and DATA_START markers", call.=FALSE)
}

idx_config <- idx_config[[1]]
idx_data <- idx_data[[1]]

manifest <- read_block_csv(lines,
                           header_line_idx = 1,
                           data_start_idx = 2,
                           data_end_idx = idx_config - 1)
manifest <- manifest %>%
  mutate(across(everything(), ~ifelse(is.character(.x), trimws(.x), .x))) %>%
  filter(!(is.na(data_column) | data_column == ""))

config_header_idx <- idx_config + 1
config_data_start <- idx_config + 2
config_data_end <- idx_data - 1
config <- read_block_csv(lines,
                         header_line_idx = config_header_idx,
                         data_start_idx = config_data_start,
                         data_end_idx = config_data_end)

data_header_idx <- idx_data + 1
data_data_start <- idx_data + 2
data_data_end <- length(lines)
wide <- read_block_csv(lines,
                       header_line_idx = data_header_idx,
                       data_start_idx = data_data_start,
                       data_end_idx = data_data_end)

if (!("GENE" %in% colnames(wide))) {
  stop("DATA_START block must have a first column named 'GENE'", call.=FALSE)
}

# ---------- extract shared parameters ----------
config$enrich_fold_change <- ifelse(tolower(config$enrich_fold_change) == "none", -Inf, as.numeric(config$enrich_fold_change))
config$enrich_p_val <- as.numeric(config$enrich_p_val)
config$normalized <- tolower(as.character(config$normalized)) %in% c("true", "t", "yes", "y", "1")

out_dir <- as.character(config$output_directory[1])
p_val_thr <- config$enrich_p_val[1]
fc_thr <- config$enrich_fold_change[1]
normalized <- config$normalized[1]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
output_location <- file.path(out_dir, "")

# ---------- build raw_counts matrix from DATA_START ----------
gene_order <- wide$GENE
wide_mat <- wide %>% select(-GENE)

for (cn in colnames(wide_mat)) {
  wide_mat[[cn]] <- suppressWarnings(as.numeric(wide_mat[[cn]]))
  wide_mat[[cn]][is.na(wide_mat[[cn]])] <- 0
}

wide_mat <- as.data.frame(wide_mat, check.names=FALSE)
rownames(wide_mat) <- gene_order

sample_names <- c()
s_num_replicates <- c()
n_num_replicates <- c()
cols_list <- list()

for (i in seq_len(nrow(config))) {
  s_names <- strsplit(as.character(config$sample_names_selected[i]), ";", fixed=TRUE)[[1]]
  n_names <- strsplit(as.character(config$sample_names_background[i]), ";", fixed=TRUE)[[1]]
  s_names <- s_names[nchar(s_names) > 0]
  n_names <- n_names[nchar(n_names) > 0]

  bait <- infer_bait_name(config$sample_names_selected[i])
  if (is.na(bait) || bait == "") next

  rename_map <- c()
  if (length(s_names) > 0) {
    for (r in seq_along(s_names)) rename_map[s_names[[r]]] <- paste0(bait, "R", r, "S")
  }
  if (length(n_names) > 0) {
    for (r in seq_along(n_names)) rename_map[n_names[[r]]] <- paste0(bait, "R", r, "N")
  }

  missing_cols <- setdiff(names(rename_map), colnames(wide_mat))
  if (length(missing_cols) > 0) {
    stop(paste0("DATA_START is missing required columns for bait '", bait, "': ",
                paste(missing_cols, collapse=", ")), call.=FALSE)
  }

  sample_names <- c(sample_names, bait)
  s_num_replicates <- c(s_num_replicates, length(s_names))
  n_num_replicates <- c(n_num_replicates, length(n_names))
  cols_list[[bait]] <- rename_map
}

keep <- !duplicated(sample_names)
sample_names <- sample_names[keep]
s_num_replicates <- s_num_replicates[keep]
n_num_replicates <- n_num_replicates[keep]

if (length(sample_names) == 0) stop("No baits inferred from CONFIG_START block", call.=FALSE)

raw_counts <- NULL
for (j in seq_along(sample_names)) {
  bait <- sample_names[[j]]
  rename_map <- cols_list[[bait]]
  s_cols <- unname(rename_map[grep("S$", names(rename_map))])
  n_cols <- unname(rename_map[grep("N$", names(rename_map))])

  orig_s <- names(rename_map)[grep("S$", names(rename_map))]
  orig_n <- names(rename_map)[grep("N$", names(rename_map))]

  bait_counts <- wide_mat[, c(orig_s, orig_n), drop=FALSE]
  colnames(bait_counts) <- c(s_cols, n_cols)

  if (is.null(raw_counts)) raw_counts <- bait_counts else raw_counts <- cbind(raw_counts, bait_counts)
}

raw_counts <- raw_counts[rowSums(raw_counts) > 0, , drop=FALSE]
raw_counts <- round(raw_counts)

if (verbose) message("converting counts to integer mode")

estimate_overdisp <- function(mat) {
  if (is.null(mat) || ncol(mat) < 2) return(NA_real_)
  mu <- rowMeans(mat)
  v <- apply(mat, 1, var)
  alpha <- (v - mu) / (mu^2)
  alpha[mu <= 0] <- NA_real_
  alpha[alpha < 0] <- 0
  mean(alpha, na.rm = TRUE)
}

n_cols_idx <- grepl("N$", colnames(raw_counts))
mean_overdisp_nonselected <- if (sum(n_cols_idx) >= 2) estimate_overdisp(raw_counts[, n_cols_idx, drop=FALSE]) else NA_real_
if (verbose && !is.na(mean_overdisp_nonselected)) {
  message(sprintf("✅ Mean overdispersion for non-selected samples: %.4f", mean_overdisp_nonselected))
}

# ---------- build colData ----------
condition <- c("S", "N")
cols <- data.frame(
  baits = rep(sample_names, s_num_replicates + n_num_replicates),
  conditions = unlist(mapply(function(srep, nrep) rep(condition, c(srep, nrep)),
                             s_num_replicates, n_num_replicates, SIMPLIFY = FALSE)),
  replication = unlist(mapply(function(srep, nrep) if (nrep > 0) c(seq_len(srep), seq_len(nrep)) else seq_len(srep),
                              s_num_replicates, n_num_replicates, SIMPLIFY = FALSE)),
  stringsAsFactors = FALSE
)
rownames(cols) <- colnames(raw_counts)
cols$conditions <- factor(cols$conditions)
cols$replication <- factor(cols$replication)
cols$group <- factor(paste0(cols$baits, cols$conditions))
saveRDS(cols, file.path(output_location, "cols.RDS"))

# ---------- normalization (v5: poscounts instead of total-read-count) ----------
# dds is built first so estimateSizeFactors() has a DESeqDataSet to work
# from. When config$normalized is TRUE (input already normalized upstream),
# every sample gets size factor 1, same as v4's behavior in that case.
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = cols, design = ~ group)

if (normalized) {
  normFactor <- rep(1, ncol(raw_counts))
  names(normFactor) <- colnames(raw_counts)
  sizeFactors(dds) <- normFactor
} else {
  dds <- estimateSizeFactors(dds, type = "poscounts")
  normFactor <- sizeFactors(dds)
}
normalized_counts <- round(raw_counts / rep(normFactor, each = nrow(raw_counts)))

saveRDS(raw_counts, file.path(output_location, "raw_counts_salmon.RDS"))
saveRDS(normalized_counts, file.path(output_location, "normalized_counts_salmon.RDS"))
saveRDS(normFactor, file.path(output_location, "normFactor_salmon.RDS"))

# ---------- run DESeq2 ----------
run_DESeq2 <- function(dds) {
  dds <- DESeq(dds, quiet = !verbose)
  dds
}

dds <- run_DESeq2(dds)

if (verbose) {
  mean_disp <- mean(dispersions(dds), na.rm = TRUE)
  if (is.finite(mean_disp)) message(sprintf("📊 Mean overdispersion values by DESeq2 group: %.4f", mean_disp))
}
dds_list <- list(dds)
saveRDS(dds_list, file.path(output_location, "dds.RDS"))

message("✅ Calculations complete. File outputs next")
if (verbose) message("✅ Deleting Temp Files")


# ---------- enrichment scoring ----------
debug <- FALSE
if (exists("opt") && !is.null(opt$debug)) debug <- isTRUE(opt$debug)

calculate_enrichment_bundle <- function(dds_list, sample_names, p_val_thr, fc_thr, output_location, debug=FALSE) {
  contrast_tables <- list()

  for (g in seq_along(dds_list)) {
    dds <- dds_list[[g]]
    samples <- unique(dds$group)
    samples <- gsub("S|N", "", samples)

    for (sample in samples) {
      if (!(paste0(sample, "S") %in% dds$group && paste0(sample, "N") %in% dds$group)) next
      contrast <- data.frame(results(dds, contrast = c("group", paste0(sample, "S"), paste0(sample, "N"))), cooksCutoff = FALSE)
      contrast <- contrast[!is.na(contrast$pvalue), , drop=FALSE]
      contrast$gene <- rownames(contrast)
      contrast$bait <- sample
      contrast_tables[[paste0(g, sample)]] <- contrast
    }
  }

  bind_contrast_tables <- dplyr::bind_rows(contrast_tables)

  everything_output <- bind_contrast_tables[, c("gene", "bait", "pvalue", "log2FoldChange")]
  if (isTRUE(debug)) {
    write.csv(everything_output, file = file.path(output_location, "everything.csv"), row.names = FALSE)
  }

  enrichment_score <- bind_contrast_tables %>%
    dplyr::filter(pvalue < p_val_thr, log2FoldChange > fc_thr)

  if (nrow(enrichment_score) > 0) {
    # Global rank across every bait's rows together (not grouped per bait) -
    # matches both v2 and the original Y2H-SCORES enrichment_score.R
    # (Software/enrichment_score.R: "enrichment_score$rank <- rank(-enrichment_score$stat)"),
    # not a per-contrast-type score. Kept unchanged here even with a third
    # contrast type (Bait2) added, since that's the validated original design.
    enrichment_score$rank <- rank(-enrichment_score$stat)
    enrichment_score$bait_enrich_score <- (max(enrichment_score$rank) - enrichment_score$rank) / max(enrichment_score$rank)
  } else {
    enrichment_score$bait_enrich_score <- numeric(0)
  }

  if (isTRUE(debug)) {
    saveRDS(enrichment_score, file.path(output_location, "enrichment_score.RDS"))
  }

  list(enrichment_score=enrichment_score, everything_output=everything_output)
}

bundle <- calculate_enrichment_bundle(dds_list, sample_names, p_val_thr, fc_thr, output_location, debug=debug)
enrichment_score_all_rel <- bundle$enrichment_score
everything_output <- bundle$everything_output

if (nrow(enrichment_score_all_rel) > 0) {
  enrichment_output <- enrichment_score_all_rel[, c("gene", "bait", "bait_enrich_score", "pvalue", "log2FoldChange")]
  colnames(enrichment_output) <- c("prey", "bait", "Enrichment_score", "pvalue", "log2FoldChange")
} else {
  enrichment_output <- data.frame(prey=character(), bait=character(), Enrichment_score=numeric(), pvalue=numeric(), log2FoldChange=numeric())
}
if (isTRUE(debug)) {
  write.csv(enrichment_output, file = file.path(output_location, "Enrichment_only_scores.csv"), row.names = FALSE)
}

if (nrow(everything_output) > 0) {
  if (nrow(enrichment_output) > 0) {
    key <- paste(enrichment_output$prey, enrichment_output$bait, sep="\t")
    val <- enrichment_output$Enrichment_score
    names(val) <- key
    everything_output$Enrichment_score <- val[paste(everything_output$gene, everything_output$bait, sep="\t")]
    everything_output$Enrichment_score[is.na(everything_output$Enrichment_score)] <- 0
  } else {
    everything_output$Enrichment_score <- 0
  }
}

# ---------- specificity contrast (new in v3) ----------
# Bait1-selected vs Bait2-selected, direct DESeq2 contrast - only run when
# both baits are actually present. This is the piece Y2H-SCORES calls the
# specificity score (Software/spec_score.R's calc_spec_score()), scaled
# down from its general N-bait combn()-based version to exactly the one
# pairwise comparison this build needs (Bait1 vs Bait2), since locking the
# design to two baits means there's only ever one pair to compare - no
# combinatorial generalization, no kde2d-based fold-change reweighting.
specificity_output <- NULL
if (all(c("Bait1", "Bait2") %in% sample_names)) {
  if (verbose) message("Running Bait1 vs Bait2 specificity contrast...")
  spec <- data.frame(results(dds, contrast = c("group", "Bait1S", "Bait2S")), cooksCutoff = FALSE)
  spec <- spec[!is.na(spec$pvalue), , drop = FALSE]
  spec$gene <- rownames(spec)
  specificity_output <- spec[, c("gene", "pvalue", "log2FoldChange")]
  colnames(specificity_output) <- c("gene", "pvalue_specificity", "log2FoldChange_specificity")
  # log2FoldChange_specificity > 0: prey more enriched under Bait1 selection
  # than Bait2 selection (i.e. more specific to Bait1). < 0: more specific
  # to Bait2. This is a direct pairwise comparison, not a combined score -
  # read alongside each bait's own enrichment columns, not in place of them.
}

# ---------- bait-vs-vector contrast(s) (v4) ----------
# Restores Vector to the specificity comparison set (see file header) - the
# same pairwise-contrast mechanism as the Bait1-vs-Bait2 block above, just
# against "VectorS" instead of the other bait. This is the comparison that
# actually answers "is this gene enriched under this bait, relative to
# background" - what pBaitN_Vec's Bayesian posterior is also estimating.
bait_vector_outputs <- list()
for (bait_id in c("Bait1", "Bait2")) {
  if (all(c(bait_id, "Vector") %in% sample_names)) {
    if (verbose) message(sprintf("Running %s vs Vector contrast...", bait_id))
    bv <- data.frame(results(dds, contrast = c("group", paste0(bait_id, "S"), "VectorS")), cooksCutoff = FALSE)
    bv <- bv[!is.na(bv$pvalue), , drop = FALSE]
    bv$gene <- rownames(bv)
    suffix <- paste0(tolower(bait_id), "_vs_vector")
    out <- bv[, c("gene", "pvalue", "log2FoldChange")]
    colnames(out) <- c("gene", paste0("pvalue_", suffix), paste0("log2FoldChange_", suffix))
    # log2FoldChange_baitN_vs_vector > 0: prey more abundant under this
    # bait's selection than under Vector's selection (enriched vs background).
    bait_vector_outputs[[bait_id]] <- out
  }
}

# ---------- ORDER/PADDING ----------
if (!is.null(gene_order) && length(gene_order) > 0 && nrow(everything_output) > 0) {
  bait_levels <- unique(everything_output$bait)
  template <- tidyr::expand_grid(gene = gene_order, bait = bait_levels)
  everything_output <- dplyr::left_join(template, everything_output, by = c("gene","bait"))
  everything_output$Enrichment_score[is.na(everything_output$Enrichment_score)] <- 0
  everything_output$gene <- factor(everything_output$gene, levels = gene_order)
  everything_output <- everything_output[order(everything_output$gene, match(everything_output$bait, bait_levels)), ]
  everything_output$gene <- as.character(everything_output$gene)
}

# ---------- LEGACY CONTRACT: write everything_combined.csv ----------
# One row per gene. Bait1/Vector columns always present (matches v2's
# Bait/Vector naming, just renamed Bait->Bait1). Bait2 and specificity
# columns only appear when Bait2 data was actually supplied - the same
# collated file with only Bait1+Vector produces the exact same columns v2
# always did, just with "_bait1" instead of "_bait" suffixes.
if (nrow(everything_output) == 0) {
  combined <- data.frame(gene=character(), pvalue_bait1=numeric(), log2FoldChange_bait1=numeric(), Enrichment_score_bait1=numeric(),
                         pvalue_vector=numeric(), log2FoldChange_vector=numeric(), Enrichment_score_vector=numeric())
  write.csv(combined, file = file.path(output_location, "everything_combined.csv"), row.names = FALSE)
} else {
  bait1_df <- everything_output[everything_output$bait == "Bait1", ]
  vector_df <- everything_output[everything_output$bait == "Vector", ]

  combined <- merge(
    bait1_df[, c("gene", "pvalue", "log2FoldChange", "Enrichment_score")],
    vector_df[, c("gene", "pvalue", "log2FoldChange", "Enrichment_score")],
    by = "gene", all = TRUE, suffixes = c("_bait1", "_vector"), sort = FALSE
  )

  if ("Bait2" %in% unique(everything_output$bait)) {
    bait2_df <- everything_output[everything_output$bait == "Bait2", ]
    bait2_df <- bait2_df[, c("gene", "pvalue", "log2FoldChange", "Enrichment_score")]
    colnames(bait2_df) <- c("gene", "pvalue_bait2", "log2FoldChange_bait2", "Enrichment_score_bait2")
    combined <- merge(combined, bait2_df, by = "gene", all = TRUE, sort = FALSE)
  }

  if (!is.null(specificity_output)) {
    combined <- merge(combined, specificity_output, by = "gene", all = TRUE, sort = FALSE)
  }

  for (bait_id in names(bait_vector_outputs)) {
    combined <- merge(combined, bait_vector_outputs[[bait_id]], by = "gene", all = TRUE, sort = FALSE)
  }

  if (!is.null(gene_order) && length(gene_order) > 0) {
    combined <- dplyr::left_join(data.frame(gene = gene_order), combined, by = "gene")
  }

  write.csv(combined, file = file.path(output_location, "everything_combined.csv"), row.names = FALSE)
}

# ---------- Cleanup ----------
if (verbose) message("✅ Deleting Temp Files")
if (!isTRUE(debug)) {
  cleanup_files <- c(
    "raw_counts_salmon.RDS",
    "normalized_counts_salmon.RDS",
    "normFactor_salmon.RDS",
    "cols.RDS",
    "dds.RDS",
    "enrichment_score.RDS",
    "everything.csv",
    "Enrichment_only_scores.csv"
  )
  for (f in cleanup_files) {
    full_path <- file.path(output_location, f)
    if (file.exists(full_path)) file.remove(full_path)
  }
}

message("✅ Job completed. Final results written to everything_combined.csv")
