#!/usr/bin/env Rscript

# Collated-input version of run_BIGenrichment4.R
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
  # lines: vector of raw lines
  # header_line_idx: index of header line within lines
  # data_start_idx: first data line index
  # data_end_idx: last data line index (inclusive)
  #
  # NOTE: collated_input.csv blocks may have trailing delimiters which can create
  # empty / NA column names. dplyr verbs fail on data frames with NA/"" names,
  # so we sanitize names here.
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
  # sample_names_selected is a single string like "Vector1S;Vector2S" or "BaitS"
  s <- strsplit(sample_names_selected, ";", fixed=TRUE)[[1]]
  s <- s[nchar(s) > 0]
  if (length(s) == 0) return(NA_character_)
  # strip trailing condition letter and trailing replicate number
  base <- s %>%
    str_replace("(S|N)$", "") %>%
    str_replace("\\d+$", "")
  u <- unique(base)
  if (length(u) == 1) return(u[[1]])
  # fallback: longest common prefix
  lcp <- Reduce(function(a,b) {
    i <- 1
    max_i <- min(nchar(a), nchar(b))
    while (i <= max_i && substr(a,i,i) == substr(b,i,i)) i <- i + 1
    substr(a,1,i-1)
  }, u)
  if (nchar(lcp) > 0) return(lcp)
  paste(u, collapse="_")
}

# predefine outputs to avoid "object not found" in edge cases
everything_output <- data.frame(gene=character(), bait=character(), pvalue=numeric(), log2FoldChange=numeric(), stringsAsFactors = FALSE)

# ---------- parse collated file ----------
lines <- lines <- readLines(collated_file, warn=FALSE)

# Markers may be written as 'CONFIG_START' followed by padding commas.
idx_config <- which(str_detect(trimws(lines), "^CONFIG_START"))
idx_data   <- which(str_detect(trimws(lines), "^DATA_START"))

if (length(idx_config) == 0 || length(idx_data) == 0) {
  stop("collated_input.csv must contain CONFIG_START and DATA_START markers", call.=FALSE)
}

idx_config <- idx_config[[1]]
idx_data <- idx_data[[1]]

# manifest block: header at line 1, data from line 2 to just before blank lines / CONFIG_START
# We'll read from line1 .. idx_config-1 and drop blank rows.
manifest <- read_block_csv(lines,
                           header_line_idx = 1,
                           data_start_idx = 2,
                           data_end_idx = idx_config - 1)
manifest <- manifest %>%
  mutate(across(everything(), ~ifelse(is.character(.x), trimws(.x), .x))) %>%
  filter(!(is.na(data_column) | data_column == ""))

# config block: header is line after CONFIG_START, data until line before DATA_START
config_header_idx <- idx_config + 1
config_data_start <- idx_config + 2
config_data_end <- idx_data - 1
config <- read_block_csv(lines,
                         header_line_idx = config_header_idx,
                         data_start_idx = config_data_start,
                         data_end_idx = config_data_end)

# data block: header is line after DATA_START, data until EOF
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
# Convert 'none' to -Inf and ensure types
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
# Preserve gene order exactly as provided in DATA_START
gene_order <- wide$GENE
wide_mat <- wide %>% select(-GENE)

# Ensure numeric
for (cn in colnames(wide_mat)) {
  wide_mat[[cn]] <- suppressWarnings(as.numeric(wide_mat[[cn]]))
  wide_mat[[cn]][is.na(wide_mat[[cn]])] <- 0
}

wide_mat <- as.data.frame(wide_mat, check.names=FALSE)
rownames(wide_mat) <- gene_order

# We will construct raw_counts with columns named like <bait>R<rep><S|N>
# based on CONFIG_START sample_names_selected/background lists.
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

  # map from original group_name columns (e.g., Vector1S) to DESeq2-style names (VectorR1S)
  rename_map <- c()
  if (length(s_names) > 0) {
    for (r in seq_along(s_names)) rename_map[s_names[[r]]] <- paste0(bait, "R", r, "S")
  }
  if (length(n_names) > 0) {
    for (r in seq_along(n_names)) rename_map[n_names[[r]]] <- paste0(bait, "R", r, "N")
  }

  # ensure all required columns exist in DATA_START
  missing_cols <- setdiff(names(rename_map), colnames(wide_mat))
  if (length(missing_cols) > 0) {
    stop(paste0("DATA_START is missing required columns for bait '", bait, "': ",
                paste(missing_cols, collapse=", ")), call.=FALSE)
  }

  # collect column metadata
  sample_names <- c(sample_names, bait)
  s_num_replicates <- c(s_num_replicates, length(s_names))
  n_num_replicates <- c(n_num_replicates, length(n_names))
  cols_list[[bait]] <- rename_map
}

# de-duplicate in case of repeats
keep <- !duplicated(sample_names)
sample_names <- sample_names[keep]
s_num_replicates <- s_num_replicates[keep]
n_num_replicates <- n_num_replicates[keep]

if (length(sample_names) == 0) stop("No baits inferred from CONFIG_START block", call.=FALSE)

# build raw_counts in the same column ordering pattern as the original script
raw_counts <- NULL
for (j in seq_along(sample_names)) {
  bait <- sample_names[[j]]
  rename_map <- cols_list[[bait]]
  # order: all S reps then all N reps (matching original)
  s_cols <- unname(rename_map[grep("S$", names(rename_map))])
  n_cols <- unname(rename_map[grep("N$", names(rename_map))])

  # But rename_map keys are original names; we need to pull from wide_mat by original names in that order
  orig_s <- names(rename_map)[grep("S$", names(rename_map))]
  orig_n <- names(rename_map)[grep("N$", names(rename_map))]

  bait_counts <- wide_mat[, c(orig_s, orig_n), drop=FALSE]
  colnames(bait_counts) <- c(s_cols, n_cols)

  if (is.null(raw_counts)) raw_counts <- bait_counts else raw_counts <- cbind(raw_counts, bait_counts)
}

raw_counts <- raw_counts[rowSums(raw_counts) > 0, , drop=FALSE]
raw_counts <- round(raw_counts)

if (verbose) message("converting counts to integer mode")

# quick NB overdispersion estimate alpha = (Var - Mean) / Mean^2
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

# ---------- normalization (match original logic) ----------
geometric.mean <- function(x, na.rm = TRUE) {
  exp(mean(log(x[x > 0]), na.rm = na.rm))
}

totCounts <- colSums(raw_counts)
if (normalized) {
  normFactor <- rep(1, ncol(raw_counts))
  names(normFactor) <- names(totCounts)
} else {
  normFactor <- totCounts / geometric.mean(totCounts[totCounts > 0], na.rm = TRUE)
}
normalized_counts <- round(raw_counts / rep(normFactor, each = nrow(raw_counts)))

saveRDS(raw_counts, file.path(output_location, "raw_counts_salmon.RDS"))
saveRDS(normalized_counts, file.path(output_location, "normalized_counts_salmon.RDS"))
saveRDS(normFactor, file.path(output_location, "normFactor_salmon.RDS"))

# ---------- build colData (match original) ----------
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

# ---------- run DESeq2 ----------
run_DESeq2 <- function(dds) {
  dds <- DESeq(dds, quiet = !verbose)
  dds
}

dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = cols, design = ~ group)
sizeFactors(dds) <- normFactor[names(normFactor) %in% rownames(cols)]
dds <- run_DESeq2(dds)

if (verbose) {
  mean_disp <- mean(dispersions(dds), na.rm = TRUE)
  if (is.finite(mean_disp)) message(sprintf("📊 Mean overdispersion values by DESeq2 group: %.4f", mean_disp))
}
dds_list <- list(dds)
saveRDS(dds_list, file.path(output_location, "dds.RDS"))

message("✅ Calculations complete. File outputs next")
if (verbose) message("✅ Deleting Temp Files")


# ---------- enrichment scoring + legacy output contract ----------
# Goal: Always produce exactly one retained output by default:
#   everything_combined.csv
# If --debug is enabled, also retain intermediates (everything.csv, Enrichment_only_scores.csv, RDS files).

# Options
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

  # Everything output (all contrasts)
  everything_output <- bind_contrast_tables[, c("gene", "bait", "pvalue", "log2FoldChange")]
  if (isTRUE(debug)) {
    write.csv(everything_output, file = file.path(output_location, "everything.csv"), row.names = FALSE)
  }

  # Enrichment score table (subset + rank-based score)
  enrichment_score <- bind_contrast_tables %>%
    dplyr::filter(pvalue < p_val_thr, log2FoldChange > fc_thr)

  if (nrow(enrichment_score) > 0) {
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

# Enrichment-only CSV (debug-retained only)
if (nrow(enrichment_score_all_rel) > 0) {
  enrichment_output <- enrichment_score_all_rel[, c("gene", "bait", "bait_enrich_score", "pvalue", "log2FoldChange")]
  colnames(enrichment_output) <- c("prey", "bait", "Enrichment_score", "pvalue", "log2FoldChange")
} else {
  enrichment_output <- data.frame(prey=character(), bait=character(), Enrichment_score=numeric(), pvalue=numeric(), log2FoldChange=numeric())
}
if (isTRUE(debug)) {
  write.csv(enrichment_output, file = file.path(output_location, "Enrichment_only_scores.csv"), row.names = FALSE)
}

# Add Enrichment_score back onto everything_output (like legacy)
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

# ORDER/PADDING: pad EVERYTHING back to input gene order (if available)
# gene_order is expected to be defined earlier from the DATA_START table.
if (!is.null(gene_order) && length(gene_order) > 0 && nrow(everything_output) > 0) {
  bait_levels <- unique(everything_output$bait)
  template <- tidyr::expand_grid(gene = gene_order, bait = bait_levels)
  everything_output <- dplyr::left_join(template, everything_output, by = c("gene","bait"))
  everything_output$Enrichment_score[is.na(everything_output$Enrichment_score)] <- 0
  everything_output$gene <- factor(everything_output$gene, levels = gene_order)
  everything_output <- everything_output[order(everything_output$gene, match(everything_output$bait, bait_levels)), ]
  everything_output$gene <- as.character(everything_output$gene)
}

# LEGACY CONTRACT: Always write everything_combined.csv
# Prefer the exact Bait/Vector merge when both are present; otherwise fall back to a wide pivot for all baits.
if (nrow(everything_output) == 0) {
  # still write an empty file with expected columns
  combined <- data.frame(gene=character(), pvalue_bait=numeric(), log2FoldChange_bait=numeric(), Enrichment_score_bait=numeric(),
                         pvalue_vector=numeric(), log2FoldChange_vector=numeric(), Enrichment_score_vector=numeric())
  write.csv(combined, file = file.path(output_location, "everything_combined.csv"), row.names = FALSE)
} else {
  baits_present <- unique(everything_output$bait)
  if (all(c("Bait","Vector") %in% baits_present)) {
    bait_df <- everything_output[everything_output$bait == "Bait", ]
    vector_df <- everything_output[everything_output$bait == "Vector", ]

    combined <- merge(
      bait_df[, c("gene", "bait", "pvalue", "log2FoldChange", "Enrichment_score")],
      vector_df[, c("gene", "bait", "pvalue", "log2FoldChange", "Enrichment_score")],
      by = "gene",
      all = TRUE,
      suffixes = c("_bait", "_vector"),
      sort = FALSE
    )

    # pad combined back to input gene order
    if (!is.null(gene_order) && length(gene_order) > 0) {
      combined <- dplyr::left_join(data.frame(gene = gene_order), combined, by = "gene")
    }

    write.csv(combined, file = file.path(output_location, "everything_combined.csv"), row.names = FALSE)
  } else {
    combined <- everything_output %>%
      tidyr::pivot_wider(names_from = bait, values_from = c(pvalue, log2FoldChange, Enrichment_score))
    if (!is.null(gene_order) && length(gene_order) > 0) {
      combined <- dplyr::left_join(data.frame(gene = gene_order), combined, by = "gene")
    }
    write.csv(combined, file = file.path(output_location, "everything_combined.csv"), row.names = FALSE)
  }
}

# Cleanup: delete intermediates unless debug
if (verbose) message("✅ Deleting Temp Files")
if (!isTRUE(debug)) {
  cleanup_files <- c(
    "raw_counts_salmon.RDS",
    "normalized_counts_salmon.RDS",
    "normFactor_salmon.RDS",
    "cols.RDS",
    "dds.RDS",
    "enrichment_score.RDS",
    "enrichment_score_all_rel.RDS",
    "random_spec_groups.RDS",
    "everything.csv",
    "Enrichment_only_scores.csv"
  )
  for (f in cleanup_files) {
    full_path <- file.path(output_location, f)
    if (file.exists(full_path)) file.remove(full_path)
  }
}

message("✅ Job completed. Final results written to everything_combined.csv")
