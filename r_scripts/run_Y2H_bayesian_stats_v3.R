#!/usr/bin/env Rscript
#
# Bayesian/MCMC statistics for Stat Maker v3, via the vendored, unmodified
# pbreheny/deepn package (see r_scripts/deepn_bayesian_src/README.md).
#
# Bait1 and Bait2 are each a single dataset; Vector is always 2 replicates
# (same design as the DESeq2 side, run_Y2H_enrichement_stats_v3.R). Bait2 is
# optional - runMCMC() itself already switches between the single-bait model
# (sModel1.jag) and two-bait model (sModel2.jag) based on whether a second
# bait is present, so no changes were needed to the vendored code to support
# this - it already matches the locked design.

suppressMessages(library(optparse))

option_list <- list(
  make_option("--vec_sel1", type="character"),
  make_option("--vec_nonsel1", type="character"),
  make_option("--vec_sel2", type="character"),
  make_option("--vec_nonsel2", type="character"),
  make_option("--bait1_sel", type="character"),
  make_option("--bait1_nonsel", type="character"),
  make_option("--bait2_sel", type="character", default=NA_character_),
  make_option("--bait2_nonsel", type="character", default=NA_character_),
  make_option("--threshold", type="integer", default=3),
  make_option("--outfile", type="character"),
  make_option("--msgfile", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressMessages(library(deepn))

has_bait2 <- !is.na(opt$bait2_sel) && !is.na(opt$bait2_nonsel)

params_file <- tempfile(fileext=".params")
lines <- c(
  sprintf("%-25s = %s", "Vector_Selected_1", opt$vec_sel1),
  sprintf("%-25s = %s", "Vector_Non-Selected_1", opt$vec_nonsel1),
  sprintf("%-25s = %s", "Vector_Selected_2", opt$vec_sel2),
  sprintf("%-25s = %s", "Vector_Non-Selected_2", opt$vec_nonsel2),
  sprintf("%-25s = %s", "Bait1_Selected", opt$bait1_sel),
  sprintf("%-25s = %s", "Bait1_Non-Selected", opt$bait1_nonsel),
  sprintf("%-25s = %d", "Threshold", opt$threshold)
)
if (has_bait2) {
  lines <- c(lines,
    sprintf("%-25s = %s", "Bait2_Selected", opt$bait2_sel),
    sprintf("%-25s = %s", "Bait2_Non-Selected", opt$bait2_nonsel))
}
writeLines(lines, params_file)

cat("Running Bayesian/MCMC analysis (analyzeDeepn, JAGS", if (has_bait2) "2-bait model" else "1-bait model", ") ...\n")
analyzeDeepn(infile=params_file, outfile=opt$outfile, msgfile=opt$msgfile, sort=0)
cat("Bayesian analysis complete. Wrote:", opt$outfile, "\n")

# Surface the overdispersion status report to stdout too, so the Python
# side's line-streaming log picks it up without a separate read.
cat("----- Overdispersion report -----\n")
cat(readLines(opt$msgfile), sep="\n")
cat("\n----------------------------------\n")
