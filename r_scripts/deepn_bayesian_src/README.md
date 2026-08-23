# Vendored copy of pbreheny/deepn

This is a vendored copy of the original JAGS-based Bayesian statistics
package for DEEPN, hardwired directly rather than installed at runtime
(`devtools::install_github`), since it's a fixed, unchanging dependency -
no need to anticipate future upstream changes.

Source: https://github.com/pbreheny/deepn
Base commit: 04d69d464f1b8f410d747fe416898ec22a8361ef (last upstream push
2020-05-04), R/ code and inst/model/sModel1.jag taken as-is from this commit.

Vendored: 2026-08-18

**Package version: 1.5-4, not 1.6-0.** `DESCRIPTION` deliberately reports
1.5-x (not the base commit's own 1.6-0), because `inst/model/sModel2.jag` -
the file that actually determines this package's statistical output - was
reverted to the 1.5-2 state, not kept at 1.6-0. Checked directly: the only
diff anywhere in `R/` between the 1.5-2 tag and this base commit is one
cosmetic plotting-argument rename in `sqrtAxis.R` (`lab=` -> `labels=`,
unrelated to the statistics and not even reachable from the CLI pipeline);
`psi` is referenced nowhere in `R/`, only inside the model file itself - so
labeling the package 1.5-x accurately describes its actual behavior, and
`packageVersion("deepn")` in R will now report a number that matches what
it actually computes. The `-2` to `-4` bump is for the second, independent
local fix below (`overdisp.R`/`applyFilter.R`) - not a claim that any
version `1.5-3` or `1.5-4` ever existed upstream.

## Deliberate exception: inst/model/sModel2.jag is NOT from the base commit

`inst/model/sModel2.jag` (the two-bait model - the one Stat Maker v3 uses
whenever Bait2 is supplied) has been rolled back to the version that existed
**before** commit `0f7afd03c4b98142072df79f7f7b5706fafaa77a` ("Version
1.6-0", 2018-07-18) - i.e. the version in commit `c0d7b59` (2016-06-15),
last touched again in that same window before 1.6-0.

Why: the 1.6-0 commit added a free per-selected-sample nuisance parameter
(`psi`) to every selected condition, including each bait's own selected
sample. Stat Maker v3's locked design gives each bait exactly ONE selected
sample (no bait replicates - only Vector has 2). With only one sample to
inform it, `psi` for that sample is perfectly confounded with the
population-level mean of that bait's `gamma` (selection effect) parameter -
a flat ridge in the posterior, not a peak. Confirmed via Gelman-Rubin
R-hat on real data: 3.75-6.96 with `psi` (should be ~1.0; >1.5 = failed
convergence on 100% of parameters), even after 13x more adapt/burnin
iterations. Removing `psi` (i.e. reverting to the pre-1.6-0 model) gives
R-hat 0.9998-1.0083 (excellent convergence) on the same data.

This also happens to be (very likely) the actual model version used to
generate the historical 2018 DEEPN_2018_RabGTPase StatMaker results, since
the commit landed after that dataset's folder date. Validated directly: a
full re-run against the original 8 gene_count_summary.csv files for that
dataset, using this reverted model, reproduces the historical output's
pBait1/pBait2 >0.9 "confident hit" groups with 100% agreement (51/51 and
39/39 genes, Jaccard 1.0) and AdjEnr1/AdjEnr2 correlate at r=0.9998/0.9996 -
vs. the 1.6-0 model's near-total disagreement (Jaccard 0.004-0.10 on the
same comparison) and failed convergence on the same data.

`inst/model/sModel1.jag` (single-bait model) was never touched by the
1.6-0 commit and needed no change.

## Deliberate exception: R/overdisp.R and R/applyFilter.R are NOT from the base commit

`overdisp()` computes "Baseline" (Vector non-selected dispersion) and
"Selected" (Vector selected dispersion) - printed as e.g. "Baseline (vector
only)" and reported to the user as if they characterize Vector alone. In
the base-commit code, both were computed from `Data$Vector` *after*
`applyFilter()` had already trimmed it down using a joint vector+bait
threshold filter (the non-selected side requires a gene to clear threshold
in Vector's own mean AND every bait's non-selected sample simultaneously;
the selected side requires clearing threshold in Vector OR any bait). That
means these two "vector only" numbers silently depended on which bait
happened to be attached to a given run - confirmed directly: the same
Vector12/Vector13 pair, paired with six different baits (Rab11, Rab17,
Rab19, Rab30, Rab31, Rab32), gave "Selected" dispersion ranging 0.46-0.56
purely from this filtering artifact, since a bait's own selected-condition
abundance was letting genes into the calculation that Vector's own selected
replicates didn't actually support well. This value is not cosmetic - it's
passed straight into the JAGS model as the shared negative-binomial
dispersion for every selected-condition sample, Vector's and both baits'
alike (`mcmc.R`'s `jData$omega`), so an artificially tight/inconsistent
value there directly affects statistical confidence for every gene.

Fixed by computing "Baseline"/"Selected" from `Data$Counts$Vector` - the
original, full, unfiltered Vector counts that `applyFilter()` already
preserves alongside its bait-filtered working copy - using a Vector-only
threshold filter (same logic as the original filter, with every bait
column removed). Verified: the same Vector12/Vector13 pair now gives an
identical 0.0040/0.6259 regardless of which bait accompanies it, across
all six datasets that use this pair. `omega[3]` ("baitEffect") is left
untouched, still computed from the bait-joint-filtered view - that value
is *supposed* to depend on the bait (it exists specifically to compare
Vector's baseline against a bait's own baseline), so its bait-dependence
is correct, not a bug. `applyFilter()` gained one line (`Data$threshold <-
thresh`) purely so the fixed `overdisp()` can reuse the same threshold
value without needing it passed in separately.

Contains:
- `R/` - the R functions (`analyzeDeepn`, `importFromDeepn`, `mcmc`/`runMCMC`,
  `psm`, `summary-psm`, etc.) - unmodified from the base commit, **except**
  `overdisp.R` and `applyFilter.R` as described above.
- `inst/model/sModel1.jag` - single-bait JAGS model - unmodified from the
  base commit.
- `inst/model/sModel2.jag` - two-bait JAGS model - reverted to pre-1.6-0
  as described above.

Wired into Stat Maker v3's Bayesian phase (`r_scripts/run_Y2H_bayesian_stats_v3.R`),
installed as a real local R package (`functions/stat_collation.py`'s
`install_deepn_package`) so `system.file(package="deepn")` model-file
lookups inside the unmodified `R/mcmc.R` resolve correctly.
