# Vendored copy of pbreheny/deepn

This is a vendored copy of the original JAGS-based Bayesian statistics
package for DEEPN, hardwired directly rather than installed at runtime
(`devtools::install_github`), since it's a fixed, unchanging dependency -
no need to anticipate future upstream changes.

Source: https://github.com/pbreheny/deepn
Base commit: 04d69d464f1b8f410d747fe416898ec22a8361ef (last upstream push
2020-05-04), R/ code and inst/model/sModel1.jag taken as-is from this commit.

Vendored: 2026-08-18

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

Contains:
- `R/` - the R functions (`analyzeDeepn`, `importFromDeepn`, `mcmc`/`runMCMC`,
  `psm`, `summary-psm`, etc.) - unmodified from the base commit.
- `inst/model/sModel1.jag` - single-bait JAGS model - unmodified from the
  base commit.
- `inst/model/sModel2.jag` - two-bait JAGS model - reverted to pre-1.6-0
  as described above.

Wired into Stat Maker v3's Bayesian phase (`r_scripts/run_Y2H_bayesian_stats_v3.R`),
installed as a real local R package (`functions/stat_collation.py`'s
`install_deepn_package`) so `system.file(package="deepn")` model-file
lookups inside the unmodified `R/mcmc.R` resolve correctly.
