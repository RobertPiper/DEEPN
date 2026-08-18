# Vendored copy of pbreheny/deepn

This is an unmodified, hardwired copy of the original JAGS-based Bayesian
statistics package for DEEPN, vendored directly rather than installed at
runtime (`devtools::install_github`), since it's a fixed, unchanging
dependency - no need to anticipate future upstream changes.

Source: https://github.com/pbreheny/deepn
Commit: 04d69d464f1b8f410d747fe416898ec22a8361ef (last upstream push 2020-05-04)
Vendored: 2026-08-18

Contains:
- `R/` - the R functions (`analyzeDeepn`, `importFromDeepn`, `mcmc`/`runMCMC`,
  `psm`, `summary-psm`, etc.)
- `inst/model/sModel1.jag` - single-bait JAGS model
- `inst/model/sModel2.jag` - two-bait JAGS model (the one Stat Maker v3's
  Bayesian phase will use, matching the DESeq2 side's Bait1+Bait2+Vector
  3-way design)

Not yet wired into Stat Maker v3 - this is prep for the MCMC/Bayesian
expansion phase (bundling JAGS itself, checking/installing rjags/runjags,
and calling `analyzeDeepn()` against this vendored package instead of a
network install).
