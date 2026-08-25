# Garbage/background gene analysis

Scripts behind the Rab GTPase HuORFeome garbage/background gene list
(2026-08-25 analysis). Run from the repo root (`python3
garbage_analysis/<script>.py`) with the project venv active; each expects
the historical SM6 output and reference files under `Dropbox
(endosomeLAB)/SM6_batch_output/` and the raw gene-count/junction data
under `/Volumes/PiperLabDataDisk/DEEPN_2018_RabGTPase`.

- `sm6_batch.py` - shared helpers (canonical dataset list, file-config
  resolution, xlsx->csv conversion) used by the other scripts.
- `rab_hit_recapitulation.py` - loads the prior manually-curated hit list
  (`TAB_RRAB_Hits.xlsx`) and all 36 SM6 datasets; also runs the
  hit-criteria sweep that validated the GUI's default criteria against
  that prior list.
- `vector_garbage_hunt.py` - the vector-alone DESeq2 consistency test (16
  independent Vector-Selected/Non-Selected replicate pairs; also contains
  the leave-one-out pseudo-bait approach, which turned out to be pure
  noise and was not used further).
- `run_bottleneck_statmaker.py` - re-runs the pooled HuORFeome bottleneck
  sample through the full Stat Maker v7 pipeline against one Vector1/
  Vector2 reference pair.
- `run_bottleneck_multi_ref.py` - the same bottleneck run repeated against
  4 more independent reference-vector pairs, to check the result isn't an
  artifact of one particular reference choice.
- `build_unified_criteria.py` - builds `Rab_unified_criteria.xlsx`, the
  revised Rab interactor table under the GUI's locked default criteria,
  with the final garbage/background union list subtracted.

Final locked thresholds: vector-only test >=7/16 independent preparations
(p-value_raw <= 0.01, positive fold), bottleneck test >=2/4 independent
reference-vector comparisons (in-frame:Forward >= 85%, raw ratio or Enr1
>= 2.5x, p-value_raw < 0.1). See `Rab_garbage_union_final.csv` and the
brief/long write-ups in Dropbox for the full method and rationale.
