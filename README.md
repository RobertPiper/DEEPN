# DEEPN

DEEPN (Dynamic Enrichment for Evaluation of Protein Networks) is a workflow
for analyzing yeast two-hybrid (Y2H) next-generation sequencing data. It maps
reads (via the companion tool [MAPster](https://github.com/RobertPiper/MAPster)),
counts reads per gene, extracts junction sequences, and statistically ranks
candidate protein-protein interactions (via the bundled Stat Maker module).

This is a Python 3 / PyQt5, native Apple Silicon port of the original
[emptyewer/DEEPN](https://github.com/emptyewer/DEEPN), described in:

> Krishnamani, V., Peterson, T.A., Piper, R.C., Stamnes, M.A. Informatic
> Analysis of Sequence Data from Batch Yeast 2-Hybrid Screens. J. Vis.
> Exp. (136), e57802, doi:10.3791/57802 (2018).

## Getting a build

Prebuilt macOS app bundles (DEEPN itself, and standalone Stat Maker) are
published as assets on this repo's [Releases](../../releases) page rather
than committed to git.

## Missing data files

Two reference gene-list files are not tracked in this repository, because
they exceed GitHub's 100MB per-file limit:

- `lists/hg38GeneList.prn` (~167MB)
- `lists/mm10GeneList.prn` (~103MB)

Download `lists_data.zip` from the [latest Release](../../releases/latest)
and unzip it into `lists/` before running from source. Prebuilt app bundles
already include these files - this only matters if you're building from
source.

## Running from source

```
python3 -m venv .venv
source .venv/bin/activate
pip install PyQt5 pyqtgraph joblib sortedcontainers xlsxwriter numpy pandas matplotlib
python3 deepn.py
```

See `setup_deepn26.py` for the py2app packaging script used to build the
distributed `.app` bundles.
