# Reproducible TranscriptMetaAnalyst analysis

This folder repeats an RNA-Seq meta-analysis made with TranscriptMetaAnalyst (TMA) 2.0.0, on the demo data bundled with TranscriptMetaAnalyst, with every setting of the run.

## Run it

From this folder:

```
Rscript setup.R      # once: installs this TMA release and checks the package versions
Rscript analysis.R   # repeats the analysis and writes the results to rerun/
```

Expected run time: about 1 minute (the analysis took 4 s in the app; saving the figures adds a little).

## Contents

| File | Content |
|---|---|
| `analysis.R` | The analysis with all settings of the run: analysis methods, plot display options and enrichment choices. Reads `data/`, writes `rerun/`. |
| `setup.R` | Installs the TMA release used and lists any package whose version differs from the original run. |
| `settings.json` | Every analysis and display setting, machine-readable. |
| `data/` | The input data as used in the analysis. |
| `results/` | The tables and figures of the original run. |
| `session_info.txt` | R and package versions of the original run. |

Comparisons: Treated vs Cont.

## Checking the re-run

`rerun/` gets the same files as `results/`. The tables (`*_meta_analysis_results.csv`, `*_DEGs.csv`, `*_heatmap_data.csv`, `*_volcano_data.csv`) should be identical when the package versions match; figures can differ in small details between systems. Enrichment results come from online services and need an internet connection; they can change when a service updates its database.
