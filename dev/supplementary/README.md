# Supplementary analyses for the manuscript

Scripts that produce the supplementary file (`output/TMA_Supplementary.docx`),
the manuscript figures (`output/figures/`) and the values of the manuscript
placeholders (`output/placeholders.md`). Run everything from the package root
(about an hour; needs an internet connection and Chrome or Chromium):

```
NOT_CRAN=true Rscript dev/supplementary/run_all.R
```

| Script | Section | Output |
|---|---|---|
| `00_setup.R` | shared paths, helpers and the batch-effect simulation helpers | |
| `01_batch_case_study.R` | S1 batch-effect case study | `S1_*`, `figures/fig3_batch.*` |
| `02_microarray_case_study.R` | S3 microarray case study (KEGG, UniProt, STRING online) | `S3_*` |
| `03_rnaseq_demo_case_study.R` | S2 RNA-Seq demo case study | `S2_*` |
| `04_performance.R` | S4 performance (machine, benchmarks, app timings) | `S4_*` |
| `05_figures.R` | manuscript figures 1 and 2 | `figures/fig1_interface.*`, `figures/fig2_results.*` |
| `06_bundle_example.R` | S5 example of the ready-to-run bundle | `S5_*` |
| `07_placeholders.R` | values of the `[[...]]` placeholders; counts the test expectations | `placeholders.md`, `placeholders.csv` |
| `TMA_Supplementary.Rmd` | the supplementary file, from the outputs above | `TMA_Supplementary.docx` |

The simulated data come from `dev/simulate_studies.R` (also used by
`dev/benchmark.R`), seed 2026. Figures are 300 dpi PNG and PDF, 17.4 cm wide
(double column). Each script writes a `.log` file to `output/`.
