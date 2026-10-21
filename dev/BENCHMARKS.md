# Pipeline benchmarks (task 2.7)

Before: the pipeline as it was before the Phase 2 refactor (`dev/bench_baseline_pipeline.R`).
After: `run_tma_pipeline()` on branch `refactor/modules-and-perf`.
Reproduce with `Rscript dev/benchmark.R` from the package root.

Machine: Linux, 12 cores, R 4.6.1, DESeq2 1.52.0, metaRNASeq 1.0.8 (2026-10-20).

## Datasets

- **Demo**: the bundled *Synechocystis* data, 3 studies × 4 samples (2 vs 2), 3,405 genes.
- **Simulated**: 6 studies × 6 samples (3 vs 3), 30,000 genes. Negative-binomial counts, 10% of genes truly DE (|log2FC| = 1.5, same direction in all studies), and a gene-wise baseline shift per study. Seed 2026.

## Results

| Dataset | Pipeline | Median time | Memory allocated¹ | Result kept in memory | Meta-analysis DEGs |
|---|---|---|---|---|---|
| Demo | before | 3.28 s | 533 MB | 6.3 MB | 373 |
| Demo | after, serial | 3.70 s | 521 MB | 2.9 MB | 373 |
| Demo | after, parallel | 2.85 s | n/a | 2.9 MB | 373 |
| Simulated | before | 23.2 s | 3.58 GB | 120.9 MB | 2,600 |
| Simulated | after, serial | 23.2 s | 3.40 GB | 44.7 MB | 2,603 |
| Simulated | after, parallel | 8.34 s | n/a | 44.7 MB | 2,603 |

¹ Allocated in the main R process during one run (`bench::mark`). It cannot be measured for the parallel run, where the work happens in worker processes.

Demo: 3 iterations; simulated: 2 iterations.

## Summary

- **Speed**: running the per-study DE in parallel (BiocParallel, one worker per study) makes the simulated 6-study analysis **2.8× faster** (23.2 s → 8.3 s). The demo has only 3 small studies, so the gain is small there (3.3 s → 2.9 s).
- **Memory**: the app now keeps standardised result tables and one normalised count matrix instead of six full `DESeqDataSet` objects. The result held per analysis shrinks by **63%** on the simulated data (121 MB → 45 MB) and by 54% on the demo (6.3 MB → 2.9 MB).
- **Same results**: on the demo data the refactored pipeline gives exactly the same 373 DEGs (checked by the test suite against a recorded reference). On the simulated data, pre-filtering with `filterByExpr` removed 418 of 30,000 genes that were low in every study. That changed the DEG calls slightly: 2,587 true and 16 false positives, against 2,582 and 18 before (3,000 genes are truly DE).
- **App responsiveness** (not part of the timings above):
  - The Demo tab no longer re-runs the whole pipeline every time it is opened. Before, every visit started a new DESeq2 run (about 3 s for the demo). Now the demo is computed once per server process and shared by all users (`bindCache`).
  - An uploaded RNA-seq analysis runs as a background task (`ExtendedTask` + `future`). The app stays responsive: the click returns in under 1 s while the analysis continues in a separate process.
