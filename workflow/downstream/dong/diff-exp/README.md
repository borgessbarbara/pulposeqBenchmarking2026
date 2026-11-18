# Differential Expression Analysis

`01_de_edger.R`: H1975 (TKI-resistant) vs HCC827 (TKI-sensitive), gene and transcript level
(Fig. 2d, e, f).

- Counts: oarfish, via `read_oarfish()` (`../../helpers/oarfish_counts.R`): transcript level
  `catchOarfish()`, gene level `catchOarfishGene()`, both divided by their read-to-transcript
  ambiguity (RTA) overdispersion, not rounded. Sequins removed; rRNA and processed pseudogenes
  dropped before filtering.
- Filter >= 10 reads in >= 3 samples on the undivided counts; TMM; edgeR glmQLFit (robust) ->
  glmTreat (|log2FC| > 0.5), FDR 0.05; display log2FC = predFC (prior.count 2). lncRNA results
  re-adjusted within lncRNAs.
- Every DE table has an `overdispersion` column.
- Outputs (`results/<level>/`): `de_TKI_resistant_vs_TKI_sensitive.csv`, `logcpm.rds`,
  `edger_input_counts.tsv` (the exact matrix edgeR received, read by `isoformic/`),
  sample Spearman correlations, volcano / heatmap / dot plots, `lncrna/`; `results/all_results.csv`.
