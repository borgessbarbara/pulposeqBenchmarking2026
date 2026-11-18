# Synthetic RNA controls (Sequins) analysis with Anaquin

`02_sequin_edger.R`: sequin ground truth for the Dong DE (Fig. 2c). H1975 + Mix A,
HCC827 + Mix B; expected log2FC = log2(Mix A / Mix B).

- Same test, filter and scoring as `01_de_edger.R` and
  `longbench/01_sequin_edger.R` (detection rule, outliers, capture, robust agreement).
- Counts: oarfish, `count_mode = "divided"` (counts / RTA overdispersion, primary, Fig. 2c) and
  `"undivided"` (sensitivity), both levels. The DE filter follows `01_de_edger.R`
  (`filter_counts = "undivided"`).
- `Rscript 02_sequin_edger.R bambu` repeats it on Bambu's counts of the same run
  (`results_bambu/`), for `dong/03_quantifier_comparison.R`.
- Outputs (`results[_bambu]/`): `<level>/<count_mode>/sequin_de.csv`, `sequin_de_results.csv`,
  `summary.csv`, `sequin_fc_outliers.csv`, `sequin_identification.csv`, `sequin_per_isoform.csv`,
  `sequin_detection_by_conc.csv`, `novel_models_on_sequins.csv`.
