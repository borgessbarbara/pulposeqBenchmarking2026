# Synthetic RNA controls (Sequins + SIRVs + ERCCs) analysis

LongBench spike-ins, per protocol (ONT cDNA, ONT DRS, PacBio), on the oarfish runs. Each script
answers one or more of four questions, with the spike-in family designed for it:

| Question | Spike-in | Script |
|---|---|---|
| 1. Identification: recovered, nothing invented? | SIRV E0 + long SIRVs, sequins | 02, 01 |
| 2. Abundance: equal amounts measured as equal? | SIRV E0 + long SIRVs; sequin capture | 02, 01 |
| 3. Dose-response | ERCC | 02 |
| 4. Differential abundance | sequins, Mix B vs Mix A | 01 |

Both scripts take an optional argument: `oarfish` (default; results in `results/sequins/` and
`results/sirv_ercc/`, read by `08_fig1_panels.R`) or `bambu` (Bambu's counts of the same runs;
results in `results/sequins_bambu/` and `results/sirv_ercc_bambu/`, read by
`dong/03_quantifier_comparison.R`). Counts come from `../../helpers/spikein_counts.R`.

## Shared definitions
- Detected = >= 3 reads in >= 2 libraries (ERCC: >= 1 read in >= 2 libraries); estimated
  reads rounded per library. `detection_grid.csv` shows SIRV F1 under other rules.
- Counts are EM-estimated reads per transcript (oarfish or Bambu), not molecules.

## `01_sequin_edger.R` - sequin ground truth (Fig. 1d)
- Contrast Mix B (SCLC-A) vs Mix A; expected log2FC from Anaquin's `RnaQuinIsoformMixture`.
- oarfish: `catchOarfish()` / `catchOarfishGene()`; `count_mode = "divided"` (counts / RTA
  overdispersion, primary) and `"undivided"` (sensitivity), both levels. Bambu: `"undivided"`.
- Filter >= 10 reads in >= 2 libraries (`filter_counts`, default on the undivided reads as the
  Dong DE); library sizes and TMM factors from the ratio-1 sequins; edgeR glmQLFit (robust) ->
  glmTreat (H0 |log2FC| <= 0.5), FDR 0.05. Truly DE = |expected log2FC| > 0.5; wrong-sign calls
  are false discoveries.
- Agreement: Spearman rho, R2 / slope with and without sequins absent from one mix, shrunken log2FC.
- Identification on `chrIS`: detection, novel models (Bambu vs curated), 70% detection limit,
  reads per library, capture efficiency vs length.
- Outputs (`results/sequins[_bambu]/`): `edger/<platform>/<level>/<count_mode>/sequin_de.csv`,
  `edger/sequin_de_results.csv`, `edger/fold_change_summary_edger.csv`,
  `edger/sequin_fc_outliers.csv`, `sequin_identification.csv`, `sequin_per_isoform.csv`,
  `sequin_detection_by_conc.csv`, `novel_models_on_sequins.csv`.

## `02_sirv_ercc.R` - SIRV-Set 4 and ERCC (Fig. 1b, c)
- SIRV identification (P / R / F1, Bambu vs curated), accuracy against 1/84 (median across
  libraries), length classes, ERCC dose-response (fit on detected ERCCs) and 70% detection limit.
- oarfish: undivided counts (nothing is tested) plus each SIRV's RTA overdispersion and the
  Spearman correlation of abundance error with it (`sirv_error_od_rho`).
- Outputs (`results/sirv_ercc[_bambu]/`): `sirv_ercc_summary.csv`, `length_bias.csv`,
  `length_class_summary.csv`, `detection_grid.csv`, `ercc_per_transcript.csv`,
  `novel_models_on_spikeins.csv`, `<platform>/sirv_e0_per_isoform.csv`,
  `<platform>/sirv_per_isoform.csv`.

## `dong/03_quantifier_comparison.R` - oarfish vs Bambu (Supplement)
- Needs 01 and 02 with both quantifiers, and `dong/02_sequin_edger.R` with both.
- Outputs (`results/quantifier_comparison/`): `comparison_summary.csv`,
  `sirv_per_isoform_paired.csv`, `open_questions.csv` (R2_18, SIRV612 / SIRV609, PacBio
  gene-level wrong signs); figure `figures/supp_quantifier_sirv.{pdf,png}`.
