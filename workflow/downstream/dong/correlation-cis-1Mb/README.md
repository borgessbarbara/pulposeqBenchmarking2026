# Correlation analysis - Novel Transcripts forming cis-acting pairs within 1Mb

`04_cis_circos.R`: pairs of DE transcripts within 1 Mb of each other (novel transcripts and their
neighbours), their expression correlation across samples, and the circos plot.

- Expression: `diff-exp/results/transcript/logcpm.rds` from `01_de_edger.R`, i.e. logCPM of the
  RTA-divided oarfish counts.
- Coordinates: `transcriptome/curated/annotation/bambu_annotated_transcriptome.gtf` and
  `novel_transcripts_validated.gtf`.
- Outputs (`results/`): `cis_pairs_all.csv`, `cis_pairs_correlated.csv`,
  `circos_cis_de_pairs.csv` (read by `isoformic/`), `heatmap_cis_1mb.csv`, circos figures.
