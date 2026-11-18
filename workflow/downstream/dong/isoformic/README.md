# Isoformic

`06_isoformic.R`: isoform-level interpretation of the Dong DE with isoformic (Mamede et al. 2025,
NAR Genom Bioinform 7:lqaf176) (Fig. 2j-m, o, p; Extended Data Figs. 5 and 6).

- DE calls and expression from `diff-exp/` (01): DE tables, `logcpm.rds`, and
  `edger_input_counts.tsv` (the RTA-divided oarfish counts edgeR received). The normalisation it
  recomputes as a check now uses those counts as they are (no rounding), so it matches 01 exactly.
- cis pairs from `correlation-cis-1Mb/` (03).
- Genomic context panels from `transcriptome/curated/annotation/annotations_final.gtf` and
  `coverage/*.bw` (pulposeq's plotgardener engine), with the reference GTF given to pulposeq.
- Outputs: `results/` (tables/, plots/, run_info.txt), `figures/` (journal panels).
