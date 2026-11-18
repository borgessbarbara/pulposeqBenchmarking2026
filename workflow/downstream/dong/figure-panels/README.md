# Fig. 2 panels (Dong)

- `08_fig2_panels.R`: Fig. 2 a-f, one file per panel (`figures/fig2_panels/`), numbers in
  `results/figure_tables/`. Needs `diff-exp/` 01, `sequins/` 02 and `read-support/` 05; sequin
  panels use `count_mode = "divided"`.
- `09_circos_coverage_biotypes.R`: genome-wide circos of RNA coverage and annotated + novel
  transcripts per strand (`figures/circos_coverage/`, `results/circos_coverage/`); reads the
  curated metadata and `coverage/*.bw`.
