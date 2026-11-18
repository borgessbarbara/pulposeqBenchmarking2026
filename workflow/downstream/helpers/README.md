# Shared helpers

`oarfish_counts.R` - `read_oarfish(run_dir, level)`, the one reader every count-based script uses.

`spikein_counts.R` - `read_counts_q(run_dir, level, quantifier, rows)`: spike-in rows from oarfish
(undivided, divided, overdispersion) or from Bambu's `counts` (no bootstraps), curated or not;
`quantifier_arg()` reads `oarfish` / `bambu` from the command line and `quant_dir()` names the
results folder (`<folder>_bambu` for Bambu). Used by the spike-in scripts.

- Transcript level: `edgeR::catchOarfish()`; gene level: `edgeR::catchOarfishGene()` with the
  transcript-to-gene table from Bambu's raw counts (`remove.version = FALSE`, so IDs match
  pulposeq's tables). Both read the per-sample `.quant`, `.meta_info.json` and `.infreps.pq`
  in `<run>/transcriptome/oarfish/`.
- Returns undivided counts, divided counts (counts / RTA overdispersion, not rounded), the
  overdispersion, its prior and, at transcript level, each transcript's gene.
- Restricted to pulposeq's curated set by default (`curated = FALSE` keeps every model).
- Optional cache (`cache_dir`), re-read when any oarfish file is newer.
- Needs edgeR >= 4.99.12 (Bioconductor devel / edgeR 5.0.0), readr, nanoparquet, jsonlite.
