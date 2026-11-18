# downstream/ - Fig. 1 (LongBench) and Fig. 2 (Dong) on the oarfish runs

These are the downstream analyses of the original (Bambu-only) runs, adapted to the pulposeq
runs in `config-runs/` (`run_oarfish: true`, `oarfish_num_bootstraps: 20`,
`min_total_counts: 5`, pulposeq `transcriptome/` output layout). The original scripts were in
`before-oarfish-integration/` (removed from this branch; see the git history).

## Where the scripts look

| What | Path (set in `longbench/00_paths.R`, `dong/00_paths.R`) |
|---|---|
| project folder (`project`) | `C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/` |
| LongBench run folders (pulposeq `--outdir`) | `<project>/ont_cdna`, `<project>/ont_drs`, `<project>/pb` |
| Dong run folder | `<project>/dong` |
| this repository (clone) | `<project>/workflow/` |
| reference files | `<repo>/references/` (e.g. the Lexogen SIRV design .xlsx) |
| sample sheets | `<repo>/downstream/<dataset>/data/` |
| GENCODE v50 GTF (Dong `06_isoformic.R` only; too large for the repo) | `ref_gtf` in `dong/00_paths.R` |
| results and figures | `<repo>/downstream/<dataset>/<analysis>/{results,figures}/` (git-ignored); assembled figures in `<project>/figures/` |

To move the project: edit `project` in both `00_paths.R`, `BASE` in `../run_all.R` and the
`source(".../00_paths.R")` line at the top of each script (the four root scripts set `project` themselves).
Scripts were renumbered in run order on 14 Nov 2026 (see "Layout and run order"); what changed
in their logic compared with the Bambu-only scripts is listed below.

## What changed compared with the original (Bambu-only) scripts

**Counts that are modelled come from oarfish, with the bootstrap replicates used as edgeR intends.**
`helpers/oarfish_counts.R` (`read_oarfish()`) is the single reader:

| Level | edgeR function | What edgeR receives |
|---|---|---|
| Transcript | `catchOarfish()` | counts / read-to-transcript ambiguity (RTA) overdispersion |
| Gene | `catchOarfishGene()` | counts / gene-level overdispersion (bootstraps summed per gene) |

- The overdispersion is estimated from each sample's `.infreps.pq` bootstraps, over every
  model oarfish quantified, then restricted to pulposeq's curated set. It is >= 1
  (1 = no ambiguity); every DE table reports it in an `overdispersion` column.
- Divided counts are not rounded: they "behave much like negative binomial counts" and go into
  the unchanged edgeR quasi-likelihood test (glmQLFit robust -> glmTreat |log2FC| > 0.5, FDR 0.05)
  (edgeR reference manual; Baldoni et al. 2024a NAR 52:e13; 2024b NARGAB lqae151; 2025 NAR gkaf1305).
- Expression filters (>= 10 reads in >= n samples) are applied to the **undivided read counts**;
  edgeR then tests the divided counts.
- The sequin scripts run both levels on **divided** (primary, figures) and **undivided** counts
  (`count_mode`), so the effect of the RTA division is visible. Their DE filter follows
  `01_de_edger.R` (`filter_counts = "undivided"`); change both together.
- SIRV/ERCC accuracy uses undivided oarfish counts (nothing is tested) and reports each
  isoform's overdispersion, with the Spearman correlation between abundance error and
  overdispersion (assignment errors come with high overdispersion; capture errors do not).
- **Quantifier comparison.** Every spike-in script also runs on Bambu's counts of the same runs
  (`Rscript <script> bambu`, or `@bambu` in `run_all.R`; results in `<folder>_bambu/`), and
  `dong/03_quantifier_comparison.R` sets oarfish and Bambu side by side
  (same libraries and transcripts; Bambu on genome alignments, oarfish on its own transcriptome mapping) for the Supplement.
- **Spike-in definitions** (identical for both quantifiers and both datasets): detected = >= 3
  reads in >= 2 libraries (ERCC: >= 1 read in >= 2 libraries), counts rounded per library;
  SIRV accuracy = median across libraries of log2(observed / expected), expected = 1/84;
  sequin outliers, capture efficiency vs length and robust fold-change agreement are written
  next to the DE tables.
- Read support (full-length / unique reads) is Bambu's, as before. The curated set is now chosen
  on oarfish counts, so transcripts kept on oarfish counts alone have no Bambu reads and fall in
  "neither".

**Paths.** Every pulposeq file is read from the `transcriptome/` layout; `00_paths.R` of each
dataset holds all of them (`run_dir`, `pp`). Results go to `<analysis>/results/` and
figures to `<analysis>/figures/`.

**Sample sheets.** `longbench/data/` holds the run sample sheets from `config-runs/`
(`*_ONTc`, `*_ONTd`, `*_PB`). `dong/data/samplesheet.csv` keeps the group labels the scripts use
(H1975 = TKI-resistant, HCC827 = TKI-sensitive) plus a `cell_line` column.

## Requirements

- **edgeR >= 4.99.12** (Bioconductor devel, released as edgeR 5.0.0), for
  `catchOarfish()` and `catchOarfishGene()`:

  ```r
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(version = "devel")
  BiocManager::install("edgeR")
  ```
- `readr`, `nanoparquet`, `jsonlite` (read by `catchOarfish`), plus the packages each script
  loads (tidyverse, Anaquin, ggrepel, pheatmap, patchwork, circlize, showtext, svglite,
  isoformic, plotgardener, rtracklayer, ...), as before.
- The run folders must contain `transcriptome/oarfish/` with the per-sample `.quant`,
  `.meta_info.json` and `.infreps.pq` files (published by pulposeq when `run_oarfish` is true).

## Layout and run order

Scripts sit directly in `downstream/longbench/` and `downstream/dong/`, numbered in the order
they run; the cross-dataset scripts (`01_`-`04_`) and `run_all.R` sit in the repository root.
Each analysis subfolder (`spike-ins/`, `diff-exp/`, ...) now holds only that analysis's
`results/` and `figures/` (and its README).

Edit `project` in `longbench/00_paths.R` and `dong/00_paths.R`, and `BASE` in `../run_all.R`.
Then `Rscript run_all.R` (from the repository root) runs everything in order, each script in a
fresh R process with its own folder as working directory. Output is printed to the console
only; no log files are written.

| Order | Script | What | Outputs in | Needs |
|---|---|---|---|---|
| longbench | `00_nanocomp.sh`, `00_nanocomp_bam.sh` | NanoComp merges (SLURM, run on the cluster; not in `run_all.R`) | `<project>/<run>/nanocomp/` | - |
| longbench | `01_sequin_edger.R` | sequin DE + identification (also `@bambu`) | `spike-ins/results/sequins{,_bambu}/` | - |
| longbench | `02_sirv_ercc.R` | SIRV E0, length bias, ERCC (also `@bambu`) | `spike-ins/results/sirv_ercc{,_bambu}/` | - |
| longbench | `03_quantifier_plots.R` | oarfish vs Bambu panels (Extended Data Fig. 4 a-h) | `spike-ins/` | 01, 02 (both quantifiers) |
| longbench | `04_ujc_report.R` | unique intron chains | `unique_intron_chains/results/` | - |
| longbench | `05_ujc_length.R` | sharing by transcript length | `unique_intron_chains/results/length/` | 04 |
| longbench | `06_read_support.R` | read-level support | `read-support/results/` | - |
| longbench | `07_supp_tables_ujc.R` | UJC supplementary tables | `unique_intron_chains/results/supp_tables/` | 04 |
| longbench | `08_fig1_panels.R` | Fig. 1 panels | `figure-panels/` | 01-06 |
| longbench | `09_assemble_figures.R` | Fig. 1 + Extended Data Fig. 4 a-h layout | `<project>/figures/` | 03, 08 |
| dong | `01_de_edger.R` | gene + transcript DE | `diff-exp/results/` | - |
| dong | `02_sequin_edger.R` | sequin DE (also `@bambu`) | `sequins/results{,_bambu}/` | - |
| dong | `03_quantifier_comparison.R` | oarfish vs Bambu, spike-ins of both datasets | `../longbench/spike-ins/` | LongBench 01-02, Dong 02 (both quantifiers) |
| dong | `04_cis_circos.R` | potential cis-acting pairs (Fig. 2q, r) | `correlation-cis-1Mb/results/` | 01 |
| dong | `05_cis_lncRNA_top10.R` | novel lncRNAs in cis pairs (Extended Data Fig. 7) | `correlation-cis-1Mb/results_lncRNA_top10/` | 01, 04 |
| dong | `06_isoformic.R` | isoform-level interpretation (Extended Data Figs. 5-6, Fig. 2 profiles and loci) | `isoformic/` | 01, 04 |
| dong | `07_read_support.R` | read-level support | `read-support/results/` | - |
| dong | `08_fig2_panels.R` | Fig. 2 panels (runs twice: before and after 10) | `figure-panels/` | 01, 02, 07; 2nd pass: 06, 10 |
| dong | `09_circos_coverage_biotypes.R` | genome-wide circos of the assembly | `figure-panels/` | - |
| dong | `10_quantifier_panels.R` | Dong quantifier panels | `figure-panels/` | 03, 08 |
| dong | `11_assemble_figures.R` | Fig. 2, Extended Data Figs. 4 (i-l) and 7 | `<project>/figures/` | 04-06, 08 (2nd pass), 10, LongBench 09 |
| root | `01_computational_resources.R` | Extended Data Fig. 1 (critical-path run time) | `downstream/resources/results/` | - |
| root | `02_transcript_structure.R` | Extended Data Fig. 3 | `downstream/transcript-structure/results/` | - |
| root | `03_sequencing_metrics.R` | Extended Data Fig. 2 (sequencing, orientation, alignment per library; tables for ST1) | `downstream/sequencing-metrics/results/` | - |
| root | `04_supplementary_tables.R` | Supplementary Tables (summaries only; per-feature tables are in the repository) | `<project>/pulposeq_supplementary_tables.xlsx` | everything above |

`read_oarfish()` caches what it reads (`spike-ins/results/oarfish_cache/`,
`diff-exp/results/oarfish_cache/`) and re-reads when the oarfish files change; delete the cache
folder to force a fresh read.
