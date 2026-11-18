## ===========================================================================
## 00_paths.R - every path used by the Dong (Fig. 2) scripts, oarfish run
## Sourced by each script in downstream/dong/<analysis>/; edit here only.
##
## Inputs  (read only): the pulposeq run folder data_dir (the --outdir of
##                      config-runs/dong/run_dong.yml, made with --run_oarfish true),
##                      downstream/dong/data/samplesheet.csv and the reference GTF
## Outputs (written)  : <analysis>/results/ and <analysis>/figures/ inside
##                      downstream/dong/. Nothing else is written.
## ===========================================================================

project    <- "C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026"   # the run folders + the clone of this repo
repo       <- file.path(project, "workflow")
base       <- repo
data_dir   <- file.path(project, "dong")                    # pulposeq --outdir of the Dong run
run_dir    <- c(dong = data_dir)                           # for the scripts shared with LongBench
ref_dir    <- file.path(repo, "references")                # reference files shipped with the repo
## GENCODE v50 GTF given to pulposeq (--annotation): too large for the repo,
## so it stays outside it; point this to your copy
ref_gtf    <- "C:/Users/mbarb/OneDrive/paper-pulposeq/refs/gencode.v50.primary_assembly.annotation.dong.sequins.gtf"

dg_dir     <- file.path(repo, "downstream/dong")
helper_dir <- file.path(dg_dir, "helpers")
shared_helper_dir <- file.path(repo, "downstream/helpers")   # oarfish_counts.R, spikein_counts.R
## columns: sample, group (TKI-resistant / TKI-sensitive), cell_line. The run's own
## samplesheet (config-runs/dong/samplesheet.csv) names the groups by cell line;
## this copy keeps the group labels the scripts use (H1975 = TKI-resistant,
## HCC827 = TKI-sensitive).
sheet      <- file.path(dg_dir, "data", "samplesheet.csv")

## pulposeq output files, relative to the run folder (pulposeq transcriptome/ layout)
pp <- list(
  counts_tx       = "transcriptome/curated/counts/BambuOutput_counts_transcript_validated.txt",
  counts_gene     = "transcriptome/curated/counts/BambuOutput_counts_gene_validated.txt",
  oarfish_tx      = "transcriptome/curated/counts/oarfish_counts_transcript_validated.txt",
  oarfish_gene    = "transcriptome/curated/counts/oarfish_counts_gene_validated.txt",
  validated_gtf   = "transcriptome/curated/annotation/BambuOutput_annotations_validated.gtf",
  final_gtf       = "transcriptome/curated/annotation/annotations_final.gtf",
  ref_tx_gtf      = "transcriptome/curated/annotation/bambu_annotated_transcriptome.gtf",
  novel_gtf       = "transcriptome/curated/annotation/novel_transcripts_validated.gtf",
  novel_meta      = "transcriptome/curated/metadata/novel_transcripts_validated_metadata.csv",
  ann_meta        = "transcriptome/curated/metadata/annotated_transcriptome_metadata.csv",
  read_support    = "transcriptome/curated/metadata/transcript_read_support.tsv",
  bambu_counts_tx = "transcriptome/bambu/BambuOutput_counts_transcript.txt",
  coverage_dir    = "coverage"
)
run_path <- function(key) file.path(data_dir, pp[[key]])

## results of each analysis (inside its own folder)
de_dir     <- file.path(dg_dir, "diff-exp/results")              # 01: <level>/de_*.csv, logcpm.rds, ...
seq_dir    <- file.path(dg_dir, "sequins/results")               # 02: sequin_de_results.csv, summary.csv
cis_dir    <- file.path(dg_dir, "correlation-cis-1Mb/results")  # 04, 05: circos_cis_de_pairs.csv, ...
iso_dir    <- file.path(dg_dir, "isoformic/results")  # 06: tables/, plots/
iso_paper  <- file.path(dg_dir, "isoformic/figures")  # 06: journal panels (Fig. 2j-m, o, p; ED5, ED6)
rs_dir     <- file.path(dg_dir, "read-support/results")          # 07
tab_dir    <- file.path(dg_dir, "figure-panels/results/figure_tables")  # 08: numbers behind the panels
panel_dir  <- file.path(dg_dir, "figure-panels/figures/fig2_panels")  # 08: one file per panel
res_dir    <- file.path(dg_dir, "figure-panels/results")         # 09 writes res_dir/circos_coverage
fig_dir    <- file.path(dg_dir, "figure-panels/figures")         # 09 writes fig_dir/circos_coverage
oarfish_cache <- file.path(de_dir, "oarfish_cache")              # read_oarfish() .rds copies
