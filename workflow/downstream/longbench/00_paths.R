## ===========================================================================
## 00_paths.R - every path used by the LongBench (Fig. 1) scripts, oarfish runs
## Sourced by each script in downstream/longbench/<analysis>/; edit here only.
##
## Inputs  (read only): the pulposeq run folders run_dir[] (the --outdir of the
##                      runs in config-runs/longbench/, made with --run_oarfish true),
##                      downstream/longbench/data/ (sample sheets) and the repository's
##                      references/ (Lexogen SIRV design)
## Outputs (written)  : <analysis>/results/ and <analysis>/figures/ inside
##                      downstream/longbench/. Nothing else is written.
## ===========================================================================

project  <- "C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026"   # the run folders + the clone of this repo
repo     <- file.path(project, "workflow")
base     <- repo                                         # downstream/ and references/ live here
lb_dir   <- file.path(repo, "downstream/longbench")

helper_dir <- file.path(lb_dir, "helpers")
data_dir   <- file.path(lb_dir, "data")                  # samplesheet_*.csv
ref_dir    <- file.path(repo, "references")              # reference files shipped with the repo
sirv_design <- file.path(ref_dir, "SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx")
shared_helper_dir <- file.path(repo, "downstream/helpers")  # oarfish_counts.R, spikein_counts.R

## pulposeq run folders (one per LongBench protocol): the --outdir of each run,
## i.e. <project>/ont_cdna, <project>/ont_drs, <project>/pb
platforms <- c("ont_cdna", "ont_drs", "pb")
run_dir   <- setNames(file.path(project, platforms), platforms)

## pulposeq output files, relative to a run folder (pulposeq transcriptome/ layout)
pp <- list(
  counts_tx       = "transcriptome/curated/counts/BambuOutput_counts_transcript_validated.txt",
  counts_gene     = "transcriptome/curated/counts/BambuOutput_counts_gene_validated.txt",
  full_length     = "transcriptome/curated/counts/BambuOutput_fullLengthCounts_transcript_validated.txt",
  unique          = "transcriptome/curated/counts/BambuOutput_uniqueCounts_transcript_validated.txt",
  oarfish_tx      = "transcriptome/curated/counts/oarfish_counts_transcript_validated.txt",
  validated_gtf   = "transcriptome/curated/annotation/BambuOutput_annotations_validated.gtf",
  novel_gtf       = "transcriptome/curated/annotation/novel_transcripts_validated.gtf",
  novel_meta      = "transcriptome/curated/metadata/novel_transcripts_validated_metadata.csv",
  ann_meta        = "transcriptome/curated/metadata/annotated_transcriptome_metadata.csv",
  read_support    = "transcriptome/curated/metadata/transcript_read_support.tsv",
  bambu_counts_tx = "transcriptome/bambu/BambuOutput_counts_transcript.txt",
  bambu_novel_gtf = "transcriptome/bambu/bambu_novel_transcripts.gtf",
  report_dir      = "transcriptome_report"
)
run_file <- function(plat, key) file.path(run_dir[[plat]], pp[[key]])

## results of each analysis (inside its own folder)
spk_dir   <- file.path(lb_dir, "spike-ins/results")                   # 01, 02
seq_dir   <- file.path(spk_dir, "sequins")                            # 01
sirv_dir  <- file.path(spk_dir, "sirv_ercc")                          # 02
ujc_dir   <- file.path(lb_dir, "unique_intron_chains/results")        # 04, 05 (length/)
supp_dir  <- file.path(ujc_dir, "supp_tables")                        # 07
rs_dir    <- file.path(lb_dir, "read-support/results")                # 06
tab_dir   <- file.path(lb_dir, "figure-panels/results/figure_tables") # 07: numbers behind the panels
panel_dir <- file.path(lb_dir, "figure-panels/figures/fig1_panels")   # 07: one file per panel
oarfish_cache <- file.path(spk_dir, "oarfish_cache")                  # read_oarfish() .rds copies
