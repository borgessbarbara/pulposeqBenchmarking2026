#!/usr/bin/env Rscript
# =============================================================================
#  06_isoformic.R - isoform-level interpretation of the Dong DE with isoformic
#                   (Fig. 2j-m, o, p; Extended Data Figs. 5 and 6)
#  isoformic: Mamede et al. 2025, NAR Genom Bioinform 7(4):lqaf176,
#  doi:10.1093/nargab/lqaf176.
#  Contrast : TKI_resistant (H1975) vs TKI_sensitive (HCC827) - ONT, Bambu (pulposeq)
#
#  Expression values are 01_de_edger.R's, computed from oarfish counts
#  divided by their read-to-transcript ambiguity (RTA) overdispersion; the
#  divided counts are not rounded (see (d) below).
#
#  Inputs (read only)
#    diff-exp/results/{gene,transcript}/de_TKI_resistant_vs_TKI_sensitive.csv,
#      diff-exp/results/{gene,transcript}/logcpm.rds,
#      diff-exp/results/{gene,transcript}/edger_input_counts.tsv   (01_de_edger.R)
#    correlation-cis-1Mb/results/circos_cis_de_pairs.csv          (04_cis_circos.R)
#    <run>/transcriptome/curated/annotation/annotations_final.gtf, <run>/coverage/*.bw,
#      data/samplesheet.csv, refs/gencode.v50 GTF (pulposeq --annotation)
#  Outputs
#    isoformic/results/
#      plots/01_profiles_panel.pdf  + _PNG/ + _SVG/   profile panel, every event
#      plots/02_profiles.pdf        + _PNG/ + _SVG/   one profile per event
#      plots/03_log2fc.pdf          + _PNG/ + _SVG/   gene + isoform log2FC per event
#      plots/04_event_panels_{PDF,PNG,SVG}/           profile + genomic context
#      plots/05_genomic_context_{PDF,PNG,SVG}/        genomic context alone
#      plots/06_cis_pairs_heatmap.pdf + _PNG/ + _SVG/ cis pairs, per-sample z-scores
#      plots/07_cis_pairs_profiles.pdf + _PNG/ + _SVG/ cis pairs, expression profiles
#      tables/ (00_normalisation_factors ... 23_cis_pairs_used), run_info.txt
#    isoformic/figures/ (journal panels at print size, section 7)
#      locus_tracks_{PDF,PNG,SVG}/<event>.*   genomic context of every event: one
#                                             coverage track per library, log2(depth + 1),
#                                             novel models shaded (EGFR = the Fig. 2 locus)
#      "Extended data figure <5-6>".rds       the ggplots of Extended Data Figs. 5 and 6; their
#                                             .pdf/.svg/.png go to <project>/figures/
#      ed_cis_pairs_heatmap.*                 cis-pair heatmap (Extended Data)
#      fig2_panels/fig2{h,j,k,l}_profile_*.*  Fig. 2j-m as printed (section 7e)
#      (the cis-pair heatmap, Fig. 2q, and the circos, Fig. 2r: 04_cis_circos.R)
#      (each .pdf, .svg, .png, .rds = the ggplot)
#  PNG at 600 dpi; PDF (cairo) and SVG are vector.
#
#  Run      : Rscript 06_isoformic.R    (or source() it from RStudio)
#
#  Deviations from the isoformic vignette (and why):
#   (a) DE calls come from YOUR edgeR `call` column (FDR), not raw p-values.
#   (b) The gene+transcript long table is built here instead of with
#       isoformic::join_DEG_DET(), which hard-codes a GENCODE biotype allow-list
#       and would drop every novel_* transcript (your lncRNA discovery set).
#   (c) log2FC bar plots use plot_gene_log2fc() below, a copy of
#       isoformic::plot_log2FC() styling with one fix: upstream decides whether
#       to use transparency by looking only at the GENE row's DE status, so when
#       the gene is DE, non-significant isoforms are drawn fully opaque.
#   (d) Expression values are the SAME "normalised counts" your edgeR script
#       (01_de_edger.R) plots: read from its
#       diff-exp/results/<level>/logcpm.rds ($norm, rows = $id). If that file is
#       missing they are recomputed exactly as that script does: the divided
#       counts it gave edgeR (edger_input_counts.tsv, not rounded), keep the
#       features it tested (its drop_biotypes + read-count filter), keep.lib.sizes = FALSE,
#       normLibSizes() (TMM), ncnt = cpm(y) * mean(eff) / 1e6. When both exist,
#       the log reports how far the recomputed values are from yours.
#       Bambu's own CPM files are never used (Bambu library size, before your
#       filtering, no TMM).
#   (e) Units are TMM-normalised counts (count scale), not TPM: one long read =
#       one molecule, so no length normalisation is needed.
#   (f) Genomic context panels use pulposeq's plotgardener engine
#       (bin/genomic_context.R) instead of isoformic::plot_genomic_context(),
#       which needs a GENCODE GFF3 + OrgDb and cannot draw novel models. Panels
#       are built from pulposeq's annotations_final.gtf and coverage/*.bw, with
#       isoforms coloured by biotype, faded when not DE, and coverage multiplied
#       by the transcript-level mean(eff) / eff - the same per-sample scaling
#       that turns raw counts into your normalised counts.
#   (g) Nothing in diff-exp/results/ is modified; everything is written to
#       cfg$out_dir and cs$paper_dir.
# =============================================================================

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")   # data_dir, de_dir, cis_dir, iso_dir, iso_paper, ref_gtf

## ---- 0a. CONFIG: events, cis pairs, panels (edit here) ------------------------
cs <- list(
  out_dir        = iso_dir,                                         # 00_paths.R
  cis_pairs_file = file.path(cis_dir, "circos_cis_de_pairs.csv"),  # 04_cis_circos.R
  cis_keep       = c("Annotated protein-coding", "Annotated lncRNA"),

  # One row per event. gene = Bambu gene symbol as in your DE tables.
  #   partner_gene  : add this gene's DE transcripts (e.g. an antisense lncRNA)
  #   extra_tx      : add these transcripts (";"-separated names or IDs)
  #   extra_label   : second label line for extra_tx
  #   novel_matched : add DE novel models gffcompare matched to the gene but
  #                   Bambu quantified under their own gene (e.g. antisense
  #                   intronic lncRNAs); they are labelled with their class
  #   panel         : include in the consolidated profile panel
  #   context       : draw genomic context + profile/context event panel
  # Every event is plotted individually (profile + log2FC) regardless of panel.
  events = tibble::tribble(
    ~gene,      ~title,                                            ~axis,                     ~partner_gene, ~extra_tx,      ~extra_label,        ~novel_matched, ~panel, ~context,
    "EGFR",     "EGFR \u00b7 AS intronic lncRNA",          "EGFR driver locus",       NA,            NA,             NA,                  TRUE,           TRUE,   TRUE,
    "NAPSA",    "NAPSA \u00b7 lineage marker loss",                      "LUAD lineage identity",   NA,            NA,             NA,                  FALSE,          TRUE,   FALSE,
    "CEACAM5",  "CEACAM5 \u00b7 lineage marker loss",                    "LUAD lineage identity",   NA,            NA,             NA,                  FALSE,          TRUE,   FALSE,
    "KRT7",     "KRT7 \u00b7 lineage marker",                            "LUAD lineage identity",   NA,            NA,             NA,                  FALSE,          FALSE,  FALSE,
    # NDUFA12, COX6B1, NDUFA6, ELOC and TPT1 are not switches with oarfish counts
    # (DE isoforms all go one way): retitled and kept out of the Extended Data
    # grid and Fig. 2
    "NDUFA12",  "NDUFA12 \u00b7 isoform loss",                           "OXPHOS isoform switching", NA,           NA,             NA,                  FALSE,          FALSE,  FALSE,
    "COX6B1",   "COX6B1 \u00b7 isoform loss",                            "OXPHOS isoform switching", NA,           NA,             NA,                  FALSE,          FALSE,  FALSE,
    "NDUFA6",   "NDUFA6 \u00b7 isoform gain",                            "OXPHOS isoform switching", NA,           NA,             NA,                  FALSE,          FALSE,  FALSE,
    "NDUFA11",  "NDUFA11 \u00b7 isoform switch",                         "OXPHOS isoform switching", NA,           NA,             NA,                  FALSE,          TRUE,   FALSE,
    "ELOC",     "ELOC \u00b7 isoform loss + AS intronic lncRNA", "Hypoxia (VHL-HIF axis)",  NA,          NA,             NA,                  TRUE,           FALSE,  TRUE,
    # BambuTx861_o is not a MIR210HG isoform (novel model ~250 kb away): it is
    # shown only with its partner, cis pair 11 (MIR210HG-201)
    "MIR210HG", "MIR210HG \u00b7 hypoxia lncRNA",                   "Hypoxia (VHL-HIF axis)",  NA,            NA,             NA,                  FALSE,          FALSE,  FALSE,
    "UBE2I",    "UBE2I / UBE2I-AS1 \u00b7 sense-antisense switch", "lncRNA regulation",      "UBE2I-AS1",   NA,             NA,                  FALSE,          TRUE,   TRUE,
    "SNHG17",   "SNHG17 \u00b7 lncRNA isoforms",                         "lncRNA regulation",       NA,            NA,             NA,                  FALSE,          TRUE,   FALSE,
    "MUC1",     "MUC1 \u00b7 isoform switch",                            "LUAD lineage identity",   NA,            NA,             NA,                  FALSE,          TRUE,   FALSE,
    # novel intronic lncRNAs (BambuTx744_i, BambuTx1164_i) up while the RAB11B-AS1
    # gene is flat; sense gene RAB11B (partner, panel below) goes down
    "RAB11B-AS1", "RAB11B-AS1 / RAB11B \u00b7 novel intronic lncRNAs", "Novel model, known locus", "RAB11B", NA,         NA,                  TRUE,           TRUE,   TRUE,
    "TPT1",     "TPT1 \u00b7 isoform loss incl. NMD",                "Translation / survival",  NA,            NA,             NA,                  FALSE,          FALSE,  FALSE,
    # lineage plasticity markers: every tested isoform
    # is shown (cs$profile_all_isoforms), kept out of the Ext. Data grid
    # (cs$panel_exclude); ASCL1 is not expressed (not tested) -> logged and skipped
    "ASCL1",    "ASCL1 \u00b7 neuroendocrine marker",                   "Lineage plasticity",      NA,            NA,             NA,                  FALSE,          FALSE,  FALSE,
    "POU2F3",   "POU2F3 \u00b7 tuft-cell marker",                       "Lineage plasticity",      NA,            NA,             NA,                  FALSE,          FALSE,  FALSE,
    "HMGA2",    "HMGA2 \u00b7 EMT marker",                              "Lineage plasticity",      NA,            NA,             NA,                  FALSE,          FALSE,  FALSE,
    # optional, EGFR amplicon neighbour (co-amplified with EGFR):
    "SEC61G", "SEC61G \u00b7 EGFR-amplicon neighbour", "EGFR driver locus", NA, NA, NA, FALSE, FALSE, FALSE,
    # ---- candidates for review: switches = DE isoforms
    # in both directions, and DE novel lncRNA / non-coding models in known loci.
    # Profiles + log2FC bars only (results/plots/02_profiles*, 03_log2fc*): kept
    # out of the Ext. Data grid by cs$review_events. Promote one by removing it
    # from cs$review_events (and set panel / context).
    "CYTOR",   "CYTOR \u00b7 lncRNA isoform switch",           "lncRNA regulation",      NA, NA, NA, FALSE, FALSE, FALSE,
    "CFLAR",   "CFLAR (c-FLIP) \u00b7 isoform switch",         "Apoptosis resistance",   NA, NA, NA, FALSE, FALSE, FALSE,
    "FHL2",    "FHL2 \u00b7 isoform switch",                   "EMT / signalling",       NA, NA, NA, FALSE, FALSE, FALSE,
    "SLC3A2",  "SLC3A2 (CD98hc) \u00b7 isoform switch",        "Nutrient transport",     NA, NA, NA, FALSE, FALSE, FALSE,
    "CD44",    "CD44 \u00b7 isoform switch",                   "EMT / stemness",         NA, NA, NA, FALSE, FALSE, FALSE,
    "NDUFS8",  "NDUFS8 \u00b7 isoform switch",                 "OXPHOS isoform switching", NA, NA, NA, FALSE, FALSE, FALSE,
    "PKM",     "PKM \u00b7 isoform switch",                    "Glycolysis",             NA, NA, NA, FALSE, FALSE, FALSE,
    "VEGFA",   "VEGFA \u00b7 isoform switch",                  "Hypoxia (VHL-HIF axis)", NA, NA, NA, FALSE, FALSE, FALSE,
    "LY6D",    "LY6D \u00b7 novel non-coding gain",            "Novel model, known locus", NA, NA, NA, FALSE, FALSE, FALSE,
    "SCIRT",   "SCIRT \u00b7 switch to novel lncRNA",          "Novel model, known locus", NA, NA, NA, FALSE, FALSE, FALSE,
    "ABL1",    "ABL1 \u00b7 switch + novel intronic lncRNAs",  "Novel model, known locus", NA, NA, NA, TRUE, FALSE, FALSE,
    "TRAF3",   "TRAF3 \u00b7 novel antisense lncRNA",          "Novel model, known locus", NA, NA, NA, TRUE, FALSE, FALSE
  ),
  review_events = c("CYTOR", "CFLAR", "FHL2", "SLC3A2", "CD44", "NDUFS8", "PKM",
                    "VEGFA", "LY6D", "SCIRT", "ABL1", "TRAF3"),
  # not switches with oarfish counts: own profiles only, out of
  # the consolidated panel / Ext. Data grid (cs$panel_all = TRUE ignores ~panel)
  retired_events = c("NDUFA12", "COX6B1", "NDUFA6", "ELOC", "TPT1"),

  profile_de_only  = TRUE,          # profiles: gene + DE isoforms (+ added transcripts)
  profile_all_isoforms = c("ASCL1", "POU2F3", "HMGA2"),   # these: every tested isoform
                                    # (POU2F3 has no DE isoform, HMGA2 only one)
  panel_exclude    = c("ASCL1", "POU2F3", "HMGA2"),       # not in 01 panel / Ext. Data grid
  max_isoforms     = 10,            # isoforms per profile / bar plot
  # partner gene (e.g. UBE2I-AS1 for UBE2I): drawn as its own profile panel
  # BELOW the gene, and both panels show the gene line + their top pair_max_iso
  # isoforms (DE first, |log2FC|, strongest of each direction always kept).
  pair_max_iso     = 5,
  min_own_iso      = 6,             # single-gene events: own isoforms never drop below this
  panel_ncol       = 2,             # consolidated panel: columns
  panel_all        = TRUE,          # profile panels (01 + Ext. Data): every event, not only panel = TRUE
  panel_cell       = c(3.7, 3.0),   # in, per profile in the panel
  panel_base       = 10,            # pt, panel text
  panel_title_wrap = 50,            # characters per title line in the panel
  inset_size       = c(5.2, 4.4),   # in, profile inside an event panel
  inset_base       = 8,            # pt, event-panel profile text

  # journal panels (section 7): drawn at print size for the Brief Communication
  run_paper        = TRUE,
  paper_dir        = iso_paper,                       # 00_paths.R
  ed_dir           = file.path(project, "figures"),  # Extended Data Figs. 5, 6 (.pdf/.svg/.png), with the other figures
  # genomic context panels (section 7a): one coverage track per library,
  # log2(depth + 1), novel models' spans shaded
  coverage_loci    = "all",                          # locus_tracks_*/: "all" events or a vector of names
  paper_locus_tracks_mm = 119.5,                     # locus_tracks_*/: tracks only, this wide
  paper_loci       = c("EGFR", "UBE2I"),  # locus_samples_*/: profile + tracks
  paper_locus      = list(width_mm  = 171.5,         # profile + locus = full Fig. 2 row (174 mm) - 2.5 mm letter
                          font_pt   = 6,             # Nature: 5-7 pt
                          max_tx    = 8,             # isoforms per locus (DE first)
                          sig_h     = 0.22,          # in, coverage track height
                          tx_box_mm = 1.6,           # exon box height
                          header = FALSE, ideogram = FALSE, legend = TRUE),
  paper_groups     = c(TKI_sensitive = "HCC827", TKI_resistant = "H1975"),
  paper_inset_mm   = 52,            # mm, DE isoform profile left of each locus (of the 171.5)
  paper_inset_min_h_mm = 40,        # mm, profile never shorter than this
  paper_font       = 7,  # pt, Extended Data panels (isoform labels 6.5)
  paper_dpi        = 600,
  paper_ed_width   = 183,           # mm, width of the profile grid (Nature double column);
                                    # a right-hand legend adds paper_legend_w
  paper_ed_ncol    = 4,             # event profiles per row (13 events -> 4 rows)
  paper_ed_row_h   = 48,            # mm per row of profiles
  paper_ed_title_wrap = 30,  # characters per profile title line (about 31 fit one panel; 42 clipped long titles)
  paper_cis_ncol   = 4,             # cis-pair profiles per row
  paper_cis_row_h  = 45,            # mm per row of cis-pair profiles
  paper_legend_pos = "right",       # profile grids without a free cell: legend "right" (top-aligned
                                    # column), "top" or "bottom"
  paper_legend_w   = 34,            # mm added to the page width for a right-hand legend
  paper_profiles_combined = TRUE,   # events + cis pairs in one Ext. Data grid
  paper_ed_all_events = TRUE,       # Ext. Data grid: every drawn event (FALSE = consolidated-panel set)
  paper_page_mm    = c(297, 420),  # Ext. Data profile pages: A3 portrait
  paper_page_ncol  = 5,             # profile columns; rows and the box follow from the page
  paper_page_nrow  = NA,  # unused: rows = what the events need
  paper_box_max_h  = 42,  # mm, box height cap
  paper_box_mm     = c(36.5, 42),   # NULL = the largest box that fits the events page;
                                    # or c(w, h) mm. Also the Fig. 2 profiles' box (7e).
                                    # Fixed at the box of the published Fig. 2 (17 Nov), so
                                    # changing the ED5 / ED6 headings does not resize Fig. 2j-m
  profiles_only    = nzchar(Sys.getenv("ISOFORMIC_PROFILES_ONLY")),  # skip the genomic-context
                                    # sections 5 and 7a (redraw the profiles only)
  paper_all_ncol   = 4,             # combined grid: 24 profiles = 4 x 6, legend below
  paper_all_row_h  = 52,            # mm per row of the combined grid (45 x 52 mm cells)
  paper_hm_tile_mm = 2.2,           # mm, width of one heatmap tile (Ext. Data cis heatmap)
  paper_hm_row_mm  = 2.6            # mm, height of one heatmap row (fits a 6 pt name)
)

## ---- 0b. ENGINE CONFIG ---------------------------------------------------------
## (the base settings; the overrides for this run are applied right below the
##  list, with modifyList)
cfg <- list(
  dong_dir     = data_dir,             # pulposeq run folder; relative paths below resolve here
  # outputs of 01_de_edger.R (read only, never modified)
  edger_dir    = de_dir,
  deg_file     = file.path(de_dir, "gene/de_TKI_resistant_vs_TKI_sensitive.csv"),
  det_file     = file.path(de_dir, "transcript/de_TKI_resistant_vs_TKI_sensitive.csv"),
  samplesheet  = sheet,               # 00_paths.R (data/samplesheet.csv)
  # the counts 01_de_edger.R gave edgeR: oarfish counts / RTA overdispersion
  tx_counts_file   = file.path(de_dir, "transcript/edger_input_counts.tsv"),
  gene_counts_file = file.path(de_dir, "gene/edger_input_counts.tsv"),
  # your plot theme; only `group_cols` is taken from it (NULL/missing = defaults)
  plot_theme   = file.path(helper_dir, "plot_theme.R"),
  out_dir      = iso_dir,

  group_col    = "group",              # "-" is converted to "_", as in your edgeR script
  control      = "TKI_sensitive",  # denominator of the contrast (HCC827)
  case         = "TKI_resistant",  # numerator (H1975)

  # log2FC = predFC(prior.count = 2), the display fold change your plots use;
  # log2FC_test = the glmTreat coefficient
  lfc_col      = "log2FC",
  alpha        = 0.05,             # FDR cut-off; identical to your `call` column
  # NOTE: your calls come from glmTreat(lfc = 0.5), i.e. the fold-change
  # threshold is already inside the test. A second |log2FC| filter would
  # double-filter, so keep 0 to reproduce your `call` column exactly.
  lfc_cut      = 0,

  genes_manual = c("EGFR", "MET", "AXL"),  # gene sets of section 10 (log only:
  n_auto_genes = 15,                       # their plots were switched off)

  # Genomic context panels (pulposeq engine) - drawn for the events (sections 5, 7)
  run_context    = TRUE,
  final_gtf      = run_path("final_gtf"),     # transcriptome/curated/annotation/annotations_final.gtf
  bigwig_dir     = run_path("coverage_dir"),  # coverage/ (*.bw, one per sample)
  # REQUIRED for genomic context: the exact GTF pulposeq was given as --annotation
  # (gffcompare class codes were computed against it)
  ref_gtf        = ref_gtf,     # 00_paths.R
  context_max_tx = 20,    # isoforms per panel (DE first, then |log2FC|)
  # reference import is cached here (NULL = R's per-user cache dir, outside OneDrive)
  ref_cache_dir  = NULL,

  group_colors   = NULL,  # NULL = group_cols from plot_theme, else blue / red

  # Figure style. "poster": every figure is drawn at (about) its final print
  # size, 6-8 in wide, with a 14 pt base font - text stays >= 11 pt when
  # printed at that size. Crowded plots are trimmed (isoforms per gene, rare
  # biotypes, pathways, overview labels). "report": the earlier large pages.
  fig_style          = "poster",
  legend_position    = "right",   # "right" or "bottom" (poster mode)
  base_size          = NULL,  # NULL = 14 (poster) / 11 (report)
  poster_max_iso     = 12,    # isoforms per log2FC bar plot (all DE kept first)
  poster_min_biotype = 10,    # biotypes with fewer DE transcripts -> "other"
  poster_ctx_font    = 1.4,   # genomic context: font scale
  poster_ctx_width   = 7,     # genomic context: panel width (in)
  # genomic context: draw validated models that were never DE-tested
  # ("| not tested"); FALSE keeps panels clean
  ctx_show_untested  = FALSE,
  # profile plots: space left of the first condition for the isoform labels,
  # in category widths (raise if long names get clipped, lower for less white)
  profile_label_room = 0.70,
  # genomic context: chromosomes kept from annotations_final.gtf (and the
  # reference): chr1-22, chrX, chrM. Drops unplaced / unlocalized scaffolds
  # (chrUn_*, *_random, KI*/GL*), chrY (female samples; PAR copies) and the
  # sequin decoy chromosome chrIS
  ctx_keep_chrom     = "^(chr)?([1-9]|1[0-9]|2[0-2]|X|M|MT)$",

  dpi                = 600,          # PNG resolution (all figures)
  ctx_coverage       = "samples",    # "samples" (one track each) or "group_mean"
  ctx_cov_scale      = "linear",     # coverage y axis: "linear" or "log2" (log2(depth + 1):
                                     # low peaks, e.g. intronic lncRNAs, stay visible)
  ctx_cov_highlight  = FALSE,        # shade the span of the novel models on the coverage tracks
  ctx_show_non_de    = TRUE,         # genomic context: draw tested but non-DE models
  ctx_label_style    = "full",       # "full", "short" (poster: no class code / Bambu id) or
                                     # "strand" (name (+/-) + class, no log2FC; journal panels)
  ctx_formats        = c("PNG", "SVG"),  # any of "PDF", "PNG", "SVG"
  # journal panels (set in section 7 below; leave NULL here)
  ctx_paper          = NULL,   # list(width_mm, font_pt, ...): context panel at print size
  ctx_group_labels   = NULL,   # condition names on coverage tracks and profile x axes,
                               # e.g. c(TKI_sensitive = "HCC827")
  lab_mm             = NULL,   # repel label size (mm); NULL = 3.6 poster / 2.4 report
  profile_compact    = FALSE,  # profiles: tighter label repel (journal panels)
  profile_label_side = "auto",       # profile labels: "auto" (less crowded side), "left", "right"
  profile_sd         = "triangles",  # per-condition SD: "triangles" (isoformic: mean +/- SD ends)
                                     # or "bar" (thin bar only)
  profile_cap_note   = TRUE    # subtitle notes "10/19 isoforms" when the isoform cap applies
)

# the values set for this run, applied over the base settings
cfg <- utils::modifyList(cfg, list(
  run_context       = TRUE,
  fig_style         = "poster",
  legend_position   = "right",
  poster_max_iso    = cs$max_isoforms,
  dpi               = 300,
  ctx_coverage      = "group_mean",   # one mean track per condition
  ctx_show_non_de   = FALSE,          # context: DE models + referenced reference only
  ctx_label_style   = "strand",       # locus labels: name (+/-) + class, no log2FC
                                      # (the profile next to the locus carries it)
  ctx_formats       = c("PDF", "PNG", "SVG")
))

## ---- 1. Packages ------------------------------------------------------------
if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_missing <- function(run_context) {
  if (!requireNamespace("isoformic", quietly = TRUE)) {
    install.packages(
      "isoformic",
      repos = c("https://luciorq.r-universe.dev", getOption("repos"))
    )
  }
  cran <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "purrr", "tibble", "svglite")
  for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("edgeR", update = FALSE, ask = FALSE)
  }
  if (run_context) {
    bioc <- c("plotgardener", "rtracklayer", "GenomicRanges", "IRanges",
              "S4Vectors", "BiocGenerics", "GenomeInfoDb", "txdbmaker")
    miss <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(miss, update = FALSE, ask = FALSE)
    }
  }
}
install_missing(cfg$run_context)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(purrr)
})

if (utils::packageVersion("isoformic") < "0.1.2") {
  stop("isoformic >= 0.1.2 is required - reinstall from https://luciorq.r-universe.dev")
}

## ---- 2. Helpers -------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

log_step <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
}

# absolute paths are kept; relative paths are resolved against dong_dir
resolve <- function(x) {
  if (grepl("^([A-Za-z]:|/|\\\\)", x)) x else file.path(cfg$dong_dir, x)
}

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

# keep only the chromosomes matching cfg$ctx_keep_chrom (logs what was dropped)
keep_chroms <- function(gr, what) {
  chr  <- as.character(GenomeInfoDb::seqnames(gr))
  keep <- grepl(cfg$ctx_keep_chrom, chr)
  if (any(!keep)) {
    dropped <- sort(table(chr[!keep]), decreasing = TRUE)
    log_step("  %s: dropped %d records on %d excluded sequences (%s%s)", what, sum(!keep),
             length(dropped), paste(head(names(dropped), 4), collapse = ", "),
             if (length(dropped) > 4) ", ..." else "")
  }
  gr <- gr[keep]
  GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- unique(chr[keep])
  gr
}

regex_escape <- function(x) gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", x, perl = TRUE)

# muffle only the expected many-to-many join warning inside isoformic internals
quiet_joins <- function(expr) {
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      if (grepl("many-to-many|multiple rows", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Rename your edgeR output to the column names isoformic expects
standardise_de <- function(df, level = c("gene", "transcript"), lfc_col) {
  level <- match.arg(level)
  need <- c("feature", "id", "biotype", lfc_col, "pvalue", "padj", "novel", "call")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop(sprintf("%s DE table is missing column(s): %s", level, paste(miss, collapse = ", ")))
  }
  out <- df |>
    mutate(.lfc = .data[[lfc_col]]) |>
    select(
      id, feature, biotype, .lfc, pvalue, padj, call, novel,
      any_of(c("category", "baseMean", "logCPM", "zero_in_one_group"))
    ) |>
    rename(log2FC = .lfc, qvalue = padj) |>
    mutate(is_de = call != "NS", novel = as.logical(novel))
  new_names <- if (level == "gene") {
    c("gene_id", "gene_name", "gene_type")
  } else {
    c("transcript_id", "transcript_name", "transcript_type")
  }
  names(out)[match(c("id", "feature", "biotype"), names(out))] <- new_names
  out
}

# Pick the Bambu file (raw counts) that covers the DE features, preferring the
# most refined pulposeq stage when several cover them equally well
find_bambu_file <- function(dir, pattern, de_ids, what) {
  cands <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  cands <- cands[!grepl("[/\\\\]work[/\\\\]", cands)]  # skip Nextflow work dirs
  if (!length(cands)) stop("No ", what, " found under ", dir, " - set its path in cfg.")
  cov <- vapply(cands, function(f) {
    ids <- readr::read_tsv(f, col_select = 1, show_col_types = FALSE, progress = FALSE)[[1]]
    max(mean(de_ids %in% ids), mean(strip_version(de_ids) %in% strip_version(ids)))
  }, numeric(1))
  prio <- ifelse(grepl("validated", cands), 1L, ifelse(grepl("filter", cands), 2L, 3L))
  ok <- which(cov >= max(cov) - 0.005)
  pick <- ok[order(prio[ok])][1]
  message("  ", what, " candidates (fraction of DE features found):")
  message(paste0("    ", sprintf("%.4f", cov), "  ", cands, collapse = "\n"))
  cands[pick]
}

# Normalised counts exactly as 01_de_edger.R computes them, for the features it
# tested: the divided counts as given (not rounded) -> DGEList (lib.size = colSums
# of the kept rows, i.e. keep.lib.sizes = FALSE) -> normLibSizes() (TMM) ->
# cpm(y) * mean(eff) / 1e6
edger_norm <- function(counts) {
  y <- edgeR::DGEList(counts = counts)
  y <- if (exists("normLibSizes", envir = asNamespace("edgeR"))) {
    edgeR::normLibSizes(y)
  } else {
    edgeR::calcNormFactors(y)   # same TMM, pre-edgeR-4 name
  }
  eff <- y$samples$lib.size * y$samples$norm.factors
  list(
    norm    = edgeR::cpm(y) * mean(eff) / 1e6,
    samples = data.frame(sample = colnames(counts),
                         lib.size = y$samples$lib.size,
                         norm.factors = y$samples$norm.factors,
                         eff.lib.size = eff,
                         rel.scale = mean(eff) / eff)
  )
}

# Prefer the normalised counts your edgeR run saved (logcpm.rds: $norm, $id);
# report how far the recomputed matrix is from them
use_saved_norm <- function(recomp, rds_path, level) {
  if (!file.exists(rds_path)) {
    log_step("  %s: %s not found - using recomputed normalised counts", level, rds_path)
    recomp$source <- "recomputed"
    return(recomp)
  }
  saved <- readRDS(rds_path)
  m <- saved$norm
  rownames(m) <- saved$id
  miss_s <- setdiff(colnames(recomp$norm), colnames(m))
  if (length(miss_s)) stop(level, ": samples missing from ", rds_path, ": ", paste(miss_s, collapse = ", "))
  common <- intersect(rownames(recomp$norm), rownames(m))
  a <- recomp$norm[common, , drop = FALSE]
  b <- m[common, colnames(a), drop = FALSE]
  rel <- max(abs(a - b) / pmax(1, abs(b)))
  log_step("  %s: using %s (%d / %d tested features); recomputed vs saved: max relative diff %.1e",
           level, basename(rds_path), length(common), nrow(recomp$norm), rel)
  if (rel > 1e-6) warning(level, ": recomputed normalisation differs from your logcpm.rds - the saved values are used.")
  if (length(common) < nrow(recomp$norm)) {
    warning(level, ": ", nrow(recomp$norm) - length(common),
            " DE features not in logcpm.rds; their recomputed values are used.")
  }
  out <- recomp$norm
  out[common, ] <- b
  recomp$norm   <- out
  recomp$source <- rds_path
  recomp
}

# Match samplesheet sample names to Bambu column names / bigWig basenames
match_samples <- function(cpm_cols, samples, what = "count columns") {
  hit <- vapply(samples, function(s) {
    if (s %in% cpm_cols) return(s)
    m <- cpm_cols[grepl(paste0("^", regex_escape(s), "([^A-Za-z0-9]|$)"), cpm_cols, perl = TRUE)]
    if (length(m) == 1) m else NA_character_
  }, character(1))
  if (anyNA(hit)) {
    stop("Could not match samplesheet sample(s) to ", what, ": ",
         paste(samples[is.na(hit)], collapse = ", "),
         "\n  ", what, " are: ", paste(cpm_cols, collapse = ", "))
  }
  if (anyDuplicated(hit)) stop("Ambiguous sample matching: ", paste(hit, collapse = ", "))
  hit  # named: sample -> count column
}

# First file matching `pattern` under `dir` (Nextflow work/ skipped), preferring
# paths that match `prefer`
find_one <- function(dir, pattern, prefer = NULL, what) {
  hits <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("[/\\\\]work[/\\\\]", hits)]
  if (!length(hits)) stop("No ", what, " found under ", dir, " - set its path in cfg.")
  if (!is.null(prefer)) {
    pref <- hits[grepl(prefer, hits)]
    if (length(pref)) hits <- pref
  }
  hits[1]
}

# isoformic tx_to_gene: transcript_id, gene_id (MUST be column 2),
# transcript_name, gene_name, transcript_type
build_tx_to_gene <- function(det, counts_tbl, deg) {
  name_by_gid <- setNames(deg$gene_name, deg$gene_id)
  gid_by_name <- setNames(deg$gene_id, deg$gene_name)
  det |>
    select(key, transcript_id, transcript_name, transcript_type, novel) |>
    left_join(distinct(counts_tbl, key, gene_id), by = "key") |>
    mutate(
      # GENCODE names: "EGFR-201" -> "EGFR" (annotated transcripts only)
      stripped  = if_else(!novel, sub("-[0-9]+$", "", transcript_name), NA_character_),
      gene_id   = coalesce(gene_id, unname(gid_by_name[stripped])),
      gene_name = coalesce(unname(name_by_gid[gene_id]), stripped)
    ) |>
    group_by(gene_id) |>
    # novel isoforms of a gene not tested at gene level inherit a sibling's name
    mutate(gene_name = if_else(
      is.na(gene_name) & !is.na(gene_id), first(gene_name[!is.na(gene_name)]), gene_name
    )) |>
    ungroup() |>
    mutate(gene_name = coalesce(gene_name, gene_id)) |>
    select(key, transcript_id, gene_id, transcript_name, gene_name, transcript_type, novel)
}

# isoformic palette + colours for pulposeq novel biotypes (paired hues:
# novel_lncRNA ~ lncRNA pink, novel_protein_coding ~ protein_coding green)
make_palette <- function(types) {
  pal <- isoformic::tx_type_palette()
  # CDS not defined: isoformic's #fc8d62 is too close to the gene line (#fb8072)
  pal["protein_coding_CDS_not_defined"] <- "#a6761d"            # ochre-brown
  extra <- c(
    novel_lncRNA         = "#ae017e",
    novel_protein_coding = "#238b45",
    novel_non_coding     = "#6a51a3",
    TEC                  = "#bdbdbd"
  )
  pal <- c(pal, extra[setdiff(names(extra), names(pal))])
  missing <- setdiff(types, names(pal))
  c(pal, setNames(rep("#d9d9d9", length(missing)), missing))
}

# Same design as isoformic::plot_log2FC(), with the alpha bug fixed
# ---- figure style ------------------------------------------------------------
is_poster <- function() identical(cfg$fig_style, "poster")
base_pt   <- function() cfg$base_size %||% (if (is_poster()) 14 else 11)
lab_mm    <- function() cfg$lab_mm %||% (if (is_poster()) 3.6 else 2.4)   # repel label text (mm)
legend_right <- function() !identical(cfg$legend_position, "bottom")
# legend keys: no underscores; long names wrapped onto two lines in poster mode
# so a right-hand legend does not squeeze the panel
pretty_type <- function(x) {
  x <- gsub("_", " ", x)
  if (is_poster()) stringr::str_wrap(x, width = 16) else x
}
# condition display names (cfg$ctx_group_labels, e.g. TKI_sensitive -> HCC827);
# the value is captured when the plot is built, not when it is drawn
group_labeller <- function() {
  gl <- cfg$ctx_group_labels
  function(x) if (is.null(gl)) x else dplyr::coalesce(unname(gl[x]), x)
}
# figure size (in): poster = final print size, report = the earlier pages
fig_size <- function(kind, n = NULL) {
  P <- is_poster()
  switch(kind,
    # poster widths leave ~2 in for a right-hand legend
    log2fc   = if (P) c(8, 5) else c(10, 6),
    profile  = if (P) c(8, 5.2) else c(10, 6),
    summary  = if (P) c(10, 5.8) else c(11, 6),
    enrich   = if (P) c(8.5, 6.8) else c(11, 8),
    overview = if (P) c(9, 6.2) else c(13, 9),
    markers  = c(if (P) 8.5 else 10, max(4, (if (P) 0.36 else 0.34) * (n %||% 10) + 1.6)),
    cis      = if (P) c(12, 6) else c(12, 5.5),
    c(8, 6)
  )
}
# DE legend with BOTH keys always visible: limits keep "No"/"Yes" even when one
# is absent from the data, and the keys get an explicit grey fill (otherwise an
# absent level is drawn with the bar's unset fill and looks empty)
de_key <- function(kind = c("bar", "point")) {
  kind <- match.arg(kind)
  guide_legend(override.aes = if (kind == "bar") {
    list(fill = "grey35", colour = "grey25")
  } else {
    list(shape = 21, fill = "grey35", colour = "grey25", size = 3)
  })
}

# multi-column legend with the title above the keys, so bottom legends fit the
# poster width (long biotype names otherwise run off the page)
leg2 <- function(ncol = if (legend_right()) 1 else 2) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    guide_legend(ncol = ncol, byrow = TRUE)
  } else {
    guide_legend(ncol = ncol, byrow = TRUE, title.position = "top")
  }
}

# poster adjustments applied to every ggplot / patchwork at export: one base
# font size for all text, legends below the panel, softer subtitles
poster_theme <- function() {
  theme(
    text            = element_text(size = base_pt(), family = theme_font),
    plot.title      = element_text(size = rel(1.15)),
    plot.subtitle   = element_text(size = rel(0.78), colour = "grey30"),
    plot.caption    = element_text(size = rel(0.62), hjust = 0),
    legend.position = if (legend_right()) "right" else "bottom",
    legend.box      = "vertical",
    legend.title    = element_text(size = rel(0.85)),
    legend.text     = element_text(size = rel(0.8)),
    legend.margin   = margin(0, 0, 0, 0)
  ) + if (utils::packageVersion("ggplot2") >= "3.5.0") theme(legend.title.position = "top") else theme()
}

plot_gene_log2fc <- function(de_long, gene, pal, contrast_lab, alpha_lab, keep_ids = character(0)) {
  d <- filter(de_long, gene_name == .env$gene)
  if (!nrow(d)) stop("gene not found in DE table")
  # poster: at most cfg$poster_max_iso isoforms - every DE isoform first, then
  # the largest |log2FC|; the subtitle says how many are shown
  n_tx  <- sum(d$feature_type != "gene")
  sub_n <- NULL
  if (is_poster() && n_tx > cfg$poster_max_iso) {
    keep_tx <- d |>
      filter(feature_type != "gene") |>
      arrange(desc(feature_id %in% keep_ids), desc(is_de == "yes"), desc(abs(log2FC))) |>
      slice_head(n = cfg$poster_max_iso) |>
      pull(feature_id)
    d <- d |> filter(feature_type == "gene" | feature_id %in% keep_tx)
    sub_n <- sprintf("%d/%d isoforms\n(rank DE and |log2FC|)",
                     cfg$poster_max_iso, n_tx)
  }
  tx_order <- d |>
    filter(feature_type != "gene") |>
    arrange(feature_type, desc(log2FC)) |>
    pull(feature_name)
  lv <- unique(c(if (any(d$feature_type == "gene")) gene, tx_order))
  d <- d |>
    mutate(
      feature_name = factor(feature_name, levels = lv),
      is_de = factor(is_de, levels = c("no", "yes"))
    )
  p <- ggplot(d, aes(x = feature_name, y = log2FC, fill = feature_type, alpha = is_de)) +
    # show.legend = TRUE: ggplot2 >= 3.5 otherwise drops the "No" key when
    # every bar is DE
    geom_col(colour = "grey25", linewidth = 0.2, show.legend = TRUE) +
    geom_hline(yintercept = 0, colour = "grey50") +
    scale_fill_manual(values = pal, labels = pretty_type, guide = leg2()) +
    scale_alpha_manual(values = c(no = 0.25, yes = 1), limits = c("no", "yes"),
                       labels = c(no = "No", yes = "Yes"), drop = FALSE,
                       guide = de_key("bar")) +
    labs(title = gene, subtitle = sub_n, x = NULL, y = contrast_lab,
         fill = "Feature type", alpha = alpha_lab) +
    theme_bw(base_size = base_pt(), base_family = theme_font)
  if (nrow(d) > 15) {
    p + coord_flip() + scale_x_discrete(limits = rev)
  } else {
    p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}

# keep        profile labels (genename) always drawn - e.g. a partner lncRNA or a
#             novel model from another Bambu gene added to this gene's profile
# de_only     drop this gene's non-DE isoforms (the gene row always stays)
# type_limits fixed transcript-type legend (identical legends merge in panels)
# partner     a second gene drawn in its own panel BELOW the gene (e.g. the
#             antisense lncRNA UBE2I-AS1 under UBE2I); each gets its gene line
#             + at most cs$pair_max_iso isoforms
#
# Isoform cap (poster / journal panels), switch-aware: DE isoforms first,
# ranked by |log2FC|, but the strongest DE isoform of EACH direction is always
# kept (two of each when 8 or more are shown), so a switch is never hidden by
# a cap filled with isoforms going one way (e.g. UBE2I-204, the only UBE2I
# isoform up, is below six down isoforms in |log2FC|).
pick_isoforms <- function(cand, n) {
  if (length(cand) <= n) return(cand)
  di  <- match(cand, det_final$transcript_name)
  de  <- coalesce(det_final$is_de[di], FALSE)
  lfc <- coalesce(det_final$log2FC[di], 0)
  ord <- order(-de, -abs(lfc))
  n_dir <- if (n >= 8) 2L else 1L
  must <- unlist(lapply(c(1, -1), function(s) head(cand[ord][(de & sign(lfc) == s)[ord]], n_dir)))
  unique(c(must, cand[ord]))[seq_len(n)]
}
plot_gene_profile <- function(profile_df, gene, pal, keep = character(0), de_only = FALSE,
                              title = gene, subtitle = NULL, type_limits = NULL,
                              partner = character(0)) {
  if (!gene %in% profile_df$parent_gene) stop("gene not in profile data")
  partner <- setdiff(intersect(partner, profile_df$parent_gene), gene)
  genes   <- c(gene, partner)
  if (de_only) {
    profile_df <- profile_df |>
      filter(!parent_gene %in% genes | transcript_type == "gene" | DE %in% "Yes" |
               genename %in% keep)
  }
  # poster: isoform cap per gene (see pick_isoforms)
  sub_n <- NULL
  if (is_poster()) {
    notes <- character(0)
    for (g in genes) {
      cand  <- unique(profile_df$genename[profile_df$parent_gene == g &
                                            profile_df$transcript_type != "gene" &
                                            !profile_df$genename %in% keep])
      k_g   <- if (g == gene) keep else character(0)
      n_cap <- if (length(partner)) cs$pair_max_iso %||% 5L else
        max(cs$min_own_iso %||% 0L, cfg$poster_max_iso - length(keep))
      if (length(cand) > n_cap) {
        keep_tx <- pick_isoforms(cand, n_cap)
        profile_df <- profile_df |>
          filter(parent_gene != g | transcript_type == "gene" | genename %in% c(keep_tx, k_g))
        notes <- c(notes, sprintf("%s%d/%d isoforms", if (length(partner)) paste0(g, ": ") else "",
                                  length(keep_tx) + length(k_g), length(cand) + length(k_g)))
      }
    }
    if (length(notes)) sub_n <- paste(notes, collapse = "; ")
  }
  if (length(partner)) {
    profile_df <- profile_df |> filter(parent_gene %in% genes) |>
      mutate(parent_gene = factor(parent_gene, levels = genes))
  }
  p <- isoformic::plot_tx_expr(genes_to_plot = genes, profile_data = profile_df)
  # gene on top, partner below (isoformic would put them side by side)
  if (length(partner)) p <- p + facet_wrap(~ parent_gene, ncol = 1, scales = "free_y")

  # isoformic labels every DE isoform at the first condition with tight repel
  # settings; replace that layer with a spread-out one placed in extra space
  # left of the first condition (falls back to isoformic's labels if the plot
  # object cannot be edited). The gene row is always labelled, DE or not.
  var   <- colnames(profile_df)[2]
  lv    <- levels(profile_df[[var]])
  # label side: "left" (first condition), "right" (last) or "auto" = the side
  # where fewer labelled values crowd together (e.g. SNHG17: 11 isoforms near
  # 0 in HCC827, spread out in H1975 -> labels on the right)
  side  <- cfg$profile_label_side %||% "auto"
  if (identical(side, "auto")) {
    lr  <- profile_df[profile_df$parent_gene %in% genes &
                        (profile_df$DE %in% "Yes" | profile_df$transcript_type %in% "gene" |
                           profile_df$genename %in% keep), , drop = FALSE]
    thr <- max(diff(range(lr$log2_mean_TPM, na.rm = TRUE)), 1) * 0.05
    crowd <- vapply(c(lv[1], lv[length(lv)]), function(l) {
      v <- sort(lr$log2_mean_TPM[lr[[var]] %in% l])
      sum(diff(v) < thr)
    }, numeric(1))
    side <- if (crowd[2] < crowd[1]) "right" else "left"
  }
  right <- identical(side, "right")
  lvl1  <- if (right) lv[length(lv)] else lv[1]        # condition the labels sit at
  room  <- cfg$profile_label_room
  p_new <- tryCatch({
    keep_l <- !vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1))
    p$layers <- p$layers[keep_l]
    p + ggrepel::geom_text_repel(
      data = function(d) d[(d$DE %in% "Yes" | d$transcript_type %in% "gene") & d[[var]] %in% lvl1, , drop = FALSE],
      mapping = aes(x = .data[[var]], y = .data$log2_mean_TPM, label = .data$genename),
      size = lab_mm() * 0.9, direction = "y",
      hjust = if (right) 0 else 1, nudge_x = if (right) 0.15 else -0.15,
      xlim = if (right) c(length(lv) + 0.08, Inf) else c(-Inf, 0.92),
      # journal panels: labels kept close to their line (short leaders)
      box.padding = if (isTRUE(cfg$profile_compact)) 0.08 else 0.25,
      point.padding = if (isTRUE(cfg$profile_compact)) 0.05 else 0.2,
      force = if (isTRUE(cfg$profile_compact)) 0.5 else 2,
      min.segment.length = 0, segment.colour = "grey60",
      segment.size = if (isTRUE(cfg$profile_compact)) 0.2 else 0.3,
      max.overlaps = Inf, seed = 1, inherit.aes = FALSE, lineheight = 0.85
    ) +
      # padding on the label side sized for the labels only (category widths)
      scale_x_discrete(expand = expansion(add = if (right) c(0.25, room) else c(room, 0.25)),
                       labels = group_labeller())
  }, error = function(e) NULL)
  if (!is.null(p_new)) p <- p_new

  # both DE keys always shown: ggplot2 >= 3.5 drops keys for absent values
  # unless the layer has show.legend = TRUE
  for (i in seq_along(p$layers)) {
    if ("alpha" %in% names(p$layers[[i]]$mapping)) p$layers[[i]]$show.legend <- TRUE
  }

  sub_txt <- paste(c(subtitle, if (!isFALSE(cfg$profile_cap_note)) sub_n), collapse = " \u00b7 ")
  if (!nzchar(sub_txt)) sub_txt <- NULL

  # suppressMessages: replacing isoformic's colour / alpha scales is intended
  suppressMessages(
    p +
      scale_colour_manual(values = pal, aesthetics = c("colour", "fill"), labels = pretty_type,
                          guide = leg2(), limits = type_limits, drop = is.null(type_limits)) +
      scale_alpha_manual(values = c(Yes = 1, No = 0.2), limits = c("Yes", "No"), drop = FALSE) +
      labs(x = NULL, y = expression(log[2](mean ~ normalized ~ count + 1)), title = title,
           subtitle = sub_txt) +
      # the facet strip only repeats the title
      theme(strip.background = element_blank(), strip.text = element_blank())
  ) |> profile_style()
}

# one profile style everywhere (poster pages, event panels, cis pairs, journal
# panels): isoformic draws thick axes and large points; scale them with the
# base font (1 = poster 14 pt). cfg$profile_point / profile_line override.
# Same look as the other figure panels (theme_bw-like, as fig2_v2): boxed
# panel, no grid, the figure font; border 0.4 at 7 pt, scaled with the base font.
# isoformic (R/profile_plot.R) draws the per-condition SD three times: a
# triangle at mean + SD (shape 24), one at mean - SD (shape 25) and an error
# bar between them. Default (cfg$profile_sd = "triangles") keeps the triangles at
# half the mean point's size, without legend keys; "bar" keeps only the bar.
profile_style <- function(p) {
  k  <- base_pt() / 14
  pt <- cfg$profile_point %||% (2.2 * k)
  ln <- cfg$profile_line  %||% (0.8 * k)
  bw <- 0.4 * base_pt() / 7
  if (!identical(cfg$profile_sd, "triangles")) {
    tri <- vapply(p$layers, function(l) {
      inherits(l$geom, "GeomPoint") && any(as.character(l$aes_params$shape) %in% c("24", "25"))
    }, logical(1))
    p$layers <- p$layers[!tri]
  }
  for (i in seq_along(p$layers)) {
    g <- p$layers[[i]]$geom
    if (inherits(g, c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel"))) {
      p$layers[[i]]$aes_params$family <- theme_font
      next
    }
    if (inherits(g, "GeomPoint")) {
      tri_i <- any(as.character(p$layers[[i]]$aes_params$shape) %in% c("24", "25"))
      p$layers[[i]]$aes_params$size <- if (tri_i) 0.5 * pt else pt            # SD ends / the mean
      if (tri_i) p$layers[[i]]$show.legend <- FALSE                            # legend: point + line
    }
    if (inherits(g, c("GeomPath", "GeomSegment"))) p$layers[[i]]$aes_params$linewidth <- ln
    if (inherits(g, c("GeomErrorbar", "GeomLinerange"))) {                      # the SD
      p$layers[[i]]$aes_params$linewidth <- 0.6 * ln
      p$layers[[i]]$aes_params$alpha     <- 0.45
      p$layers[[i]]$show.legend          <- FALSE
    }
  }
  p + theme(text = element_text(family = theme_font, face = "plain"),
            plot.title = element_text(face = "plain"), legend.title = element_text(face = "plain"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = bw),
            panel.background = element_rect(fill = "white", colour = NA),
            panel.grid = element_blank(),
            axis.line  = element_blank(),
            axis.ticks = element_line(linewidth = 0.75 * bw, colour = "black"),
            axis.ticks.length = unit(1.2 * k, "mm"))
}

safe_plot <- function(expr, label) {
  tryCatch(expr, error = function(e) {
    message("    ! skipped ", label, ": ", conditionMessage(e))
    NULL
  })
}

# file-name-safe version of a label
slug_file <- function(x) substr(gsub("(^-|-$)", "", gsub("[^A-Za-z0-9]+", "-", x)), 1, 80)

# Every figure is written three ways:
#   <name>.pdf          all pages in one file
#   <name>_PNG/NN_<page>.png   one 300-dpi PNG per page
#   <name>_SVG/NN_<page>.svg   one SVG per page (svglite)
# page_names label the per-page files (gene names etc.); defaults to page01...
save_pages <- function(plots, file, width = 10, height = 6, page_names = NULL, kind = NULL, n = NULL) {
  if (is.null(page_names)) page_names <- sprintf("page%02d", seq_along(plots))
  keep <- !vapply(plots, is.null, logical(1))
  plots <- plots[keep]
  page_names <- page_names[keep]
  if (!length(plots)) return(invisible(NULL))
  if (!is.null(kind)) {
    sz <- fig_size(kind, n)
    width  <- sz[1]
    height <- sz[2]
  }
  if (is_poster()) {
    plots <- lapply(plots, function(p) {
      if (inherits(p, "patchwork")) p & poster_theme() else p + poster_theme()
    })
  }

  # cairo_pdf embeds unicode (rho, arrows, middle dots) that pdf() cannot
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(file, width = width, height = height, onefile = TRUE)
  } else {
    grDevices::pdf(file, width = width, height = height, onefile = TRUE)
  }
  for (p in plots) {
    tryCatch(suppressMessages(print(p)),
             error = function(e) message("    ! page failed: ", conditionMessage(e)))
  }
  grDevices::dev.off()

  base    <- tools::file_path_sans_ext(basename(file))
  png_dir <- file.path(dirname(file), paste0(base, "_PNG"))
  svg_dir <- file.path(dirname(file), paste0(base, "_SVG"))
  dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(svg_dir, showWarnings = FALSE, recursive = TRUE)
  stems <- sprintf("%02d_%s", seq_along(plots), slug_file(page_names))
  for (i in seq_along(plots)) {
    tryCatch(suppressMessages({
      ggsave(file.path(png_dir, paste0(stems[i], ".png")), plots[[i]],
             width = width, height = height, dpi = cfg$dpi, bg = "white")
      ggsave_svg(file.path(svg_dir, paste0(stems[i], ".svg")), plots[[i]],
                 width = width, height = height, bg = "white")
    }), error = function(e) message("    ! PNG/SVG export failed for ", stems[i], ": ", conditionMessage(e)))
  }
  log_step("  wrote %s (%d pages) + %s_PNG/ + %s_SVG/", basename(file), length(plots), base, base)
  invisible(file)
}

# Import the reference GTF once and cache it (keyed by file size + mtime), so
# reruns skip the multi-minute parse. Only the attributes used here are kept.
load_reference <- function(path) {
  cache_dir <- cfg$ref_cache_dir %||% tools::R_user_dir("isoformic_dong", which = "cache")
  info  <- file.info(path)
  key   <- sprintf("%s_%s_%s.rds", tools::file_path_sans_ext(basename(path)),
                   format(info$size, scientific = FALSE), format(as.numeric(info$mtime), scientific = FALSE))
  cache <- file.path(cache_dir, key)
  if (file.exists(cache)) {
    log_step("  reference: cached import %s", cache)
    return(readRDS(cache))
  }
  log_step("  reading reference annotation (first run only: a few minutes, several GB RAM): %s", path)
  gr <- rtracklayer::import(path, feature.type = c("transcript", "exon", "CDS"))
  keep_cols <- intersect(c("type", "gene_id", "gene_name", "gene_type", "gene_biotype",
                           "transcript_id", "transcript_name", "transcript_type",
                           "transcript_biotype", "phase"),
                         colnames(S4Vectors::mcols(gr)))
  S4Vectors::mcols(gr) <- S4Vectors::mcols(gr)[, keep_cols, drop = FALSE]
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    saveRDS(gr, cache)
    log_step("  reference cached to %s", cache)
  }, error = function(e) log_step("  (reference cache not written: %s)", conditionMessage(e)))
  gr
}

## ---- 3. Output folders ------------------------------------------------------
out_root <- resolve(cfg$out_dir)
out_t <- function(x) file.path(out_root, "tables", x)
out_p <- function(x) file.path(out_root, "plots", x)
dir.create(file.path(out_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "plots"), recursive = TRUE, showWarnings = FALSE)
log_step("isoformic %s | output: %s", utils::packageVersion("isoformic"), out_root)

# your condition colours (group_cols) from helpers/plot_theme.R, read-only.
# Sourced into its own environment; any theme it sets also applies here.
theme_group_cols <- NULL
if (!is.null(cfg$plot_theme) && file.exists(resolve(cfg$plot_theme))) {
  theme_env <- new.env()
  tryCatch({
    sys.source(resolve(cfg$plot_theme), envir = theme_env)
    if (exists("group_cols", envir = theme_env, inherits = FALSE)) {
      theme_group_cols <- get("group_cols", envir = theme_env)
      log_step("  group colours from %s", cfg$plot_theme)
    }
  }, error = function(e) stop("plot_theme.R not loaded: ", conditionMessage(e), call. = FALSE))
}
# the figure font of plot_theme.R (Source Sans 3) for every plot here; no
# fallback: plot_theme.R must load (it registers the font for showtext)
if (!exists("theme_env") || !exists("font_family", envir = theme_env, inherits = FALSE))
  stop("helpers/plot_theme.R not found (cfg$plot_theme): it sets the figure font Source Sans 3")
theme_font <- get("font_family", envir = theme_env)
## ---- editable SVG text ------------------------------------------------------
## plot_theme.R switches showtext on, and showtext draws every letter as a glyph
## outline on EVERY device, SVG included. These helpers switch showtext off while
## an SVG is written, so svglite writes real <text> elements in the figure font:
## text boxes you can edit in Inkscape / Illustrator / Affinity / Vectornator.
##   svg_font        font name written into the SVG. It must be installed on the
##                   computer that edits the file (Google Fonts "Source Sans 3",
##                   static TTFs, "Install for all users"), or the editor swaps it.
##   fix_text_size   FALSE: no textLength attribute, so edited text is not
##                   stretched or squeezed to the old width.
## PNG/PDF are unchanged (showtext stays on for them).
if (!exists("svg_editable_helpers", inherits = TRUE)) {
  svg_editable_helpers <- TRUE
  svg_font <- get0("svg_font", ifnotfound = "Source Sans 3")
  ## R family names used in the plots -> the real font name
  svg_aliases <- function() {
    fams <- unique(c("sans", get0("font_family", ifnotfound = NULL), get0("theme_font", ifnotfound = NULL)))
    fams <- setdiff(fams, c("serif", "mono", "symbol", svg_font))
    stats::setNames(rep(list(svg_font), length(fams)), fams)
  }
  ## is the font known to systemfonts (the text measurer svglite uses)? checked once
  svg_font_ok <- local({
    ok <- NULL
    function() {
      if (!is.null(ok)) return(ok)
      have <- tryCatch(svg_font %in% systemfonts::system_fonts()$family, error = function(e) FALSE)
      if (!have && "require_font" %in% getNamespaceExports("systemfonts")) {
        ## systemfonts >= 1.2: fetch from Google Fonts for this R session (metrics only)
        have <- isTRUE(tryCatch({ systemfonts::require_font(svg_font); TRUE },
                                error = function(e) FALSE))
      }
      if (!have) message("  ! '", svg_font, "' is not installed: SVG text is still editable, ",
                         "but measured with a fallback font. Install it and restart R.")
      ok <<- have
      ok
    }
  })
  svg_device <- function(filename, width, height, bg = "white", ...) {
    svg_font_ok()
    svglite::svglite(filename, width = width, height = height, bg = bg,
                     system_fonts = svg_aliases(), fix_text_size = FALSE)
  }
  showtext_pause  <- function() {
    on <- "showtext" %in% loadedNamespaces()
    if (on) showtext::showtext_auto(FALSE)
    on
  }
  showtext_resume <- function(was_on) if (isTRUE(was_on)) showtext::showtext_auto(TRUE)
  ## ggsave() for SVG with editable text
  ggsave_svg <- function(filename, plot, width, height, ...) {
    st <- showtext_pause(); on.exit(showtext_resume(st))
    ggplot2::ggsave(filename, plot, width = width, height = height, device = svg_device, ...)
  }
  ## grid / base graphics: st <- svg_open(file, w, h); draw ...; svg_close(st)
  svg_open  <- function(filename, width, height) {
    st <- showtext_pause(); svg_device(filename, width, height); invisible(st)
  }
  svg_close <- function(st) { invisible(grDevices::dev.off()); showtext_resume(st) }
}
## ---------------------------------------------------------------------------

## ---- 4. Read + standardise DE tables ---------------------------------------
log_step("Reading DE tables")
deg <- standardise_de(read_csv(resolve(cfg$deg_file), show_col_types = FALSE), "gene", cfg$lfc_col)
det <- standardise_de(read_csv(resolve(cfg$det_file), show_col_types = FALSE), "transcript", cfg$lfc_col)
log_step("  genes: %d tested, %d DE | transcripts: %d tested, %d DE (%d novel DE)",
         nrow(deg), sum(deg$is_de), nrow(det), sum(det$is_de), sum(det$is_de & det$novel))

agree <- function(d) 100 * mean(d$is_de == (d$qvalue < cfg$alpha & abs(d$log2FC) > cfg$lfc_cut))
log_step("  edgeR `call` vs (FDR < %.2g & |log2FC| > %.2g): genes %.2f%%, transcripts %.2f%% agree",
         cfg$alpha, cfg$lfc_cut, agree(deg), agree(det))

## ---- 5. Samplesheet ---------------------------------------------------------
log_step("Reading samplesheet")
sheet <- read_csv(resolve(cfg$samplesheet), show_col_types = FALSE)
if (!all(c("sample", cfg$group_col) %in% names(sheet))) {
  stop("samplesheet needs columns 'sample' and '", cfg$group_col, "'")
}
# TKI-resistant -> TKI_resistant, exactly as make_coldata() in your edgeR script
sheet[[cfg$group_col]] <- gsub("-", "_", sheet[[cfg$group_col]])
groups_found <- unique(sheet[[cfg$group_col]])
if (!all(c(cfg$control, cfg$case) %in% groups_found)) {
  stop("cfg$control/cfg$case not found in samplesheet. Groups present: ",
       paste(groups_found, collapse = ", "))
}
sheet <- sheet |>
  filter(.data[[cfg$group_col]] %in% c(cfg$control, cfg$case)) |>
  transmute(sample = as.character(sample), group = .data[[cfg$group_col]]) |>
  as.data.frame()
log_step("  %s", paste(sprintf("%s=%d", names(table(sheet$group)), table(sheet$group)), collapse = ", "))

## ---- 6. Raw Bambu counts (the edgeR input) + transcript-to-gene map --------
tx_counts_path <- if (is.null(cfg$tx_counts_file)) {
  log_step("Searching for Bambu raw transcript counts under %s", cfg$dong_dir)
  find_bambu_file(cfg$dong_dir, "^BambuOutput_counts_transcript.*\\.txt$",
                  det$transcript_id, "BambuOutput_counts_transcript*.txt")
} else {
  resolve(cfg$tx_counts_file)
}
log_step("Transcript counts: %s", tx_counts_path)

cnt_tx   <- read_tsv(tx_counts_path, show_col_types = FALSE, progress = FALSE)
tx_col   <- intersect(c("TXNAME", "transcript_id"), names(cnt_tx))[1]
gene_col <- intersect(c("GENEID", "gene_id"), names(cnt_tx))[1]
if (is.na(tx_col) || is.na(gene_col)) {
  stop("Counts file needs TXNAME + GENEID columns; found: ", paste(head(names(cnt_tx), 5), collapse = ", "))
}
samples <- sheet$sample
smap <- match_samples(setdiff(names(cnt_tx), c(tx_col, gene_col)), samples, "transcript count columns")
cnt_tx <- cnt_tx |> select(tx = all_of(tx_col), gid = all_of(gene_col), all_of(smap))

# version-robust ID keys (only strip versions if that improves matching)
tx_key <- if (mean(det$transcript_id %in% cnt_tx$tx) <
              mean(strip_version(det$transcript_id) %in% strip_version(cnt_tx$tx))) strip_version else identity
g_key  <- if (mean(cnt_tx$gid %in% deg$gene_id) <
              mean(strip_version(cnt_tx$gid) %in% strip_version(deg$gene_id))) strip_version else identity
cnt_tx$key <- tx_key(cnt_tx$tx)
det$key    <- tx_key(det$transcript_id)
cnt_tx <- distinct(cnt_tx, key, .keep_all = TRUE)
gid_to_de <- setNames(deg$gene_id, g_key(deg$gene_id))
cnt_tx$gene_id <- coalesce(unname(gid_to_de[g_key(cnt_tx$gid)]), cnt_tx$gid)

cov_tx <- mean(det$key %in% cnt_tx$key)
log_step("  DE transcripts found in counts file: %.2f%%", 100 * cov_tx)
if (cov_tx < 0.99) {
  warning("Not every DE transcript is in the counts file - is this the matrix you gave edgeR? ",
          "Set cfg$tx_counts_file explicitly.")
}

tx2g <- build_tx_to_gene(det, cnt_tx, deg)
n_unmapped <- sum(is.na(tx2g$gene_id))
if (n_unmapped) {
  warning(n_unmapped, " DE transcripts have no parent gene; they are excluded from gene-linked analyses.")
}
tx2g_iso <- tx2g |> filter(!is.na(gene_id)) |> select(-key)  # isoformic-format tx_to_gene
write_csv(tx2g_iso, out_t("01_tx_to_gene.csv"))

# ---- 6b. Gene counts --------------------------------------------------------
gene_counts_path <- if (!is.null(cfg$gene_counts_file)) {
  resolve(cfg$gene_counts_file)
} else {
  gc <- list.files(dirname(tx_counts_path), "^BambuOutput_counts_gene.*\\.txt$", full.names = TRUE)
  stage <- sub("^BambuOutput_counts_transcript", "", basename(tx_counts_path))
  same  <- gc[sub("^BambuOutput_counts_gene", "", basename(gc)) == stage]
  if (length(same)) same[1] else if (length(gc)) gc[1] else NA_character_
}
if (!is.na(gene_counts_path)) {
  log_step("Gene counts: %s", gene_counts_path)
  cnt_gene <- read_tsv(gene_counts_path, show_col_types = FALSE, progress = FALSE)
  gid_col  <- intersect(c("GENEID", "gene_id"), names(cnt_gene))[1]
  if (is.na(gid_col)) gid_col <- names(cnt_gene)[1]
  gmap <- match_samples(setdiff(names(cnt_gene), c(gid_col, tx_col)), samples, "gene count columns")
  cnt_gene <- cnt_gene |> select(gid = all_of(gid_col), all_of(gmap))
} else {
  log_step("Gene counts: no BambuOutput_counts_gene*.txt next to the transcript file - summing transcript counts per gene")
  cnt_gene <- cnt_tx |> group_by(gid) |> summarise(across(all_of(samples), sum), .groups = "drop")
}
cnt_gene$gene_id <- coalesce(unname(gid_to_de[g_key(cnt_gene$gid)]), cnt_gene$gid)
cov_gene <- mean(deg$gene_id %in% cnt_gene$gene_id)
log_step("  DE genes found in gene counts: %.2f%%", 100 * cov_gene)
if (cov_gene < 0.99) {
  warning("Not every DE gene is in the gene counts - set cfg$gene_counts_file to the matrix you gave edgeR.")
}

# ---- 6c. Normalised counts: your saved edgeR values (or an exact recompute) --
tested_tx <- cnt_tx |> inner_join(select(det, key, transcript_id), by = "key")
m_tx <- as.matrix(tested_tx[, samples])
rownames(m_tx) <- tested_tx$transcript_id

tested_gene <- cnt_gene |> filter(gene_id %in% deg$gene_id) |> distinct(gene_id, .keep_all = TRUE)
m_gene <- as.matrix(tested_gene[, samples])
rownames(m_gene) <- tested_gene$gene_id

log_step("Normalised counts")
norm_tx   <- use_saved_norm(edger_norm(m_tx),
                            file.path(resolve(cfg$edger_dir), "transcript", "logcpm.rds"), "transcript")
norm_gene <- use_saved_norm(edger_norm(m_gene),
                            file.path(resolve(cfg$edger_dir), "gene", "logcpm.rds"), "gene")

norm_tab <- bind_rows(transcript = norm_tx$samples, gene = norm_gene$samples, .id = "level")
write_csv(norm_tab, out_t("00_normalisation_factors.csv"))
log_step("  TMM norm.factors - should equal y$samples$norm.factors in your edgeR run:")
log_step("    transcript: %s", paste(sprintf("%s=%.3f", norm_tx$samples$sample,
                                                norm_tx$samples$norm.factors), collapse = ", "))
log_step("    gene      : %s", paste(sprintf("%s=%.3f", norm_gene$samples$sample,
                                                norm_gene$samples$norm.factors), collapse = ", "))

## ---- 7. DET "final" table (isoformic core table) ----------------------------
log_step("Building DET final table")
sig_gene_ids <- deg$gene_id[deg$is_de]
det_final <- det |>
  select(-key) |>
  left_join(select(tx2g_iso, transcript_id, gene_id, gene_name), by = "transcript_id")

# isoformic::is_deg_sig() -> DEG_sig = YES/NO (parent gene is DE)
det_final <- isoformic::is_deg_sig(
  deg_sig_vector = det_final$transcript_name[det_final$gene_id %in% sig_gene_ids],
  det_table = det_final
)

det_final <- det_final |>
  left_join(
    deg |> select(gene_id, gene_log2FC = log2FC, gene_qvalue = qvalue, gene_is_de = is_de),
    by = "gene_id"
  ) |>
  mutate(
    gene_tested = !is.na(gene_is_de),
    det_class = case_when(
      !is_de       ~ "not DET",
      !gene_tested ~ "DET, gene not tested",
      gene_is_de   ~ "DET & DEG",
      TRUE         ~ "DET only (gene not DE)"
    ),
    direction_vs_gene = case_when(
      det_class == "DET & DEG" & sign(log2FC) == sign(gene_log2FC) ~ "concordant",
      det_class == "DET & DEG" ~ "discordant",
      TRUE ~ NA_character_
    )
  )
write_csv(det_final, out_t("02_DET_final.csv"))

# vignette definition: DE transcript whose gene is NOT a significant DEG
dets_not_degs <- det_final |> filter(is_de, DEG_sig == "NO") |> arrange(qvalue)
write_csv(dets_not_degs, out_t("04_DETs_not_DEGs.csv"))

## ---- 8. Gene + transcript long table (join_DEG_DET-compatible) -------------
de_long <- bind_rows(
  deg |> transmute(
    feature_id = gene_id, feature_name = gene_name, gene_id, gene_name,
    feature_type = "gene", log2FC, pvalue, qvalue,
    is_de = if_else(is_de, "yes", "no"), novel
  ),
  det_final |> filter(!is.na(gene_id)) |> transmute(
    feature_id = transcript_id, feature_name = transcript_name, gene_id, gene_name,
    feature_type = transcript_type, log2FC, pvalue, qvalue,
    is_de = if_else(is_de, "yes", "no"), novel
  )
)
write_csv(de_long, out_t("03_DEG_DET_long.csv"))

# transparency: what isoformic::join_DEG_DET() would have kept
iso_join <- tryCatch(
  isoformic::join_DEG_DET(
    deg_table = deg |> select(gene_id, gene_name, log2FC, pvalue),
    det_table = det_final |> select(transcript_id, transcript_name, transcript_type, log2FC, pvalue),
    logfc_cut = cfg$lfc_cut, pval_cut = cfg$alpha
  ),
  error = function(e) NULL
)
if (!is.null(iso_join)) {
  kept <- sum(iso_join$feature_type != "gene")
  log_step("  join_DEG_DET() would keep %d / %d transcripts (drops %d, incl. all %d novel)",
           kept, nrow(det_final), nrow(det_final) - kept, sum(det_final$novel))
}

## ---- 10. Gene sets to plot --------------------------------------------------
switch_genes <- det_final |>
  filter(is_de, !is.na(gene_id)) |>
  group_by(gene_id, gene_name) |>
  summarise(
    n_up = sum(log2FC > 0), n_down = sum(log2FC < 0),
    lfc_max = max(log2FC), lfc_min = min(log2FC),
    types = paste(sort(unique(transcript_type)), collapse = ";"),
    .groups = "drop"
  ) |>
  filter(n_up > 0, n_down > 0) |>
  mutate(switch_span = lfc_max - lfc_min) |>
  left_join(deg |> select(gene_id, gene_log2FC = log2FC, gene_qvalue = qvalue, gene_is_de = is_de),
            by = "gene_id") |>
  arrange(desc(switch_span))
write_csv(switch_genes, out_t("05_isoform_switch_genes.csv"))

det_only_genes <- det_final |>
  filter(det_class == "DET only (gene not DE)") |>
  group_by(gene_id, gene_name) |>
  summarise(n_det = n(), best_q = min(qvalue), .groups = "drop") |>
  arrange(best_q, desc(n_det))

novel_in_annot <- det_final |>
  filter(novel, is_de, !is.na(gene_id), !startsWith(gene_id, "BambuGene")) |>
  group_by(gene_id, gene_name) |>
  summarise(
    n_novel_de = n(),
    novel_types = paste(sort(unique(transcript_type)), collapse = ";"),
    best_q = min(qvalue), .groups = "drop"
  ) |>
  arrange(best_q)
write_csv(novel_in_annot, out_t("06_novel_DE_in_annotated_genes.csv"))

missing_manual <- setdiff(cfg$genes_manual, de_long$gene_name)
if (length(missing_manual)) log_step("  manual genes not in DE tables: %s", paste(missing_manual, collapse = ", "))

n_auto <- cfg$n_auto_genes
plot_sets <- list(
  manual             = intersect(cfg$genes_manual, de_long$gene_name),
  isoform_switch     = head(switch_genes$gene_name, n_auto),
  det_only           = head(det_only_genes$gene_name, n_auto),
  novel_in_annotated = head(novel_in_annot$gene_name, n_auto)
)
log_step("  plot sets: %s", paste(sprintf("%s=%d", names(plot_sets), lengths(plot_sets)), collapse = ", "))

## ---- 11. Profile data (isoformic::prepare_profile_data) --------------------
log_step("Preparing expression profile data (edgeR TMM-normalised counts)")
cpm_tx <- tibble::as_tibble(norm_tx$norm, rownames = "transcript_id") |>
  inner_join(select(tx2g_iso, transcript_id, transcript_name), by = "transcript_id") |>
  select(transcript_id, transcript_name, all_of(samples))
cpm_gene <- tibble::as_tibble(norm_gene$norm, rownames = "gene_id") |>
  filter(gene_id %in% tx2g_iso$gene_id)

profile_df <- quiet_joins(isoformic::prepare_profile_data(
  txi_gene             = cpm_gene,
  txi_transcript       = cpm_tx,
  sample_metadata      = sheet,
  tx_to_gene           = tx2g_iso,
  de_result_gene       = deg |> select(gene_name, log2FC, qvalue),
  de_result_transcript = det_final |> select(transcript_name, log2FC, qvalue),
  var                  = "group",
  var_levels           = c(cfg$control, cfg$case),
  gene_col             = "gene_name",
  tx_col               = "transcript_name",
  pvalue_cutoff        = cfg$alpha,
  lfc_cutoff           = cfg$lfc_cut,
  use_fdr              = TRUE
))

## ---- 12. Plot settings (log2FC bars + expression profiles) -------------------
pal <- make_palette(unique(de_long$feature_type))
contrast_lab <- sprintf("log2FC (%s / %s)", cfg$case, cfg$control)
alpha_lab <- sprintf("DE (FDR < %.2g)", cfg$alpha)

## ---- 14. Genomic context panels (adapted from pulposeq bin/genomic_context.R) --
# Per gene: title + region, chromosome bar, annotated exon footprint, every isoform
# of the gene (plus novel models gffcompare classified against it, e.g. antisense
# or intronic lncRNAs carrying their own BambuGene id), then one coverage track per
# sample grouped by condition, a coordinate axis and a biotype legend.
#
# Kept from pulposeq: plotTranscripts is called once PER ISOFORM (transcriptFilter)
# at a row packed here, with the label drawn separately and clamped inside the
# track; the TxDb is built from a small sub-GTF whose transcript_id IS the label;
# no OrgDb (gene.id.column = display.column = "GENEID"); the chromosome bar uses
# the bigWig header instead of AnnotationHub, so nothing needs the network.
# For DE: labels carry the isoform's log2FC (* = DE), non-DE isoforms are
# drawn in a lightened colour, and coverage is scaled per million assigned reads
# and coloured by condition.
if (cfg$run_context) {
  log_step("Genomic context panels")
  bare_id <- strip_version
  ctx_dir <- out_p("genomic_context")
  dir.create(ctx_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- 14a. pulposeq inputs -------------------------------------------------
  final_gtf <- if (is.null(cfg$final_gtf)) {
    find_one(cfg$dong_dir, "^annotations_final\\.gtf$", "curated", "annotations_final.gtf")
  } else {
    resolve(cfg$final_gtf)
  }
  bw_dir <- if (is.null(cfg$bigwig_dir)) {
    dirname(find_one(cfg$dong_dir, "\\.bw$", "[/\\\\]coverage[/\\\\]", "bigWig (.bw) files"))
  } else {
    resolve(cfg$bigwig_dir)
  }
  bw_all   <- list.files(bw_dir, pattern = "\\.bw$", full.names = TRUE)
  bw_names <- sub("\\.bw$", "", basename(bw_all))
  bw_map   <- match_samples(bw_names, sheet$sample, "bigWig files")
  sample_order <- sheet$sample[order(match(sheet$group, c(cfg$control, cfg$case)), sheet$sample)]
  bw_files <- setNames(bw_all[match(bw_map, bw_names)], names(bw_map))[sample_order]
  sample_group <- setNames(sheet$group, sheet$sample)[sample_order]
  grp_col <- cfg$group_colors %||% theme_group_cols %||%
    setNames(c("#2c7fb8", "#c0392b"), c(cfg$control, cfg$case))
  log_step("  GTF: %s", final_gtf)
  log_step("  bigWigs (%d): %s", length(bw_files), bw_dir)

  # bigWig health check. rtracklayer reads bigWigs through the UCSC C library,
  # whose only error is "UCSC library operation failed" - so test every file up
  # front and say WHY: unreadable file (e.g. OneDrive online-only placeholder,
  # truncated copy) vs chromosome naming (chr7 vs 7) vs window past the end.
  # Three ways to hand a bigWig to the UCSC library, tried in order on the
  # first file; the first that works is used for every file and every import:
  #   direct   full path "C:/.../x.bw"
  #   relative bare file name from inside its folder - the kent library opens
  #            files through its URL cache layer, and a drive letter "C:" can
  #            be taken for a URL scheme
  #   copy     local copy in R's temp folder (forces OneDrive to materialise
  #            online-only files), then relative
  bw_local  <- file.path(tempdir(), "isoformic_bw")
  bw_access <- list(
    direct   = function(f, FUN) FUN(f),
    relative = function(f, FUN) withr::with_dir(dirname(f), FUN(basename(f))),
    copy     = function(f, FUN) {
      dir.create(bw_local, showWarnings = FALSE, recursive = TRUE)
      lf <- file.path(bw_local, basename(f))
      if (!file.exists(lf) || file.info(lf)$size != file.info(f)$size) {
        if (!file.copy(f, lf, overwrite = TRUE)) stop("could not copy ", f, " to ", bw_local)
      }
      withr::with_dir(bw_local, FUN(basename(f)))
    }
  )
  bw_seqinfo <- function(path) GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(path))
  bw_mode <- NULL
  for (m in names(bw_access)) {
    ok <- tryCatch({
      bw_access[[m]](bw_files[[1]], bw_seqinfo)
      TRUE
    }, error = function(e) {
      log_step("  bigWig access '%s' failed: %s", m, conditionMessage(e))
      FALSE
    })
    if (ok) {
      bw_mode <- m
      log_step("  bigWig access: '%s'%s", m,
               if (m == "copy") sprintf(" (local copies in %s)", bw_local) else "")
      break
    }
  }
  # read a bigWig with the access mode that worked
  bw_read <- function(f, FUN) {
    if (is.null(bw_mode)) stop("no bigWig access mode works on this system")
    bw_access[[bw_mode]](f, FUN)
  }

  bw_info <- lapply(bw_files, function(f) {
    tryCatch({
      si <- bw_read(f, bw_seqinfo)
      list(ok = TRUE, lens = GenomeInfoDb::seqlengths(si), size = file.info(f)$size)
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e), size = file.info(f)$size))
  })
  bw_ok <- vapply(bw_info, `[[`, logical(1), "ok")
  for (s in names(bw_files)[!bw_ok]) {
    log_step("  ! bigWig cannot be opened: %s (%s bytes): %s", bw_files[[s]],
             format(bw_info[[s]]$size, big.mark = ","), bw_info[[s]]$msg)
  }
  if (!any(bw_ok)) {
    log_step("  ! no bigWig readable by any access mode - panels are drawn WITHOUT coverage.")
    log_step("    Run the check in the chat reply (rtracklayer version / file header) and send the output.")
  }
  bw_lens <- if (any(bw_ok)) bw_info[[which(bw_ok)[1]]]$lens else NULL
  if (!is.null(bw_lens)) {
    log_step("  bigWig chromosomes: %d (e.g. %s)", length(bw_lens), paste(head(names(bw_lens), 3), collapse = ", "))
  }
  # GTF chromosome name -> bigWig chromosome name (handles chr7 vs 7, chrM vs MT)
  bw_chrom <- function(chrom) {
    if (is.null(bw_lens) || chrom %in% names(bw_lens)) return(chrom)
    alt <- if (startsWith(chrom, "chr")) sub("^chr", "", chrom) else paste0("chr", chrom)
    if (chrom == "chrM") alt <- "MT"
    if (chrom == "MT") alt <- "chrM"
    if (alt %in% names(bw_lens)) alt else NA_character_
  }
  bw_state <- new.env()
  bw_state$warned <- FALSE

  # Coverage scaling: the bigWigs are raw depth (pulposeq bam_coverage.R). Each
  # sample is multiplied by mean(eff) / eff from the transcript-level TMM, the
  # same factor that turns your raw counts into normalised counts.
  scale_f   <- setNames(norm_tx$samples$rel.scale, norm_tx$samples$sample)[names(bw_files)]
  depth_lab <- "read depth x TMM scaling (as edgeR normalized counts)"

  # ---- 14b. Transcript models from annotations_final.gtf -------------------
  gtf <- rtracklayer::import(final_gtf, feature.type = c("transcript", "exon", "CDS"))
  gtf <- keep_chroms(gtf, "annotations_final.gtf")
  md  <- S4Vectors::mcols(gtf)
  tx_gr <- gtf[as.character(md$type) == "transcript"]
  tmd   <- S4Vectors::mcols(tx_gr)
  gcol  <- function(n) {
    if (n %in% colnames(tmd)) as.character(tmd[[n]]) else rep(NA_character_, length(tx_gr))
  }
  # known rows use GENCODE/Ensembl attribute names, novel rows pulposeq_* ones
  first_filled <- function(...) {
    vs <- list(...)
    out <- vs[[1]]
    for (v in vs[-1]) {
      fill <- (is.na(out) | !nzchar(out)) & !is.na(v) & nzchar(v)
      out[fill] <- v[fill]
    }
    out
  }
  gtx <- data.frame(
    transcript_id          = gcol("transcript_id"),
    gene_id                = gcol("gene_id"),
    gtf_name               = gcol("transcript_name"),
    gtf_biotype            = first_filled(gcol("pulposeq_transcript_biotype"),
                                          gcol("transcript_biotype"), gcol("transcript_type")),
    status                 = gcol("transcript_status"),
    classification         = gcol("classification"),
    compared_gene_id       = gcol("compared_gene_id"),
    # the reference transcript gffcompare matched each novel model against
    compared_transcript_id = gcol("compared_transcript_id"),
    chrom                  = as.character(GenomeInfoDb::seqnames(tx_gr)),
    strand                 = as.character(BiocGenerics::strand(tx_gr)),
    start                  = BiocGenerics::start(tx_gr),
    end                    = BiocGenerics::end(tx_gr),
    stringsAsFactors       = FALSE
  )
  gtx <- gtx[!is.na(gtx$transcript_id) & !duplicated(gtx$transcript_id), , drop = FALSE]
  gtx$source <- "final"

  # ---- 14c. Reference annotation (the GTF pulposeq was given as --annotation) --
  # Required: gffcompare assigned every class code against the FULL annotation,
  # including isoforms not detected here. Without it the gene row omits their
  # exons and a panel can contradict its own labels (pulposeq genomic_context.R).
  if (is.null(cfg$ref_gtf) || !file.exists(resolve(cfg$ref_gtf))) {
    stop("cfg$ref_gtf not found (", cfg$ref_gtf %||% "NULL", "). Genomic context needs the ",
         "reference GTF pulposeq used as --annotation (ref_gtf in 00_paths.R).")
  }
  ref_all <- keep_chroms(load_reference(resolve(cfg$ref_gtf)), "reference GTF")
  rmd    <- S4Vectors::mcols(ref_all)
  r_type <- as.character(rmd$type)
  r_gid  <- bare_id(as.character(rmd$gene_id))
  r_tid  <- bare_id(as.character(rmd$transcript_id))
  gtype_col <- intersect(c("gene_type", "gene_biotype"), colnames(rmd))[1]
  is_rt  <- r_type == "transcript"
  ref_gene_type <- setNames(as.character(rmd[[gtype_col]])[is_rt], r_gid[is_rt])
  ref_gene_type <- ref_gene_type[!duplicated(names(ref_gene_type))]
  ref_gene_name <- setNames(as.character(rmd$gene_name)[is_rt], r_gid[is_rt])
  ref_gene_name <- ref_gene_name[!duplicated(names(ref_gene_name))]

  # keep what can be drawn: every reference exon (gene rows for ANY gene, so
  # other scripts can draw their own gene lists) + every row of the reference
  # transcripts that a novel model was classified against
  want_ref_tx <- gtx$compared_transcript_id[gtx$status %in% "novel"]
  want_ref_tx <- unique(bare_id(want_ref_tx[!is.na(want_ref_tx) & nzchar(want_ref_tx)]))
  keep <- r_type == "exon" | (!is.na(r_tid) & r_tid %in% want_ref_tx)
  ref_sub <- ref_all[keep]
  rm(ref_all, rmd, r_type, r_gid, r_tid, keep, is_rt, is_pc)
  invisible(gc())
  rsm <- S4Vectors::mcols(ref_sub)

  # gene row: the gene's full annotated exon footprint. Index-based split: a
  # plain split() on a GRanges falls back to base R and returns a list, which
  # reduce() cannot take.
  # (reduced per gene on demand in gene_structure(); no split() on GRanges)
  rex   <- ref_sub[as.character(rsm$type) == "exon"]
  rex_g <- bare_id(as.character(S4Vectors::mcols(rex)$gene_id))

  # referenced annotated transcripts with no record in the validated set
  rt  <- ref_sub[as.character(rsm$type) == "transcript"]
  rtm <- S4Vectors::mcols(rt)
  rcol <- function(n) if (n %in% colnames(rtm)) as.character(rtm[[n]]) else rep(NA_character_, length(rt))
  ref_tx <- data.frame(
    transcript_id          = rcol("transcript_id"),
    gene_id                = rcol("gene_id"),
    gtf_name               = rcol("transcript_name"),
    gtf_biotype            = first_filled(rcol("transcript_type"), rcol("transcript_biotype")),
    status                 = rep("annotated", length(rt)),
    classification         = NA_character_,
    compared_gene_id       = NA_character_,
    compared_transcript_id = NA_character_,
    chrom                  = as.character(GenomeInfoDb::seqnames(rt)),
    strand                 = as.character(BiocGenerics::strand(rt)),
    start                  = BiocGenerics::start(rt),
    end                    = BiocGenerics::end(rt),
    stringsAsFactors       = FALSE
  )
  ref_tx <- ref_tx[!is.na(ref_tx$transcript_id) &
                     !(bare_id(ref_tx$transcript_id) %in% bare_id(gtx$transcript_id)), , drop = FALSE]
  ref_tx$source <- rep("reference", nrow(ref_tx))
  log_step("  reference: %d genes with exon footprints, %d referenced transcripts not detected",
           length(unique(rex_g)), nrow(ref_tx))
  gtx <- rbind(gtx, ref_tx)

  # ---- 14d. DE results, labels and colours ------------------------------------
  di <- match(bare_id(gtx$transcript_id), bare_id(det_final$transcript_id))
  gtx$tested  <- !is.na(di)
  gtx$name    <- coalesce(det_final$transcript_name[di], gtx$gtf_name, gtx$transcript_id)
  gtx$biotype <- coalesce(det_final$transcript_type[di], gtx$gtf_biotype, "biotype not set")
  gtx$log2FC  <- det_final$log2FC[di]
  gtx$is_de   <- coalesce(det_final$is_de[di], FALSE)
  # gene membership = Bambu's assignment (GENEID of the counts edgeR used)
  gtx$gene_id <- ifelse(gtx$tested, coalesce(det_final$gene_id[di], gtx$gene_id), gtx$gene_id)

  short_lab <- identical(cfg$ctx_label_style, "short")
  cls_short <- sub("\\s*\\([^)]*\\)\\s*$", "", gtx$classification)   # "antisense intronic (i)" -> "antisense intronic"
  cls_txt <- ifelse(gtx$status %in% "novel" & !is.na(gtx$classification) &
                      nzchar(gtx$classification),
                    paste0(" | ", if (short_lab) cls_short else gtx$classification), "")
  gtx$label <- ifelse(
    gtx$status %in% "annotated", paste0(gtx$name, " | not detected"),
    ifelse(gtx$tested,
           sprintf("%s | %+.2f%s", gtx$name, gtx$log2FC, ifelse(gtx$is_de, " *", "")),
           paste0(gtx$name, " | not tested"))
  )
  # novel models whose Bambu gene differs from gffcompare's match are drawn in
  # the matched gene's panel for context only - the label names the Bambu gene
  # they were quantified and tested under
  bambu_lab <- coalesce(deg$gene_name[match(bare_id(gtx$gene_id), bare_id(deg$gene_id))], gtx$gene_id)
  other_gene <- gtx$status %in% "novel" & !is.na(gtx$compared_gene_id) &
    bare_id(gtx$gene_id) != bare_id(gtx$compared_gene_id)
  if (!short_lab) cls_txt <- paste0(cls_txt, ifelse(other_gene, paste0(" | Bambu: ", bambu_lab), ""))
  # one label per transcript, unique across the whole TxDb (pulposeq: uniquify the
  # transcript -> label mapping, never the GTF rows, or exons lose their parent)
  gtx$label <- make.unique(paste0(gtx$label, cls_txt), sep = " #")

  # same colours as the log2FC plots; non-DE, untested and undetected models
  # lightened (opaque blend with white - the Windows png device may lack alpha)
  lighten <- function(col, f = 0.35) {
    m <- grDevices::col2rgb(col) / 255
    grDevices::rgb(t(m * f + (1 - f)))
  }
  biotype_col <- function(bt) {
    out <- unname(pal[bt])
    out[is.na(out)] <- "#969696"
    out
  }
  gtx$colour <- ifelse(gtx$is_de, biotype_col(gtx$biotype), lighten(biotype_col(gtx$biotype)))

  gene_structure <- function(gid) {
    k <- bare_id(gid)
    hit <- which(rex_g == k)
    if (length(hit)) return(GenomicRanges::reduce(rex[hit]))
    # novel Bambu gene: no reference entry, use its own validated exons
    ex <- gtf[as.character(md$type) == "exon"]
    GenomicRanges::reduce(ex[bare_id(as.character(S4Vectors::mcols(ex)$gene_id)) %in% k])
  }

  # ---- 14e. Panels -----------------------------------------------------------
  # What a panel contains (pulposeq rule):
  #   detected   every validated model of the gene (up to context_max_tx, DE first)
  #   mine       novel models gffcompare classified against the gene but that
  #              Bambu assigned elsewhere (e.g. antisense / intronic lncRNAs with
  #              their own BambuGene id). Drawn for locus context only, labelled
  #              "| Bambu: <gene>", never counted as this gene's isoforms
  #   forced     the lncRNA of an lncRNA-known gene pair (always drawn)
  #   referenced undetected annotated transcripts that one of the drawn novel
  #              models was classified against, so each class code can be checked
  #              against the picture; always drawn, not counted in the cap
  build_panel <- function(gid, set, title, force_ids = character(0), note = "") {
    is_ann <- gtx$status %in% "annotated"
    own    <- bare_id(gtx$gene_id) %in% bare_id(gid) & !is_ann
    mine   <- gtx$status %in% "novel" & bare_id(gtx$compared_gene_id) %in% bare_id(gid)
    forced <- gtx[bare_id(gtx$transcript_id) %in% bare_id(force_ids) & !is_ann, , drop = FALSE]
    # validated models that were not DE-tested (filtered by expression) are
    # left out unless cfg$ctx_show_untested - they only add clutter
    shown  <- gtx$tested | isTRUE(cfg$ctx_show_untested)
    if (!isTRUE(cfg$ctx_show_non_de)) shown <- shown & gtx$is_de
    mem    <- gtx[(own | mine) & shown & !(gtx$transcript_id %in% forced$transcript_id), , drop = FALSE]
    if (!nrow(mem) && !nrow(forced)) return(NULL)
    main_chr <- if (nrow(forced)) forced$chrom[1] else names(which.max(table(mem$chrom)))
    mem <- mem[mem$chrom == main_chr, , drop = FALSE]
    mem <- mem[order(-mem$is_de, -mem$tested, -abs(coalesce(mem$log2FC, 0))), , drop = FALSE]
    n_all <- nrow(mem)
    mem <- rbind(forced, head(mem, cfg$context_max_tx))

    st <- gene_structure(gid)
    st <- st[as.character(GenomeInfoDb::seqnames(st)) == main_chr]
    s <- min(c(mem$start, BiocGenerics::start(st)))
    e <- max(c(mem$end, BiocGenerics::end(st)))
    pad <- ceiling((e - s) * 0.05)
    win_s <- max(1, s - pad)
    win_e <- e + pad

    ref_ids <- mem$compared_transcript_id[mem$status %in% "novel"]
    ref_ids <- unique(bare_id(ref_ids[!is.na(ref_ids) & nzchar(ref_ids)]))
    refd <- gtx[is_ann & bare_id(gtx$transcript_id) %in% ref_ids & gtx$chrom == main_chr &
                  gtx$end >= win_s & gtx$start <= win_e, , drop = FALSE]

    list(set = set, gene = title, gene_id = gid, note = note,
         chrom = main_chr, win_s = win_s, win_e = win_e,
         structure = st, tx = rbind(mem, refd),
         n_hidden = max(0L, n_all - cfg$context_max_tx))
  }
}

if (cfg$run_context) {
  # ---- 14f. TxDb for a set of panels (validated + referenced reference rows) ---
  lab_map <- setNames(gtx$label, gtx$transcript_id)
  minimal_gtf <- function(gr, ids) {
    m   <- S4Vectors::mcols(gr)
    tid <- as.character(m$transcript_id)
    gr  <- gr[!is.na(tid) & tid %in% ids]
    m   <- S4Vectors::mcols(gr)
    out <- GenomicRanges::granges(gr)
    cols <- S4Vectors::DataFrame(
      type          = as.character(m$type),
      gene_id       = as.character(m$gene_id),
      transcript_id = unname(lab_map[as.character(m$transcript_id)])
    )
    if ("phase" %in% colnames(m)) cols$phase <- m$phase
    S4Vectors::mcols(out) <- cols
    out
  }
  # sub-GTF of exactly the drawn models (also loadable in IGV) -> TxDb -> assembly
  make_assembly <- function(panels, regions_gtf) {
    all_drawn <- do.call(rbind, lapply(panels, `[[`, "tx"))
    ids_final <- unique(all_drawn$transcript_id[all_drawn$source == "final"])
    ids_ref   <- unique(all_drawn$transcript_id[all_drawn$source == "reference"])
    sub_gr <- c(minimal_gtf(gtf, ids_final), minimal_gtf(ref_sub, ids_ref))
    dir.create(dirname(regions_gtf), showWarnings = FALSE, recursive = TRUE)
    rtracklayer::export(sub_gr, regions_gtf, format = "gtf")
    txdb <- suppressWarnings(suppressMessages(
      if (requireNamespace("txdbmaker", quietly = TRUE)) {
        txdbmaker::makeTxDbFromGFF(regions_gtf, format = "gtf")
      } else {
        GenomicFeatures::makeTxDbFromGFF(regions_gtf, format = "gtf")
      }
    ))
    plotgardener::assembly(
      Genome = "dong_pulposeq", TxDb = txdb, OrgDb = NULL,
      gene.id.column = "GENEID", display.column = "GENEID"
    )
  }
  chrom_lens <- bw_lens   # from the bigWig header (health check above)

  # ---- 14g. Drawing (pulposeq geometry) ----------------------------------------
  # Geometry (inches) and font sizes (pt), computed for EVERY panel, so a script
  # that source()s this one can switch cfg$ctx_paper on for a single export
  # (isoformic_consolidated.R, section 7):
  #   ctx_paper = NULL   poster / report pages exactly as before (poster: fonts
  #                      x cfg$poster_ctx_font, every text band scaled with them)
  #   ctx_paper = list(width_mm, font_pt = 6, header = FALSE, ideogram = FALSE,
  #                    legend = TRUE, sig_h, tx_box_mm)
  #                      a journal panel at its final print size: fixed width,
  #                      absolute font sizes (Nature: 5-7 pt, never below 5) and,
  #                      by default, no title or chromosome bands - the figure
  #                      legend gives gene, region and FDR (the panel table
  #                      written by render_panels() has the window coordinates)
  ctx_geom <- function() {
    pp <- cfg$ctx_paper
    if (is.null(pp)) {
      FS <- if (is_poster()) cfg$poster_ctx_font else 1
      g <- list(
        width = if (is_poster()) cfg$poster_ctx_width else 8,
        margin_l = 0.4, margin_r = 0.4, margin_t = 0.4, margin_b = 0.4,
        header = TRUE, ideogram = TRUE, legend = TRUE, wrap_legend = TRUE,
        f_title = 11 * FS, f_sub = 8 * FS, f_small = 7 * FS, f_lab = 8 * FS, f_tx = 6 * FS,
        title_h = 0.52 * FS, sub1_y = 0.18 * FS, sub2_y = 0.32 * FS,
        ideo_h = 0.17, ideo_lab = 0.22 * FS, zoom_h = 0.40, ideo_lwd = 1.2, zoom_lwd = 0.8,
        gene_h = 0.26, sig_h = 0.55, sig_lab_h = 0.17 * FS, gap = 0.12,
        label_h = 0.45 * FS, no_sig_h = 0.3,
        tx_label_h = 0.13 * FS, tx_row_gap = 0.08 * FS, tx_box_mm = 4, min_tx_h = 0.4,
        lgnd_char = 0.048 * FS, lgnd_key = 0.30, lgnd_h = 0.18 * FS, lgnd_step = 0.20 * FS,
        lgnd_pad = 0.14, inset_pad = 0.1, label_col = NULL, inline_legend = FALSE
      )
    } else {
      pt <- pp$font_pt %||% 6
      ln <- function(f, k = 1.35) k * f / 72             # height of one text line (in)
      ft <- max(5, pt - 0.5)                             # isoform labels
      g <- list(
        width = (pp$width_mm %||% 82.5) / 25.4,
        margin_l = 0.06 + 4 * pt * 0.0075,               # room for the "gene" label
        margin_r = 0.03, margin_t = 0.03, margin_b = 0.02,
        header = isTRUE(pp$header), ideogram = isTRUE(pp$ideogram),
        legend = !isFALSE(pp$legend), wrap_legend = FALSE,
        f_title = pt + 1, f_sub = pt, f_small = pt, f_lab = pt, f_tx = ft,
        title_h = ln(pt + 1) + 2 * ln(pt) + 0.03, sub1_y = ln(pt + 1), sub2_y = ln(pt + 1) + ln(pt),
        ideo_h = 0.07, ideo_lab = ln(pt), zoom_h = 0.12, ideo_lwd = 0.6, zoom_lwd = 0.4,
        gene_h = 0.08, sig_h = pp$sig_h %||% 0.26, sig_lab_h = ln(pt), gap = 0.04,
        label_h = 2.4 * pt / 72, no_sig_h = ln(pt),
        tx_label_h = ln(ft, 1.25), tx_row_gap = 0.02, tx_box_mm = pp$tx_box_mm %||% 1.6,
        min_tx_h = 0.1,
        lgnd_char = pt * 0.0069, lgnd_key = 0.14, lgnd_h = ln(pt), lgnd_step = ln(pt, 1.5),
        lgnd_pad = 0.03, lgnd_box = 0.75 * ln(pt), inset_pad = 0.01,
        label_col = NULL,                                # labels keep the biotype colour
        inline_legend = TRUE
      )
    }
    g$tx_space_h <- 0.2
    g$tx_stroke  <- 0.1
    g$tx_model_h <- (g$tx_box_mm / 25.4) * (1 + g$tx_space_h)
    g$tx_row_h   <- g$tx_label_h + g$tx_model_h + g$tx_row_gap
    g
  }

  # text drawn above each model. ctx_label_style = "strand" (journal panels):
  # name, strand and, for novel models, the class - no log2FC, which the
  # profile beside the locus shows. Other styles: the TxDb label (gtx$label).
  # The TxDb id is never changed, only the text on the page.
  ctx_disp_label <- function(ptx) {
    if (!identical(cfg$ctx_label_style, "strand")) return(ptx$label)
    st  <- ifelse(ptx$strand %in% c("+", "-"), sprintf(" (%s)", ptx$strand), "")
    cls <- ifelse(ptx$status %in% "novel" & !is.na(ptx$classification) & nzchar(ptx$classification),
                  paste0(" \u00b7 ", sub("\\s*\\([^)]*\\)\\s*$", "", ptx$classification)), "")
    nd  <- ifelse(ptx$status %in% "annotated", " \u00b7 not detected", "")
    paste0(ptx$name, st, cls, nd)
  }

  # one-line biotype legend (journal panels): key + text placed by text width,
  # wrapping only when the track is full (plotLegend spaces items evenly and
  # misaligns long names)
  inline_legend_layout <- function(bts, track_w, G) {
    txt <- gsub("_", " ", bts)
    w   <- G$lgnd_box + 0.03 + nchar(txt) * G$f_small * 0.0072 + 0.10
    x   <- numeric(length(w)); row <- integer(length(w)); cx <- 0; r <- 1L
    for (k in seq_along(w)) {
      if (cx > 0 && cx + w[k] - 0.10 > track_w) { r <- r + 1L; cx <- 0 }
      x[k] <- cx; row[k] <- r; cx <- cx + w[k]
    }
    list(bt = bts, txt = txt, x = x, row = row, n_rows = if (length(w)) max(row) else 0L)
  }

  slug <- function(x) substr(gsub("(^-|-$)", "", gsub("[^A-Za-z0-9]+", "-", x)), 1, 60)

  # plotText in the figure font (plotgardener passes it on to grid::gpar)
  pg_text <- function(...) plotgardener::plotText(..., fontfamily = theme_font)

  # plotgardener's internal plotChromGenomeLabel hard-codes fontface
  # "bold" for the chromosome name; patch that one assignment to "plain" once
  # (no bold in figures). If the patch fails the run goes on with bold names.
  local({
    ok <- tryCatch({
      f <- utils::getFromNamespace("plotChromGenomeLabel", "plotgardener")
      b <- deparse(body(f), width.cutoff = 500L)
      hit <- grepl('object$gp$fontface <- "bold"', b, fixed = TRUE)
      if (any(hit)) {
        b <- gsub('object$gp$fontface <- "bold"', 'object$gp$fontface <- "plain"', b, fixed = TRUE)
        body(f) <- parse(text = b, keep.source = FALSE)[[1]]
        utils::assignInNamespace("plotChromGenomeLabel", f, "plotgardener")
      }
      sum(hit)
    }, error = function(e) { log_step("    chromosome label left bold: %s", conditionMessage(e)); 0L })
    log_step("    plotgardener chromosome label: %d bold assignment(s) set to plain", ok)
  })

  # rows packed so no model or label overlaps; also sizes the track exactly
  pack_rows <- function(starts, ends, labels, win_s, win_e, track_w, fontsize) {
    n <- length(starts)
    if (!n) return(integer(0))
    span    <- max(1, win_e - win_s)
    char_in <- fontsize * 0.0075
    lab_bp  <- nchar(labels) * char_in / track_w * span
    vs  <- pmax(starts, win_s)
    ve  <- pmin(ends, win_e)
    mid <- (vs + ve) / 2
    l   <- pmin(vs, mid - lab_bp / 2)
    r   <- pmax(ve, mid + lab_bp / 2)
    pad <- 0.02 * span
    row_end <- numeric(0)
    rows    <- integer(n)
    for (i in order(l)) {
      slot <- which(row_end + pad < l[i])
      k <- if (length(slot)) slot[1] else length(row_end) + 1L
      row_end[k] <- r[i]
      rows[i]    <- k
    }
    rows
  }

  # chromosome bar + zoom lines, from the bigWig header (no AnnotationHub)
  draw_ideogram <- function(chrom, win_s, win_e, chrom_len, x, y, w, h,
                            panel_x0, panel_x1, zoom_y1,
                            fontsize = 8, lwd = 1.2, zoom_lwd = 0.8) {
    plotgardener::plotRect(x = x, y = y, width = w, height = h, just = c("left", "top"),
                           default.units = "inches", fill = "#F0F0F0", linecolor = "#BDBDBD")
    f0  <- max(0, min(1, win_s / chrom_len))
    f1  <- max(0, min(1, win_e / chrom_len))
    rx0 <- x + f0 * w
    rx1 <- max(x + f1 * w, rx0 + 0.03)
    plotgardener::plotRect(x = rx0, y = y, width = rx1 - rx0, height = h,
                           just = c("left", "top"), default.units = "inches",
                           fill = "#C6E48B", linecolor = NA)
    plotgardener::plotSegments(x0 = (rx0 + rx1) / 2, y0 = y - 0.05,
                               x1 = (rx0 + rx1) / 2, y1 = y + h + 0.05,
                               default.units = "inches", linecolor = "#D62728", lwd = lwd)
    pg_text(label = sprintf("Chromosome %s", sub("^chr", "", chrom)),
                           x = x + w, y = y + h + 0.07, just = c("right", "top"),
                           fontsize = fontsize, fontcolor = "grey45", default.units = "inches")
    plotgardener::plotSegments(x0 = rx0, y0 = y + h + 0.03, x1 = panel_x0, y1 = zoom_y1,
                               default.units = "inches", linecolor = "#CCCCCC", lty = 2, lwd = zoom_lwd)
    plotgardener::plotSegments(x0 = rx1, y0 = y + h + 0.03, x1 = panel_x1, y1 = zoom_y1,
                               default.units = "inches", linecolor = "#CCCCCC", lty = 2, lwd = zoom_lwd)
  }

  # returns the showtext state to give back to showtext_resume() after dev.off()
  open_dev <- function(file, width, height) {
    st <- FALSE
    if (grepl("\\.pdf$", file)) {
      if (isTRUE(capabilities("cairo"))) {
        grDevices::cairo_pdf(file, width = width, height = height)
      } else {
        grDevices::pdf(file, width = width, height = height)
      }
    } else if (grepl("\\.svg$", file)) {
      if (requireNamespace("svglite", quietly = TRUE)) {
        st <- showtext_pause()                   # editable <text>, not glyph outlines
        svg_device(file, width = width, height = height)
      } else {
        stop("svglite is needed for SVG with editable text (install.packages(\"svglite\")); base svg() would turn the text into glyph outlines")
      }
    } else if (isTRUE(capabilities("cairo"))) {
      grDevices::png(file, width = width, height = height, units = "in", res = cfg$dpi, type = "cairo")
    } else {
      grDevices::png(file, width = width, height = height, units = "in", res = cfg$dpi)
    }
    invisible(st)
  }

  # inset: optional list(plot = <ggplot>, w = <in>, h = <in>) drawn at the left
  # of the page (e.g. the gene's expression profile); the tracks move right
  draw_context_panel <- function(p, out_png, ctx_assembly, inset = NULL) {
    G <- ctx_geom()                              # poster / report / paper geometry
    # coverage: bigWig chromosome naming, window clamped to the chromosome end,
    # one failing file never takes the whole panel down
    bc   <- bw_chrom(p$chrom)
    clen <- if (!is.null(bw_lens) && !is.na(bc)) as.numeric(bw_lens[[bc]]) else Inf
    sig <- list()
    sig_err <- character(0)
    for (s in names(bw_files)[bw_ok]) {
      g <- tryCatch({
        if (is.na(bc)) stop("chromosome ", p$chrom, " is not in the bigWig")
        w_s <- max(1, p$win_s)
        w_e <- min(p$win_e, clen)
        x <- bw_read(bw_files[[s]], function(path) {
          rtracklayer::import.bw(path, which = GenomicRanges::GRanges(bc, IRanges::IRanges(w_s, w_e)))
        })
        if (!length(x)) {
          GenomicRanges::GRanges(p$chrom, IRanges::IRanges(w_s, w_e), score = 0)
        } else {
          GenomicRanges::GRanges(p$chrom, IRanges::ranges(x), score = x$score * scale_f[[s]])
        }
      }, error = function(e) {
        sig_err <<- c(sig_err, sprintf("%s: %s", s, conditionMessage(e)))
        NULL
      })
      if (!is.null(g)) sig[[s]] <- g
    }
    if (length(sig_err) && !bw_state$warned) {
      log_step("  ! coverage read failed (first case, %s %s:%d-%d): %s", p$gene, p$chrom,
               p$win_s, p$win_e, sig_err[1])
      bw_state$warned <- TRUE
    }
    # one track per sample, or the mean of each condition's samples
    track_group <- sample_group[names(sig)]
    track_lab   <- if (is.null(cfg$ctx_group_labels)) sprintf("%s  (%s)", names(sig), track_group) else
      names(sig)                                 # journal panels: the sample name says the line
    if (identical(cfg$ctx_coverage, "group_mean") && length(sig)) {
      w_s <- max(1, p$win_s)
      w_e <- min(p$win_e, clen)
      to_rle <- function(g) {
        cv <- IRanges::coverage(g, weight = "score")[[p$chrom]]
        if (length(cv) < w_e) cv <- c(cv, S4Vectors::Rle(0, w_e - length(cv)))
        cv[w_s:w_e]
      }
      grp_of <- sample_group[names(sig)]
      n_grp  <- table(grp_of)
      sig <- lapply(split(names(sig), grp_of), function(ss) {
        m <- Reduce(`+`, lapply(sig[ss], to_rle)) / length(ss)
        l <- S4Vectors::runLength(m)
        e <- w_s - 1 + cumsum(l)
        GenomicRanges::GRanges(p$chrom, IRanges::IRanges(e - l + 1, e),
                               score = as.numeric(S4Vectors::runValue(m)))
      })
      sig <- sig[intersect(c(cfg$control, cfg$case), names(sig))]
      track_group <- setNames(names(sig), names(sig))
      glab <- if (is.null(cfg$ctx_group_labels)) names(sig) else
        dplyr::coalesce(unname(cfg$ctx_group_labels[names(sig)]), names(sig))
      track_lab   <- sprintf("%s  (mean of %d)", glab, as.integer(n_grp[names(sig)]))
    }

    # optional log2(depth + 1) axis (track labels say so)
    log_cov <- identical(cfg$ctx_cov_scale, "log2")
    if (log_cov) sig <- lapply(sig, function(g) { g$score <- log2(g$score + 1); g })
    # shared y range across tracks so heights compare directly
    ymax <- if (length(sig)) {
      max(c(vapply(sig, function(g) max(g$score, na.rm = TRUE), numeric(1)), 1e-6))
    } else 1

    ptx <- p$tx
    disp <- ctx_disp_label(ptx)                  # text on the page (TxDb id stays ptx$label)
    g_row <- deg[match(bare_id(p$gene_id), bare_id(deg$gene_id)), ]
    g_txt <- if (nrow(g_row) && !is.na(g_row$log2FC)) {
      sprintf("gene log2FC %+.2f, FDR %.2g%s", g_row$log2FC, g_row$qvalue,
              if (isTRUE(g_row$is_de)) " (DE)" else " (not DE)")
    } else "gene not tested"
    sub1 <- if (is_poster()) {
      paste0(g_txt, if (nzchar(p$note %||% "")) paste0(" | ", p$note) else "")
    } else {
      sprintf("%s | plot set: %s%s", g_txt, p$set,
              if (nzchar(p$note %||% "")) paste0(" | ", p$note) else "")
    }
    cov_txt <- if (length(sig)) {
      paste0(if (is_poster()) "TMM-scaled depth" else depth_lab,
             if (identical(cfg$ctx_coverage, "group_mean")) ", condition mean" else ", per library",
             if (identical(cfg$ctx_cov_scale, "log2")) ", log2(depth + 1)" else "",
             if (isTRUE(cfg$ctx_cov_highlight) && any(ptx$status %in% "novel")) "; shaded: novel models" else "")
    } else "NOT DRAWN (bigWig unreadable, see log)"
    n_nd <- sum(ptx$status %in% "annotated")
    sub2 <- if (is_poster()) {
      sprintf("%d models (%d DE%s)%s%s | coverage: %s",   # (15 Nov c) "* FDR" dropped: no asterisks drawn
              nrow(ptx), sum(ptx$is_de), if (n_nd) sprintf(", %d not detected", n_nd) else "",
              if (p$n_hidden > 0) sprintf(", +%d hidden", p$n_hidden) else "",
              if (isTRUE(cfg$ctx_show_non_de)) "" else ", non-DE isoforms hidden",
              cov_txt)
    } else {
      sprintf("%d models drawn (%d DE, %d annotated not detected)%s | lightened = not DE | coverage: %s",
              nrow(ptx), sum(ptx$is_de), sum(ptx$status %in% "annotated"),
              if (p$n_hidden > 0) sprintf(", %d more not shown", p$n_hidden) else "",
              cov_txt)
    }

    width <- G$width; track_w <- width - G$margin_l - G$margin_r
    ix <- if (is.null(inset)) 0 else inset$w    # inset width: tracks start at xl
    xl <- G$margin_l + ix
    title_h <- if (G$header) G$title_h else 0
    ideo_h <- G$ideo_h; ideo_lab <- G$ideo_lab; zoom_h <- G$zoom_h
    gene_h <- G$gene_h; sig_h <- G$sig_h; sig_lab_h <- G$sig_lab_h; gap <- G$gap; label_h <- G$label_h

    tx_rows <- pack_rows(ptx$start, ptx$end, disp, p$win_s, p$win_e, track_w, G$f_tx)
    n_rows  <- if (length(tx_rows)) max(tx_rows) else 1L
    tx_h    <- max(G$min_tx_h, n_rows * G$tx_row_h)

    chrom_len <- if (!is.null(chrom_lens) && !is.na(bc) && bc %in% names(chrom_lens)) {
      as.numeric(chrom_lens[[bc]])
    } else NA_real_
    show_ideo <- G$ideogram && is.finite(chrom_len) && chrom_len > 0
    head_h    <- if (show_ideo) ideo_h + ideo_lab + zoom_h else 0

    bts <- sort(unique(ptx$biotype))
    per_row <- max(1L, floor(track_w / (max(nchar(bts)) * G$lgnd_char + G$lgnd_key)))
    lgnd_rows <- if (G$legend) unname(split(bts, ceiling(seq_along(bts) / per_row))) else list()
    lgnd_h <- if (G$legend) length(lgnd_rows) * G$lgnd_step + G$lgnd_pad else 0
    lg_inline <- if (G$legend && isTRUE(G$inline_legend)) inline_legend_layout(bts, track_w, G) else NULL
    if (!is.null(lg_inline)) {
      lgnd_rows <- list()                        # drawn as the inline legend below
      lgnd_h    <- lg_inline$n_rows * G$lgnd_step + G$lgnd_pad
    }

    height <- G$margin_t + G$margin_b + title_h + head_h + gene_h + tx_h + gap * 3 +
      length(sig) * (sig_lab_h + sig_h + gap) + (if (!length(sig)) G$no_sig_h else 0) +
      label_h + lgnd_h
    if (!is.null(inset)) {
      # inset height NA = as tall as the tracks (journal panels: profile beside the locus)
      if (is.na(inset$h %||% NA_real_)) {
        inset$h <- max(height - G$margin_t - G$margin_b - title_h, inset$h_min %||% 0)
      }
      height <- max(height, G$margin_t + G$margin_b + title_h + inset$h)
    }

    st_dev <- open_dev(out_png, width + ix, height)
    on.exit({ grDevices::dev.off(); showtext_resume(st_dev) }, add = TRUE)
    plotgardener::pageCreate(width = width + ix, height = height, default.units = "inches",
                             showGuides = FALSE)
    pars <- plotgardener::pgParams(
      chrom = p$chrom, chromstart = p$win_s, chromend = p$win_e,
      assembly = ctx_assembly, x = xl, width = track_w, default.units = "inches"
    )
    if (!is.null(inset)) {
      plotgardener::plotGG(plot = inset$plot, x = G$inset_pad, y = G$margin_t + title_h,
                           width = inset$w - 2 * G$inset_pad, height = inset$h,
                           just = c("left", "top"), default.units = "inches")
    }

    if (G$header) {
      pg_text(
        label = sprintf("%s  |  %s:%s-%s", p$gene, p$chrom,
                        format(p$win_s, big.mark = ","), format(p$win_e, big.mark = ",")),
        x = G$margin_l, y = G$margin_t, just = c("left", "top"),
        fontsize = G$f_title, default.units = "inches"
      )
      pg_text(label = sub1, x = G$margin_l, y = G$margin_t + G$sub1_y,
                             just = c("left", "top"), fontsize = G$f_sub, fontcolor = "grey35",
                             default.units = "inches")
      pg_text(label = sub2, x = G$margin_l, y = G$margin_t + G$sub2_y,
                             just = c("left", "top"), fontsize = G$f_small, fontcolor = "grey45",
                             default.units = "inches")
    }

    y <- G$margin_t + title_h
    if (show_ideo) {
      draw_ideogram(p$chrom, p$win_s, p$win_e, chrom_len,
                    x = xl, y = y, w = track_w, h = ideo_h,
                    panel_x0 = xl, panel_x1 = xl + track_w,
                    zoom_y1 = y + ideo_h + ideo_lab + zoom_h - 0.04,
                    fontsize = G$f_lab, lwd = G$ideo_lwd, zoom_lwd = G$zoom_lwd)
      y <- y + ideo_h + ideo_lab + zoom_h
    }

    if (length(p$structure)) {
      plotgardener::plotRanges(params = pars, data = p$structure, collapse = TRUE,
                               fill = "#4D4D4D", linecolor = NA, y = y, height = gene_h)
      pg_text(label = "gene", x = xl - 0.06, y = y + gene_h / 2,
                             just = c("right", "center"), fontsize = G$f_small,
                             fontcolor = "grey35", default.units = "inches")
    }
    y <- y + gene_h + gap

    # isoform track: one plotTranscripts call per model (see pulposeq notes)
    gx <- function(pos) {
      xl + (min(max(pos, p$win_s), p$win_e) - p$win_s) / max(1, p$win_e - p$win_s) * track_w
    }
    for (i in seq_len(nrow(ptx))) {
      lb <- ptx$label[i]
      ry <- y + (tx_rows[i] - 1L) * G$tx_row_h
      plotgardener::plotTranscripts(
        params = pars, y = ry + G$tx_label_h, height = G$tx_model_h,
        labels = NULL,
        boxHeight = grid::unit(G$tx_box_mm, "mm"), spaceHeight = G$tx_space_h,
        stroke = G$tx_stroke, limitLabel = FALSE,
        fill = c("#969696", "#969696"), colorbyStrand = FALSE,
        transcriptFilter = lb,
        transcriptHighlights = data.frame(transcript = lb, color = ptx$colour[i],
                                          stringsAsFactors = FALSE)
      )
      lw <- nchar(disp[i]) * G$f_tx * 0.0075
      lx <- min(max((gx(ptx$start[i]) + gx(ptx$end[i])) / 2, xl + lw / 2),
                xl + track_w - lw / 2)
      pg_text(
        label = disp[i], x = lx, y = ry, just = c("center", "top"),
        fontsize = G$f_tx,
        fontface = "plain",                      # no bold in figures; non-DE models are grey
        fontcolor = if (ptx$is_de[i]) (G$label_col %||% biotype_col(ptx$biotype[i])) else "grey45",
        default.units = "inches"
      )
    }
    y <- y + tx_h + gap

    # coverage: control samples first, then case; coloured by condition
    if (!length(sig)) {
      pg_text(label = "coverage not drawn: bigWig could not be read (see log)",
                             x = xl, y = y, just = c("left", "top"), fontsize = G$f_lab,
                             fontcolor = "#b2182b", default.units = "inches")
      y <- y + G$no_sig_h
    }
    # optional shading of the novel models' spans across the coverage tracks
    y_cov0 <- y
    for (k in seq_along(sig)) {
      s  <- names(sig)[k]
      gc <- unname(grp_col[track_group[[s]]])
      if (is.na(gc)) gc <- "#3B6FB6"
      if (isTRUE(cfg$ctx_cov_highlight) && k == 1L) {
        nv <- p$tx[p$tx$status %in% "novel", , drop = FALSE]
        cov_h <- length(sig) * (sig_lab_h + sig_h + gap) - gap
        for (r in seq_len(nrow(nv))) {
          x0 <- xl + (max(nv$start[r], p$win_s) - p$win_s) / (p$win_e - p$win_s) * track_w
          x1 <- xl + (min(nv$end[r], p$win_e) - p$win_s) / (p$win_e - p$win_s) * track_w
          plotgardener::plotRect(x = x0, y = y_cov0, width = max(x1 - x0, 0.01), height = cov_h,
                                 just = c("left", "top"), default.units = "inches",
                                 fill = "#fff2b3", linecolor = NA)
        }
      }
      pg_text(
        label = sprintf(if (log_cov) "[0 - %s] log2(depth + 1)" else "[0 - %s]",
                        format(signif(ymax, 3), scientific = FALSE, trim = TRUE)),
        x = xl, y = y, just = c("left", "top"),
        fontsize = G$f_small, fontcolor = "grey35", default.units = "inches"
      )
      pg_text(
        label = track_lab[k], x = xl + track_w, y = y,
        just = c("right", "top"), fontsize = G$f_lab, fontcolor = gc, default.units = "inches"
      )
      plotgardener::plotSignal(
        params = pars, data = sig[[s]], y = y + sig_lab_h, height = sig_h,
        range = c(0, ymax), linecolor = gc, fill = gc
      )
      y <- y + sig_lab_h + sig_h + gap
    }

    plotgardener::plotGenomeLabel(
      chrom = p$chrom, chromstart = p$win_s, chromend = p$win_e,
      assembly = ctx_assembly, x = xl, y = y, length = track_w,
      default.units = "inches", scale = "bp", fontsize = G$f_lab, fontfamily = theme_font
    )

    ly <- y + label_h
    for (r in lgnd_rows) {
      plotgardener::plotLegend(
        legend = if (G$wrap_legend) pretty_type(r) else gsub("_", " ", r),
        fill = biotype_col(r), border = FALSE, orientation = "h",
        fontsize = G$f_small, fontfamily = theme_font,
        x = xl, y = ly, width = track_w * length(r) / per_row,
        height = G$lgnd_h, just = c("left", "top"), default.units = "inches"
      )
      ly <- ly + G$lgnd_step
    }
    if (!is.null(lg_inline)) {
      for (k in seq_along(lg_inline$txt)) {
        yk <- y + label_h + (lg_inline$row[k] - 1L) * G$lgnd_step
        xk <- xl + lg_inline$x[k]
        plotgardener::plotRect(x = xk, y = yk, width = G$lgnd_box, height = G$lgnd_box,
                               just = c("left", "top"), default.units = "inches",
                               fill = biotype_col(lg_inline$bt[k]), linecolor = NA)
        pg_text(label = lg_inline$txt[k], x = xk + G$lgnd_box + 0.03,
                               y = yk + G$lgnd_box / 2, just = c("left", "center"),
                               fontsize = G$f_small, default.units = "inches")
      }
    }
    plotgardener::pageGuideHide()
    basename(out_png)
  }

  # Draw panels -> <plot_dir>/genomic_context_<set>_PNG/ and _SVG/, one file per
  # panel; writes a table of what was drawn and returns it
  # folder: output folder stem (default genomic_context_<set>); insets: list
  # parallel to panels of list(plot, w, h) or NULL (see draw_context_panel)
  render_panels <- function(panels, plot_dir, table_file, regions_gtf,
                            folder = NULL, insets = NULL, formats = cfg$ctx_formats) {
    keep_p <- !vapply(panels, is.null, logical(1))
    panels <- panels[keep_p]
    if (!is.null(insets)) insets <- insets[keep_p]
    if (!length(panels)) return(NULL)
    assembly <- make_assembly(panels, regions_gtf)
    out <- lapply(seq_along(panels), function(i) {
      p    <- panels[[i]]
      stem <- if (is.null(folder)) sprintf("%s_%s", p$set, slug(p$gene)) else slug(p$gene)
      fig  <- NA_character_
      for (fmt in formats) {
        fdir <- file.path(plot_dir, sprintf("%s_%s", folder %||% paste0("genomic_context_", p$set), fmt))
        dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
        f   <- file.path(fdir, paste0(stem, ".", tolower(fmt)))
        res <- tryCatch(
          suppressMessages(draw_context_panel(p, f, assembly,
                                              inset = if (is.null(insets)) NULL else insets[[i]])),
          error = function(e) {
            message("    ! context panel (", fmt, ") failed for ", p$gene, ": ", conditionMessage(e))
            NA_character_
          }
        )
        if (fmt == formats[1]) fig <- res
      }
      data.frame(set = p$set, gene_name = p$gene, gene_id = p$gene_id, chrom = p$chrom,
                 win_s = p$win_s, win_e = p$win_e, n_drawn = nrow(p$tx),
                 n_de = sum(p$tx$is_de), n_novel = sum(p$tx$status %in% "novel"),
                 n_not_detected = sum(p$tx$status %in% "annotated"),
                 n_hidden = p$n_hidden, figure = fig, stringsAsFactors = FALSE)
    })
    out <- dplyr::bind_rows(out)
    write_csv(out, table_file)
    log_step("  drew %d / %d genomic context panels -> %s/%s_<format>/",
             sum(!is.na(out$figure)), nrow(out), plot_dir, folder %||% "genomic_context_<set>")
    out
  }
}

## ---- 15. Run info -----------------------------------------------------------
writeLines(
  c(
    paste("Run:", format(Sys.time())),
    paste("Transcript counts:", tx_counts_path),
    paste("Gene counts:", if (is.na(gene_counts_path)) "summed from transcript counts" else gene_counts_path),
    paste("Normalised counts (transcript):", norm_tx$source),
    paste("Normalised counts (gene):", norm_gene$source),
    "", "Config:", capture.output(str(cfg)),
    "", capture.output(sessionInfo())
  ),
  file.path(out_root, "run_info.txt")
)

log_step("Done.")
log_step("  DETs whose gene is not DE : %d", nrow(dets_not_degs))
log_step("  isoform-switch genes      : %d", nrow(switch_genes))
log_step("  annotated genes with DE novel isoforms: %d", nrow(novel_in_annot))
log_step("  results in %s", out_root)

## =============================================================================
## Curated events, cis pairs and journal panels (isoformic_consolidated2.R)
## =============================================================================
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
log_step("06_isoformic: %d events, output %s", nrow(cs$events), out_root)

## ---- 1. Helpers -------------------------------------------------------------
plot_dir <- file.path(out_root, "plots")
pout     <- function(x) file.path(plot_dir, x)
bare     <- strip_version

# temporarily change the base font size used by poster_theme() / base_pt()
with_base <- function(size, expr) {
  old <- cfg$base_size
  cfg$base_size <<- size
  on.exit(cfg$base_size <<- old, add = TRUE)
  force(expr)
}

# "antisense intronic (i)" -> "antisense intronic"
short_class <- function(x) {
  x <- sub("\\s*\\([^)]*\\)\\s*$", "", x)
  ifelse(is.na(x) | !nzchar(x), "novel", x)
}

# gene symbol for a Bambu gene id; drops the "_ENSG..." suffix edgeR added to
# duplicated symbols (e.g. LSP1P5_ENSG00000288905.3 -> LSP1P5)
gene_label <- function(gid) {
  g <- deg$gene_name[match(bare(gid), bare(deg$gene_id))]
  g <- coalesce(g, unname(ref_gene_name[bare(gid)]), gid)
  sub("_ENS[GT][0-9]+(\\.[0-9]+)?$", "", g)
}

## ---- 2. Resolve events --------------------------------------------------------
# extras = transcripts added to the gene's plots: partner-gene DE transcripts,
# DE novel models gffcompare matched to the gene, explicit extra_tx
resolve_event <- function(r) {
  gene <- r$gene
  if (!gene %in% de_long$gene_name) {
    log_step("  ! %s: not in the DE tables - skipped", gene)
    return(NULL)
  }
  gid <- de_long$gene_id[match(gene, de_long$gene_name)]
  ex  <- tibble::tibble(transcript_id = character(), transcript_name = character(),
                        label = character(), source = character())

  # partner gene: drawn in its own profile panel below the gene (merged into
  # the gene's panel, its isoforms would fill the isoform cap). Its top DE
  # isoforms (same
  # switch-aware rule, cs$pair_max_iso) go to the log2FC bars and the locus.
  partners <- character(0)
  pt_ids   <- character(0)
  if (!is.na(r$partner_gene)) {
    pt <- det_final |>
      filter(is_de, gene_name %in% r$partner_gene |
               sub("-[0-9]+$", "", transcript_name) %in% r$partner_gene)
    if (!nrow(pt)) log_step("  ! %s: no DE transcript of partner %s", gene, r$partner_gene)
    if (!r$partner_gene %in% profile_df$parent_gene) {
      log_step("  ! %s: partner %s not in the profile data - not drawn", gene, r$partner_gene)
    } else partners <- r$partner_gene
    top <- pick_isoforms(pt$transcript_name, cs$pair_max_iso %||% 5L)
    pt_ids <- pt$transcript_id[match(top, pt$transcript_name)]
    log_step("  %s: partner %s below it (%d of %d DE isoforms: %s)", gene, r$partner_gene,
             length(top), nrow(pt), paste(top, collapse = ", "))
  }

  if (isTRUE(r$novel_matched)) {
    nv <- gtx |>
      filter(status %in% "novel", tested, is_de,
             bare(compared_gene_id) %in% bare(gid), !bare(gene_id) %in% bare(gid))
    di <- match(bare(nv$transcript_id), bare(det_final$transcript_id))
    ex <- bind_rows(ex, tibble::tibble(
      transcript_id = det_final$transcript_id[di], transcript_name = nv$name,
      label = nv$name,
      source = "novel, gffcompare-matched"))
  }

  if (!is.na(r$extra_tx)) {
    ids <- trimws(strsplit(r$extra_tx, ";", fixed = TRUE)[[1]])
    i   <- coalesce(match(ids, det_final$transcript_name),
                    match(bare(ids), bare(det_final$transcript_id)))
    if (anyNA(i)) log_step("  ! %s: extra transcript(s) not found: %s", gene, paste(ids[is.na(i)], collapse = ", "))
    i   <- i[!is.na(i)]
    nm  <- det_final$transcript_name[i]
    ex  <- bind_rows(ex, tibble::tibble(
      transcript_id = det_final$transcript_id[i], transcript_name = nm,
      label = if (!is.na(r$extra_label)) sprintf("%s\n(%s)", nm, r$extra_label) else nm,
      source = "extra"))
  }

  ex <- ex |> filter(!is.na(transcript_id)) |> distinct(transcript_id, .keep_all = TRUE)

  # an "isoform switch" needs DE isoforms in BOTH directions; events chosen on
  # the Bambu counts may no longer be switches with oarfish counts
  own  <- det_final |> filter(is_de, gene_name == gene)
  n_up <- sum(own$log2FC > 0, na.rm = TRUE); n_dn <- sum(own$log2FC < 0, na.rm = TRUE)
  if (grepl("switch", r$title, ignore.case = TRUE) && (n_up == 0 || n_dn == 0))
    log_step("  ! %s: title says 'switch' but its DE isoforms all go one way (%d up, %d down)",
             gene, n_up, n_dn)
  list(row = r, gene = gene, gene_id = gid, extras = ex, partners = partners,
       partner_ids = pt_ids, n_up = n_up, n_down = n_dn)
}

log_step("Resolving events")
evs <- lapply(seq_len(nrow(cs$events)), function(i) resolve_event(cs$events[i, ]))
evs <- evs[!vapply(evs, is.null, logical(1))]
names(evs) <- vapply(evs, `[[`, character(1), "gene")
for (e in evs) {
  log_step("  %-9s + %s", e$gene,
           if (nrow(e$extras)) paste(e$extras$transcript_name, collapse = ", ") else "-")
}
bind_rows(lapply(evs, function(e) {
  tibble::tibble(gene = e$gene, gene_id = e$gene_id, title = e$row$title, axis = e$row$axis,
                 panel = e$row$panel, context = e$row$context,
                 de_isoforms_up = e$n_up, de_isoforms_down = e$n_down,
                 switch = e$n_up > 0 & e$n_down > 0,
                 added_transcripts = paste(e$extras$transcript_name, collapse = ";"),
                 added_sources = paste(e$extras$source, collapse = ";"))
})) |> write_csv(out_t("20_events.csv"))

## ---- 3. Profiles --------------------------------------------------------------
# the gene's profile rows + the extras (relabelled, moved under the gene)
event_profile_df <- function(ev) {
  base <- profile_df |> filter(parent_gene == ev$gene, !genename %in% ev$extras$transcript_name)
  add  <- profile_df |> filter(genename %in% ev$extras$transcript_name, transcript_type != "gene")
  add$parent_gene <- ev$gene
  add$genename    <- ev$extras$label[match(add$genename, ev$extras$transcript_name)]
  # partner gene (e.g. UBE2I-AS1): all its rows, kept under its own name, so
  # plot_gene_profile() draws it as a second panel below the gene
  pg <- profile_df |> filter(parent_gene %in% (ev$partners %||% character(0)),
                             parent_gene != ev$gene)
  bind_rows(base, add, pg)
}

event_profile <- function(ev, title = ev$row$title, type_limits = NULL, subtitle = ev$row$axis) {
  safe_plot(
    plot_gene_profile(event_profile_df(ev), ev$gene, pal, keep = ev$extras$label,
                      de_only = cs$profile_de_only && !ev$gene %in% cs$profile_all_isoforms,
                      title = title, subtitle = subtitle,
                      type_limits = type_limits, partner = ev$partners %||% character(0)),
    ev$gene
  )
}

log_step("Profiles")
# 3a consolidated panel: one shared legend (fixed transcript-type limits make
# the per-plot legends identical, so patchwork merges them)
panel_evs <- if (isTRUE(cs$panel_all)) evs else evs[vapply(evs, function(e) isTRUE(e$row$panel), logical(1))]
panel_evs <- panel_evs[!names(panel_evs) %in% c(cs$panel_exclude, cs$review_events, cs$retired_events)]
if (length(panel_evs)) {
  wrap_t <- function(ev) stringr::str_wrap(ev$row$title, cs$panel_title_wrap)
  first  <- lapply(panel_evs, function(ev) event_profile(ev, title = wrap_t(ev)))
  types_of <- function(ev, p) {
    t <- tryCatch(as.character(p$data$transcript_type), error = function(e) NULL)
    if (is.null(t)) t <- as.character(event_profile_df(ev)$transcript_type)
    t
  }
  ok     <- !vapply(first, is.null, logical(1))
  types  <- unique(unlist(Map(types_of, panel_evs[ok], first[ok])))
  types  <- c(intersect("gene", types), sort(setdiff(types, "gene")))
  profs  <- lapply(panel_evs, function(ev) event_profile(ev, title = wrap_t(ev), type_limits = types))
  profs  <- profs[!vapply(profs, is.null, logical(1))]
  # single shared legend: only the first profile keeps its guides (the fixed
  # type limits make it list every transcript type in the panel); patchwork
  # collects that one legend to the right of the whole panel
  if (length(profs) > 1) {
    profs[-1] <- lapply(profs[-1], function(p) p + guides(colour = "none", fill = "none", alpha = "none"))
  }
  lay    <- if (utils::packageVersion("patchwork") >= "1.2.0") {
    patchwork::plot_layout(ncol = cs$panel_ncol, guides = "collect", axis_titles = "collect")
  } else {
    patchwork::plot_layout(ncol = cs$panel_ncol, guides = "collect")
  }
  panel  <- patchwork::wrap_plots(profs) + lay
  n_row <- ceiling(length(profs) / cs$panel_ncol)
  with_base(cs$panel_base, save_pages(
    list(panel), pout("01_profiles_panel.pdf"),
    width = cs$panel_ncol * cs$panel_cell[1] + 2.0, height = n_row * cs$panel_cell[2] + 0.3,
    page_names = "profiles_panel"
  ))
}

# 3b one profile per event
save_pages(lapply(evs, event_profile), pout("02_profiles.pdf"), kind = "profile",
           page_names = names(evs))

## ---- 4. log2FC bars -----------------------------------------------------------
event_bars <- function(ev) {
  ids <- c(ev$extras$transcript_id, ev$partner_ids)
  d <- bind_rows(
    de_long |> filter(gene_name == ev$gene, !feature_id %in% ids),
    de_long |> filter(feature_id %in% ids, feature_type != "gene") |> mutate(gene_name = ev$gene)
  )
  n_tx <- sum(d$feature_type != "gene")
  note <- if (n_tx > cfg$poster_max_iso) sprintf("%d of %d isoforms shown", cfg$poster_max_iso, n_tx)
  safe_plot(
    plot_gene_log2fc(d, ev$gene, pal, contrast_lab, alpha_lab, keep_ids = ids) +
      labs(title = ev$row$title, subtitle = paste(c(ev$row$axis, note), collapse = " \u00b7 ")),
    ev$gene
  )
}
log_step("log2FC bars")
save_pages(lapply(evs, event_bars), pout("03_log2fc.pdf"), kind = "log2fc", page_names = names(evs))

## ---- 5. Genomic context + event panels --------------------------------------
ctx_evs <- evs[vapply(evs, function(e) isTRUE(e$row$context), logical(1))]
if (isTRUE(cs$profiles_only)) ctx_evs <- list()      # profiles-only rerun
if (length(ctx_evs)) {
  log_step("Genomic context (%d events)", length(ctx_evs))
  ctx_panels <- lapply(ctx_evs, function(ev) {
    build_panel(ev$gene_id, "event", ev$row$title, force_ids = c(ev$extras$transcript_id, ev$partner_ids),
                note = ev$row$axis)
  })
  # 5a profile (left) + genomic context (right) on one page
  insets <- lapply(ctx_evs, function(ev) {
    p <- event_profile(ev)
    if (is.null(p)) return(NULL)
    list(plot = with_base(cs$inset_base, p + poster_theme() + labs(title = NULL, subtitle = NULL)),
         w = cs$inset_size[1], h = cs$inset_size[2])
  })
  render_panels(ctx_panels, plot_dir, out_t("21_event_panels.csv"),
                pout("04_event_panels_regions.gtf"), folder = "04_event_panels", insets = insets)
  # 5b genomic context alone
  render_panels(ctx_panels, plot_dir, out_t("22_genomic_context.csv"),
                pout("05_genomic_context_regions.gtf"), folder = "05_genomic_context")
}

## ---- 6. Cis pairs ---------------------------------------------------------------
log_step("Cis pairs")
cis_raw <- read_csv(resolve(cs$cis_pairs_file), show_col_types = FALSE)
find_tx <- function(x) {
  coalesce(match(x, det_final$transcript_name), match(bare(x), bare(det_final$transcript_id)))
}
cis <- cis_raw |>
  mutate(pair_no = row_number()) |>               # numbering of your pair table
  filter(partner_category %in% cs$cis_keep, !as.logical(partner_novel)) |>
  mutate(
    n_i = find_tx(novel), p_i = find_tx(partner),
    novel_id       = det_final$transcript_id[n_i], partner_id     = det_final$transcript_id[p_i],
    novel_type     = det_final$transcript_type[n_i], partner_type = det_final$transcript_type[p_i],
    novel_log2FC   = det_final$log2FC[n_i],        partner_log2FC = det_final$log2FC[p_i],
    novel_is_de    = det_final$is_de[n_i],         partner_is_de  = det_final$is_de[p_i],
    # Bambu's gene assignment, from the DE tables (not from the pair file)
    novel_gene     = det_final$gene_id[n_i],       partner_gene_b = det_final$gene_id[p_i],
    dist_txt = ifelse(distance_bp == 0, "overlapping",
                      paste(prettyNum(round(distance_bp / 1000), big.mark = ","), "kb"))
  )
miss <- cis |> filter(is.na(novel_id) | is.na(partner_id))
if (nrow(miss)) log_step("  ! not found in the DE tables (skipped): %s", paste(miss$pair_lab, collapse = "; "))
cis <- cis |> filter(!is.na(novel_id), !is.na(partner_id))
log_step("  %d pairs (partner annotated protein-coding or lncRNA)", nrow(cis))
write_csv(cis |> select(-n_i, -p_i), out_t("23_cis_pairs_used.csv"))

if (nrow(cis)) {
  # 6a per-sample z-score strip: rows = novel / partner, columns = samples
  # the biotype palette of make_palette() above
  cat_cols <- c("Annotated protein-coding" = "#A6D854", "Annotated lncRNA" = "#E78AC3",
                "Annotated other" = "#C9C9C9", "Novel protein-coding" = "#238B45",
                "Novel lncRNA" = "#AE017E", "Novel non-coding" = "#6A51A3")
  hm_cols  <- c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2", "#F2CBB7", "#EE8468", "#B40426")
  samp  <- sheet$sample[order(match(sheet$group, c(cfg$control, cfg$case)), sheet$sample)]
  sgrp  <- sheet$group[match(samp, sheet$sample)]
  xpos  <- seq_along(samp) + ifelse(sgrp == cfg$case, 0.3, 0)   # small gap between lines
  nm    <- norm_tx$norm
  zrow  <- function(id) {
    v <- log2(nm[id, samp] + 1)
    s <- stats::sd(v)
    if (!is.finite(s) || s == 0) v * 0 else (v - mean(v)) / s
  }
  rows <- purrr::map_dfr(seq_len(nrow(cis)), function(i) {
    r <- cis[i, ]
    tibble::tibble(pair_no = r$pair_no, name = c(r$novel, r$partner),
                   id = c(r$novel_id, r$partner_id),
                   category = c(r$novel_category, r$partner_category),
                   y = -(i - 1) * 2.6 - c(0, 1))
  })
  tiles <- purrr::map_dfr(seq_len(nrow(rows)), function(j) {
    tibble::tibble(x = xpos, y = rows$y[j], z = as.numeric(zrow(rows$id[j])))
  })
  mids <- rows |> group_by(pair_no) |> summarise(y = mean(y), .groups = "drop") |>
    left_join(cis |> select(pair_no, rho, novel_chr, dist_txt), by = "pair_no")
  lim  <- max(abs(tiles$z), na.rm = TRUE)
  xr   <- max(xpos) + 0.8
  ctl  <- xpos[sgrp == cfg$control]
  cse  <- xpos[sgrp == cfg$case]

  # ts scales every text and key size (1 = poster); section 7 draws it at 6 pt
  cis_heatmap <- function(ts = 1) {
    ggplot() +
      geom_tile(data = tiles, aes(x, y, fill = z), width = 0.95, height = 0.92) +
      geom_point(data = rows, aes(x = 0.25, y = y, colour = category), shape = 15, size = 3.2 * ts) +
      geom_text(data = rows, aes(x = -0.15, y = y, label = name), hjust = 1, size = 3.2 * ts) +
      geom_text(data = mids, aes(x = -4.0, y = y, label = sprintf("%02d", pair_no)),
                size = 3.6 * ts) +
      geom_text(data = mids, aes(x = xr, y = y + 0.28, label = sprintf("\u03c1 = %+.2f", rho)),
                hjust = 0, size = 3.2 * ts) +
      geom_text(data = mids, aes(x = xr, y = y - 0.3, label = sprintf("%s \u00b7 %s", novel_chr, dist_txt)),
                hjust = 0, size = max(2.9 * ts, 5 / ggplot2::.pt), colour = "grey40") +
      annotate("text", x = mean(ctl), y = 1.25, label = sprintf("HCC827 (%s)", cfg$control),
               colour = unname(grp_col[cfg$control]), size = 3.3 * ts) +
      annotate("text", x = mean(cse), y = 1.25, label = sprintf("H1975 (%s)", cfg$case),
               colour = unname(grp_col[cfg$case]), size = 3.3 * ts) +
      scale_fill_gradientn(colours = hm_cols, limits = c(-lim, lim), name = "row z-score") +
      scale_colour_manual(values = cat_cols, name = NULL) +
      scale_x_continuous(limits = c(-4.7, xr + 2.6), expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(add = c(0.8, 0.4))) +
      coord_cartesian(clip = "off") +
      labs(title = "Top cis pairs: novel transcript vs annotated partner",
           subtitle = "row z-score of log2(normalized count + 1) per sample",
           caption = "With 3 vs 3 samples |\u03c1| = 1 whenever both transcripts separate the two lines.") +
      theme_void(base_size = 11 * ts, base_family = theme_font) +
      theme(legend.position = "right", plot.margin = margin(8, 8, 8, 8) * ts)
  }
  p_hm <- cis_heatmap()

  # Ext. Data version: the poster layout (pair number, names + biotype square,
  # per-sample tiles, rho above "chr . distance"), only the tiles narrower.
  # Built in millimetres (1 data unit = 1 mm, coord_fixed) so text and tiles
  # keep their print size; biotype squares are exactly as tall as the tiles.
  # Returns the plot and the page size (mm).
  # (the Fig. 2q heatmap is drawn by 04_cis_circos.R) Extra arguments: scl scales every fixed mm
  # distance, case_gap = mm between the two lines' tiles (before scl), side_pt =
  # size of the rho / location text, fill_name = colour-bar title, line_col =
  # colour of the line names (NULL = group colours), cat_legend = biotype key
  cis_heatmap_paper <- function(pt = 6, tile_w = 2.2, row_h = 2.6, pair_gap = 1.6, scl = 1,
                                case_gap = 1.0, side_pt = pt, fill_name = "row z-score",
                                line_col = NULL, cat_legend = TRUE) {
    txt  <- pt / ggplot2::.pt                         # geom_text size for pt
    stxt <- side_pt / ggplot2::.pt
    lcol <- function(g) if (is.null(line_col)) unname(grp_col[g]) else line_col
    chw  <- pt * 0.17                                 # ~ one character at pt (mm)
    n    <- nrow(cis)
    rr   <- rows |> mutate(k = rep(seq_len(n), each = 2), pos = rep(c(0, 1), n),
                           y = -((k - 1) * (2 * row_h + pair_gap) + pos * row_h))
    pr   <- rr |> filter(pos == 0) |>
      transmute(k, y, pair_no = cis$pair_no, rho = cis$rho,
                loc = sprintf("%s \u00b7 %s", cis$novel_chr, cis$dist_txt))
    num  <- sprintf("%d", pr$pair_no)                   # as in "cis-acting pair 1"
    x_num <- max(nchar(num)) * chw                     # right edge of the pair numbers
    x_nm  <- x_num + 2.5 * scl + max(nchar(rr$name)) * chw   # right edge of the names
    cat_w <- 1.6 * scl
    x_cat <- x_nm + 1.0 * scl + cat_w / 2                # biotype square centre
    x0    <- x_cat + cat_w / 2 + 1.0 * scl               # left edge of the first sample tile
    xs    <- x0 + (seq_along(samp) - 0.5) * tile_w + ifelse(sgrp == cfg$case, case_gap * scl, 0)
    x_txt <- max(xs) + tile_w / 2 + 1.5 * scl
    x_end <- x_txt + max(nchar(c(sprintf("\u03c1 = %+.2f", pr$rho), pr$loc))) * side_pt * 0.17 + 1 * scl
    y_top <- row_h * 1.6                               # room for the line names
    y_bot <- min(rr$y) - row_h / 2 - 0.3 * scl
    tl <- data.frame(x = rep(xs, times = nrow(rr)), y = rep(rr$y, each = length(samp)), z = tiles$z)
    p <- ggplot() +
      geom_tile(data = tl, aes(x, y, fill = z), width = 0.94 * tile_w, height = 0.92 * row_h) +
      geom_tile(data = rr, aes(x = x_cat, y = y), width = cat_w, height = 0.92 * row_h,
                fill = unname(cat_cols[rr$category])) +                     # same height as the tiles
      geom_point(data = rr, aes(x = x_cat, y = y, colour = category), alpha = 0, size = 0) +  # legend only
      geom_text(data = pr, aes(x = x_num, y = y - row_h / 2, label = num), hjust = 1,
                size = txt, family = theme_font) +
      geom_text(data = rr, aes(x = x_nm, y = y, label = name), hjust = 1, size = txt, family = theme_font) +
      geom_text(data = pr, aes(x = x_txt, y = y, label = sprintf("\u03c1 = %+.2f", rho)), hjust = 0,
                size = stxt, family = theme_font) +
      geom_text(data = pr, aes(x = x_txt, y = y - row_h, label = loc), hjust = 0, size = stxt,
                colour = "grey40", family = theme_font) +
      annotate("text", x = mean(xs[sgrp == cfg$control]), y = row_h * 1.05, label = "HCC827",
               size = txt, family = theme_font, colour = lcol(cfg$control)) +
      annotate("text", x = mean(xs[sgrp == cfg$case]), y = row_h * 1.05, label = "H1975",
               size = txt, family = theme_font, colour = lcol(cfg$case)) +
      scale_fill_gradientn(colours = hm_cols, limits = c(-lim, lim), name = fill_name) +
      scale_colour_manual(values = cat_cols, name = NULL) +
      guides(fill = guide_colourbar(order = 1, barwidth = unit(2 * scl, "mm"), barheight = unit(14 * scl, "mm")),
             colour = if (cat_legend) guide_legend(order = 2, override.aes = list(alpha = 1, size = 2.2, shape = 15))
                      else "none") +
      coord_fixed(ratio = 1, xlim = c(0, x_end), ylim = c(y_bot, y_top), expand = FALSE, clip = "off") +
      theme_void(base_size = pt) +
      theme(text = element_text(family = theme_font, face = "plain"),
            legend.text = element_text(size = pt), legend.title = element_text(size = pt),
            legend.key.size = unit(2.5 * scl, "mm"), legend.position = "right",
            plot.margin = margin(1, 1, 1, 1, "mm"))
    # height = the data range + margins, so 1 data unit is exactly 1 mm; the
    # width adds the legend
    list(plot = p, w = ceiling(x_end + 40 * scl), h = ceiling(y_top - y_bot + 2))
  }
  with_base(11, save_pages(list(p_hm), pout("06_cis_pairs_heatmap.pdf"),
                           width = 8.5, height = nrow(cis) * 2.6 * 0.2 + 1.8,
                           page_names = "cis_pairs_heatmap"))

  # 6b per pair: expression profile of the novel transcript, its annotated
  # partner and the partner's gene - the same isoformic profile as the events
  # (replaces the log2FC bars; the profile carries both expression and change)
  cis_profile_df <- function(r) {
    pn <- det_final$transcript_name[r$p_i]            # names as in profile_df
    nn <- det_final$transcript_name[r$n_i]
    pg <- profile_df$parent_gene[match(pn, profile_df$genename)]
    if (is.na(pg)) stop("partner ", pn, " not in the profile data")
    base <- profile_df |> filter(parent_gene == pg, transcript_type == "gene" | genename == pn)
    add  <- profile_df |> filter(genename == nn, transcript_type != "gene")
    if (!nrow(add)) stop("novel transcript ", nn, " not in the profile data")
    # drawn under the partner's gene; the "_ENSG..." suffix edgeR added to
    # duplicated symbols is dropped from the shown gene name (LSP1P5_ENSG... -> LSP1P5)
    gl <- sub("_ENS[GT][0-9]+(\\.[0-9]+)?$", "", pg)
    df <- bind_rows(base, add)
    df$genename[df$transcript_type == "gene"] <- gl
    df$parent_gene <- gl
    list(df = df, gene = gl, keep = c(pn, nn))
  }
  ## ed = TRUE (Extended Data Fig. 6): title = the partner gene; subtitle = the
  ## pair number (as in Fig. 2q), then chromosome, distance and rho
  cis_profile <- function(i, type_limits = NULL, title_wrap = NULL, ed = FALSE) {
    r <- cis[i, ]
    safe_plot({
      x   <- cis_profile_df(r)
      ttl <- if (ed) x$gene else sprintf("cis-acting pair %d", r$pair_no)   # transcripts are named on the lines
      if (!is.null(title_wrap)) ttl <- stringr::str_wrap(ttl, title_wrap)
      loc <- sprintf("%s \u00b7 %s \u00b7 \u03c1 = %+.2f", r$novel_chr, r$dist_txt, r$rho)
      plot_gene_profile(x$df, x$gene, pal, keep = x$keep, de_only = FALSE, title = ttl,
                        subtitle = if (ed) sprintf("cis-acting pair %d\n%s", r$pair_no, loc) else loc,
                        type_limits = type_limits)
    }, r$pair_lab)
  }
  # every transcript type across the pairs: one shared legend in panels
  cis_types <- unique(unlist(lapply(seq_len(nrow(cis)), function(i) {
    tryCatch(as.character(cis_profile_df(cis[i, ])$df$transcript_type), error = function(e) NULL)
  })))
  cis_types <- c(intersect("gene", cis_types), sort(setdiff(cis_types, "gene")))
  save_pages(lapply(seq_len(nrow(cis)), cis_profile), pout("07_cis_pairs_profiles.pdf"),
             kind = "profile",
             page_names = sprintf("%02d_%s_%s", cis$pair_no, cis$novel, cis$partner))
}

log_step("Consolidated figures done: %s", plot_dir)

## ---- 7. Journal panels (final print size; Nature: text 5-7 pt) --------------
# Drawn at the size they are printed, so nothing is rescaled at submission:
#   Loci         DE isoform profile (left) + genomic context of cs$paper_loci,
#                together 171.5 mm wide (the 174 mm figure row minus 2.5 mm
#                for the letter); height follows the tracks
#                and the profile fills it. Locus labels: name (+/-) + class.
#   Ext. Data    profiles of every event, profiles of every cis pair and the
#                cis-pair heatmap, 180 mm wide; .rds = the ggplot, for a composite
# Output: isoformic/figures/
#   locus_{PDF,PNG,SVG}/<gene>.*  + locus_panels.csv (windows) + locus_regions.gtf
#   ed_profiles_events_cis.* (or ed_isoform_profiles.* + ed_cis_pairs_profiles.*
#   with paper_profiles_combined = FALSE), ed_cis_pairs_heatmap.*
#   (each .pdf, .svg, .png, .rds)
if (isTRUE(cs$run_paper)) {
  log_step("Journal panels")
  paper_dir <- resolve(cs$paper_dir)
  dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
  ppt <- cs$paper_font

  # replace several cfg values while expr runs (restored on exit, also on error)
  with_cfg <- function(vals, expr) {
    old <- lapply(names(vals), function(k) cfg[[k]])
    names(old) <- names(vals)
    for (k in names(vals)) cfg[[k]] <<- vals[[k]]
    on.exit(for (k in names(old)) cfg[[k]] <<- old[[k]], add = TRUE)
    force(expr)
  }
  # settings shared by every journal profile: 6 pt, repel labels ppt - 0.5 pt,
  # HCC827 / H1975 on the x axis (Ext. Data keeps the poster isoform cap)
  paper_cfg <- list(base_size = ppt, legend_position = "bottom",
                    lab_mm = (ppt - 0.5) / ggplot2::.pt / 0.9,
                    ctx_group_labels = cs$paper_groups,
                    profile_compact = TRUE,      # labels close to their lines
                    profile_cap_note = FALSE,    # no "10/19 isoforms" in subtitles
                    profile_label_room = 1.5)    # label space left of HCC827 (category widths)
  # every text element set explicitly: theme_bw's rel(0.8) would put axis text
  # below 5 pt at a 6 pt base
  paper_theme <- function(pt, axes = TRUE, strips = TRUE) {
    th <- theme(
      text            = element_text(size = pt, family = theme_font, face = "plain"),
      plot.title      = element_text(size = pt + 1, face = "plain"),
      plot.subtitle   = element_text(size = pt, colour = "grey30"),
      plot.caption    = element_text(size = pt, colour = "grey30", hjust = 0),
      strip.text      = element_text(size = pt),
      legend.text     = element_text(size = pt),
      legend.title    = element_text(size = pt),
      legend.key.size = unit(3, "mm"),
      legend.margin   = margin(0, 0, 0, 0)
    )
    if (!strips) th <- th + theme(strip.text = element_blank(), strip.background = element_blank())
    if (!axes) return(th)
    th + theme(axis.text = element_text(size = pt), axis.title = element_text(size = pt),
               # boxed panel as the other figure panels (border 0.4 at 7 pt)
               panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4 * pt / 7),
               axis.line = element_blank(),
               axis.ticks = element_line(linewidth = 0.3 * pt / 7, colour = "black"),
               axis.ticks.length = unit(0.8, "mm"))
  }
  # isoformic draws poster-sized points and lines; set them for print size
  paper_geoms <- function(p, point = 1.2, line = 0.55) {
    for (i in seq_along(p$layers)) {
      g <- p$layers[[i]]$geom
      if (inherits(g, c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel"))) next
      if (inherits(g, "GeomPoint")) p$layers[[i]]$aes_params$size <- point
      if (inherits(g, c("GeomPath", "GeomSegment"))) p$layers[[i]]$aes_params$linewidth <- line
      if (inherits(g, c("GeomErrorbar", "GeomLinerange"))) p$layers[[i]]$aes_params$linewidth <- 0.6 * line
    }
    p
  }
  # n profiles in ncol columns, one shared legend at the bottom
  # the shared legend goes into the empty cells of the last row when there are
  # any (guide_area), otherwise below the grid; title above the keys
  profile_grid <- function(pl, ncol) {
    if (length(pl) > 1) {
      pl[-1] <- lapply(pl[-1], function(p) p + guides(colour = "none", fill = "none", alpha = "none"))
    }
    spare <- ceiling(length(pl) / ncol) * ncol - length(pl)
    if (spare > 0) {
      pl <- c(pl, list(patchwork::guide_area()))
    } else {                                   # horizontal legend: types in 2 rows, DE in 1
      horiz <- (cs$paper_legend_pos %||% "right") %in% c("top", "bottom")
      pl[[1]] <- pl[[1]] + if (horiz) {
        guides(colour = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE),
               alpha = guide_legend(nrow = 1))
      } else {
        guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1), alpha = guide_legend(ncol = 1))
      }
    }
    lay <- if (utils::packageVersion("patchwork") >= "1.2.0") {
      patchwork::plot_layout(ncol = ncol, guides = "collect", axis_titles = "collect")
    } else {
      patchwork::plot_layout(ncol = ncol, guides = "collect")
    }
    in_cell <- spare > 0
    pos     <- if (in_cell) "right" else cs$paper_legend_pos %||% "right"
    horiz   <- pos %in% c("top", "bottom")
    title_pos <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
      theme(legend.title.position = if (horiz) "left" else "top")
    } else theme()
    (patchwork::wrap_plots(pl) + lay) & paper_theme(ppt, strips = FALSE) &   # title names the gene
      theme(legend.position  = pos,
            legend.direction = if (horiz) "horizontal" else "vertical",
            legend.box       = if (horiz) "horizontal" else "vertical",     # right: one column
            legend.box.just  = if (horiz) "top" else "left",
            legend.justification = if (horiz) "center" else if (in_cell) c(0, 0.5) else c(0, 1),  # right: top
            legend.spacing.y = unit(2, "mm"),
            legend.spacing.x = unit(1.5, "mm"),
            legend.key.size  = unit(2.5, "mm")) & title_pos
  }
  # PDF (cairo, keeps unicode) + SVG + PNG + the ggplot itself (.rds)
  save_paper <- function(p, stem, w_mm, h_mm, rds_stem = stem) {  # rds_stem: .rds kept apart
    w <- w_mm / 25.4; h <- h_mm / 25.4
    dev_pdf <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"   # pdf() has no filename arg
    ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = dev_pdf)
    ggsave_svg(paste0(stem, ".svg"), p, width = w, height = h, bg = "white")   # editable text
    ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = cs$paper_dpi, bg = "white")
    saveRDS(p, paste0(rds_stem, ".rds"))
    log_step("  wrote %s.{pdf,svg,png,rds} (%.0f x %.0f mm)", basename(stem), w_mm, h_mm)
  }
  # showtext draws raster text at its own dpi: match the PNG resolution
  st_dpi <- function(d) if ("showtext" %in% loadedNamespaces()) showtext::showtext_opts(dpi = d)
  st_dpi(cs$paper_dpi)

  tryCatch({
    # 7a. Loci - genomic context with one coverage track per library, log2(depth
    # + 1) axis, the novel models' spans shaded:
    #   locus_tracks_*/  every event (cs$coverage_loci = "all"): tracks only,
    #                    cs$paper_locus_tracks_mm wide
    #   locus_samples_*/ cs$paper_loci: the same with the DE isoform profile on the left
    cov_names <- if (identical(cs$coverage_loci, "all")) names(evs) else intersect(cs$coverage_loci, names(evs))
    if (isTRUE(cs$profiles_only)) cov_names <- character(0)   # profiles-only rerun
    miss <- setdiff(cs$paper_loci, names(evs))
    if (length(miss)) log_step("  ! paper loci not among the events (skipped): %s", paste(miss, collapse = ", "))
    if (length(cov_names)) {
      ploc <- cs$paper_locus
      ploc$width_mm <- ploc$width_mm - cs$paper_inset_mm        # tracks = total - profile
      with_cfg(modifyList(paper_cfg, list(ctx_paper = ploc, ctx_label_style = "strand",
                                 poster_max_iso = 50,    # profile = every model the locus draws
                                 lab_mm = 5.5 / ggplot2::.pt / 0.9, # isoform names 5.5 pt
                                 profile_label_room = 1.2,          # ~13 mm for the names
                                 context_max_tx = ploc$max_tx %||% cfg$context_max_tx,
                                 dpi = cs$paper_dpi,
                                 ctx_coverage = "samples", ctx_cov_scale = "log2",
                                 ctx_cov_highlight = TRUE)), {
        cov_panels <- lapply(evs[cov_names], function(ev) tryCatch(
          build_panel(ev$gene_id, "paper", ev$gene, force_ids = c(ev$extras$transcript_id, ev$partner_ids),
                      note = ev$row$axis),
          error = function(e) { log_step("  ! locus %s: %s", ev$gene, conditionMessage(e)); NULL }))
        cov_panels <- cov_panels[!vapply(cov_panels, is.null, logical(1))]
        log_step("  coverage panels (per library): %d events", length(cov_panels))
        with_cfg(list(ctx_paper = modifyList(ploc, list(sig_h = 0.15, width_mm = cs$paper_locus_tracks_mm,
                                                         ## pipeline-style header: "GENE | chr:start-end",
                                                         ## two grey lines, ideogram + zoom lines
                                                         header = TRUE, ideogram = TRUE))), {
          render_panels(cov_panels, paper_dir, file.path(paper_dir, "locus_tracks_panels.csv"),
                        file.path(paper_dir, "locus_tracks_regions.gtf"), folder = "locus_tracks",
                        insets = NULL, formats = c("PDF", "PNG", "SVG"))
        })
        ## the profiles are their own panels (no profile inset next to the locus)
      })
    }

    # 7b/c. Ext. Data - profiles of every event (cs$panel_all) and every cis
    # pair: one shared legend (transcript types of both sets). With
    # cs$paper_profiles_combined = TRUE they fill one grid (events first, then
    # cis pairs; 24 profiles = 4 x 6, legend below); FALSE = two figures.
    # EVERY drawn event profile (main, retired, lineage markers, review
    # candidates), not only the consolidated-panel set; their transcript types
    # join the shared legend
    ed_evs <- if (isTRUE(cs$paper_ed_all_events)) evs else panel_evs
    ed_types <- unlist(lapply(ed_evs, function(ev)
      tryCatch(as.character(event_profile_df(ev)$transcript_type), error = function(e) NULL)))
    all_types <- unique(c(types, ed_types, if (exists("cis_types")) cis_types))
    all_types <- c(intersect("gene", all_types), sort(setdiff(all_types, "gene")))
    eds <- if (length(ed_evs)) with_cfg(paper_cfg, {
      ## title = the gene(s) (text before " · " of the event title); subtitle =
      ## the description, then the biological axis, one line each
      ed_title <- function(ev) stringr::str_wrap(sub(" · .*$", "", ev$row$title), cs$paper_ed_title_wrap)
      ed_sub   <- function(ev) {
        d <- if (grepl(" · ", ev$row$title)) sub("^.*? · ", "", ev$row$title, perl = TRUE) else NULL
        paste(c(d, ev$row$axis), collapse = "\n")
      }
      x <- lapply(ed_evs, function(ev) event_profile(ev, title = ed_title(ev), subtitle = ed_sub(ev),
                                                     type_limits = all_types))
      x <- x[!vapply(x, is.null, logical(1))]
      lapply(x, paper_geoms, point = 0.8, line = 0.35)
    }) else list()
    cps <- if (exists("cis_profile") && nrow(cis)) with_cfg(paper_cfg, {
      x <- lapply(seq_len(nrow(cis)), function(i) {
        cis_profile(i, type_limits = all_types, title_wrap = cs$paper_ed_title_wrap, ed = TRUE)
      })
      x <- x[!vapply(x, is.null, logical(1))]
      lapply(x, paper_geoms, point = 0.8, line = 0.35)
    }) else list()
    # page height: rows x row height (+ a bottom strip when no cell is free for the legend)
    grid_h <- function(n, ncol, row_h) {
      spare <- ceiling(n / ncol) * ncol - n
      side  <- spare > 0 || identical(cs$paper_legend_pos %||% "right", "right")
      ceiling(n / ncol) * row_h + if (side) 4 else 18
    }
    grid_w <- function(n, ncol) {
      spare <- ceiling(n / ncol) * ncol - n
      cs$paper_ed_width +
        if (spare == 0 && identical(cs$paper_legend_pos %||% "right", "right")) cs$paper_legend_w else 0
    }
    ## one standard plot box for EVERY profile (each facet of the
    ## two-gene profiles, UBE2I / UBE2I-AS1 and RAB11B-AS1 / RAB11B, too), set by
    ## what fits on an A3 portrait page: the event profiles on one page
    ## ("Extended data figure 5": cs$paper_page_ncol columns, a two-gene profile
    ## takes two rows of its column, its second box on the next row's grid line,
    ## the shared legend in the free cells of the last row) and the cis-pair
    ## profiles on their own page ("Extended data figure 6", same box). The box
    ## (cs$paper_box_mm, written to paper_box_mm.rds) is also used by the Fig. 2
    ## profiles (7e). Boxes line up in rows and columns; one y-axis title per page.
    if (length(c(eds, cps))) {
      ## rooms around the box(es) of a ggplot, mm (left / right / top / bottom)
      grob_rooms <- function(g) {
        pl <- g$layout[grepl("^panel", g$layout$name), ]
        w  <- grid::convertWidth(g$widths, "mm", valueOnly = TRUE)
        h  <- grid::convertHeight(g$heights, "mm", valueOnly = TRUE)
        c(L = sum(w[seq_len(min(pl$l) - 1)]), R = sum(w[-seq_len(max(pl$r))]),
          T = sum(h[seq_len(min(pl$t) - 1)]), B = sum(h[-seq_len(max(pl$b))]),
          W = sum(w), H = sum(h))
      }
      n_facets <- function(p) length(unique(ggplot2::ggplot_build(p)$layout$layout$PANEL))
      cell_theme <- function(p) p + paper_theme(ppt, strips = FALSE) +
        theme(legend.position = "none", plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")) + labs(y = NULL)
      ## a two-gene profile shares ONE box: each facet half its height
      half_gap <- 1.5
      fix_box <- function(p, box, nf = 1) p +
        ggh4x::force_panelsizes(rows = unit(if (nf > 1) (box[2] - (nf - 1) * half_gap) / nf else box[2], "mm"),
                                cols = unit(box[1], "mm")) +
        theme(panel.spacing = unit(half_gap, "mm"))
      ## the shared legend (every transcript type: fixed limits), 2 columns of keys
      ## one column, as the Fig. 2 profile legend (7e): transcript types on top,
      ## each name on one line (no poster wrapping here: the column is free),
      ## then "Differentially expressed", then the "mean ± s.d." key of the
      ## triangles; sentence-case titles, more room between keys and blocks
      key_gap <- 0.9; block_gap <- 2.5                 # mm, between keys / between blocks
      legend_grob <- function(p) {
        guide_of <- function(g) {
          i <- which(grepl("guide-box", g$layout$name) & !vapply(g$grobs, inherits, logical(1), "zeroGrob"))
          g$grobs[[i[1]]]
        }
        leg_theme <- paper_theme(ppt, strips = FALSE) +
          theme(legend.position = "right", legend.direction = "vertical",
                legend.box = "vertical", legend.box.just = "left",
                legend.key.size = unit(2.5, "mm"),
                legend.key.spacing.y = unit(key_gap, "mm"),
                legend.spacing.y = unit(block_gap, "mm"))
        ## the type labels are a function of cfg$fig_style, evaluated here
        l1 <- with_cfg(list(fig_style = "paper"), guide_of(ggplotGrob(p + leg_theme +
                          labs(colour = "Transcript type", fill = "Transcript type",
                               alpha = "Differentially expressed") +
                          guides(colour = guide_legend(ncol = 1, order = 1), fill = guide_legend(ncol = 1, order = 1),
                                 alpha = guide_legend(ncol = 1, order = 2)))))
        l2 <- guide_of(ggplotGrob(
          ggplot(data.frame(x = 1, y = 1, k = "mean ± s.d."), aes(x, y, shape = k)) +
            geom_point(colour = "grey35", fill = "grey35", size = 1.1) +   # ED profile points are 0.8
            scale_shape_manual(values = 24, name = NULL) + leg_theme))
        w <- vapply(list(l1, l2), function(l) grid::convertWidth(sum(l$widths), "mm", valueOnly = TRUE), 0)
        h <- vapply(list(l1, l2), function(l) grid::convertHeight(sum(l$heights), "mm", valueOnly = TRUE), 0)
        l1 <- gtable::gtable_add_cols(l1, unit(max(w) - w[1], "mm"))
        l2 <- gtable::gtable_add_cols(l2, unit(max(w) - w[2], "mm"))
        gtable::gtable_col("ed_profile_legend", list(l1, grid::nullGrob(), l2), width = unit(max(w), "mm"),
                           heights = unit(c(h[1], block_gap, h[2]), "mm"))
      }
      ## grid placement: row-major, a two-facet profile takes the cell below too
      place <- function(nf, ncol, reserve = 0) {
        occ <- matrix(FALSE, 0, ncol); pos <- matrix(NA_integer_, length(nf), 2)
        k <- 1; r <- 1
        while (k <= length(nf)) {
          while (nrow(occ) < r + 1) occ <- rbind(occ, rep(FALSE, ncol))
          for (c in seq_len(ncol)) {
            if (k > length(nf)) break
            if (occ[r, c]) next
            if (nf[k] > 1) {
              while (nrow(occ) < r + nf[k] - 1) occ <- rbind(occ, rep(FALSE, ncol))
              if (any(occ[r:(r + nf[k] - 1), c])) next
              occ[r:(r + nf[k] - 1), c] <- TRUE
            } else occ[r, c] <- TRUE
            pos[k, ] <- c(r, c); k <- k + 1
          }
          r <- r + 1
        }
        used <- max(which(rowSums(occ) > 0))
        occ <- occ[seq_len(used), , drop = FALSE]
        ## free cells for the legend at the end of the last row; else one more row
        free_last <- sum(!occ[used, ])
        if (free_last < reserve) { occ <- rbind(occ, rep(FALSE, ncol)); used <- used + 1 }
        list(pos = pos, nrow = used, occ = occ)
      }
      mar <- 5; gx <- 2; gy <- 2; ylab_w <- 4
      page <- cs$paper_page_mm                        # A4 portrait
      ylab <- expression(log[2] * "(mean normalised count + 1)")
      ## measure every profile once (box size does not change the rooms)
      grDevices::cairo_pdf(tempfile(fileext = ".pdf"))
      ## figure data behind every profile, kept for sharing
      source(file.path(repo, "downstream/helpers/source_data.R"))
      save_source_data(c(setNames(eds, paste0("ED5_", names(eds))),                 # Extended Data Fig. 5
                         setNames(cps, sprintf("ED6_pair_%02d", seq_along(cps)))),     # Extended Data Fig. 6
                       file.path(paper_dir, "source_data_profiles.rds"))
      prep <- function(pl) lapply(pl, cell_theme)
      eds_c <- prep(eds); cps_c <- prep(cps)
      nf_e  <- vapply(eds_c, n_facets, integer(1)); nf_c <- vapply(cps_c, n_facets, integer(1))
      rooms <- sapply(seq_along(c(eds_c, cps_c)), function(i) grob_rooms(ggplotGrob(fix_box(c(eds_c, cps_c)[[i]], c(30, 30), c(nf_e, nf_c)[i])))[1:4])
      rownames(rooms) <- c("L", "R", "T", "B")
      leg   <- legend_grob(c(eds, cps)[[1]])
      leg_w <- grid::convertWidth(grid::grobWidth(leg), "mm", valueOnly = TRUE)
      leg_h <- grid::convertHeight(grid::grobHeight(leg), "mm", valueOnly = TRUE)
      grDevices::dev.off()
      rm_ <- apply(rooms, 1, max) + 0.3
      ## box that fits the events page: the largest w x h (not forced square)
      ## every profile takes one cell (two-gene loci split their box); the
      ## column count is the one that gives the largest box (min of w and h)
      fit <- function(ncol) {
        nr <- place(rep(1, length(nf_e)), ncol, reserve = 0)$nrow  # legend on the right, not in a row
        c(ncol = ncol, nrow = nr,
          bw = unname((page[1] - 2 * mar - ylab_w - ncol * gx - leg_w) / ncol - rm_["L"] - rm_["R"]),
          bh = unname((page[2] - 2 * mar - (nr - 1) * gy) / nr - rm_["T"] - rm_["B"]))
      }
      fits <- t(sapply(3:8, fit))
      best <- fits[which.max(pmin(fits[, "bw"], fits[, "bh"])), ]
      ncol <- best[["ncol"]]; lay_e <- list(nrow = best[["nrow"]]); bw <- best[["bw"]]; bh <- best[["bh"]]
      box <- unname(floor(2 * c(bw, bh)) / 2)        # mm, w x h: the largest that fits, rounded down to 0.5 mm
      if (!is.null(cs$paper_box_max_h)) box[2] <- min(box[2], cs$paper_box_max_h)
      if (!is.null(cs$paper_box_mm)) box <- cs$paper_box_mm
      saveRDS(box, file.path(paper_dir, "paper_box_mm.rds"))
      prof_box <- box                                 # for the Fig. 2 profiles (7e)
      log_step("  standard profile box: %.1f x %.1f mm (fits %d x %d cells on the page)",
               box[1], box[2], ncol, lay_e$nrow)
      cw <- rm_["L"] + box[1] + rm_["R"]; ch <- rm_["T"] + box[2] + rm_["B"]
      ## one page: profiles on the grid, the legend in a column on the right
      profile_page <- function(pl, nf, stem, legend = TRUE) {
        lay  <- place(rep(1, length(nf)), ncol, reserve = 0)
        gw   <- ylab_w + ncol * cw + (ncol - 1) * gx + if (legend) gx + leg_w else 0
        gh   <- lay$nrow * ch + (lay$nrow - 1) * gy
        pw   <- page[1]; ph <- min(page[2], gh + 2 * mar)        # height trimmed to the content
        x0   <- (pw - gw) / 2 + ylab_w; y0 <- (ph - gh) / 2   # grid centred on the page
        kids <- list(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)))
        grDevices::cairo_pdf(tempfile(fileext = ".pdf"))
        for (k in seq_along(pl)) {
          g  <- ggplotGrob(fix_box(pl[[k]], box, nf[k]))
          r  <- grob_rooms(g)
          cx <- x0 + (lay$pos[k, 2] - 1) * (cw + gx); cy <- y0 + (lay$pos[k, 1] - 1) * (ch + gy)
          vp <- grid::viewport(x = unit(cx + rm_["L"] - r["L"], "mm"),
                               y = unit(ph - (cy + rm_["T"] - r["T"]), "mm"),
                               width = unit(r["W"], "mm"), height = unit(r["H"], "mm"),
                               just = c("left", "top"))
          kids[[length(kids) + 1]] <- grid::grobTree(g, vp = vp)
        }
        grDevices::dev.off()
        if (legend) {                                  # right of the grid, level with the first box
          lx <- x0 + ncol * (cw + gx); ly <- y0 + rm_["T"]
          kids[[length(kids) + 1]] <- grid::grobTree(leg, vp = grid::viewport(
            x = unit(lx, "mm"), y = unit(ph - ly, "mm"), width = unit(leg_w, "mm"),
            height = unit(leg_h, "mm"), just = c("left", "top")))
        }
        kids[[length(kids) + 1]] <- grid::textGrob(ylab, x = unit(x0 - ylab_w + 1.5, "mm"),
          y = unit(ph - y0 - gh / 2, "mm"), rot = 90,
          gp = grid::gpar(fontsize = ppt, fontfamily = theme_font))
        list(grob = grid::gTree(children = do.call(grid::gList, kids)), w = pw, h = ph)
      }
      ## "Extended data figure 5": every event; "... 6": every cis pair
      ## .pdf/.svg/.png in <project>/figures with the other figures; .rds stays here
      ed_dir <- cs$ed_dir %||% paper_dir
      dir.create(ed_dir, recursive = TRUE, showWarnings = FALSE)
      if (length(eds_c)) {
        pg <- profile_page(eds_c, nf_e, "events")
        save_paper(pg$grob, file.path(ed_dir, "Extended data figure 5"), pg$w, pg$h,
                   rds_stem = file.path(paper_dir, "Extended data figure 5"))
      }
      if (length(cps_c)) {
        pg <- profile_page(cps_c, nf_c, "cis")
        save_paper(pg$grob, file.path(ed_dir, "Extended data figure 6"), pg$w, pg$h,
                   rds_stem = file.path(paper_dir, "Extended data figure 6"))
      }
      unlink(file.path(paper_dir, paste0(rep(c("Extended data figure 2", "Extended data figure 3", "Extended data figure 4"), each = 4), rep(c(".pdf", ".svg", ".png", ".rds"), 3))))  # files of the earlier numbering (this folder only)
      log_step("  %d event profiles (ED5) + %d cis-pair profiles (ED6), box %.1f mm",
               length(eds_c), length(cps_c), box[1])
    }

    # 7d. Ext. Data - cis pairs heatmap: small square tiles, page sized to fit
    if (exists("cis_heatmap_paper") && nrow(cis)) {
      hp <- cis_heatmap_paper(ppt, tile_w = cs$paper_hm_tile_mm, row_h = cs$paper_hm_row_mm)
      save_paper(hp$plot, file.path(paper_dir, "ed_cis_pairs_heatmap"), hp$w, hp$h)
    }

    # 7e. Fig. 2 panels as printed in the final Fig. 2 (A3 page): j, k, l, m =
    # single DE isoform profiles (titles, subtitles and axis text 8 pt, isoform
    # labels 5.5 pt), with one shared legend file (fig2_profile_legend) for the row.
    # Size = the panel's extent on the page + 1 mm.
    # (Fig. 2o, p, the loci, are locus_*/; Fig. 2r, the circos, and Fig. 2q, the
    # cis-pair heatmap: 04_cis_circos.R. File names keep the letters fig2h-l.)
    fig2_dir <- file.path(paper_dir, "fig2_panels")
    dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)
    fig2_prof <- tibble::tribble(
      ~gene,   ~file,                       ~title,                                       ~subtitle,                    ~w,   ~h,
      # sizes = extent on the A3 Figure 2 page + 1 mm
      "EGFR",  "fig2h_profile_EGFR",        "EGFR \u00b7 AS intronic lncRNA",             "Driver locus",               48.7, 46.9,
      # title shortened: "UBE2I / UBE2I-AS1 \u00b7 AS lncRNA" ran past the 48.7 mm panel
      "UBE2I", "fig2j_profile_UBE2I",       "UBE2I / UBE2I-AS1",                          "Sense-antisense switch",     48.7, 95.0,   # two stacked profiles (UBE2I above UBE2I-AS1)
      "RAB11B-AS1", "fig2k_profile_RAB11B-AS1", "RAB11B-AS1 / RAB11B",                "Novel intronic lncRNAs",     60.1, 95.0,   # two stacked profiles, as j
      "MUC1",  "fig2l_profile_MUC1",        "MUC1 \u00b7 isoform switch",                 "LUAD lineage identity",      48.7, 54.9,
      # LY6D in Fig. 2 (MUC1 still written)
      "LY6D",  "fig2m_profile_LY6D",        "LY6D \u00b7 novel non-coding gain",          "Novel model, known locus",   48.7, 54.9,
      # lineage plasticity markers: not in the printed Fig. 2, same style and size as l
      "ASCL1", "marker_profile_ASCL1",      "ASCL1 \u00b7 neuroendocrine marker",         "Lineage plasticity",         48.7, 54.9,
      "POU2F3","marker_profile_POU2F3",     "POU2F3 \u00b7 tuft-cell marker",             "Lineage plasticity",         48.7, 54.9,
      "HMGA2", "marker_profile_HMGA2",      "HMGA2 \u00b7 EMT marker",                    "Lineage plasticity",         48.7, 54.9
    )
    ## the Extended Data label size (paper_cfg: paper_font - 0.5 pt);
    ## titles, subtitles and axes 8 pt as the rest of Fig. 2. The rooms of each
    ## file go to fig2_profile_rooms.rds (06: volcano of row 3; 09: alignment)
    prof_rooms <- list()
    ## Fig. 2 row 3 (helpers/fig2_rows.R: fig2_row3_profiles): one shared legend after
    ## the four profiles (transcript type, DE, triangles = mean +/- s.d.); their boxes
    ## are narrowed so that the legend and one more gutter fit in the same row width
    row3_files  <- c("fig2h_profile_EGFR", "fig2j_profile_UBE2I", "fig2k_profile_RAB11B-AS1",
                     "fig2m_profile_LY6D")
    fig2_gut_mm <- 4                                       # the Fig. 2 gutter (helpers/fig2_rows.R)
    fig2_cfg    <- modifyList(paper_cfg, list(lab_mm = (ppt - 1.5) / ggplot2::.pt / 0.9))  # isoform names 1 pt smaller
    row3_genes  <- fig2_prof$gene[match(row3_files, fig2_prof$file)]
    leg_mm <- c(0, 0); box_cut <- 0
    with_cfg(fig2_cfg, {
      p3 <- Filter(Negate(is.null), lapply(row3_genes, function(g) if (!is.null(evs[[g]])) event_profile(evs[[g]])))
      if (length(p3)) {
        types <- unique(unlist(lapply(p3, function(p)
          ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$get_limits())))
        types <- names(pal)[names(pal) %in% types]       # the palette's order
        guide_box <- function(g) {
          k <- which(grepl("guide-box", g$layout$name) & !vapply(g$grobs, inherits, logical(1), "zeroGrob"))
          g$grobs[[k[1]]]
        }
        grDevices::cairo_pdf(tempfile(fileext = ".pdf"))  # text metrics for the measurements
        l1 <- guide_box(ggplotGrob(
          event_profile(evs[[row3_genes[1]]], type_limits = types) + paper_theme(8, strips = FALSE) +
            theme(legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
                  legend.box.just = "left", legend.key.size = unit(3, "mm")) +
            labs(colour = "Transcript type", fill = "Transcript type", alpha = "Differentially expressed") +
            guides(colour = guide_legend(ncol = 1, order = 1), fill = guide_legend(ncol = 1, order = 1),
                   alpha = guide_legend(ncol = 1, order = 2))))
        l2 <- guide_box(ggplotGrob(
          ggplot(data.frame(x = 1, y = 1, k = "mean ± s.d."), aes(x, y, shape = k)) +
            geom_point(colour = "grey35", fill = "grey35", size = 1.6) +
            scale_shape_manual(values = 24, name = NULL) + paper_theme(8) +
            theme(legend.position = "right", legend.key.size = unit(3, "mm"))))
        w  <- vapply(list(l1, l2), function(l) grid::convertWidth(sum(l$widths), "mm", valueOnly = TRUE), 0)
        h  <- vapply(list(l1, l2), function(l) grid::convertHeight(sum(l$heights), "mm", valueOnly = TRUE), 0)
        grDevices::dev.off()
        l2 <- gtable::gtable_add_cols(l2, unit(max(w) - w[2], "mm"))   # left-aligned under the first block
        l1 <- gtable::gtable_add_cols(l1, unit(max(w) - w[1], "mm"))
        leg <- gtable::gtable_col("fig2_profile_legend", list(l1, l2), width = unit(max(w), "mm"),
                                  heights = unit(h, "mm"))
        leg_mm  <- c(max(w) + 1, sum(h) + 1)
        box_cut <- (leg_mm[1] + fig2_gut_mm) / length(row3_files)
        save_paper(leg, file.path(fig2_dir, "fig2_profile_legend"), leg_mm[1], leg_mm[2])
        prof_rooms[["fig2_profile_legend"]] <- c(W = leg_mm[1], H = leg_mm[2], T = 0, B = 0, box_h = NA)
        log_step("  Fig. 2 row-3 legend %.1f x %.1f mm; row-3 profile boxes %.1f mm narrower",
                 leg_mm[1], leg_mm[2], box_cut)
      }
    })
    with_cfg(fig2_cfg, {
      for (i in seq_len(nrow(fig2_prof))) {
        r  <- fig2_prof[i, ]
        ev <- evs[[r$gene]]
        if (is.null(ev)) { log_step("  ! Fig. 2 profile: %s is not an event - skipped", r$gene); next }
        p <- event_profile(ev, title = r$title)
        if (is.null(p)) next
        p <- paper_geoms(p, point = 0.8, line = 0.35) + paper_theme(ppt, strips = FALSE) +
          theme(plot.title    = element_text(size = 8, face = "plain"),  # 8 pt as the rest of Fig. 2
                plot.subtitle = element_text(size = 8, colour = "grey30"),
                axis.text     = element_text(size = 8, colour = "black"),
                axis.title    = element_text(size = 8),
                legend.position = "none") +
          labs(subtitle = r$subtitle)
        ## the standard profile box (each facet of a two-gene profile
        ## too); the file is the plot's natural size around it (r$w, r$h unused)
        bx <- if (exists("prof_box")) prof_box else cs$paper_box_mm %||% c(30, 30)
        if (r$file %in% row3_files) bx <- c(bx[1] - box_cut, bx[2])   # room for the row-3 legend
        nf2 <- length(unique(ggplot2::ggplot_build(p)$layout$layout$PANEL))   # two-gene loci: half box each
        p  <- p + ggh4x::force_panelsizes(rows = unit(if (nf2 > 1) (bx[2] - (nf2 - 1) * 1.5) / nf2 else bx[2], "mm"),
                                          cols = unit(bx[1], "mm")) +
          theme(panel.spacing = unit(1.5, "mm"))
        grDevices::cairo_pdf(tempfile(fileext = ".pdf")); gg <- ggplotGrob(p)   # pdf(NULL) has no CID fonts
        wh <- c(sum(grid::convertWidth(gg$widths, "mm", valueOnly = TRUE)),
                sum(grid::convertHeight(gg$heights, "mm", valueOnly = TRUE)))
        pl <- gg$layout[grepl("^panel", gg$layout$name), ]
        hh <- grid::convertHeight(gg$heights, "mm", valueOnly = TRUE)
        prof_rooms[[r$file]] <- c(W = wh[1], H = wh[2], T = sum(hh[seq_len(min(pl$t) - 1)]),
                                  B = sum(hh[-seq_len(max(pl$b))]), box_h = bx[2])
        grDevices::dev.off()
        save_paper(p, file.path(fig2_dir, r$file), wh[1], wh[2])   # titles kept (no _no_titles copy)
      }
    })
    saveRDS(prof_rooms, file.path(fig2_dir, "fig2_profile_rooms.rds"))
    # Fig. 2q (cis-pair heatmap) is drawn by 04_cis_circos.R in the printed
    # layout (paper_cis_heatmap.*).
  }, finally = st_dpi(300))                        # back to the plot_theme.R default
  log_step("Journal panels done: %s", paper_dir)
}
