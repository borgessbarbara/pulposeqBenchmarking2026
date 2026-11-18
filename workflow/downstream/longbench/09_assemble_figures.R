## ===========================================================================
## 09_assemble_figures.R - put the single panels together on their pages
##
##   Fig. 1   340 x 485 mm: row 1 left empty for the pulposeq workflow (panel
##            a, drawn in the vector editor); b-k from 08_fig1_panels.R
##   Extended Data Fig. 4 a-h, quantifier comparison: 340 x 347 mm (Fig. 1's
##            width); a-h from 03_quantifier_plots.R. Only the layout
##            is saved here; dong/11 adds the Dong row (i-l) and writes the page
## Both: 5 mm margins, 4 mm gutters, no empty areas; the panels WITHOUT titles
## (what the subtitles said goes in the figure legends); panel letters 12 pt,
## plain, in the gutter above each panel's top-left corner.
##
## Needs (run first): 08_fig1_panels.R, 03_quantifier_plots.R
## Outputs (in <project>/figures/, next to the vector-editor files, which are
## never touched):
##   fig1_assembled.svg / .png                  (SVG: nested panels, editable text)
##   spike-ins/results/quantifier_comparison/ed1_longbench_layout.rds  (ED4 a-h, read by dong/11)
## Slot sizes here must match the sizes the panel scripts save; a mismatch is
## reported as a warning.
## Needs: magick
## ===========================================================================

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
source(file.path(helper_dir, "compose_figure.R"))   # compose_figure, layout_rows, check_page
out_dir <- file.path(project, "figures")
dir.create(out_dir, showWarnings = FALSE)

## ------------------------------------------------------------- Fig. 1 -----
## row 1 (137)  a 330 (empty: the workflow)
## row 2 (70)   b 165 | c 161
## row 3 (72)   d 100 | e 226
## row 4 (100)  f 117 | g 105 | h 100 x 48 above i 100 x 48
## row 5 (80)   j 160 | k 166
fig1_page <- c(w = 340, h = 485)
fig1 <- layout_rows(list(
  list(h = 137, cells = list(list("a", NA, 330))),
  list(h = 70,  cells = list(list("b", "fig1b_spikein_identification", 165),
                             list("c", "fig1c_length_bias", 161))),
  list(h = 72,  cells = list(list("d", "fig1d_ercc", 100),
                             list("e", "fig1e_sequin_isoform", 226))),
  list(h = 100, cells = list(list("f", "fig1f_quantifier", 117),
                             list("g", "fig1g_read_support", 105),
                             list(w = 100, stack = list(list("h", "fig1h_novel_length", 48),
                                                        list("i", "fig1i_uic_sharing", 48))))),
  list(h = 80,  cells = list(list("j", "fig1j_uic_sharing_by_length", 160),
                             list("k", "fig1k_transcript_ends", 166)))))
check_page(fig1, fig1_page)
print(fig1, row.names = FALSE)
compose_figure(fig1, fig1_page, paste0(panel_dir, "_no_titles"), file.path(out_dir, "fig1_assembled"))

## ------------------------------------- supplementary: quantifiers -------
## row 1 (85)  a 150 | b 176
## row 2 (80)  c 216 | d 110
## row 3 (80)  e 150 | f 176
## row 4 (80)  g 200 | h 126
supp_page <- c(w = 340, h = 347)
supp <- layout_rows(list(
  list(h = 85, cells = list(list("a", "supp_q_a1_spikein_metrics", 150),
                            list("b", "supp_q_a3_spikein_quant_metrics", 176))),
  list(h = 80, cells = list(list("c", "supp_q_c2_human_counts", 216),
                            list("d", "supp_q_c1_jaccard", 110))),
  list(h = 80, cells = list(list("e", "supp_q_b1_divided_vs_undivided", 150),
                            list("f", "supp_q_b2_overdispersion", 176))),
  list(h = 80, cells = list(list("g", "supp_q_b3_overdispersion_isoforms", 200),
                            list("h", "supp_q_c3_disagreement", 126)))))
check_page(supp, supp_page)
## these are ED4 a-h; the Dong panels
## (i-l) go below them. dong/11_assemble_figures.R writes
## "Extended data figure 4" from this layout, so rerun dong/11 after this script.
saveRDS(list(lay = supp, page = supp_page, dir = file.path(lb_dir, "spike-ins/figures/quantifier_no_titles")),
        file.path(lb_dir, "spike-ins/results/quantifier_comparison/ed1_longbench_layout.rds"))
message("Done: figures in ", out_dir, " (Extended Data Fig. 4 a-h layout saved for dong/11)")
