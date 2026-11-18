## ===========================================================================
## 08_fig2_panels.R - Fig. 2 a, b, h, i, n: the Dong panels drawn with ggplot2
## H1975 (TKI-resistant) vs HCC827 (TKI-sensitive), ONT PCR-cDNA, Bambu full
## counts, 3 independent cultures per line (Dong et al. 2023 Nat Methods).
##
## Panel  object  content
##   a     pA      read-level support across biotypes
##   b     pB      novel transcripts by biotype and GffCompare class code
##   -     pID     sequins, Bambu vs oarfish (same run): identification recall / F1
##                 and isoform DE sensitivity / precision / F1 / fold-change R2
##                 (undivided counts), values printed; validation note in the subtitle
##   -     pC      sequin ground truth, isoform level: observed vs expected log2FC,
##                 DE precision / sensitivity / F1 top-left, sequins tested +
##                 R2 / slope bottom-right
##   h     pD      sample-sample Spearman correlation, isoform level (edgeR logCPM)
##   i     pE      volcano, gene level, Up / Down totals; marker genes labelled
##                 (only those tested AND DE)
##   n     pF      DE isoforms of the marker genes (Bambu gene assignment), rows
##                 grouped by gene and ranked by log2FC; z-score of logCPM per
##                 isoform; last column = isoform DE call
##   j-m           DE isoform profiles: 06_isoformic.R, section 7e
##                 (isoformic/figures/fig2_panels/)
##   o, p          loci: 06_isoformic.R (isoformic/figures/locus_*/)
##   r             cis-acting pairs circos: 04_cis_circos.R (paper_cis_circos.*)
##   q             cis-pair heatmap: 04_cis_circos.R (paper_cis_heatmap.*)
##   c-g           quantifier and sequin panels: 10_quantifier_panels.R; pID and
##                 pC (files fig2c_*, fig2d_*) are not in the assembled Fig. 2
##   ED Figs. 5, 6 every event + cis-pair profile: 06_isoformic.R
##
## Layout = the final Fig. 2 (A3 page, 297 x 420 mm):
## every panel is saved at the size it has on that page, with the font sizes,
## legend positions, key sizes and panel-border widths measured there
## (titles / subtitles 10 pt, legend and tick text 8 pt; see each panel), and the
## orientation as printed (vertical heatmap in n). No bold text. Theme base
## size 6.
##
## Needs (run first): 01_de_edger.R, 02_sequin_edger.R (oarfish and bambu: panel c),
##   07_read_support.R
## Outputs
##   figures/fig2_panels/fig2{a..g}_*.{pdf,svg,png}   (svg = editable text)
##   results/figure_tables/supp_novel_biotype_classcode_dong.csv   (Fig. 2b numbers)
##   results/figure_tables/sequin_identification_dong.csv           (pID numbers)
## Needs: tidyverse, ggrepel, RColorBrewer, svglite, showtext, systemfonts
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")
source(file.path(helper_dir, "plot_theme.R"))       # font_family, base_theme, group_cols, npg_cols
source(file.path(helper_dir, "figure_helpers.R"))   # save_panel, keep_level, read_support_plot
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------- settings ---
alpha    <- 0.05                  # as in 01_de_edger.R / 02_sequin_edger.R
lfc      <- 0.5
line_lab <- c(TKI_sensitive = "HCC827", TKI_resistant = "H1975")
lv_lab   <- c(gene = "gene", transcript = "isoform")
markers  <- c(
  # EGFR driver locus
  "EGFR",
  # LUAD lineage identity
  "NAPSA", "CEACAM5", "KRT7", "NKX2-1", "SFTPB", "MUC1", "SLC2A1", "TOP2A", "TPX2",
  "MALAT1", "CYTOR", "MIR4435-2HG", "DARS1-AS1", "HAGLR",
  # OXPHOS isoform switching (NDUFA11; NDUFA12 / COX6B1 / NDUFA6 are not
  # switches with oarfish counts)
  "NDUFA11",
  # hypoxia (VHL-HIF axis)
  "ELOC", "MIR210HG",
  # lncRNA regulation
  "UBE2I", "UBE2I-AS1", "SNHG17",
  # novel intronic lncRNAs in a known locus (Fig. 2l); sense gene RAB11B
  "RAB11B-AS1", "RAB11B",
  # translation / survival
  "TPT1",
  # cis-pair partner (with BambuTx946_k)
  "LSP1P5"
)  # panels i (labels) and n (heatmap)
hm_only_de <- TRUE  # panel n: TRUE = DE isoforms only; FALSE = every tested isoform
max_iso_per_gene <- Inf  # panel n: top isoforms per gene by padj (Inf = all)
z_lim      <- 2  # panel n: colour scale limits (row z-score)
y_cap      <- 9.5  # pC: |observed log2FC| drawn at most this far
cor_level  <- "transcript"  # panel h: "gene" or "transcript" (isoform)

dir_cols <- c(Down = "#1170ff", Up = "darkred")
## DE heatmap colours (as in the DE scripts): Moreland cool-warm
hm_cols  <- colorRampPalette(c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2",
                               "#F2CBB7", "#EE8468", "#B40426"), space = "Lab")(100)
cor_cols <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(50)
cc_cols  <- c("junction match (j)"      = npg_cols[["red"]],  # as Fig. 1
              "contains reference (k)"  = npg_cols[["salmon"]],
              "retained intron (m, n)"  = npg_cols[["brown"]],
              "sense overlap (i, o, y)" = npg_cols[["beige"]],
              "antisense (i, x)"        = npg_cols[["greyblue"]],
              "intergenic (u)"          = npg_cols[["mint"]])

## ------------------------------------------------ figure-wide styling -----
## base theme: line widths / spacing of the original panels (base size 6);
## every text size is set explicitly from the final figure
gg35 <- utils::packageVersion("ggplot2") >= "3.5.0"
theme_set(theme_bw(base_family = font_family, base_size = 6) + base_theme +
  theme(text            = element_text(face = "plain", lineheight = 1.2),   # no bold; line spacing as in the figure
        plot.title      = element_text(size = 10, face = "plain"),
        plot.subtitle   = element_text(size = 10, colour = "grey30", lineheight = 1),
        plot.title.position = "plot",                            # titles flush with the plot's left edge
        strip.text      = element_text(size = 8, face = "plain", margin = margin(1.5, 1.5, 1.5, 1.5)),
        axis.text       = element_text(size = 8),
        axis.title      = element_text(size = 10, face = "plain"),
        legend.text     = element_text(size = 8),
        legend.title    = element_text(size = 8, face = "plain"),
        legend.key.size = unit(4, "mm"),
        legend.position = "right",
        legend.margin   = margin(0, 0, 0, 0),
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        plot.margin     = margin(0.5, 0.5, 0.5, 0.5, "mm")))
update_geom_defaults("text", list(family = font_family, size = 8 / .pt))
pt2lw  <- function(pt) pt / (.pt * 0.75)                          # line width in pt -> ggplot linewidth
border <- function(pt) theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = pt2lw(pt)))
## legend beside the panel: vertical position (0 = bottom, 1 = top), as placed in the figure
leg_right <- function(y) theme(legend.position = "right", legend.justification = c(0, y))

## small boxed labels with a leader line (7.5 pt, as in the figure)
lab_repel <- function(data, mapping, ...) {
  geom_label_repel(data = data, mapping = mapping, size = 7 / .pt, family = font_family,  # 7 pt
                   colour = "black", fill = "white", label.size = 0.15,
                   label.r = unit(0.05, "lines"), label.padding = unit(0.08, "lines"),
                   box.padding = 0.3, point.padding = 0.15, min.segment.length = 0,
                   max.overlaps = Inf, segment.colour = "black", segment.size = 0.25,
                   show.legend = FALSE, seed = 1, ...)
}

## ================================================================ a =====
## read-level support (07_read_support.R; same design as Fig. 1g)
rs <- read_csv(file.path(rs_dir, "read_support_all_runs_long.csv"), show_col_types = FALSE)
pA <- read_support_plot(rs, label_fs = 6, min_lab = 7, n_sep = "\n") +   # 7: no clash at the narrower width
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.title.x = element_text(size = 8), axis.text.y = element_text(size = 8, colour = "grey20"),  # 8 pt (rows 1-2)
        axis.ticks.y = element_line()) +
  border(0.5) + leg_right(0.69)

## ================================================================ b =====
## novel transcripts by biotype x class code (validated pulposeq metadata;
## the class-code grouping of Fig. 1)
cc_group <- function(cls) case_when(
  grepl("\\(j\\)", cls)                   ~ "junction match (j)",
  grepl("\\(k\\)", cls)                   ~ "contains reference (k)",
  grepl("retained intron", cls)           ~ "retained intron (m, n)",
  grepl("^antisense", cls)                ~ "antisense (i, x)",
  grepl("\\(u\\)", cls)                   ~ "intergenic (u)",
  TRUE                                    ~ "sense overlap (i, o, y)")
bt_lab <- c(novel_lncRNA = "lncRNA", novel_protein_coding = "protein-coding",
            novel_non_coding = "non-coding")

dn <- read_csv(run_path("novel_meta"),
               show_col_types = FALSE)
dnov <- dn %>%
  transmute(biotype = factor(pulposeq_transcript_biotype, levels = names(bt_lab)),
            cls     = factor(cc_group(classification), levels = names(cc_cols)))
if (anyNA(dnov$biotype))
  message("  ! ", sum(is.na(dnov$biotype)), " novel transcripts with an unexpected biotype: ",
          paste(unique(dn$pulposeq_transcript_biotype[is.na(dnov$biotype)]), collapse = ", "))

dnov %>%
  dplyr::count(biotype, cls, .drop = FALSE) %>%
  mutate(`% of novel transcripts` = 100 * n / sum(n)) %>%
  group_by(biotype) %>%
  mutate(`% within biotype` = 100 * n / na_if(sum(n), 0)) %>%   # NA when a biotype has no novel transcripts
  ungroup() %>%
  transmute(Dataset = "Dong", Protocol = "Dong ONT cDNA", Biotype = bt_lab[as.character(biotype)],
            `Class code (GffCompare)` = cls, `Novel transcripts` = n,
            `% of novel transcripts`, `% within biotype`) %>%
  write_csv(file.path(tab_dir, "supp_novel_biotype_classcode_dong.csv"))

## share of each biotype's novel transcripts (0-100%), not counts; no totals
## above the bars and no y title (the totals are in a)
pB <- ggplot(dplyr::count(dnov, biotype, cls), aes(biotype, n, fill = cls)) +
  geom_col(width = 0.75, position = "fill") +
  scale_fill_manual(values = cc_cols, drop = FALSE) +
  scale_x_discrete(labels = bt_lab) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = NULL, fill = "GffCompare class code", title = "Novel transcripts") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  border(0.5) + leg_right(1)

## ============================================================== a + b ===
## a and b as one line plot in the style of Fig. 1g:
## one facet per biotype (annotated + novel of the same biotype), x = reference
## isoforms then the GffCompare classes of the novel ones, y = % of the
## isoforms in each read-support class (the a colours); the isoforms of each x
## are printed in one row above 100 %.
## Input: read-support/results/read_support_per_transcript_dong.csv (07).
k_lab <- function(x) ifelse(x >= 9950, paste0(round(x / 1000), "k"),   # 9,950+ -> "10k", not "10.0k"
                     ifelse(x >= 1000, sprintf("%.1fk", x / 1000), as.character(x)))   # 1.7k, 227
ab_facet <- c("annotated protein-coding" = "protein-coding", "predicted protein-coding" = "protein-coding",
              "annotated lncRNA" = "lncRNA", "novel lncRNA candidate" = "lncRNA",
              "other biotypes" = "other / non-coding", "predicted non-coding" = "other / non-coding")
rsp <- read_csv(file.path(rs_dir, "read_support_per_transcript_dong.csv"), show_col_types = FALSE) %>%
  left_join(dn %>% transmute(id = qry_id, classification), by = "id") %>%
  mutate(facet = factor(ab_facet[biotype], levels = unique(ab_facet)),
         cls   = factor(if_else(status == "annotated", "reference", cc_group(classification)),
                        levels = c("reference", names(cc_cols))),
         support = factor(support, levels = sup_lv))
if (anyNA(rsp$facet) || anyNA(rsp$cls))
  message("  ! panel a+b: ", sum(is.na(rsp$facet) | is.na(rsp$cls)), " transcripts without a biotype / class")
ab <- rsp %>% filter(!is.na(facet), !is.na(cls)) %>%
  dplyr::count(facet, cls, support, .drop = FALSE) %>%
  group_by(facet, cls) %>% mutate(total = sum(n), pct = 100 * n / total) %>% ungroup() %>%
  filter(total > 0)
write_csv(ab, file.path(tab_dir, "read_support_by_class_dong.csv"))
ab_n <- ab %>% group_by(facet, cls) %>%
  summarise(lab = k_lab(first(total)), .groups = "drop")          # one row above 100 %
## the read-support bars of a (Fig. 1g design) and the novel-transcript bars of
## b side by side in the 171 x 80 slot of panel a. The class table above is
## written for the Supplementary Tables.
pAB <- patchwork::wrap_plots(patchwork::free(pA), pB, nrow = 1, widths = c(1, 0.24))   # free(): a whole (labels, legend) ~105 mm; b bars ~25 mm

## ---- spike-in identification panel (pID; the design of Fig. 1b, longbench/08_fig1_panels.R) ----
## One row per spike-in set x protocol: expected (isoforms in the mix), detected
## (TP), fp_raw / fp_val = novel models on that set's sequences, Bambu output /
## after pulposeq validation (every spike-in is annotated, so any novel model on
## one is a false positive). Recall = TP / expected; precision = TP / (TP + FP).
id_metrics <- function(id) {
  id %>% mutate(recall        = detected / expected,
                precision     = detected / (detected + fp_val),
                f1            = 2 * precision * recall / (precision + recall),
                precision_raw = detected / (detected + fp_raw),
                f1_raw        = 2 * precision_raw * recall / (precision_raw + recall))
}
## Bambu vs oarfish on the Dong sequins: one dataset, one protocol, one
##  spike-in set, so the panel compares the two quantifiers of the same run, as
##  supp a1 does for LongBench. Rows = metrics, grouped: sequin identification
##  (recall = detected / in the mix, F1 with precision after validation) and
##  isoform-level DE (sensitivity, precision, F1, fold-change R2; undivided
##  counts for both quantifiers, so only the quantifier differs; the primary,
##  divided-count analysis is Fig. 2c). Two dots per row, Bambu light green
##  and oarfish blue, values printed right of each dot;
##  x-axis spans the data only (rounded down to 0.05).
## Fig. 2's quantifier colours. ED4 b (longbench/03_quantifier_plots.R) no longer
## uses them: there they repeated the ONT cDNA / ONT DRS colours of other panels
q_cols_dong <- c(Bambu = "#9bcf63", oarfish = "#2b8fd6")
## The two dots share the metric's row (Bambu just above, oarfish just below),
## and the values are a two-column table right of the axis, in the legend's
## order (Bambu, oarfish) and colours.
q_dots <- function(qd, label_pt = 7) {
  x_lo  <- floor((min(qd$value) - 0.02) * 20) / 20             # 0.05 below the lowest value
  rng   <- 1 - x_lo
  x_tab <- 1 + c(Bambu = 0.10, oarfish = 0.30) * rng           # value columns right of 1
  qd <- mutate(qd, yo = as.numeric(metric) + if_else(quantifier == "Bambu", -0.13, 0.13),
               xt = unname(x_tab[as.character(quantifier)]))
  ## values in a darker shade of the dot colour (the light green is too pale for text)
  lab_cols <- vapply(q_cols_dong, function(cl) colorRampPalette(c(cl, "black"))(10)[3], "")
  ggplot(qd, aes(value, yo, colour = quantifier)) +
    geom_point(size = 1.6) +
    geom_text(aes(x = xt, y = as.numeric(metric), label = sprintf("%.2f", value)),
              colour = lab_cols[as.character(qd$quantifier)],
              hjust = 0, size = label_pt / .pt, family = font_family) +
    ## group titles on top; panel heights follow the number of metrics
    facet_wrap(~ group, ncol = 1, scales = "free_y", space = "free_y") +
    scale_y_continuous(breaks = seq_along(levels(qd$metric)), labels = function(b) levels(qd$metric)[b],
                       trans = "reverse", expand = expansion(add = 0.5)) +
    scale_x_continuous(limits = c(x_lo, 1 + 0.48 * rng),
                       breaks = seq(x_lo, 1, if (rng > 0.3) 0.1 else 0.05),
                       labels = function(b) sub("\\.?0+$", "", sprintf("%.2f", b)),
                       expand = expansion(add = c(0.03 * rng, 0))) +
    scale_colour_manual(values = q_cols_dong, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1.5, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, colour = "grey20", margin = margin(1, 0, 0.8, 0)),
          legend.position = "top", legend.justification = c(0, 1), legend.location = "plot",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
          legend.key.size = unit(3, "mm"))
}
## what validation did to the false (novel) models on the set, for the subtitle
id_note <- function(id, what = "sequin") {
  d <- id_metrics(id)
  if (d$fp_val[1] == 0 && d$fp_raw[1] == 0)
    return(sprintf("No false %s models; precision = 1", what))
  if (d$fp_val[1] == 0)
    return(sprintf("Validation removed the only false %s models (%d); precision = 1",
                   what, as.integer(d$fp_raw[1])))
  sprintf("False %s models: %d after validation (%d before); precision = %.2f",
          what, as.integer(d$fp_val[1]), as.integer(d$fp_raw[1]), d$precision[1])
}

## =============================================================== pID =====
## Bambu vs oarfish on the Dong sequins (02_sequin_edger.R, run once per
## quantifier: sequins/results/ = oarfish, sequins/results_bambu/ = Bambu).
## Identification (sequin_identification.csv): of the 164 sequin isoforms in the
## mix, detected = >= 3 reads in >= 2 libraries; false positive = a novel model
## on a sequin (the same models for both quantifiers: validation does not depend
## on the quantifier). DE (summary.csv): isoform level, undivided counts.
q_seq_dir <- c(Bambu = paste0(seq_dir, "_bambu"), oarfish = seq_dir)   # as quant_dir() (helpers/spikein_counts.R)
qid <- imap_dfr(q_seq_dir, function(d, q) {
  s <- read_csv(file.path(d, "sequin_identification.csv"), show_col_types = FALSE)
  id_metrics(tibble(quantifier = q, set = "Sequins", protocol = "ONT cDNA",
                    detected = s$detected[1], expected = s$annotated[1],
                    fp_raw = s$novel_raw[1], fp_val = s$novel_validated[1]))
})
qde <- imap_dfr(q_seq_dir, function(d, q)
  read_csv(file.path(d, "summary.csv"), show_col_types = FALSE) %>%
    filter(level == "transcript", count_mode == "undivided") %>%
    transmute(quantifier = q, sensitivity, precision, f1, r2))
write_csv(left_join(qid, qde, by = "quantifier", suffix = c("_id", "_de")),
          file.path(tab_dir, "sequin_identification_dong.csv"))
## row labels (top to bottom); "F1 " = identification F1, "F1" = DE F1
## identification precision is shown too (after validation; = 1 when
## validation removed every false sequin model - id_note() goes to the log and
## the figure legend)
m_lv <- c(recall = "Recall", precision_id = "Precision ", f1_id = "F1 ",
          sensitivity = "Sensitivity", precision = "Precision",
          f1_de = "F1", fc_r2 = "Fold-change R\u00b2")
message("Fig. 2c, oarfish: ", id_note(filter(qid, quantifier == "oarfish")))
qd <- bind_rows(
  qid %>% transmute(quantifier, group = "Identification", recall, precision_id = precision, f1_id = f1) %>%
    pivot_longer(-c(quantifier, group), names_to = "metric"),
  qde %>% transmute(quantifier, group = "Isoform DE", sensitivity, precision, f1_de = f1, fc_r2 = r2) %>%
    pivot_longer(-c(quantifier, group), names_to = "metric")) %>%
  mutate(metric = factor(m_lv[metric], levels = m_lv),
         quantifier = factor(quantifier, levels = names(q_cols_dong)))
## grouped horizontal bar plot, axis from 0 (q_dots() kept above)
q_bars <- function(qd, label_pt = 6.5) {
  qd <- mutate(qd, metric = factor(metric, levels = rev(levels(metric))),   # first metric on top
               quantifier = factor(quantifier, levels = rev(names(q_cols_dong))))   # Bambu bar above oarfish
  dodge <- position_dodge(width = 0.82)
  ggplot(qd, aes(value, metric, fill = quantifier)) +
    geom_col(position = dodge, width = 0.78, orientation = "y") +
    geom_text(aes(label = sprintf("%.2f", value)), position = dodge, hjust = -0.15,
              size = label_pt / .pt, family = font_family, colour = "grey15") +
    facet_wrap(~ group, ncol = 1, scales = "free_y", space = "free_y") +
    scale_x_continuous(limits = c(0, 1.16), breaks = seq(0, 1, 0.25),
                       labels = function(b) sub("\\.?0+$", "", sprintf("%.2f", b)),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = 0.55)) +
    scale_fill_manual(values = q_cols_dong, breaks = names(q_cols_dong), name = NULL) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1.5, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, colour = "grey20", margin = margin(1, 0, 0.8, 0)),
          legend.position = "top", legend.justification = c(0, 1), legend.location = "plot",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
          legend.key.size = unit(3, "mm"))
}
pID <- q_bars(qd) +
  labs(title = "Sequins: Bambu vs oarfish",
       subtitle = "Same run, isoform level\nDE on undivided counts") +
  border(0.5)

## ================================================================ pC =====
## sequin DE ground truth, isoform level (02_sequin_edger.R: same test as the DE)
out_lev  <- c("true positive", "true negative", "false positive", "false negative", "wrong sign")
out_cols <- c("true positive" = "grey15", "true negative" = "grey75", "false positive" = "#D55E00",
              "false negative" = "#E69F00", "wrong sign" = "#CC79A7")
out_labs <- c(out_lev[-5], "wrong sign\n(opposite DE direction)")

## primary analysis: counts divided by their RTA overdispersion (count_mode "divided")
sq <- read_csv(file.path(seq_dir, "sequin_de_results.csv"), show_col_types = FALSE) %>%
  filter(level == "transcript", count_mode == "divided",
         !is.na(pvalue), is.finite(obs_lfc)) %>%   # as summarise_calls()
  mutate(outcome = factor(outcome, levels = out_lev),
         capped  = abs(obs_lfc) > y_cap,
         y       = pmax(pmin(obs_lfc, y_cap), -y_cap))
## top-left: DE calls among the TESTED sequins (the points); bottom-right: how
## many sequins could be tested (>= 10 reads in >= 3 libraries), then R2 / slope;
## that corner is empty by design (expected down, observed strongly up).
## Sequins detected / in the mix: panel c.
sq_s <- read_csv(file.path(seq_dir, "summary.csv"), show_col_types = FALSE) %>%
  filter(level == "transcript", count_mode == "divided") %>%
  mutate(lab_top = sprintf("precision = %.2f\nsensitivity = %.2f\nF1 = %.2f",
                           precision, sensitivity, f1),
         lab_bot = sprintf("%d tested\nR\u00b2 = %.2f\nslope = %.2f",
                           as.integer(tested), r2, slope))

pC <- ggplot(sq, aes(exp_lfc, y)) +
  annotate("rect", xmin = -lfc, xmax = lfc, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_point(aes(colour = outcome, shape = capped), size = 1.0, alpha = 0.9, show.legend = TRUE,
             position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  geom_text(data = sq_s, aes(x = -Inf, y = Inf, label = lab_top),
            hjust = -0.05, vjust = 1.1, lineheight = 1.46, size = 7 / .pt) +
  geom_text(data = sq_s, aes(x = Inf, y = -Inf, label = lab_bot),
            hjust = 1.05, vjust = -0.3, lineheight = 1.46, size = 7 / .pt) +
  scale_colour_manual(values = out_cols, labels = out_labs, drop = FALSE, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), breaks = "TRUE",
                     labels = sprintf("|observed log2FC| > %.1f", y_cap), name = NULL) +
  ## bottom room: the tested + R2 / slope block sits below the data
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.25))) +
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = 2)),
         shape  = guide_legend(order = 2, override.aes = list(size = 2, shape = 2))) +
  labs(x = "Expected log2FC\n(Mix A / Mix B)", y = "Observed log2FC (edgeR)",
       title = "Sequin ground truth",
       ## glmTreat TESTS against |log2FC| <= 0.5; it is not a fold-change filter
       ## three lines, so the second line stays inside the panel edge
       subtitle = paste0("Mix A (H1975) vs Mix B (HCC827)\nedgeR glmTreat (H0 |log2FC| <= ", lfc,
                         "), FDR < 0.05\nshaded: not truly DE")) +
  theme(legend.spacing.y = unit(0, "mm")) + leg_right(1)

## ================================================================ h =====
## sample-sample Spearman correlation (01_de_edger.R: cor(logCPM, method =
## "spearman") over the isoforms it tested), samples ordered by line
rho_m <- as.matrix(read.csv(file.path(de_dir, cor_level, "sample_spearman.csv"),
                            row.names = 1, check.names = FALSE))
ord <- c(sort(grep(paste0("^", line_lab[["TKI_sensitive"]]), colnames(rho_m), value = TRUE)),
         sort(grep(paste0("^", line_lab[["TKI_resistant"]]), colnames(rho_m), value = TRUE)))
rho_l <- as_tibble(rho_m, rownames = "s1") %>%
  pivot_longer(-s1, names_to = "s2", values_to = "rho") %>%
  mutate(s1 = factor(s1, levels = rev(ord)), s2 = factor(s2, levels = ord))

pD <- ggplot(rho_l, aes(s2, s1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradientn(colours = cor_cols, name = "Spearman \u03c1") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ## no coord_equal: the Fig. 2 file has a fixed square box (helpers/fig2_rows.R)
  labs(x = NULL, y = NULL, title = "Sample correlation",
       subtitle = paste0(str_to_sentence(lv_lab[[cor_level]]), "-level expression (logCPM)")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),  # 8 pt (rows 1-2)
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.key.height = unit(4, "mm")) +                    # colour bar 5 x 4 = 20 mm: inside the
                                                                # ~28 mm square of row 2
  leg_right(1)

## ================================================================ i =====
## volcano, gene level; only markers that were tested and are DE are labelled
res <- read_csv(file.path(de_dir, "all_results.csv"), show_col_types = FALSE)

vol <- res %>%
  filter(level == "gene") %>%
  mutate(symbol = sub("_ENSG.*$", "", feature),               # "LSP1P5_ENSG..." -> "LSP1P5"
         call   = factor(call, levels = c("Up", "Down", "NS")),
         y      = -log10(pmax(padj, 1e-300))) %>%
  arrange(dplyr::desc(call))                                 # NS drawn first
mk_l <- filter(vol, symbol %in% markers, call != "NS")
not_tested <- setdiff(markers, vol$symbol)
not_de     <- setdiff(filter(vol, symbol %in% markers, call == "NS")$symbol, mk_l$symbol)
if (length(not_tested))
  message("marker genes not tested at gene level (not labelled): ", paste(not_tested, collapse = ", "))
if (length(not_de))
  message("marker genes not DE at gene level (not labelled): ", paste(not_de, collapse = ", "))
n_up <- sum(vol$call == "Up"); n_dn <- sum(vol$call == "Down")

pE <- ggplot(vol, aes(log2FC, y, colour = call)) +
  geom_point(data = filter(vol, call == "NS"), size = 0.7, alpha = 0.4, stroke = 0) +
  geom_point(data = filter(vol, call != "NS"), size = 1.0, alpha = 0.6, stroke = 0) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(alpha), linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_point(data = mk_l, shape = 21, fill = NA, colour = "black", size = 1.6, stroke = 0.5) +
  ## the gene labels stay below the data's top; the
  ## Down / Up counts sit in the 12% headroom of the y scale (CEACAM5 ran into "Down: 4,272")
  lab_repel(mk_l, aes(label = symbol), ylim = c(NA, max(vol$y))) +
  annotate("text", x = -Inf, y = Inf, label = paste("Down:", scales::comma(n_dn)),
           hjust = -0.15, vjust = 1.6, colour = dir_cols[["Down"]], fontface = "plain",
           size = 8 / .pt) +  # 8 pt (Fig. 2 text)
  annotate("text", x = Inf, y = Inf, label = paste("Up:", scales::comma(n_up)),
           hjust = 1.15, vjust = 1.6, colour = dir_cols[["Up"]], fontface = "plain",
           size = 8 / .pt) +
  scale_colour_manual(values = c(dir_cols[c("Up", "Down")], NS = "grey80"),
                      labels = c(Up = "upregulated", Down = "downregulated", NS = "not significant"),
                      name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "log2FC, H1975 / HCC827 (predFC)", y = "-log10(FDR)",
       title = "Differential Gene Expression",
       subtitle = paste0("FDR ", alpha, ", |log2FC| > ", lfc)) +
  leg_right(0.66)

## ================================================================ n =====
## DE isoforms of the marker genes. Gene = Bambu's gene assignment (GENEID of
## the transcript count table, i.e. what the DE used), named via the
## reference metadata. Rows = isoforms grouped by gene (marker order) and
## ranked by log2FC; columns = samples + the isoform DE call.
iso_val <- read_tsv(run_path("counts_tx"),
                    show_col_types = FALSE, col_select = c(TXNAME, GENEID))
xt   <- readRDS(file.path(de_dir, "transcript", "logcpm.rds"))
grp  <- setNames(as.character(xt$group), colnames(xt$mat))
smp  <- c(names(grp)[grp == "TKI_sensitive"], names(grp)[grp == "TKI_resistant"])
refm <- read_csv(run_path("ann_meta"),
                 col_select = c(gene_id_version, gene_name), show_col_types = FALSE) %>%
  distinct(gene_id_version, .keep_all = TRUE)

iso <- res %>%
  filter(level == "transcript", feature %in% rownames(xt$mat)) %>%
  left_join(dplyr::select(iso_val, id = TXNAME, gene_id = GENEID), by = "id") %>%
  left_join(refm, by = c(gene_id = "gene_id_version")) %>%
  filter(gene_name %in% markers, !hm_only_de | call != "NS") %>%
  mutate(gene = factor(gene_name, levels = intersect(markers, gene_name))) %>%
  group_by(gene) %>%
  { if (is.finite(max_iso_per_gene))
      slice_min(., padj, n = max_iso_per_gene, with_ties = FALSE) else . } %>%
  ungroup() %>%
  arrange(gene, dplyr::desc(log2FC_test))
message(sprintf("panel g: %d isoforms from %d marker genes%s", nrow(iso), n_distinct(iso$gene),
                if (hm_only_de) " (DE isoforms only)" else ""))
no_iso <- setdiff(markers, as.character(iso$gene))
if (length(no_iso))
  message("marker genes with no ", if (hm_only_de) "DE " else "", "isoform in panel g: ",
          paste(no_iso, collapse = ", "))

## aligned with the cis-pair heatmap of 04_cis_circos.R, so the two stack one
## above the other. The panel is drawn on a millimetre grid
##  shared with that heatmap (HM_FIG there; KEEP THE TWO IN STEP): gene names
##  right-aligned at x_name on the left, the DE call in the square column
##  (where the cis heatmap has its biotype square), HCC827_1-3 | H1975_1-3
##  tiles at the same x, line names on top, the colour bar top right and the
##  Up / Down key under it. Rows: isoforms grouped by gene (marker order),
##  ranked by log2FC, a gap between genes; names spread apart where blocks are
##  closer than one name height, with a leader line to their block.
hm_geo <- list(width = 82.7, x_tile1 = 32.82, x_tile2 = 41.96, tile_dx = 2.94, tile_w = 2.78,
               x_cat = 29.94, cat_w = 2.14, x_name = 28.45, y_head = 1.54, y_tile0 = 3.66,
               bar_x = 72.25, bar_w = 2.73, bar_top = 13.75, bar_h = 18.9,
               x_bar_title = 71.7, y_bar_title = c(6.15, 9.55), x_bar_lab = 76.5)
hm_h   <- 215.8                                                 # panel height, mm (A3 Fig. 2)
hm_pt  <- 8                                                     # names, heads, key (as the cis heatmap)
smp_hm <- c(sort(names(grp)[grp == "TKI_sensitive"]), sort(names(grp)[grp == "TKI_resistant"]))
z <- t(scale(t(xt$mat[iso$feature, smp_hm, drop = FALSE])))
n_ref <- sum(grp[smp_hm] == "TKI_sensitive")
x_left <- setNames(with(hm_geo, c(x_tile1 + tile_dx * (seq_len(n_ref) - 1),
                                  x_tile2 + tile_dx * (seq_len(length(smp_hm) - n_ref) - 1))), smp_hm)
hm_gap <- 0.85                                                  # rows of white between genes
rows <- iso %>%
  transmute(feature, gene, call = factor(call, levels = c("Up", "Down", "NS")),
            g_no = cumsum(gene != lag(gene, default = first(gene))),
            k    = row_number() - 1 + hm_gap * g_no)            # row offset from the top, in rows
row_mm <- (hm_h - hm_geo$y_tile0 - 1) / (max(rows$k) + 1)       # row height, mm
rows <- rows %>% mutate(top = hm_geo$y_tile0 + k * row_mm, mid = top + row_mm / 2)
hm <- as_tibble(z, rownames = "feature") %>%
  pivot_longer(-feature, names_to = "sample", values_to = "z") %>%
  left_join(dplyr::select(rows, feature, top), by = "feature") %>%
  mutate(x0 = unname(x_left[sample]), z = pmax(pmin(z, z_lim), -z_lim))
## gene names at the centre of their block, spread apart (top to bottom) where
## closer than one name height; a leader line joins a moved name to its block
lab_h <- hm_pt * 25.4 / 72 * 1.15                               # one name height, mm
spread_labels <- function(y, d) {                               # y: top to bottom (increasing)
  for (it in seq_len(500)) {
    moved <- FALSE
    for (i in seq_along(y)[-1]) {
      ov <- d - (y[i] - y[i - 1])
      if (ov > 1e-6) { y[i - 1] <- y[i - 1] - ov / 2; y[i] <- y[i] + ov / 2; moved <- TRUE }
    }
    if (!moved) break
  }
  y
}
g_lab <- rows %>% group_by(gene) %>%
  summarise(y = (min(top) + max(top) + row_mm) / 2, .groups = "drop") %>%
  arrange(y) %>%
  mutate(y_lab = spread_labels(y, lab_h), moved = abs(y_lab - y) > 0.4)
fill_of <- scales::gradient_n_pal(hm_cols)
z_col   <- function(v) fill_of((v + z_lim) / (2 * z_lim))
hm$col <- z_col(hm$z)
bar_y   <- function(v) with(hm_geo, bar_top + bar_h * (z_lim - v) / (2 * z_lim))
bar     <- tibble(v0 = seq(-z_lim, z_lim, length.out = 101)) %>%
  mutate(v1 = lead(v0)) %>% filter(!is.na(v1)) %>% mutate(col = z_col((v0 + v1) / 2))
tk      <- seq(-z_lim, z_lim, by = 1)
key_y   <- hm_geo$bar_top + hm_geo$bar_h + 10 + c(0, 4.2)       # Up / Down key rows (title "DE" above)
hm_txt  <- function(...) annotate("text", family = font_family, size = hm_pt / .pt, ...)
pF <- ggplot() +
  geom_rect(data = hm, aes(xmin = x0, xmax = x0 + hm_geo$tile_w, ymin = -top, ymax = -(top + row_mm),
                           fill = col), colour = NA) +
  scale_fill_identity() +
  geom_rect(data = rows, aes(xmin = hm_geo$x_cat, xmax = hm_geo$x_cat + hm_geo$cat_w,
                             ymin = -top, ymax = -(top + row_mm),
                             colour = call, fill = c(dir_cols, NS = "grey80")[as.character(call)]),
            linewidth = 0) +
  scale_colour_manual(values = c(dir_cols, NS = "grey80"), guide = "none") +
  geom_segment(data = filter(g_lab, moved),
               aes(x = hm_geo$x_name + 0.25, y = -y_lab, xend = hm_geo$x_cat - 0.25, yend = -y),
               colour = "grey55", linewidth = 0.2) +
  geom_text(data = g_lab, aes(hm_geo$x_name, -y_lab, label = gene), hjust = 1,
            size = hm_pt / .pt, colour = "black", family = font_family) +
  ## line names over their tiles, "DE" over the call column
  hm_txt(x = mean(range(x_left[seq_len(n_ref)]) + c(0, hm_geo$tile_w)), y = -hm_geo$y_head,
         label = line_lab[["TKI_sensitive"]]) +
  hm_txt(x = mean(range(x_left[-seq_len(n_ref)]) + c(0, hm_geo$tile_w)), y = -hm_geo$y_head,
         label = line_lab[["TKI_resistant"]]) +
  ## colour bar (top right) and the Up / Down key under it
  geom_rect(data = bar, aes(xmin = hm_geo$bar_x, xmax = hm_geo$bar_x + hm_geo$bar_w,
                            ymin = -bar_y(v0), ymax = -bar_y(v1), fill = col), colour = NA) +
  annotate("segment", x = hm_geo$bar_x, xend = hm_geo$bar_x + 0.2 * hm_geo$bar_w,
           y = -bar_y(tk), yend = -bar_y(tk), colour = "white", linewidth = 0.25) +
  annotate("segment", x = hm_geo$bar_x + 0.8 * hm_geo$bar_w, xend = hm_geo$bar_x + hm_geo$bar_w,
           y = -bar_y(tk), yend = -bar_y(tk), colour = "white", linewidth = 0.25) +
  hm_txt(x = hm_geo$x_bar_lab, y = -bar_y(tk), label = formatC(tk, format = "f", digits = 0), hjust = 0) +
  hm_txt(x = hm_geo$x_bar_title, y = -hm_geo$y_bar_title, label = c("z-score", "(logCPM)"), hjust = 0) +
  annotate("rect", xmin = hm_geo$bar_x, xmax = hm_geo$bar_x + hm_geo$bar_w,
           ymin = -(key_y - 1.2), ymax = -(key_y + 1.2), fill = dir_cols[c("Up", "Down")]) +
  hm_txt(x = hm_geo$x_bar_title, y = -(key_y[1] - 4.6), label = "DTE", hjust = 0) +  # DTE key title
  hm_txt(x = hm_geo$x_bar_lab, y = -key_y, label = c("up in H1975", "down in H1975"), hjust = 0) +
  coord_cartesian(xlim = c(0, hm_geo$width), ylim = c(-hm_h, 0), expand = FALSE, clip = "off") +
  labs(title = "Differential Transcript Expression",
       subtitle = paste(c(if (is.finite(max_iso_per_gene)) paste("Top", max_iso_per_gene),
                          if (hm_only_de) "DE" else "Tested",
                          "isoforms per locus, ranked by log2FC"), collapse = " ")) +
  theme_void(base_family = font_family) +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 10, colour = "grey30"),
        plot.title.position = "plot", plot.margin = margin(0, 0, 0, 0), legend.position = "none")

## ===================================================== single panels =====
## one file per panel; size (width, height, mm) = the panel's extent on the
## final Fig. 2 page (A3) + the 0.5 mm plot margin on each side
## pID, pC (fig2c_*, fig2d_*) are drawn but not used in the assembled Fig. 2;
## ab = a and b side by side (171 x 80); e (Fig. 2h), f (Fig. 2i) and g
## (Fig. 2n): their slots in the assembled Fig. 2 (11_assemble_figures.R);
## g 82.7 wide, the width of the cis-pair heatmap it aligns with (no-title
## version: exactly the shared millimetre grid; the titled one is 13 mm taller)
panels <- list(
  fig2a_read_support          = list(p = pA,  size = c(100.0,  69.8)),  # as in the a slot
  fig2b_novel_biotype         = list(p = pB,  size = c( 60.9,  74.2)),
  fig2ab_read_support_class   = list(p = pAB, size = c(171.0,  80.0)),
  fig2c_sequin_identification = list(p = pID, size = c( 62.0,  80.0)),
  fig2d_sequin_isoform        = list(p = pC,  size = c( 89.0,  80.0)),
  fig2e_sample_correlation    = list(p = pD,  size = c( 61.6,  70.0)),
  fig2f_volcano               = list(p = pE,  size = c(130.0,  80.0)),
  fig2g_isoform_heatmap       = list(p = pF,  size = c(hm_geo$width, hm_h))
)
## figure data behind every panel, kept for sharing
source(file.path(repo, "downstream/helpers/source_data.R"))
save_source_data(panels, file.path(tab_dir, "source_data_fig2.rds"))
title_h <- c(fig2g_isoform_heatmap = 13)                   # extra height for title + subtitle
for (nm in names(panels))
  save_panel(panels[[nm]]$p, nm, panels[[nm]]$size + c(0, if (is.na(title_h[nm])) 0 else title_h[[nm]]), panel_dir)
## the same panels without title and subtitle (same file names and sizes), in
## figures/fig2_panels_no_titles/ (as longbench/08_fig1_panels.R for Fig. 1)
panel_dir_nt <- paste0(panel_dir, "_no_titles")
no_titles <- function(p) {
  if (inherits(p, "patchwork")) p & labs(title = NULL, subtitle = NULL)  # every sub-plot
  else p + labs(title = NULL, subtitle = NULL)
}
for (nm in names(panels))
  save_panel(no_titles(panels[[nm]]$p), nm, panels[[nm]]$size, panel_dir_nt)

## Fig. 2 rows 1-2 aligned: a, b and the sample correlation
## (row 2, h) are saved again in fig2_panels_no_titles/ with the box sizes and
## top / bottom room shared with 10's panels (helpers/fig2_rows.R; ggh4x).
## These overwrite the three no-title files above; the titled ones are unchanged.
source(file.path(helper_dir, "fig2_rows.R"))
## a's legend on the right of its bars (one column), so a keeps fig2_a_w of
## bars (helpers/fig2_rows.R); a legend on top would push every row-1 box down
f2_06 <- list(fig2a_read_support       = no_titles(pA) + guides(fill = guide_legend(ncol = 1)) +
                theme(legend.position = "right", legend.justification = c(0, 0.5),
                      legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0)),
              fig2b_novel_biotype      = no_titles(pB),
              fig2e_sample_correlation = no_titles(pD))
fig2_update_rooms(f2_06)
lay <- fig2_layout()
for (nm in names(f2_06)) fig2_save(f2_06[[nm]], nm, lay, panel_dir_nt)
message(sprintf("Fig. 2 rows 1-2: boxes %.1f mm tall (row 2 square); row 1: a's bars %.1f mm, c and d boxes %.1f mm wide", lay$box2, lay$a_w, lay$w1))

## row 3 = volcano | the four profiles of 06 (section 7e,
## fig2_profile_rooms.rds). The volcano gets the profiles' file height, its box
## the profiles' box height with top / bottom edges on theirs, and the width
## the profiles leave on the page. 8 pt text (Fig. 2).
pr_file <- file.path(iso_paper, "fig2_panels", "fig2_profile_rooms.rds")
if (!file.exists(pr_file)) {
  message("  ! no ", pr_file, " (run 04 first): the volcano keeps its old size")
} else {
  lg  <- readRDS(pr_file)[["fig2_profile_legend"]]   # 06's shared legend after the profiles (NULL: none)
  pr  <- readRDS(pr_file)[fig2_row3_profiles]
  if (any(sapply(pr, is.null))) stop("fig2_profile_rooms.rds lacks a row-3 profile (rerun 04): ",
                                     paste(fig2_row3_profiles[sapply(pr, is.null)], collapse = ", "))
  Tp  <- max(sapply(pr, `[[`, "T")); Bp <- max(sapply(pr, `[[`, "B")); bh <- pr[[1]][["box_h"]]
  W3  <- fig2_W - 2 * fig2_mar - fig2_gut * length(pr) - sum(sapply(pr, `[[`, "W")) -
    (if (is.null(lg)) 0 else lg[["W"]] + fig2_gut)
  pv  <- no_titles(pE) + fig2_text
  fig2_measure_open(); rv <- fig2_room(pv); fig2_measure_close()  # the SVG device's metrics
  bw  <- W3 - rv[["L"]] - rv[["R"]] - 2 * fig2_slack
  pv  <- fig2_fix(pv, bw, bh) +
    theme(plot.margin = margin(Tp - rv[["T"]], fig2_slack, Bp - rv[["B"]], fig2_slack, "mm"))
  save_panel(pv, "fig2f_volcano", c(W3, Tp + bh + Bp), panel_dir_nt)
  message(sprintf("Fig. 2 row 3: volcano %.1f x %.1f mm (box %.1f x %.1f)", W3, Tp + bh + Bp, bw, bh))
}
message("Done: ", length(panels), " panels in ", panel_dir, " and ", panel_dir_nt)
