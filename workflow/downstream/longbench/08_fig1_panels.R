## ===========================================================================
## 08_fig1_panels.R - Fig. 1 b-k: the LongBench panels used in the paper
##
## Panel  object  content
##   a     -       pulposeq workflow schematic (drawn in the vector editor; not here)
##   b     pID     spike-in identification: recall (detected / expected, counts at
##                 the points) and F1 of each spike-in set (SIRV E0, long SIRV, ERCC,
##                 sequins) per protocol, one line per set;
##                 note under the legend: false spike-in models removed by validation
##                 (false positive = novel model on a spike-in; precision / F1 in
##                 spikein_identification_by_set.csv)
##   c     pB      length bias, equimolar SIRVs; corner text: % of E0 isoforms
##                 whose median observed / expected abundance is within 2-fold
##   d     pD      ERCC dose response; R2, LOD (70 % detection)
##   e     pC      sequin ground truth, isoform level: observed vs expected
##                 log2FC per protocol, points coloured by outcome (truly DE =
##                 |expected log2FC| > 0.5; called = FDR < 0.05; wrong-sign calls
##                 count as false discoveries); DE precision / sensitivity / F1
##                 top-left, sequins tested + R2 / slope bottom-right
##   f     pQ      where the quantifier matters: log2(oarfish / Bambu reads) per
##                 spike-in transcript against its length (same runs; from the
##                 03_quantifier_plots.R table; E0 isoforms > 8-fold labelled)
##   g     pE      read-level support across biotypes: one facet per biotype, x =
##                 protocol, one line per read-support class; isoforms per protocol printed
##   h     pF      novel transcript length (density), medians labelled
##   i     pG      UIC sharing per protocol (multi-exon novel transcripts)
##   j     pH      UIC sharing by length (n printed at each point)
##   k     pI      transcript ends within an identical UIC (ECDF of |d5| / |d3|
##                 per protocol pair); note = annotation agreement among shared UICs
##
## Layout = the final Fig. 1 (A3 page): every panel is saved at the size it has
## on that page, with the font sizes,
## legend positions, key sizes and panel-border widths measured there
## (panel titles 10 pt, subtitles 9 pt, axis / legend text 8 pt, strips 9-9.5 pt;
## see `spec` below). No bold text. Theme base size 6.
##
## UIC = unique intron chain (LRGASP, Pardo-Palacios et al. 2024): chromosome +
## strand + the exact coordinates of every intron; ends are free, so protocols
## sharing a UIC can differ in transcript ends. Mono-exon models have no UIC.
## The inputs keep their original names (unique_intron_chains/results/, column UJC_ID).
##
## Needs (run first): 01_sequin_edger.R, 02_sirv_ercc.R, 04_ujc_report.R,
##   05_ujc_length.R, 06_read_support.R, 03_quantifier_plots.R (panel f);
##   Anaquin (sequin totals)
## Outputs
##   figures/fig1_panels/fig1{b..k}_*.{pdf,svg,png}       (svg = editable text)
##   figures/fig1_panels_no_titles/                        (the same, no title / subtitle)
##   figure-panels/results/figure_tables/
##     spikeins_detected.csv                 detected / in the mix, every spike-in type
##     table_spikeins_detected.csv           the same as "n / total (%)"
##     table_sequin_metrics.csv              sequin DE metrics, gene and isoform level
##     table_spikein_identification.csv      novel models on spike-ins, Bambu vs validated
##     spikein_identification_by_set.csv     panel b: detected, expected, FP, recall,
##                                           precision, F1 per spike-in set x protocol
##     table_ujc_annotation_agreement.csv    annotation agreement among shared UICs
##     novel_class_codes.csv                 novel transcripts by biotype x class code
##     supp_novel_biotype_classcode_longbench.csv  (the same, supplementary-table layout)
##     table_read_support.csv                panel g: n (%) per protocol x biotype x read-support class
## Needs: tidyverse, Anaquin, svglite, showtext, systemfonts
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
source(file.path(helper_dir, "plot_theme.R"))       # font_family, base_theme, protocol_cols, npg_cols, plat_labels
source(file.path(helper_dir, "figure_helpers.R"))   # save_panel; read-support labels / colours (bt_disp, sup_*)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

prot_levels <- c("ONT cDNA", "ONT DRS", "PacBio")
prot        <- c(ONTc = "ONT cDNA", ONTd = "ONT DRS", PB = "PacBio")

## ------------------------------------------------ figure-wide styling -----
## base theme: line widths / spacing of the original panels (base size 6);
## every text size is set explicitly from the final figure
gg35 <- utils::packageVersion("ggplot2") >= "3.5.0"
theme_set(theme_bw(base_family = font_family, base_size = 6) + base_theme +
  theme(text            = element_text(face = "plain", lineheight = 1.2),   # no bold; line spacing as in the figure
        plot.title      = element_text(size = 10, face = "plain"),
        plot.subtitle   = element_text(size = 9, colour = "grey30"),
        plot.title.position = "plot",                            # titles flush with the plot's left edge
        strip.text      = element_text(size = 9, face = "plain", margin = margin(1.5, 1.5, 1.5, 1.5)),
        axis.text       = element_text(size = 8),
        axis.title      = element_text(size = 8, face = "plain"),
        legend.text     = element_text(size = 8),
        legend.title    = element_text(size = 8, face = "plain"),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        plot.margin     = margin(0.5, 0.5, 0.5, 0.5, "mm")))
update_geom_defaults("text", list(family = font_family, size = 8 / .pt))
pt2lw  <- function(pt) pt / (.pt * 0.75)                          # line width in pt -> ggplot linewidth
border <- function(pt) theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = pt2lw(pt)))
## legend beside the panel: vertical position (0 = bottom, 1 = top), as placed in the figure
leg_right <- function(y) theme(legend.position = "right", legend.justification = c(0, y))
## legends on top: keys in a row, title above them, flush with the plot's left edge
leg_top <- theme(legend.position = "top", legend.justification = "left",
                 legend.direction = "horizontal", legend.box = "vertical",
                 legend.box.just = "left", legend.margin = margin(0, 0, 0, 0),
                 legend.box.spacing = unit(2, "mm")) +
  if (gg35) theme(legend.title.position = "top", legend.location = "plot") else theme()
## text placed at an exact distance (mm) from the top-left / top-right corner of
## a single panel: one line per element of `labels`, `pitch` mm apart
corner_text <- function(labels, colours, size_pt, top_mm, pitch_mm, left = TRUE, inset_mm = 1.5) {
  lapply(seq_along(labels), function(i) {
    annotation_custom(grid::textGrob(
      labels[i],
      x = if (left) unit(inset_mm, "mm") else unit(1, "npc") - unit(inset_mm, "mm"),
      y = unit(1, "npc") - unit(top_mm + (i - 1) * pitch_mm, "mm"),
      hjust = if (left) 0 else 1, vjust = 1,
      gp = grid::gpar(col = colours[i], fontsize = size_pt, fontfamily = font_family)))
  })
}

## ------------------------------------------------------------- colours ----
set_cols   <- c("SIRV E0" = "grey35", "long SIRV" = npg_cols[["red"]])
## every spike-in set (panels b and f)
q_spk_cols <- c(set_cols, "ERCC" = npg_cols[["greyblue"]], "sequin" = npg_cols[["mint"]])
cc_cols    <- c("junction match (j)"      = npg_cols[["red"]],
                "contains reference (k)"  = npg_cols[["salmon"]],
                "retained intron (m, n)"  = npg_cols[["brown"]],
                "sense overlap (i, o, y)" = npg_cols[["beige"]],
                "antisense (i, x)"        = npg_cols[["greyblue"]],
                "intergenic (u)"          = npg_cols[["mint"]])
share_cols <- c("all 3 protocols" = "grey15", "1 other protocol" = "grey55",
                "this protocol only" = "grey85")

## =================================================== spike-in data =======
## primary analysis only: counts divided by their RTA overdispersion (01 also
## writes count_mode "undivided" at gene level, as a sensitivity check)
fc_sum <- read_csv(file.path(seq_dir, "edger/fold_change_summary_edger.csv"), show_col_types = FALSE) %>%
  filter(count_mode == "divided") %>%
  mutate(protocol = factor(plat_labels[platform], levels = prot_levels))
se  <- read_csv(file.path(sirv_dir, "sirv_ercc_summary.csv"), show_col_types = FALSE) %>%
  mutate(protocol = factor(protocol, levels = prot_levels))
lb  <- read_csv(file.path(sirv_dir, "length_bias.csv"), show_col_types = FALSE) %>%
  mutate(protocol = factor(protocol, levels = prot_levels))
er  <- read_csv(file.path(sirv_dir, "ercc_per_transcript.csv"), show_col_types = FALSE) %>%
  mutate(protocol = factor(protocol, levels = prot_levels))
sid <- read_csv(file.path(seq_dir, "sequin_identification.csv"), show_col_types = FALSE) %>%
  mutate(protocol = factor(protocol, levels = prot_levels))
e0  <- map_dfr(platforms, function(p)
  read_csv(file.path(sirv_dir, p, "sirv_e0_per_isoform.csv"), show_col_types = FALSE)) %>%
  mutate(protocol = factor(protocol, levels = prot_levels))

## ------------------------------------- spike-ins detected (tables) ------
## sequins in the mix (Anaquin) and SIRVs: "detected" = >= 3 reads in >= 2 libraries
## (01, column `detected`; 02); ERCC: >= 1 read in >= 2 libraries (02; ERCCs are
## single-exon and mostly low input, so the stricter rule would only measure depth)
data("RnaQuinIsoformMixture", package = "Anaquin")
sq_iso <- as.character(RnaQuinIsoformMixture$Name)
sq_tot <- c(gene = n_distinct(sub("_[0-9]+$", "", sq_iso)), transcript = length(sq_iso))

det <- bind_rows(
  fc_sum %>% transmute(protocol, n = detected, total = unname(sq_tot[level]),
                       type = if_else(level == "gene", "Sequin genes", "Sequin isoforms")),
  se %>% transmute(protocol, n = e0_detected, type = "SIRV E0") %>%
    left_join(dplyr::count(e0, protocol, name = "total"), by = "protocol"),
  se %>% transmute(protocol, n = long_detected, type = "long SIRV") %>%
    left_join(dplyr::count(filter(lb, set == "long SIRV"), protocol, name = "total"),
              by = "protocol"),
  se %>% transmute(protocol, n = ercc_detected, total = ercc_total, type = "ERCC")) %>%
  mutate(type = factor(type, levels = rev(c("Sequin genes", "Sequin isoforms",
                                            "SIRV E0", "long SIRV", "ERCC"))))
write_csv(det, file.path(tab_dir, "spikeins_detected.csv"))

det %>%
  mutate(type = fct_rev(type),
         cell = sprintf("%d / %d (%.0f%%)", n, total, 100 * n / total)) %>%
  arrange(type, protocol) %>%
  dplyr::select(`Spike-in` = type, protocol, cell) %>%
  pivot_wider(names_from = protocol, values_from = cell) %>%
  write_csv(file.path(tab_dir, "table_spikeins_detected.csv"))

fc_sum %>%
  transmute(Tool = "edgeR", Level = if_else(level == "gene", "gene", "isoform"), Protocol = protocol,
            Detected = detected, Tested = tested,
            `FC R2` = round(r2, 2), `FC slope` = round(slope, 2),
            `Truly DE` = truly_de, Called = n_called, TP = true_pos, FP = false_pos,
            `Wrong sign` = wrong_sign, FN = false_neg,
            AUC = round(auc, 2), Sensitivity = round(sensitivity, 2),
            Precision = round(precision, 2), FDR = round(fdr, 2), F1 = round(f1, 2)) %>%
  arrange(Level, Protocol) %>%
  write_csv(file.path(tab_dir, "table_sequin_metrics.csv"))

## spike-in identification, Bambu vs pulposeq-validated: every spike-in
## transcript is annotated, so any novel model on a spike-in sequence is a false
## discovery. Recall = annotated spike-ins detected (>= 3 reads in >= 2 libraries); precision =
## detected / (detected + novel models). ERCCs: novel models only (their
## detection is set by concentration, panel d).
f2 <- function(x) sprintf("%.2f", x)
bind_rows(
  sid %>% transmute(`Spike-in` = "Sequin isoforms", Protocol = protocol,
                    `In annotation` = annotated, Detected = detected, Recall = f2(recall),
                    `Novel models (Bambu)` = novel_raw, `Novel models (validated)` = novel_validated,
                    `Precision (Bambu)` = f2(precision_raw), `Precision (validated)` = f2(precision_validated),
                    `F1 (Bambu)` = f2(f1_raw), `F1 (validated)` = f2(f1_validated)),
  se %>% transmute(`Spike-in` = "SIRVs (E0 + long)", Protocol = protocol,
                   `In annotation` = sirv_tp + sirv_fn, Detected = sirv_tp, Recall = f2(sirv_recall),
                   `Novel models (Bambu)` = novel_sirv_raw, `Novel models (validated)` = novel_sirv_validated,
                   `Precision (Bambu)` = f2(sirv_precision_raw), `Precision (validated)` = f2(sirv_precision),
                   `F1 (Bambu)` = f2(sirv_f1_raw), `F1 (validated)` = f2(sirv_f1)),
  se %>% transmute(`Spike-in` = "ERCC", Protocol = protocol,
                   `In annotation` = ercc_total, Detected = ercc_detected, Recall = "-",
                   `Novel models (Bambu)` = novel_ercc_raw, `Novel models (validated)` = novel_ercc_validated,
                   `Precision (Bambu)` = "-", `Precision (validated)` = "-",
                   `F1 (Bambu)` = "-", `F1 (validated)` = "-")) %>%
  mutate(`Spike-in` = factor(`Spike-in`, levels = c("Sequin isoforms", "SIRVs (E0 + long)", "ERCC"))) %>%
  arrange(`Spike-in`, Protocol) %>%
  write_csv(file.path(tab_dir, "table_spikein_identification.csv"))

## ---- spike-in identification panel (Fig. 1b) ----
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
## line plot: two facets, recall (detected / expected, the counts printed at each
## point) and F1 (value printed at each point); y-axis from 0.25; x = protocol
## (spike-in depth rises left to right), one line per set in the Fig. 1f colours.
## Precision is not drawn: after validation it is 1 in every set; the false
## models validation removed are in a note under the legend. Sensitivity = recall
## here, and accuracy is undefined (no true negatives: every annotated spike-in
## is in the mix).
id_lines <- function(id, set_cols, prot_levels) {
  d <- id_metrics(id) %>%
    mutate(set = factor(if_else(set == "Sequins", "sequin", set), levels = names(set_cols)),
           protocol = factor(protocol, levels = prot_levels))
  m_lab <- c(Recall = "Recall (detected / expected)", F1 = "F1-score")
  long <- bind_rows(
    d %>% transmute(set, protocol, metric = "Recall", value = recall,
                    label = sprintf("%d/%d", as.integer(detected), as.integer(expected))),
    d %>% transmute(set, protocol, metric = "F1", value = f1, label = sprintf("%.2f", f1))) %>%
    mutate(metric = factor(metric, levels = names(m_lab)))
  ## count labels in a darker shade of the set colour (mint is too pale for text)
  lab_cols <- vapply(set_cols, function(cl) colorRampPalette(c(cl, "black"))(10)[3], "")
  ## fixed label places, clear of the lines: right of the last protocol's points
  ## (lines only arrive from the left); elsewhere long SIRV above, E0 below, and of
  ## ERCC / sequin the higher above, the lower below
  last_p <- tail(prot_levels, 1)
  lab_d  <- long %>%
    group_by(metric, protocol) %>%
    mutate(above = case_when(set == "long SIRV" ~ TRUE, set == "SIRV E0" ~ FALSE,
                             TRUE ~ value >= max(value[set %in% c("ERCC", "sequin")]))) %>%
    ungroup() %>%
    mutate(last  = protocol == last_p,
           lx    = as.numeric(protocol) + if_else(last, 0.1, 0),
           ly    = case_when(last ~ value, above ~ value + 0.035, TRUE ~ value - 0.035),
           hjust = if_else(last, 0, 0.5),
           vjust = case_when(last ~ 0.5, above ~ 0, TRUE ~ 1))
  ## the note is only true if validation left no false model
  stopifnot("validation left false spike-in models: rewrite the panel b note" = all(d$fp_val == 0))
  nb <- function(x) gsub(" ", "\u00a0", x)       # no line break inside a set / protocol name
  fp_txt <- d %>% filter(fp_raw > 0) %>% arrange(set, protocol) %>%
    group_by(set) %>%
    summarise(t = paste(sprintf("%d\u00a0in %s", as.integer(fp_raw), nb(protocol)), collapse = ", "),
              .groups = "drop") %>%
    transmute(t = paste0(nb(set), ": ", t)) %>% pull(t) %>% paste(collapse = "; ")
  note <- paste(strwrap(sprintf("Validation removed the only false spike-in models (%s); precision = 1 in every set",
                                fp_txt), width = 32), collapse = "\n")
  ggplot(long, aes(protocol, value, colour = set, group = set)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    ## detected / expected at each recall point, the value at each F1 point
    geom_text(data = lab_d, aes(x = lx, y = ly, label = label, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, colour = lab_cols[as.character(lab_d$set)],
              size = 7 / .pt, family = font_family) +
    scale_colour_manual(values = set_cols, name = "Spike-in set") +
    facet_wrap(~ metric, nrow = 1, labeller = labeller(metric = m_lab)) +
    scale_x_discrete(expand = expansion(add = c(0.4, 0.7))) +       # room for the counts on the right
    ## from 0.25 (every value is >= 0.40), headroom above 1 for the long SIRV
    ## labels; the axis stops at 1
    scale_y_continuous(limits = c(0.25, 1.08), breaks = seq(0.25, 1, 0.25),
                       labels = c("0.25", "0.5", "0.75", "1"),
                       expand = expansion(mult = c(0.02, 0.01))) +
    ## note under the legend, as in panel k (a legend entry, so it follows the legend)
    guides(colour = guide_legend(order = 1),
           custom = guide_custom(grid::textGrob(note, x = 0, hjust = 0,
                                                gp = grid::gpar(fontsize = 8, col = "grey30",
                                                                fontfamily = font_family, lineheight = 1.19)),
                                 order = 2)) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
          panel.spacing.x = unit(3, "mm"),
          legend.justification = c(0, 0.5), legend.box.just = "left",
          legend.spacing.y = unit(4, "mm"))
}

## ------------------------------------------- b: spike-in identification ---
## expected / detected from `det` (detection rules above); false positives =
## novel models on the spike-in sequences (02: novel_models_on_spikeins.csv for
## SIRV / ERCC, Bambu output "raw" and after pulposeq validation; 01:
## sequin_identification.csv for sequins). E0 SIRVs sit on contigs SIRV1-SIRV7,
## each long SIRV on its own contig (SIRV4001, ...).
nm_sp <- read_csv(file.path(sirv_dir, "novel_models_on_spikeins.csv"),
                  col_types = cols(.default = "c")) %>%
  mutate(set = case_when(family == "ERCC" ~ "ERCC",
                         grepl("^SIRV[1-7]$", seqname) ~ "SIRV E0",
                         TRUE ~ "long SIRV"))
fp_of <- function(p, s, st) sum(nm_sp$protocol == p & nm_sp$set == s & nm_sp$stage == st)
id_sets <- c("SIRV E0", "long SIRV", "ERCC", "Sequins")
idf <- det %>%
  mutate(type = as.character(type), protocol = as.character(protocol)) %>%
  filter(type %in% c("SIRV E0", "long SIRV", "ERCC", "Sequin isoforms")) %>%
  transmute(set = if_else(type == "Sequin isoforms", "Sequins", type), protocol,
            detected = n, expected = total)
sid_c <- mutate(sid, protocol = as.character(protocol))
idf$fp_raw <- mapply(function(p, s) if (s == "Sequins") sid_c$novel_raw[sid_c$protocol == p][1]
                     else fp_of(p, s, "raw"), idf$protocol, idf$set)
idf$fp_val <- mapply(function(p, s) if (s == "Sequins") sid_c$novel_validated[sid_c$protocol == p][1]
                     else fp_of(p, s, "validated"), idf$protocol, idf$set)
id_metrics(idf) %>%
  mutate(set = factor(set, levels = id_sets), protocol = factor(protocol, levels = prot_levels)) %>%
  arrange(set, protocol) %>%
  write_csv(file.path(tab_dir, "spikein_identification_by_set.csv"))

pID <- id_lines(idf, q_spk_cols, prot_levels) +
  labs(title = "Spike-in identification",
       subtitle = "Detected: >= 3 reads in >= 2 libraries (ERCC: >= 1 read)") +
  border(0.25)

## ------------------------------------------------- c: length bias -------
e0_acc <- e0 %>% group_by(protocol) %>%
  summarise(within2 = mean(abs(median) <= 1), .groups = "drop")
## identification numbers (E0 / long SIRV detected, F1) are in panel b
b_lab <- e0_acc %>%
  mutate(lab = sprintf("%.0f%% within 2\u00d7", 100 * within2),
         ## empty corner per protocol: bottom-left for ONT cDNA and ONT DRS (long
         ## SIRVs sit low on the right), bottom-right for PacBio
         left  = protocol != "PacBio",
         x     = if_else(left, min(lb$len), max(lb$len)),   # data range: -Inf fails on a log axis
         hjust = if_else(left, 0, 1))

pB <- ggplot(lb, aes(len, log2_ratio)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black",
              fill = "grey85", linewidth = 0.5) +
  geom_point(aes(colour = set), size = 1.5) +
  geom_text(data = b_lab, aes(x = x, y = -Inf, label = lab, hjust = hjust),
            vjust = -0.25, lineheight = 1.46, size = 7 / .pt, inherit.aes = FALSE) +
  scale_x_log10(breaks = c(300, 1000, 3000, 10000), labels = c("0.3", "1", "3", "10")) +
  scale_colour_manual(values = set_cols) +
  facet_wrap(~ protocol, nrow = 1) +
  labs(x = "Transcript length (kb)", y = "log2(observed / expected)", colour = "SIRV set",
       title = "Length bias",
       subtitle = paste0("Equimolar SIRVs; 0 = expected abundance;\nwithin 2\u00d7 = E0 isoforms with median ",
                         "abundance within 2-fold of expected")) +
  theme(strip.text = element_text(size = 9.5)) + border(0.25) + leg_right(0.78)

## ---------------------------------------------------------- d: ERCC -------
## ERCCs detected / in the mix: panel b
c_lab <- se %>% arrange(protocol) %>%
  mutate(label = sprintf("%s \u00b7 R\u00b2 = %.2f, LOD = %.1f",
                         protocol, ercc_r2, ercc_lod70_amol_ul))

## plot and fit the same ERCCs: those passing the detection rule (02, `detected_rule`)
pD <- ggplot(filter(er, is.finite(log2_cpm), detected_rule), aes(log2(conc), log2_cpm, colour = protocol)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
  geom_point(size = 1.3, alpha = 0.8) +
  ## one line per protocol, 8 pt, 2.5 mm below the panel top, 3.9 mm apart
  corner_text(c_lab$label, unname(protocol_cols[as.character(c_lab$protocol)]),
              size_pt = 8, top_mm = 2.5, pitch_mm = 3.9) +
  scale_colour_manual(values = protocol_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.4))) +
  labs(x = "ERCC input (log2 amol/\u00b5l)", y = "log2 CPM (within ERCC)", colour = NULL,
       title = "ERCC dose response",
       subtitle = "R\u00b2; LOD (amol/\u00b5l)") +
  border(0.25) + leg_right(0.85)

## ------------------------------------------------- e: sequin DE -------
## sequin ground truth, isoform level (edgeR, 01_sequin_edger.R)
out_lev  <- c("true positive", "true negative", "false positive", "false negative", "wrong sign")
out_cols <- c("true positive" = "grey15", "true negative" = "grey75", "false positive" = "#D55E00",
              "false negative" = "#E69F00", "wrong sign" = "#CC79A7")  # as Fig. 2c
out_labs <- c(out_lev[-5], "wrong sign (opposite DE direction)")
y_cap    <- 9.5                    # |observed log2FC| drawn at most this far (triangles beyond)
sq_lfc   <- 0.5                    # truth threshold of 01 (shaded band)

fc <- map_dfr(platforms, function(p)
  read_csv(file.path(seq_dir, "edger", p, "transcript", "divided", "sequin_de.csv"), show_col_types = FALSE) %>%
    dplyr::select(id, exp_lfc, obs_lfc, pvalue, outcome) %>%
    mutate(platform = p)) %>%
  filter(!is.na(pvalue), is.finite(obs_lfc)) %>%                  # as summarise_calls()
  mutate(protocol = factor(plat_labels[platform], levels = prot_levels),
         outcome  = factor(outcome, levels = out_lev),
         capped   = abs(obs_lfc) > y_cap,
         y        = pmax(pmin(obs_lfc, y_cap), -y_cap))

## top-left: DE calls among the TESTED sequins (the points). Bottom-right: how
## many sequins could be tested (>= 10 reads in >= 2 libraries), then R2 / slope;
## that corner is empty by design (expected up, observed strongly down).
## Sequins detected / in the mix: panel b.
d_lab <- fc_sum %>%
  filter(level == "transcript") %>%
  mutate(lab_top = sprintf("precision = %.2f\nsensitivity = %.2f\nF1 = %.2f",
                           precision, sensitivity, f1),
         lab_bot = sprintf("%d tested\nR\u00b2 = %.2f\nslope = %.2f",
                           as.integer(tested), r2, slope))

pC <- ggplot(fc, aes(exp_lfc, y)) +
  annotate("rect", xmin = -sq_lfc, xmax = sq_lfc, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_point(aes(colour = outcome, shape = capped), size = 1.5, alpha = 0.9,
             show.legend = TRUE,
             position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  ## x at the data range with hjust 0 / 1: every line aligns flush
  geom_text(data = d_lab, aes(x = min(fc$exp_lfc), y = Inf, label = lab_top),
            hjust = 0, vjust = 1.1, lineheight = 1.46, size = 7 / .pt) +
  geom_text(data = d_lab, aes(x = max(fc$exp_lfc), y = -Inf, label = lab_bot),
            hjust = 1, vjust = -0.3, lineheight = 1.46, size = 7 / .pt) +
  scale_colour_manual(values = out_cols, labels = out_labs, drop = FALSE, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), breaks = "TRUE",
                     labels = sprintf("|observed log2FC| > %.1f", y_cap), name = NULL) +
  ## bottom room: the tested + R2 / slope block sits below the data
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.3))) +
  facet_wrap(~ protocol, nrow = 1) +
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = 2)),
         shape  = guide_legend(order = 2, override.aes = list(size = 2, shape = 2))) +
  labs(x = "Expected log2FC (Mix B / Mix A)", y = "Observed log2FC (edgeR)",
       title = "Sequin ground truth",
       subtitle = paste0("Mix B (SCLC-A) vs Mix A; edgeR glmTreat (H0 |log2FC| <= ", sq_lfc,
                         "), FDR < 0.05; shaded: not truly DE")) +
  theme(strip.text = element_text(size = 9.5), legend.spacing.y = unit(0, "mm")) +
  border(0.85) + leg_right(0.53)

## ------------------------------------------ f: quantifier (Bambu vs oarfish) ---
## the same runs quantified twice (same libraries and transcripts; Bambu on genome alignments, oarfish on its own transcriptome mapping):
## log2((oarfish + 1) / (Bambu + 1)) of each spike-in transcript's mean reads over
## libraries, against its length; transcripts with >= 1 mean read in either
## quantifier (03_quantifier_plots.R, per_transcript_quantifiers.csv.gz; the
## supplementary version is supp_q_a2). Labelled: E0 isoforms > 8-fold apart.
q_tx <- read_csv(file.path(spk_dir, "quantifier_comparison", "per_transcript_quantifiers.csv.gz"),
                 show_col_types = FALSE) %>%
  filter(category %in% names(q_spk_cols), pmax(bambu_mean, oarfish_mean) >= 1) %>%
  mutate(set      = factor(category, levels = names(q_spk_cols)),
         protocol = factor(protocol, levels = prot_levels))
q_lab <- filter(q_tx, set == "SIRV E0", abs(lr) > 3)

pQ <- ggplot(q_tx, aes(len, lr, colour = set)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 1.4, alpha = 0.8) +
  repel_labels(q_lab, aes(label = id), size = 7 / .pt, box.padding = 0.3,
               direction = "y", nudge_x = 0.35) +
  scale_colour_manual(values = q_spk_cols, name = NULL) +
  scale_x_log10(breaks = c(300, 1000, 3000, 10000), labels = c("0.3", "1", "3", "10")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(~ protocol, nrow = 1) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Transcript length (kb)", y = "log2(oarfish / Bambu reads)",
       title = "Where the quantifier matters",
       subtitle = "Each point is one spike-in transcript (mean over\nlibraries, +1); 0 = same reads") +
  theme(strip.text = element_text(size = 9.5)) + border(0.25) + leg_right(0.5)

## ---------------------------------------------------------------- g -------
## read-level support (06_read_support.R)
prot_map <- c("LongBench ONT cDNA" = "ONT cDNA", "LongBench ONT DRS" = "ONT DRS",
              "LongBench PacBio" = "PacBio")
rs <- read_csv(file.path(rs_dir, "read_support_all_runs_long.csv"), show_col_types = FALSE) %>%
  filter(run %in% names(prot_map)) %>%
  mutate(protocol = factor(prot_map[run], levels = prot_levels))

## exact numbers and % per protocol x biotype x support class: a table
rs %>%
  transmute(Protocol = protocol, Biotype = biotype, `Isoforms (total)` = total,
            Support = support, cell = sprintf("%d (%.1f%%)", n, pct)) %>%
  pivot_wider(names_from = Support, values_from = cell) %>%
  arrange(Protocol) %>%
  write_csv(file.path(tab_dir, "table_read_support.csv"))

## line plot in the style of panel b (the bar version is drawn for Fig. 2a by
## read_support_plot()). One facet per
## biotype - annotated in the left column, novel in the right - x = protocol,
## y = % of the biotype's isoforms, one line per read-support class (the bar
## colours). The biotype's isoforms in each protocol are printed above that
## protocol's highest dot (below it for the novel biotypes, near 100%).
k_lab <- function(x) ifelse(x >= 9950, paste0(round(x / 1000), "k"),   # 9,950+ -> "10k", not "10.0k"
                     ifelse(x >= 1000, sprintf("%.1fk", x / 1000), as.character(x)))   # 1.7k, 227
rs2 <- rs %>%
  filter(biotype %in% names(bt_disp)) %>%
  mutate(bt = factor(bt_disp[biotype], levels = bt_disp),
         support = factor(support, levels = sup_lv))
## counts above the highest dot; below it when that dot is near 100% (novel
## biotypes), where the facet title would crowd it
bt_n <- rs2 %>% group_by(bt, protocol) %>%
  summarise(top = max(pct), lab = k_lab(first(total)), .groups = "drop") %>%
  mutate(vj = if_else(top > 85, 1.9, -0.9))
pE <- ggplot(rs2, aes(protocol, pct, colour = support, group = support)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  geom_text(data = bt_n, aes(protocol, top, label = lab, vjust = vj), inherit.aes = FALSE,
            size = 7 / .pt, family = font_family, colour = "grey30") +
  facet_wrap(~ bt, ncol = 2, dir = "v") +
  scale_colour_manual(values = sup_cols, labels = sup_lab, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100), labels = c("0", "50%", "100%"),
                     expand = expansion(mult = 0.06)) +
  labs(x = NULL, y = "Share of isoforms", title = "Read-level support across biotypes",
       subtitle = "Numbers: isoforms of the biotype in each protocol") +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        panel.spacing = unit(2, "mm"), legend.key.spacing.y = unit(1.5, "mm")) +
  border(0.25) + leg_right(0.5)

## ============================================ novel transcript data ======
cc_group <- function(cls) case_when(
  grepl("\\(j\\)", cls)                   ~ "junction match (j)",
  grepl("\\(k\\)", cls)                   ~ "contains reference (k)",
  grepl("retained intron", cls)           ~ "retained intron (m, n)",
  grepl("^antisense", cls)                ~ "antisense (i, x)",
  grepl("\\(u\\)", cls)                   ~ "intergenic (u)",
  TRUE                                    ~ "sense overlap (i, o, y)")

runs_k <- c(ONTc = "ont_cdna", ONTd = "ont_drs", PB = "pb")
nov <- imap_dfr(runs_k, function(r, k)
  read_csv(file.path(run_dir[[r]], pp$novel_meta),
           show_col_types = FALSE) %>%
    transmute(protocol = prot[[k]], len, biotype = pulposeq_transcript_biotype,
              cls = cc_group(classification))) %>%
  mutate(protocol = factor(protocol, levels = prot_levels),
         biotype  = factor(biotype, levels = c("novel_lncRNA", "novel_protein_coding", "novel_non_coding")),
         cls      = factor(cls, levels = names(cc_cols)))

## class-code tables (a supplementary table)
nov %>% dplyr::count(protocol, biotype, cls) %>%
  group_by(protocol) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup() %>%
  write_csv(file.path(tab_dir, "novel_class_codes.csv"))
bt_lab <- c(novel_lncRNA = "lncRNA", novel_protein_coding = "protein-coding",
            novel_non_coding = "non-coding")
nov %>%
  dplyr::count(protocol, biotype, cls, .drop = FALSE) %>%
  group_by(protocol) %>% mutate(`% of novel transcripts` = 100 * n / sum(n)) %>%
  group_by(protocol, biotype) %>%
  mutate(`% within biotype` = 100 * n / na_if(sum(n), 0)) %>%   # NA when a biotype has no novel transcripts
  ungroup() %>%
  transmute(Dataset = "LongBench", Protocol = protocol, Biotype = bt_lab[as.character(biotype)],
            `Class code (GffCompare)` = cls, `Novel transcripts` = n,
            `% of novel transcripts`, `% within biotype`) %>%
  write_csv(file.path(tab_dir, "supp_novel_biotype_classcode_longbench.csv"))

## ---------------------------------------------------------------- h -------
f_lab <- nov %>% group_by(protocol) %>% summarise(m = median(len), .groups = "drop") %>%
  mutate(label = sprintf("%s = %.1f kb", protocol, m / 1000))

pF <- ggplot(nov, aes(len, colour = protocol)) +
  geom_density(linewidth = 0.5, show.legend = FALSE) +
  ## medians top-right, 8 pt, 3 mm below the panel top, 4.35 mm apart
  corner_text(f_lab$label, unname(protocol_cols[as.character(f_lab$protocol)]),
              size_pt = 8, top_mm = 3, pitch_mm = 4.35, left = FALSE) +
  scale_colour_manual(values = protocol_cols) +
  scale_x_log10(breaks = c(500, 1000, 3000, 10000), labels = c("0.5", "1", "3", "10")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(x = "Spliced length (kb)", y = "Density", title = "Novel transcripts' length", subtitle = "Medians") +
  border(0.25)

## ---------------------------------------------------------------- i -------
u <- read_tsv(file.path(ujc_dir, "ujc_table.tsv"), show_col_types = FALSE) %>%
  filter(structure == "multi-exon") %>%
  mutate(across(starts_with("present_"), ~ as.logical(as.character(.x))))

nG <- map_dfr(names(prot), function(k) {
  ux <- u[u[[paste0("present_", k)]] %in% TRUE, ]
  tibble(protocol = prot[[k]], n_protocols = ux$n_protocols)
}) %>%
  mutate(share = case_when(n_protocols == 3 ~ "all 3 protocols",
                           n_protocols == 2 ~ "1 other protocol",
                           TRUE             ~ "this protocol only"),
         share    = factor(share, levels = rev(names(share_cols))),   # "all 3" at 0 %
         protocol = factor(protocol, levels = prot_levels)) %>%
  dplyr::count(protocol, share) %>%
  group_by(protocol) %>% mutate(pct = 100 * n / sum(n), n_tot = sum(n)) %>% ungroup()
g_lab <- nG %>% dplyr::distinct(protocol, n_tot) %>%
  { setNames(sprintf("%s\n(%s)", .$protocol, formatC(.$n_tot, big.mark = ",", format = "d")), as.character(.$protocol)) }

pG <- ggplot(nG, aes(pct, protocol, fill = share)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.0f%%", pct), colour = share == "all 3 protocols"),
            position = position_stack(vjust = 0.5), size = 7 / .pt, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_manual(values = share_cols, breaks = names(share_cols), name = NULL) +
  scale_y_discrete(limits = rev(prot_levels), labels = g_lab) +           # ONT cDNA on top
  scale_x_continuous(breaks = c(0, 50, 100), labels = c("0", "50%", "100%"), expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Share of isoforms", y = NULL,
       title = "UIC sharing", subtitle = "UIC = Unique Intron Chain") +
  leg_right(0.5) + border(0.21) +
  theme(axis.text.y = element_text(size = 7, colour = "black"),
        legend.key.size = unit(3.9, "mm"))

## ---------------------------------------------------------------- j -------
bin_levels <- c("<0.5 kb", "0.5-1 kb", "1-2 kb", "2-5 kb", ">5 kb")
cb <- read_csv(file.path(ujc_dir, "length/corroboration_by_length.csv"), show_col_types = FALSE) %>%
  filter(n >= 10) %>%
  mutate(bin = factor(bin, levels = bin_levels),
         from = factor(from, levels = prot_levels),
         found_in = factor(found_in, levels = prot_levels))

## j: n printed at each point, in a darker protocol shade.
## Label side: if the two points of a bin are < 8 % apart, the higher one's label
## goes above and the lower one's below; otherwise the side (above / below) that
## lies furthest, in angle, from the point's own line segments
## (approx. mm per unit: 11.8 per bin on x, 0.78 per % on y, at 160 x 80 mm)
j_lab_cols <- vapply(protocol_cols, function(cl) colorRampPalette(c(cl, "black"))(10)[3], "")
ang_gap <- function(a, b) abs(((a - b + 180) %% 360) - 180)   # NA-safe below
seg_ang <- function(dx, dy) atan2(dy * 0.78, dx * 11.8) * 180 / pi
j_lab <- cb %>%
  mutate(xi = as.integer(bin)) %>%
  group_by(from, found_in) %>% arrange(xi, .by_group = TRUE) %>%
  mutate(a_prev = seg_ang(lag(xi) - xi, lag(pct_found) - pct_found),
         a_next = seg_ang(lead(xi) - xi, lead(pct_found) - pct_found)) %>%
  ungroup() %>%
  mutate(up_gap   = pmin(ang_gap(90, a_prev), ang_gap(90, a_next), na.rm = TRUE),
         down_gap = pmin(ang_gap(-90, a_prev), ang_gap(-90, a_next), na.rm = TRUE)) %>%
  group_by(from, bin) %>%
  mutate(close = n() > 1 && diff(range(pct_found)) < 8,
         above = if_else(close, pct_found >= max(pct_found), up_gap >= down_gap)) %>%
  ungroup() %>%
  mutate(ly = if_else(above, pct_found + 3.5, pct_found - 3.5),
         vj = if_else(above, 0, 1))
j_note <- "Numbers at the points =\ntranscripts in the bin (n)"

pH <- ggplot(cb, aes(bin, pct_found, colour = found_in, group = found_in)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  geom_text(data = j_lab, aes(x = bin, y = ly, label = formatC(n, big.mark = ",", format = "d"), vjust = vj),
            inherit.aes = FALSE, colour = j_lab_cols[as.character(j_lab$found_in)],
            size = 7 / .pt, family = font_family) +
  scale_colour_manual(values = protocol_cols) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +            # room below 0 for the labels
  scale_x_discrete(labels = function(x) sub(" kb", "", x)) +
  facet_wrap(~ from, nrow = 1) +
  guides(colour = guide_legend(order = 1)) +
  labs(x = "Length in the source protocol (kb)", y = "% with UIC also in",
       colour = "Nº of isoforms sharing\nUICs within",
       title = "UIC sharing by length",
       subtitle = "Facet = source protocol; point size = n;\nbins with n < 10 omitted") +
  leg_right(0.5) + border(0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ---------------------------------------------------------------- k -------
## annotation agreement among shared UICs (table + note in the panel)
cj <- read_tsv(file.path(ujc_dir, "summary_disagreements.tsv"), show_col_types = FALSE) %>%
  filter(category == "ALL_shared") %>%
  mutate(attribute = recode(attribute, pulposeq_transcript_biotype = "biotype",
                            class_code = "class code", host_gene = "host gene"),
         pct_conc = 100 - pct_discordant)
cj %>%
  transmute(Attribute = attribute, `Shared UICs compared` = n_compared,
            Concordant = n_concordant, `% concordant` = round(pct_conc, 1)) %>%
  write_csv(file.path(tab_dir, "table_ujc_annotation_agreement.csv"))
conc <- setNames(cj$pct_conc, cj$attribute)
i_note <- sprintf("%.1f%% shared UICs agree\non biotype, %.1f%% on class\ncode and %.1f%% on host gene",
                  conc[["biotype"]], conc[["class code"]], conc[["host gene"]])

## ECDF of |d5| / |d3| per protocol pair (numbers = protocol pairs)
pair_cols <- c(shared_all  = "grey15",            "ONTd+ONTc" = npg_cols[["red"]],
               "ONTd+PB"   = npg_cols[["brown"]], "ONTc+PB"   = npg_cols[["greyblue"]])
pair_name <- c(shared_all  = "all 3 protocols",   "ONTd+ONTc" = "ONT DRS + cDNA",
               "ONTd+PB"   = "ONT DRS + PacBio",  "ONTc+PB"   = "ONT cDNA + PacBio")
pd <- read_tsv(file.path(ujc_dir, "ujc_pair_deltas.tsv"), show_col_types = FALSE) %>%
  filter(category %in% names(pair_cols))
i_n   <- dplyr::count(pd, category, name = "n_pair")
i_lab <- with(i_n, setNames(sprintf("%s (%s)", pair_name[category], formatC(n_pair, big.mark = ",", format = "d")), category))

ends <- bind_rows(transmute(pd, category, end = "5' end |d5|", dist = abs_d5),
                  transmute(pd, category, end = "3' end |d3|", dist = abs_d3)) %>%
  filter(!is.na(dist)) %>%
  mutate(end = factor(end, levels = c("5' end |d5|", "3' end |d3|")),
         category = factor(category, levels = names(pair_cols)))

brk <- c(0, 10, 100, 1000)                                # plotted at dist + 1 on a log axis
pI <- ggplot(ends, aes(dist + 1, colour = category)) +
  geom_vline(xintercept = c(50, 500) + 1, colour = "grey85", linewidth = 0.3) +
  stat_ecdf(geom = "step", linewidth = 0.5, pad = FALSE) +
  scale_colour_manual(values = pair_cols, labels = i_lab) +
  scale_x_log10(breaks = brk + 1, labels = c("0", "10", "100", "1k")) +
  scale_y_continuous(labels = function(y) paste0(round(100 * y), "%"),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  facet_wrap(~ end, nrow = 1) +
  labs(x = "Distance between the two protocols' ends (bp)",
       y = "Cumulative % of protocol pairs", colour = "Protocol pair",
       title = "Transcript ends within an identical UIC",
       subtitle = "Share of protocol pairs whose ends lie within a distance;\ngrey lines = 50 and 500 bp") +
  ## agreement note under the legend (8 pt, grey30, line height 1.19): a legend
  ## entry, so it follows the legend at any panel size
  guides(colour = guide_legend(order = 1),
         custom = guide_custom(grid::textGrob(i_note, x = 0, hjust = 0,
                                              gp = grid::gpar(fontsize = 8, col = "grey30",
                                                              fontfamily = font_family, lineheight = 1.19)),
                               order = 2)) +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        panel.spacing.x = unit(3, "mm"),
        legend.justification = c(0, 1),
        legend.box.just = "left",
        legend.spacing.y = unit(4, "mm")) +
  border(0.25)

## ======================================================= save panels ======
## one file per panel; size (width, height, mm) = the panel's extent on the
## final Fig. 1 page (340 x 485 mm, 5 mm margins, 4 mm gutters; assembled by
## 09_assemble_figures.R):
##   row 1 (137)   a 330 (workflow, full width; drawn in the vector editor)
##   row 2 (70)    b 165 | c 161
##   row 3 (72)    d 100 | e 226
##   row 4 (100)   f 117 | g 105 | h 100 x 48 above i 100 x 48
##   row 5 (80)    j 160 | k 166
panels <- list(
  fig1b_spikein_identification = list(p = pID, size = c(165.0,  70.0)),
  fig1c_length_bias            = list(p = pB,  size = c(161.0,  70.0)),
  fig1d_ercc                   = list(p = pD,  size = c(100.0,  72.0)),
  fig1e_sequin_isoform         = list(p = pC,  size = c(226.0,  72.0)),
  fig1f_quantifier             = list(p = pQ,  size = c(117.0, 100.0)),
  fig1g_read_support           = list(p = pE,  size = c(105.0, 100.0)),
  fig1h_novel_length           = list(p = pF,  size = c(100.0,  48.0)),
  fig1i_uic_sharing            = list(p = pG,  size = c(100.0,  48.0)),
  fig1j_uic_sharing_by_length  = list(p = pH,  size = c(160.0,  80.0)),
  fig1k_transcript_ends        = list(p = pI,  size = c(166.0,  80.0))
)
for (nm in names(panels)) save_panel(panels[[nm]]$p, nm, panels[[nm]]$size, panel_dir)
## figure data behind every panel, kept for sharing
source(file.path(repo, "downstream/helpers/source_data.R"))
save_source_data(panels, file.path(repo, "downstream/longbench/figure-panels/results/source_data_fig1.rds"))
## the same panels without title and subtitle (same file names and sizes), in
## figures/fig1_panels_no_titles/
panel_dir_nt <- paste0(panel_dir, "_no_titles")
for (nm in names(panels))
  save_panel(panels[[nm]]$p + labs(title = NULL, subtitle = NULL), nm, panels[[nm]]$size, panel_dir_nt)
message("Done: ", length(panels), " panels in ", panel_dir, "; tables in ", tab_dir)
