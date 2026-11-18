## ===========================================================================
## 10_quantifier_panels.R - Dong: Bambu vs oarfish, and the read-to-transcript
##                          ambiguity (RTA) overdispersion, on the SAME pulposeq
##                          run (same libraries and transcripts; Bambu on genome alignments, oarfish on its own transcriptome mapping)
##                          -> Extended Data figure panels + tables
## The Dong version of longbench/03_quantifier_plots.R. Dong has
## sequins only (no SIRV / ERCC), so the spike-in panels (LongBench a, b) are
## replaced by sequin panels (a, b); c-h are the LongBench panels for one run.
##
##   a  sequins: mean reads per isoform, Bambu vs oarfish (colour = RTA overdispersion)
##   b  sequin isoform DE: observed vs expected log2FC, Bambu | oarfish
##      (edgeR, undivided counts for both, so the RTA division is not part of the difference)
##   c  human transcripts: mean counts, Bambu vs oarfish (2-D histogram)
##   d  detection agreement (Jaccard) by detection threshold (1-100 reads)
##   e  sequins: observed log2FC with divided vs undivided oarfish counts
##   f  RTA overdispersion by transcript category (expressed transcripts)
##   g  RTA overdispersion of human transcripts by isoforms per gene in the index
##   h  disagreement (> 2-fold, Bambu vs oarfish) by RTA overdispersion
##
## Inputs  (read only): the Dong run (curated counts, metadata, oarfish .quant,
##          Bambu table), RTA overdispersion from read_oarfish(curated = FALSE)
##          (cache in diff-exp/results/oarfish_cache/), sequins/results{,_bambu}/
## Outputs
##   figures/ed_quantifier_panels/ed_q_{a..h}_*.{pdf,svg,png}
##   figures/ed_quantifier_panels_no_titles/   (the same, no title / subtitle)
##   results/quantifier_comparison/  per_transcript_quantifiers.csv.gz,
##     jaccard_per_sample.csv, sequin_quantifiers.csv, overdispersion_*.csv,
##     human_agreement.csv, disagreement_by_overdispersion.csv
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")
source(file.path(helper_dir, "plot_theme.R"))              # font_family, base_theme, npg_cols
source(file.path(helper_dir, "figure_helpers.R"))          # save_panel
source(file.path(shared_helper_dir, "oarfish_counts.R"))   # read_oarfish, strip_bam
source(file.path(shared_helper_dir, "spikein_counts.R"))   # quant_dir

outroot  <- file.path(res_dir, "quantifier_comparison")
figdir   <- file.path(fig_dir, "ed_quantifier_panels")
figdir_nt <- paste0(figdir, "_no_titles")
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_family = font_family, base_size = 8) + base_theme +
  theme(plot.title          = element_text(size = 10, face = "plain"),
        plot.subtitle       = element_text(size = 8, colour = "grey30"),
        plot.title.position = "plot",
        strip.text          = element_text(size = 8.5, face = "plain"),
        strip.background    = element_rect(fill = "grey94"),
        panel.grid.minor    = element_blank(),
        legend.position     = "right",
        legend.key.size     = unit(3, "mm"),
        plot.margin         = margin(1, 1, 1, 1, "mm")))
update_geom_defaults("text", list(family = font_family, size = 7 / .pt))

run_col <- c("ONT cDNA" = "#27408B")                      # protocol_cols[["ONT cDNA"]]
cat_levels <- c("protein-coding", "lncRNA", "other annotated",
                "novel protein-coding", "novel lncRNA", "novel non-coding", "sequin")
## the biotype palette of every figure (as make_palette in 06_isoformic.R);
## "other" in light grey
cat_cols <- c("protein-coding" = "#A6D854", "lncRNA" = "#E78AC3", "other annotated" = "#C9C9C9",
              "novel protein-coding" = "#238B45", "novel lncRNA" = "#AE017E",
              "novel non-coding" = "#6A51A3", "sequin" = "#91D1C2")   # LongBench cat_cols + spk_cols
spike_re <- "^(R[12]_|SIRV|ERCC)"
od_lv    <- c("1-1.05", "1.05-1.25", "1.25-1.5", "1.5-2", "2-3", "> 3")
od_cut   <- function(x) cut(x, c(0, 1.05, 1.25, 1.5, 2, 3, Inf), labels = od_lv)

## ===================================================== per-transcript data ==
cat_of <- function(id, gene_bt, novel_bt) case_when(
  grepl("^R[12]_", id)                    ~ "sequin",
  novel_bt == "novel_lncRNA"              ~ "novel lncRNA",
  novel_bt == "novel_protein_coding"      ~ "novel protein-coding",
  novel_bt == "novel_non_coding"          ~ "novel non-coding",
  gene_bt  == "protein_coding"            ~ "protein-coding",
  gene_bt  == "lncRNA"                    ~ "lncRNA",
  TRUE                                    ~ "other annotated")

read_mat <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(x[, -(1:2)])
  rownames(m) <- x$TXNAME
  colnames(m) <- strip_bam(colnames(m))
  list(m = m, gene = setNames(x$GENEID, x$TXNAME))
}

o <- read_mat(run_path("oarfish_tx"))
b <- read_mat(run_path("counts_tx"))
ids <- rownames(o$m); smp <- colnames(o$m)
stopifnot(all(ids %in% rownames(b$m)), all(smp %in% colnames(b$m)))
om <- o$m; bm <- b$m[ids, smp, drop = FALSE]

of <- read_oarfish(data_dir, "transcript", curated = FALSE, cache_dir = oarfish_cache)
od <- unname(of$overdispersion[ids])

q1  <- list.files(file.path(data_dir, "transcriptome", "oarfish"), "[.]quant$", full.names = TRUE)[1]
qq  <- read_tsv(q1, show_col_types = FALSE, col_select = c(tname, len))
len <- setNames(qq$len, qq$tname)
full  <- read_tsv(run_path("bambu_counts_tx"), col_select = 1:2, show_col_types = FALSE)
n_iso <- table(full[[2]])
ann <- read_csv(run_path("ann_meta"), show_col_types = FALSE, col_select = c(transcript_id_version, gene_biotype))
nov <- read_csv(run_path("novel_meta"), show_col_types = FALSE, col_select = c(qry_id, pulposeq_transcript_biotype))

tx <- tibble(id             = ids,
             gene           = unname(o$gene[ids]),
             len            = unname(len[ids]),
             n_iso          = as.integer(n_iso[unname(o$gene[ids])]),
             gene_bt        = ann$gene_biotype[match(ids, ann$transcript_id_version)],
             novel_bt       = nov$pulposeq_transcript_biotype[match(ids, nov$qry_id)],
             overdispersion = od,
             bambu_mean     = rowMeans(bm),
             oarfish_mean   = rowMeans(om),
             libs_ge10      = rowSums(om >= 10)) %>%
  mutate(category = factor(cat_of(id, gene_bt, novel_bt), levels = cat_levels),
         lr       = log2((oarfish_mean + 1) / (bambu_mean + 1)))          # oarfish vs Bambu
write_csv(tx, file.path(outroot, "per_transcript_quantifiers.csv.gz"))

hum_rows <- !grepl(spike_re, ids)
jac <- map_dfr(c(1, 2, 5, 10, 20, 50, 100), function(th) tibble(          # detection thresholds (reads)
  threshold = th, sample = smp,
  n_bambu   = map_dbl(smp, ~ sum(bm[hum_rows, .x] >= th)),
  n_oarfish = map_dbl(smp, ~ sum(om[hum_rows, .x] >= th)),
  jaccard   = map_dbl(smp, function(s) {
    a <- bm[hum_rows, s] >= th; z <- om[hum_rows, s] >= th
    sum(a & z) / sum(a | z)
  })))
write_csv(jac, file.path(outroot, "jaccard_per_sample.csv"))

## ================================================================ A: sequins ==
## a: mean reads per sequin isoform (every isoform in the mix with reads in either)
data("RnaQuinIsoformMixture", package = "Anaquin")
seq_ids <- as.character(RnaQuinIsoformMixture$Name)
sq <- tx %>% filter(category == "sequin") %>%
  transmute(id, len, overdispersion, bambu_mean, oarfish_mean, lr) %>%
  filter(id %in% seq_ids, pmax(bambu_mean, oarfish_mean) >= 1)
write_csv(sq, file.path(outroot, "sequin_quantifiers.csv"))
a_lab <- sprintf("Spearman = %.3f\nmedian log2 ratio = %.2f\n%d sequin isoforms",   # ratio = oarfish / Bambu
                 cor(sq$bambu_mean, sq$oarfish_mean, method = "spearman"), median(sq$lr), nrow(sq))
p_a <- ggplot(sq, aes(log10(bambu_mean + 1), log10(oarfish_mean + 1))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = overdispersion), size = 1.1) +
  ## bottom right: the points run along the diagonal, so this corner is empty
  annotate("text", x = Inf, y = -Inf, label = a_lab, hjust = 1.03, vjust = -0.3, lineheight = 1.2) +
  scale_colour_viridis_c(name = "RTA\noverdispersion", option = "C", end = 0.9) +
  # (no coord_equal: the 45 mm box is fixed and square, see the save section)
  labs(x = "Bambu, log10(mean reads + 1)", y = "oarfish, log10(mean reads + 1)",
       title = "Sequins: reads per isoform, Bambu vs oarfish",
       subtitle = "Mean over libraries; dashed: identity")

## a: the Fig. 1f design for the Dong sequins (the scatter above is also
## written, as ed_q_a0_sequin_counts):
## log2((oarfish + 1) / (Bambu + 1)) of each sequin isoform's mean reads against
## its length; labelled: sequins > 2-fold apart
k_lab <- function(x) ifelse(x >= 9950, paste0(round(x / 1000), "k"),   # 9,950+ -> "10k", not "10.0k"
                     ifelse(x >= 1000, sprintf("%.1fk", x / 1000), as.character(x)))   # 1.7k, 227
sq_lab <- filter(sq, abs(lr) > 1)
p_a_len <- ggplot(sq, aes(len, lr)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(colour = cat_cols[["sequin"]], size = 1.4, alpha = 0.8) +
  ## labels kept above the lowest point, so they stay
  ## out of the note band below the data (R2_24_1 ran into the median note)
  ggrepel::geom_text_repel(data = sq_lab, aes(label = id), size = 7 / .pt, family = font_family,
                           box.padding = 0.3, min.segment.length = 0, segment.size = 0.2,
                           ylim = c(min(sq$lr) - 0.15, NA), max.overlaps = Inf, seed = 1) +
  ## note bottom left, two lines: the top left
  ## clashed with the labelled sequins
  annotate("text", x = min(sq$len), y = -Inf, hjust = -0.05, vjust = -0.4, size = 7 / .pt, lineheight = 1.1,
           label = sprintf("median = %.2f\n%d isoforms", median(sq$lr), nrow(sq))) +
  scale_x_log10(breaks = c(300, 1000, 3000), labels = c("0.3", "1", "3")) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.15))) +            # room at the bottom for the note
  labs(x = "Transcript length (kb)", y = "log2(oarfish /\nBambu reads)",
       title = "Sequins: where the quantifier matters",
       subtitle = "Each point is one sequin isoform (mean over\nlibraries, +1); 0 = same reads")

## b: sequin isoform DE, observed vs expected log2FC, both quantifiers (undivided)
qn <- c(oarfish = "oarfish", bambu = "Bambu")
de_q <- map_dfr(names(qn), ~ read_csv(file.path(quant_dir(seq_dir, .x), "sequin_de_results.csv"),
                                      show_col_types = FALSE)) %>%
  filter(level == "transcript", count_mode == "undivided", is.finite(obs_lfc)) %>%
  mutate(quantifier = factor(qn[quantifier], levels = c("Bambu", "oarfish")),
         outcome = factor(outcome, levels = c("true positive", "true negative", "false positive",
                                              "false negative", "wrong sign")))  # order as Fig. 1e
sum_q <- map_dfr(names(qn), ~ read_csv(file.path(quant_dir(seq_dir, .x), "summary.csv"),
                                       show_col_types = FALSE)) %>%
  filter(level == "transcript")
## each quantifier's sequin DE metrics are printed in its facet: precision / sensitivity / F1 top left,
## tested / R2 / slope bottom right (undivided counts, as the points); the
## oarfish facet adds the F1 of the divided counts (the analysis of Fig. 2)
f1_div <- sum_q %>% filter(quantifier == "oarfish", count_mode == "divided") %>% pull(f1)
sum_q <- sum_q %>% filter(count_mode == "undivided") %>%
  mutate(lab = sprintf("precision = %.2f\nsensitivity = %.2f\nF1 = %.2f%s", precision, sensitivity, f1,
                       ifelse(quantifier == "oarfish", sprintf("\nF1 divided = %.2f", f1_div), "")),
         lab_br = sprintf("%d tested\nR\u00b2 = %.2f\nslope = %.2f", tested, r2, slope),
         quantifier = factor(qn[quantifier], levels = c("Bambu", "oarfish")))
## one outcome palette for every sequin DE panel: the colours of Fig. 1e
## (longbench/08_fig1_panels.R)
out_cols <- c("true positive" = "grey15", "true negative" = "grey75", "false positive" = "#D55E00",
              "false negative" = "#E69F00", "wrong sign" = "#CC79A7")
p_b <- ggplot(de_q, aes(exp_lfc, obs_lfc)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = outcome), size = 1, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  geom_text(data = sum_q, aes(x = min(de_q$exp_lfc), y = Inf, label = lab), hjust = 0, vjust = 1.2,
            lineheight = 1.05, inherit.aes = FALSE) +
  geom_text(data = sum_q, aes(x = max(de_q$exp_lfc), y = -Inf, label = lab_br), hjust = 1, vjust = -0.25,
            lineheight = 1.2, inherit.aes = FALSE) +
  scale_colour_manual(values = out_cols, name = NULL, drop = FALSE) +
  ## room above the data for the 4-line note top left
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  facet_wrap(~ quantifier, nrow = 1) +
  labs(x = "Expected log2FC (H1975 / HCC827)", y = "Observed log2FC",
       title = "Sequin isoform DE, Bambu vs oarfish",
       subtitle = "edgeR on undivided counts for both; dashed: identity")

## ===================================================== B: the RTA division ==
de_o <- read_csv(file.path(seq_dir, "sequin_de_results.csv"), show_col_types = FALSE) %>%
  filter(level == "transcript")
e_dat <- inner_join(
  de_o %>% filter(count_mode == "divided") %>%
    transmute(id, divided = obs_lfc, od = overdispersion, out_div = outcome),
  de_o %>% filter(count_mode == "undivided") %>% transmute(id, undivided = obs_lfc, out_undiv = outcome),
  by = "id") %>%
  filter(is.finite(divided), is.finite(undivided)) %>%
  mutate(call = if_else(out_div == out_undiv, "same call", "call changed"))
f1o <- read_csv(file.path(seq_dir, "summary.csv"), show_col_types = FALSE) %>% filter(level == "transcript")
e_lab <- sprintf("F1 divided = %.3f\nF1 undivided = %.3f",            # as LongBench ED4 e
                 f1o$f1[f1o$count_mode == "divided"], f1o$f1[f1o$count_mode == "undivided"])
p_e <- ggplot(e_dat, aes(undivided, divided)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = od, shape = call), size = 1.1) +
  annotate("text", x = -Inf, y = Inf, label = e_lab, hjust = -0.05, vjust = 1.2, lineheight = 1.2) +
  scale_colour_viridis_c(name = "RTA\noverdispersion", option = "C", end = 0.9) +
  scale_shape_manual(values = c("same call" = 16, "call changed" = 4), name = NULL, drop = FALSE) +
  labs(x = "Observed log2FC, undivided counts", y = "Observed log2FC, divided counts",
       title = "Sequins: the RTA division changes almost nothing",
       subtitle = "Isoform level; each point is one tested sequin; dashed: identity")

## f: overdispersion by category, expressed transcripts (>= 10 reads in >= 2 libraries)
od_x <- tx %>% filter(libs_ge10 >= 2, is.finite(overdispersion))
od_sum <- od_x %>% group_by(category) %>%
  summarise(n = n(), median = median(overdispersion),
            q25 = quantile(overdispersion, 0.25), q75 = quantile(overdispersion, 0.75),
            q90 = quantile(overdispersion, 0.90),
            pct_gt_1.5 = 100 * mean(overdispersion > 1.5), pct_gt_2 = 100 * mean(overdispersion > 2),
            .groups = "drop")
write_csv(od_sum, file.path(outroot, "overdispersion_by_category.csv"))
print(od_sum, width = Inf)
od_share <- od_x %>% mutate(od_bin = od_cut(overdispersion)) %>%
  dplyr::count(category, od_bin, .drop = FALSE) %>%
  group_by(category) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  filter(is.finite(pct)) %>% mutate(spike = category == "sequin")
write_csv(od_share, file.path(outroot, "overdispersion_share_by_category.csv"))
p_f <- ggplot(od_share, aes(od_bin, pct, colour = category, group = category)) +
  geom_line(aes(linetype = spike), linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = cat_cols, breaks = cat_levels, name = NULL) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = ifelse(cat_levels == "sequin", "22", "solid")))) +
  labs(x = "RTA overdispersion", y = "% of the category's\ntranscripts",   # 2 lines: the box is small
       title = "RTA overdispersion by transcript category",
       subtitle = "Expressed transcripts (>= 10 reads in >= 2 libraries); dashed = sequins") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

## g: human transcripts, by isoforms per gene in the whole index
iso_lv <- c("1", "2-3", "4-9", "10-19", ">= 20")
od_iso_x <- od_x %>% filter(!grepl(spike_re, id), !is.na(n_iso)) %>%
  mutate(iso_bin = cut(n_iso, c(0, 1, 3, 9, 19, Inf), labels = iso_lv), od_bin = od_cut(overdispersion))
write_csv(od_iso_x %>% group_by(iso_bin) %>%
            summarise(n = n(), median = median(overdispersion), q25 = quantile(overdispersion, 0.25),
                      q75 = quantile(overdispersion, 0.75), .groups = "drop"),
          file.path(outroot, "overdispersion_by_isoforms.csv"))
od_iso_share <- od_iso_x %>% dplyr::count(iso_bin, od_bin, .drop = FALSE) %>%
  group_by(iso_bin) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% filter(is.finite(pct))
write_csv(od_iso_share, file.path(outroot, "overdispersion_share_by_isoforms.csv"))
iso_cols <- setNames(c("#B8DE29", "#55C667", "#1F968B", "#39568C", "#440154"), iso_lv)
p_g <- ggplot(od_iso_share, aes(od_bin, pct, colour = iso_bin, group = iso_bin)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = iso_cols, name = "Isoforms per gene") +
  labs(x = "RTA overdispersion", y = "% of the group's\ntranscripts",
       title = "Human transcripts: ambiguity grows with isoforms per gene",
       subtitle = "Expressed reference and novel transcripts; isoforms per gene in the oarfish index") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

## ====================================== C: human transcripts, Bambu vs oarfish ==
hum <- tx %>% filter(!grepl(spike_re, id))
any1 <- pmax(hum$bambu_mean, hum$oarfish_mean) >= 1
ge10 <- pmax(hum$bambu_mean, hum$oarfish_mean) >= 10
agree <- tibble(n_any = sum(any1),
                spearman_any = cor(hum$bambu_mean[any1], hum$oarfish_mean[any1], method = "spearman"),
                n_ge10 = sum(ge10),
                pct_within_2fold_ge10 = 100 * mean(abs(hum$lr[ge10]) <= 1),
                median_abs_lr_ge10 = median(abs(hum$lr[ge10])))
write_csv(agree, file.path(outroot, "human_agreement.csv"))
print(agree, width = Inf)

c_lab <- sprintf("Spearman = %.2f\n%.0f%% within 2-fold\n(>= 10 reads)", agree$spearman_any,
                 agree$pct_within_2fold_ge10)
p_c <- ggplot(filter(hum, pmax(bambu_mean, oarfish_mean) >= 1),
              aes(log10(bambu_mean + 1), log10(oarfish_mean + 1))) +
  geom_bin_2d(bins = 70) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey30", linewidth = 0.3) +
  annotate("text", x = -Inf, y = Inf, label = c_lab, hjust = -0.05, vjust = 1.2, lineheight = 1.2) +
  scale_fill_viridis_c(trans = "log10", name = "Transcripts", option = "D") +
  ## headroom above the data for the three note lines
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.5))) +
  # (no coord_equal: the 45 mm box is fixed and square, see the save section)
  labs(x = "Bambu,\nlog10(mean reads + 1)", y = "oarfish,\nlog10(mean reads + 1)",
       title = "Human transcripts: Bambu vs oarfish counts",
       subtitle = "Curated set, mean over libraries; dashed: identity")

## d: detection agreement; libraries (light) and their median
jac_med <- jac %>% group_by(threshold) %>% summarise(jaccard = median(jaccard), .groups = "drop")
p_d <- ggplot(jac, aes(factor(threshold), jaccard)) +
  geom_line(aes(group = sample), colour = run_col, alpha = 0.25, linewidth = 0.3) +
  geom_line(data = jac_med, aes(group = 1), colour = run_col, linewidth = 0.5) +
  geom_point(data = jac_med, colour = run_col, size = 1.2) +
  labs(x = "Detected at (>= reads)", y = "Jaccard (Bambu, oarfish)",
       title = "Detection agreement",
       subtitle = "Human transcripts of the curated set; libraries (light) and their median")

## h: disagreement by overdispersion
dis <- hum %>% filter(ge10, is.finite(overdispersion)) %>%
  mutate(od_bin = od_cut(overdispersion)) %>%
  group_by(od_bin) %>%
  summarise(n = n(), median_abs_lr = median(abs(lr)), pct_gt_2fold = 100 * mean(abs(lr) > 1), .groups = "drop")
write_csv(dis, file.path(outroot, "disagreement_by_overdispersion.csv"))
print(dis)
p_h <- ggplot(dis, aes(od_bin, pct_gt_2fold, group = 1)) +
  geom_line(colour = run_col, linewidth = 0.5) +
  ## plain dots, n printed above each, as Fig. 1b, g, j
  geom_point(colour = run_col, size = 2) +
  geom_text(aes(label = k_lab(n), hjust = if_else(as.integer(od_bin) == 1, 0.5, 1.15),
                vjust = if_else(as.integer(od_bin) == 1, -1.2, -0.5)), size = 7 / .pt, family = font_family,   # upper left: the line rises to the right
            colour = colorRampPalette(c(run_col, "black"))(10)[3]) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = "RTA overdispersion", y = "% with > 2-fold difference\n(Bambu vs oarfish)",
       title = "Quantifiers disagree where reads are ambiguous",
       subtitle = "Human transcripts with >= 10 reads (either quantifier);\nnumbers = transcripts in the bin") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

## ================================================================ save ======
## Every plot box (each facet of b too) is the same 50 x 50 mm square. The
## Extended Data page
## (11_assemble_figures.R) packs the panels row by row at their own width:
##   ED4 i-l: sequin counts, divided vs undivided, Jaccard, disagreement
##   Fig. 2:  sequin length ratio, sequin DE, human counts, overdispersion
## Within a row every panel gets the row's largest top / bottom room, so the
## boxes of a row share their top and bottom edges; left / right room is each
## panel's own (no padding to the widest legend of the page).
box <- 50; b_gap <- 3
box_fig2 <- 50  # Fig. 2 copies
fix_box <- function(p, b = box) p + ggh4x::force_panelsizes(rows = unit(b, "mm"), cols = unit(b, "mm"))
## room around the box(es), mm: left / right / top / bottom
room <- function(p) {
  g  <- ggplotGrob(p)
  pl <- g$layout[grepl("^panel", g$layout$name), ]
  w  <- grid::convertWidth(g$widths, "mm", valueOnly = TRUE)
  h  <- grid::convertHeight(g$heights, "mm", valueOnly = TRUE)
  c(L = sum(w[seq_len(min(pl$l) - 1)]), R = sum(w[-seq_len(max(pl$r))]),
    T = sum(h[seq_len(min(pl$t) - 1)]), B = sum(h[-seq_len(max(pl$b))]))
}
plots <- list(
  ed_q_a_sequin_length           = p_a_len,
  ed_q_a0_sequin_counts          = p_a,
  ed_q_b_sequin_de               = p_b + theme(panel.spacing.x = unit(b_gap, "mm")),
  ed_q_c_human_counts            = p_c,
  ed_q_d_jaccard                 = p_d,
  ed_q_e_divided_vs_undivided    = p_e,
  ed_q_f_overdispersion          = p_f,
  ed_q_g_overdispersion_isoforms = p_g,
  ed_q_h_disagreement            = p_h)
## the sequin length ratio, sequin DE Bambu | oarfish, human counts and
## overdispersion by category and by isoforms per gene go to Fig. 2 (one aligned
## set; helpers/fig2_rows.R: rows 1-2 aligned with 08's panels, row 2 = five
## equal squares, 8 pt text). The four other Dong panels are ED4 i-l, one
## row under the LongBench a-h (11 assembles the page): shared top / bottom
## room, box size chosen so the row fills the 340 mm page like the LongBench rows
ed_rows <- list(c("ed_q_a0_sequin_counts", "ed_q_e_divided_vs_undivided", "ed_q_d_jaccard",
                  "ed_q_h_disagreement"))
source(file.path(helper_dir, "fig2_rows.R"))           # fig2_rows, fig2_text, fig2_layout, fig2_save
f2_q <- intersect(unlist(fig2_rows), names(plots))     # c, d, e, f, g of Fig. 2
stopifnot(setequal(c(unlist(ed_rows), f2_q), names(plots)))
inner_w <- function(nm, b = box) if (nm == "ed_q_b_sequin_de") 2 * b + b_gap else b
fig2_measure_open()  # text metrics of the SVG device
nt   <- lapply(plots, function(p) fix_box(p + labs(title = NULL, subtitle = NULL)))
rm_  <- sapply(nt, room)
fig2_measure_close()
slack <- 0.5
ed_W  <- fig2_W - 2 * fig2_mar                         # the ED4 page's usable width (340 - 2 x 5 mm)
box_ed <- (ed_W - fig2_gut * (length(ed_rows[[1]]) - 1) -
             sum(rm_["L", ed_rows[[1]]] + rm_["R", ed_rows[[1]]] + 2 * slack)) / length(ed_rows[[1]])
message(sprintf("ED4 i-l (Dong): boxes %.1f mm", box_ed))
sizes <- list()
for (rw in ed_rows) {
  bx <- box_ed; dir_nt <- figdir_nt
  tb <- apply(rm_[c("T", "B"), rw, drop = FALSE], 1, max) + slack   # the row's common top / bottom room  (room() includes the 1 mm theme plot.margin, which the new margin replaces)
  for (nm in rw) {
    r  <- rm_[, nm]
    sz <- c(w = unname(r["L"] + inner_w(nm, bx) + r["R"] + 2 * slack), h = unname(tb["T"] + bx + tb["B"]))
    p  <- fix_box(plots[[nm]] + labs(title = NULL, subtitle = NULL), bx) + theme(plot.margin = margin(1 + tb["T"] - r["T"], 1 + slack, 1 + tb["B"] - r["B"],
                                                1 + slack, "mm"))
    save_panel(p, nm, sz, dir_nt)
    sizes[[nm]] <- sz
    ## titled version: the same box, its own natural size
    pt <- fix_box(plots[[nm]])
    grDevices::pdf(NULL); rt <- room(pt); grDevices::dev.off()
    save_panel(pt, nm, c(rt["L"] + inner_w(nm) + rt["R"], rt["T"] + box + rt["B"]), figdir)
  }
}
## Fig. 2 rows 1-2: no-title files in figures/fig2_quantifier_panels/
## at the shared layout (helpers/fig2_rows.R); titled versions (50 mm box, natural
## size) in figures/ed_quantifier_panels/, also 8 pt
f2_plots <- lapply(setNames(f2_q, f2_q), function(nm) plots[[nm]] + labs(title = NULL, subtitle = NULL) + fig2_text)
fig2_update_rooms(f2_plots)
lay <- fig2_layout()
for (nm in f2_q) {
  fig2_save(f2_plots[[nm]], nm, lay, file.path(fig_dir, "fig2_quantifier_panels"))
  sizes[[nm]] <- lay$size[[nm]]
  pt <- fix_box(plots[[nm]] + fig2_text)
  grDevices::pdf(NULL); rt <- room(pt); grDevices::dev.off()
  save_panel(pt, nm, c(rt["L"] + inner_w(nm) + rt["R"], rt["T"] + box + rt["B"]), figdir)
}
message(sprintf("Fig. 2 rows 1-2: row 2 boxes %.1f mm; a's bars %.1f mm wide (rerun 06 if it ran before this)",
                lay$box2, lay$a_w))

## figure data behind every panel, kept for sharing
source(file.path(repo, "downstream/helpers/source_data.R"))
save_source_data(plots, file.path(outroot, "source_data_quantifiers.rds"))
saveRDS(list(rows = ed_rows, fig2 = f2_q, sizes = sizes), file.path(outroot, "ed_panel_sizes_mm.rds"))   # read by 11
unlink(file.path(outroot, "ed_cell_size_mm.rds"))  # remove a file of an earlier one-cell-size layout
message("Done: ", length(plots), " panels in ", figdir, "; tables in ", outroot)
