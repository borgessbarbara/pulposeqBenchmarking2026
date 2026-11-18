## ===========================================================================
## 03_quantifier_plots.R - Bambu vs oarfish, and the read-to-transcript
##                         ambiguity (RTA) overdispersion, on the SAME pulposeq
##                         runs (same alignments, annotation and libraries)
##                         -> supplementary figure panels + tables
##
## A  where the quantifier matters (spike-ins)
##    a1  recall / F1 / R2 per spike-in set (and sequin DE), Bambu vs oarfish
##    a3  the other spike-in metrics (E0 within 2-fold, long SIRV vs E0, ERCC
##        reads per library, ERCC LOD70), Bambu vs oarfish
##    (a2, oarfish / Bambu reads vs length, is Fig. 1f: drawn by 08_fig1_panels.R)
## B  the RTA division
##    b1  sequins: observed log2FC with divided vs undivided oarfish counts
##    b2  RTA overdispersion by transcript category (human reference, novel,
##        spike-ins), expressed transcripts only
##    b3  RTA overdispersion of human transcripts by isoforms per gene in the index
## C  Bambu vs oarfish on human transcripts (curated set: same rows in both)
##    c1  detection agreement (Jaccard) by detection threshold (1-100 reads)
##    c2  mean counts, Bambu vs oarfish (2-D histogram)
##    c3  disagreement |log2(oarfish / Bambu)| by RTA overdispersion
##
## Inputs  (read only)
##   pulposeq runs: transcriptome/curated/counts/{BambuOutput,oarfish}_counts_transcript_validated.txt,
##                  transcriptome/curated/metadata/{annotated_transcriptome,novel_transcripts_validated}_metadata.csv,
##                  transcriptome/bambu/BambuOutput_counts_transcript.txt (isoforms per gene, whole index),
##                  transcriptome/oarfish/*.quant (lengths)
##   RTA overdispersion: read_oarfish(curated = FALSE) - the .rds caches 01 wrote in
##                  spike-ins/results/oarfish_cache/ (else catchOarfish() reads the bootstraps)
##   01 / 02 outputs for both quantifiers: spike-ins/results/{sirv_ercc,sequins}{,_bambu}/
## Outputs
##   spike-ins/figures/quantifier/supp_q_{a1,a3,b1,b2,b3,c1,c2,c3}_*.{pdf,svg,png}
##   spike-ins/figures/quantifier_no_titles/   (the same, no title / subtitle)
##   spike-ins/results/quantifier_comparison/
##     q_metrics_spikeins.csv            a3: every summary metric, both quantifiers
##     q_metrics_by_set.csv              a1: recall / precision / F1 / R2 (/ accuracy) per set
##     overdispersion_by_category.csv    b2: n, median, quartiles, % > 1.5 / > 2
##     overdispersion_by_isoforms.csv    b3
##     human_agreement.csv               c2/c3: Spearman, % within 2-fold, median |log2 ratio|
##     jaccard_per_sample.csv            c1
##     per_transcript_quantifiers.csv.gz one row per curated transcript x protocol
## Needs: tidyverse, edgeR >= 4.99.12 only if the caches are missing, svglite, showtext
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
source(file.path(helper_dir, "plot_theme.R"))              # font_family, base_theme, protocol_cols, npg_cols, plat_labels
source(file.path(helper_dir, "figure_helpers.R"))          # save_panel
source(file.path(shared_helper_dir, "oarfish_counts.R"))   # read_oarfish, strip_bam
source(file.path(shared_helper_dir, "spikein_counts.R"))   # quant_dir

prot_levels <- c("ONT cDNA", "ONT DRS", "PacBio")
outroot <- file.path(spk_dir, "quantifier_comparison")
figdir  <- file.path(lb_dir, "spike-ins/figures/quantifier")
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir,  recursive = TRUE, showWarnings = FALSE)

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

q_cols   <- c(Bambu = "grey55", oarfish = npg_cols[["red"]])
spk_cols <- c("SIRV E0" = "grey30", "long SIRV" = npg_cols[["red"]],
              "ERCC" = npg_cols[["greyblue"]], "sequin" = npg_cols[["mint"]])
cat_levels <- c("protein-coding", "lncRNA", "other annotated",
                "novel protein-coding", "novel lncRNA", "novel non-coding",
                "sequin", "SIRV E0", "long SIRV", "ERCC")
spike_re <- "^(R[12]_|SIRV|ERCC)"

## ===================================================== per-transcript data ==
cat_of <- function(id, gene_bt, novel_bt) case_when(
  grepl("^R[12]_", id)                    ~ "sequin",
  grepl("^SIRV[0-9]{3}$", id)             ~ "SIRV E0",
  grepl("^SIRV[0-9]{4}", id)              ~ "long SIRV",
  grepl("^ERCC", id)                      ~ "ERCC",
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

per_run <- map(platforms, function(p) {
  message("Reading ", p, " ...")
  rd <- run_dir[[p]]
  o  <- read_mat(run_file(p, "oarfish_tx"))
  b  <- read_mat(run_file(p, "counts_tx"))
  ids <- rownames(o$m); smp <- colnames(o$m)
  stopifnot(all(ids %in% rownames(b$m)), all(smp %in% colnames(b$m)))
  om <- o$m; bm <- b$m[ids, smp, drop = FALSE]

  ## RTA overdispersion of every model oarfish quantified (cached by 01), then
  ## the curated rows
  of <- read_oarfish(rd, "transcript", curated = FALSE, cache_dir = oarfish_cache)
  od <- unname(of$overdispersion[ids])

  ## lengths (oarfish .quant), isoforms per gene in the whole index (Bambu table)
  q1  <- list.files(file.path(rd, "transcriptome", "oarfish"), "\\.quant$", full.names = TRUE)[1]
  qq  <- read_tsv(q1, show_col_types = FALSE, col_select = c(tname, len))
  len <- setNames(qq$len, qq$tname)
  full  <- read_tsv(run_file(p, "bambu_counts_tx"), col_select = 1:2, show_col_types = FALSE)
  n_iso <- table(full[[2]])
  ann <- read_csv(run_file(p, "ann_meta"), show_col_types = FALSE,
                  col_select = c(transcript_id_version, gene_biotype))
  nov <- read_csv(run_file(p, "novel_meta"), show_col_types = FALSE,
                  col_select = c(qry_id, pulposeq_transcript_biotype))

  tab <- tibble(protocol       = plat_labels[[p]],
                id             = ids,
                gene           = unname(o$gene[ids]),
                len            = unname(len[ids]),
                n_iso          = as.integer(n_iso[unname(o$gene[ids])]),
                gene_bt        = ann$gene_biotype[match(ids, ann$transcript_id_version)],
                novel_bt       = nov$pulposeq_transcript_biotype[match(ids, nov$qry_id)],
                overdispersion = od,
                bambu_mean     = rowMeans(bm),
                oarfish_mean   = rowMeans(om),
                libs_ge10      = rowSums(om >= 10)) %>%
    mutate(category = factor(cat_of(id, gene_bt, novel_bt), levels = cat_levels))

  ## per-sample detection agreement, human transcripts only
  hum <- !grepl(spike_re, ids)
  jac <- map_dfr(c(1, 2, 5, 10, 20, 50, 100), function(th) tibble(   # detection thresholds (reads)
    protocol = plat_labels[[p]], threshold = th, sample = smp,
    n_bambu  = map_dbl(smp, ~ sum(bm[hum, .x] >= th)),
    n_oarfish = map_dbl(smp, ~ sum(om[hum, .x] >= th)),
    jaccard  = map_dbl(smp, function(s) {
      a <- bm[hum, s] >= th; z <- om[hum, s] >= th
      sum(a & z) / sum(a | z)
    })))
  list(tab = tab, jac = jac)
})

tx <- map_dfr(per_run, "tab") %>%
  mutate(protocol = factor(protocol, levels = prot_levels),
         lr       = log2((oarfish_mean + 1) / (bambu_mean + 1)))   # oarfish vs Bambu
jac <- map_dfr(per_run, "jac") %>% mutate(protocol = factor(protocol, levels = prot_levels))
write_csv(tx,  file.path(outroot, "per_transcript_quantifiers.csv.gz"))
write_csv(jac, file.path(outroot, "jaccard_per_sample.csv"))
message("Categories without a biotype (check): ",
        sum(tx$category == "other annotated" & is.na(tx$gene_bt) & !grepl("^Bambu", tx$id)))

## ============================================ A: where the quantifier matters ==
rd_both <- function(dir, file) bind_rows(
  read_csv(file.path(quant_dir(dir, "oarfish"), file), show_col_types = FALSE),
  read_csv(file.path(quant_dir(dir, "bambu"),   file), show_col_types = FALSE))

se_q <- rd_both(sirv_dir, "sirv_ercc_summary.csv")
fc_q <- rd_both(seq_dir, "edger/fold_change_summary_edger.csv")

metric_lab <- c(sirv_f1       = "SIRV F1",
                e0_within2    = "E0 within 2-fold (%)",
                long_vs_e0    = "Long SIRV vs E0 (log2)",
                ercc_reads    = "ERCC reads per library\n(log10)",
                ercc_r2       = "ERCC R\u00b2",
                ercc_lod      = "ERCC LOD70, amol/\u00b5l\n(log10)",
                seq_f1        = "Sequin isoform DE F1")
qm <- bind_rows(
  se_q %>% transmute(quantifier, protocol,
                     sirv_f1    = sirv_f1,
                     e0_within2 = 100 * e0_within_2fold,
                     long_vs_e0 = long_vs_e0_log2,
                     ercc_reads = log10(ercc_reads_per_library),
                     ercc_r2    = ercc_r2,
                     ercc_lod   = log10(ercc_lod70_amol_ul)) %>%
    pivot_longer(-c(quantifier, protocol), names_to = "metric"),
  ## undivided counts for both, so only the quantifier differs
  fc_q %>% filter(level == "transcript", count_mode == "undivided") %>%
    transmute(quantifier, protocol = plat_labels[platform], metric = "seq_f1", value = f1)) %>%
  mutate(quantifier = recode(quantifier, bambu = "Bambu"),
         quantifier = factor(quantifier, levels = c("Bambu", "oarfish")),
         protocol   = factor(protocol, levels = prot_levels),
         metric     = factor(metric_lab[metric], levels = metric_lab))
write_csv(qm, file.path(outroot, "q_metrics_spikeins.csv"))

## a1: recall, F1 and R2 per spike-in set, both quantifiers in one plot:
## x = protocol, colour = set, Bambu dashed / oarfish solid (Bambu an open
## point, oarfish a filled one on top, so ties show as a filled dot in a ring).
##   SIRV E0, long SIRV, ERCC, sequin: identification, as Fig. 1b - recall =
##     detected / in the mix; F1 with precision = detected / (detected + novel
##     models on the set after validation). No R2 (ERCC: the dose-response R2).
##   sequin DE: isoform-level DE (edgeR, undivided counts for both quantifiers,
##     so only the quantifier differs) - recall = sensitivity, F1, R2 of observed
##     vs expected log2FC.
## Sensitivity = recall; accuracy is not drawn (identification has no true
## negatives; for sequin DE it is in q_metrics_by_set.csv).
set_lv   <- c("SIRV E0", "long SIRV", "ERCC", "sequin", "sequin DE")
set_cols <- c(spk_cols, "sequin DE" = npg_cols[["green"]])
qn <- c(oarfish = "oarfish", bambu = "Bambu")
n_in_mix <- map_dfr(names(qn), function(q) bind_rows(
  map_dfr(platforms, ~ read_csv(file.path(quant_dir(sirv_dir, q), .x, "sirv_e0_per_isoform.csv"),
                                show_col_types = FALSE)) %>%
    dplyr::count(protocol, name = "expected") %>% mutate(set = "SIRV E0"),
  read_csv(file.path(quant_dir(sirv_dir, q), "length_bias.csv"), show_col_types = FALSE) %>%
    filter(set == "long SIRV") %>% dplyr::count(protocol, name = "expected") %>% mutate(set = "long SIRV")) %>%
  mutate(quantifier = qn[[q]]))
fp_val <- map_dfr(names(qn), function(q)
  read_csv(file.path(quant_dir(sirv_dir, q), "novel_models_on_spikeins.csv"), col_types = cols(.default = "c")) %>%
    filter(stage == "validated") %>%
    mutate(set = case_when(family == "ERCC" ~ "ERCC", grepl("^SIRV[1-7]$", seqname) ~ "SIRV E0",
                           TRUE ~ "long SIRV"), quantifier = qn[[q]]) %>%
    dplyr::count(quantifier, protocol, set, name = "fp"))
sid_q <- rd_both(seq_dir, "sequin_identification.csv")
id_q <- bind_rows(
  se_q %>% transmute(quantifier = qn[quantifier], protocol, set = "SIRV E0", detected = e0_detected),
  se_q %>% transmute(quantifier = qn[quantifier], protocol, set = "long SIRV", detected = long_detected),
  se_q %>% transmute(quantifier = qn[quantifier], protocol, set = "ERCC", detected = ercc_detected,
                     expected = ercc_total),
  sid_q %>% transmute(quantifier = qn[quantifier], protocol, set = "sequin", detected,
                      expected = annotated, fp = novel_validated)) %>%
  rows_update(n_in_mix, by = c("quantifier", "protocol", "set"), unmatched = "ignore") %>%
  rows_update(fp_val, by = c("quantifier", "protocol", "set"), unmatched = "ignore") %>%
  mutate(fp = replace_na(fp, 0), recall = detected / expected,
         precision = detected / (detected + fp), f1 = 2 * precision * recall / (precision + recall),
         r2 = NA_real_)
id_q$r2[id_q$set == "ERCC"] <- se_q$ercc_r2[match(
  paste(id_q$quantifier, id_q$protocol)[id_q$set == "ERCC"], paste(qn[se_q$quantifier], se_q$protocol))]
de_q <- fc_q %>% filter(level == "transcript", count_mode == "undivided") %>%
  transmute(quantifier = qn[quantifier], protocol = plat_labels[platform], set = "sequin DE",
            recall = sensitivity, precision, f1, r2, accuracy = (true_pos + true_neg) / tested)
qs <- bind_rows(id_q, de_q) %>%
  mutate(quantifier = factor(quantifier, levels = c("Bambu", "oarfish")),
         protocol   = factor(protocol, levels = prot_levels),
         set        = factor(set, levels = set_lv)) %>%
  arrange(set, protocol, quantifier)
## the table also carries the oarfish DIVIDED sequin DE rows (the primary analysis,
## Fig. 1e) next to the undivided ones drawn here, with a count_mode column.
## Identification and ERCC rows use per-library rounded undivided counts.
de_div <- fc_q %>% filter(level == "transcript", count_mode == "divided", quantifier == "oarfish") %>%
  transmute(quantifier = "oarfish", protocol = plat_labels[platform], set = "sequin DE",
            recall = sensitivity, precision, f1, r2, accuracy = (true_pos + true_neg) / tested)
qs_out <- bind_rows(qs %>% mutate(count_mode = "undivided"),
                    de_div %>% mutate(count_mode = "divided (primary)",
                                      quantifier = factor(quantifier, levels = c("Bambu", "oarfish")),
                                      protocol   = factor(protocol, levels = prot_levels),
                                      set        = factor(set, levels = set_lv))) %>%
  relocate(count_mode, .after = set) %>%
  arrange(set, protocol, quantifier, desc(count_mode))
write_csv(qs_out, file.path(outroot, "q_metrics_by_set.csv"))
qs_long <- qs %>%
  pivot_longer(c(recall, f1, r2), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric, levels = c("recall", "f1", "r2"),
                         labels = c("Recall", "F1", "R\u00b2")))
p_a1 <- ggplot(arrange(qs_long, quantifier),                 # Bambu first: oarfish drawn on top
               aes(protocol, value, colour = set, group = interaction(set, quantifier))) +
  geom_line(aes(linetype = quantifier), linewidth = 0.5) +
  geom_point(aes(shape = quantifier, size = quantifier), fill = "white", stroke = 0.6) +
  scale_colour_manual(values = set_cols, breaks = set_lv, name = NULL) +
  ## one quantifier legend: line type + point (Bambu dashed, open; oarfish solid, filled)
  scale_linetype_manual(values = c(Bambu = "22", oarfish = "solid"), name = NULL) +
  scale_shape_manual(values = c(Bambu = 21, oarfish = 16), name = NULL) +
  scale_size_manual(values = c(Bambu = 2.2, oarfish = 1.3), name = NULL) +
  ## from 0.25, as Fig. 1b (every value is >= 0.32)
  scale_y_continuous(limits = c(0.25, 1), breaks = seq(0.25, 1, 0.25),
                     labels = c("0.25", "0.5", "0.75", "1")) +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid", shape = 16, size = 1.3)),
         linetype = guide_legend(order = 2, override.aes = list(colour = "grey30"),
                                 keywidth = unit(8, "mm")),
         shape = guide_legend(order = 2), size = guide_legend(order = 2)) +
  facet_wrap(~ metric, nrow = 1) +
  labs(x = NULL, y = NULL, title = "Spike-in metrics, Bambu vs oarfish (same runs)",
       subtitle = paste0("Identification of each set (precision after validation) and sequin ",
                         "isoform DE (undivided counts for both)")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),   # as Fig. 1b
        panel.spacing.x = unit(3, "mm"),
        legend.spacing.y = unit(3, "mm"))

## a3: the other spike-in metrics (x = protocol, one line per quantifier, one
## facet per metric with its own y-axis). Quantifiers drawn as in a1 (Bambu
## dashed with an open point, oarfish solid with a filled point on top), in one
## neutral colour: green and blue would repeat the ONT cDNA and ONT DRS colours
## of c1 and c3.
p_a3 <- ggplot(filter(qm, !metric %in% metric_lab[c("sirv_f1", "ercc_r2", "seq_f1")]) %>%
                 arrange(quantifier),                        # Bambu first: oarfish drawn on top
               aes(protocol, value, group = quantifier)) +
  geom_line(aes(linetype = quantifier), linewidth = 0.5, colour = "grey25") +
  geom_point(aes(shape = quantifier, size = quantifier), colour = "grey25", fill = "white", stroke = 0.6) +
  scale_linetype_manual(values = c(Bambu = "22", oarfish = "solid"), name = NULL) +
  scale_shape_manual(values = c(Bambu = 21, oarfish = 16), name = NULL) +
  scale_size_manual(values = c(Bambu = 2.2, oarfish = 1.3), name = NULL) +
  guides(linetype = guide_legend(keywidth = unit(8, "mm"))) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "Spike-in quantification, Bambu vs oarfish (same runs)",
       subtitle = paste0("Undivided counts for both quantifiers; lower is better for the ERCC LOD70 ",
                         "and closer to 0 for long SIRV vs E0")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),   # as Fig. 1b
        panel.spacing.x = unit(3, "mm"),
        axis.text.x = element_text(angle = 35, hjust = 1))

## ===================================================== B: the RTA division ==
de_res <- read_csv(file.path(seq_dir, "edger/sequin_de_results.csv"), show_col_types = FALSE) %>%
  filter(level == "transcript")
b1 <- inner_join(
  de_res %>% filter(count_mode == "divided") %>%
    transmute(platform, id, divided = obs_lfc, od = overdispersion, out_div = outcome),
  de_res %>% filter(count_mode == "undivided") %>%
    transmute(platform, id, undivided = obs_lfc, out_undiv = outcome),
  by = c("platform", "id")) %>%
  filter(is.finite(divided), is.finite(undivided)) %>%
  mutate(protocol = factor(plat_labels[platform], levels = prot_levels),
         call     = if_else(out_div == out_undiv, "same call", "call changed"))
b1_lab <- fc_q %>%
  filter(quantifier == "oarfish", level == "transcript") %>%
  dplyr::select(platform, count_mode, f1) %>%
  pivot_wider(names_from = count_mode, values_from = f1) %>%
  mutate(protocol = factor(plat_labels[platform], levels = prot_levels),
         lab = sprintf("F1 divided = %.3f\nF1 undivided = %.3f", divided, undivided))

p_b1 <- ggplot(b1, aes(undivided, divided)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = od, shape = call), size = 1.1) +
  geom_text(data = b1_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2,
            lineheight = 1.2, inherit.aes = FALSE) +
  scale_colour_viridis_c(name = "RTA\noverdispersion", option = "C", end = 0.9) +
  scale_shape_manual(values = c("same call" = 16, "call changed" = 4), name = NULL) +
  facet_wrap(~ protocol, nrow = 1) +
  labs(x = "Observed log2FC, undivided counts", y = "Observed log2FC, divided counts",
       title = "Sequins: the RTA division changes almost nothing",
       subtitle = "Isoform level; each point is one tested sequin; dashed: identity")

## b2: overdispersion by category, expressed transcripts (>= 10 reads in >= 2 libraries)
od_x <- tx %>% filter(libs_ge10 >= 2, is.finite(overdispersion))
od_sum <- od_x %>%
  group_by(protocol, category) %>%
  summarise(n = n(), median = median(overdispersion),
            q25 = quantile(overdispersion, 0.25), q75 = quantile(overdispersion, 0.75),
            q90 = quantile(overdispersion, 0.90),
            pct_gt_1.5 = 100 * mean(overdispersion > 1.5),
            pct_gt_2   = 100 * mean(overdispersion > 2), .groups = "drop")
write_csv(od_sum, file.path(outroot, "overdispersion_by_category.csv"))
print(od_sum, n = Inf, width = Inf)

od_brk <- c(1, 1.25, 1.5, 2, 3, 5, 10)
## line plot: x = the overdispersion bins of c3, y = % of each category's
## transcripts in the bin; one line per category, one facet per protocol;
## spike-ins dashed (spk_cols)
od_bin_lv <- c("1-1.05", "1.05-1.25", "1.25-1.5", "1.5-2", "2-3", "> 3")
cat_cols_b2 <- c("protein-coding"       = unname(cat_cols[["Annotated protein-coding"]]),
                 "lncRNA"               = unname(cat_cols[["Annotated lncRNA"]]),
                 "other annotated"      = unname(cat_cols[["Annotated other"]]),
                 "novel protein-coding" = unname(cat_cols[["Novel protein-coding"]]),
                 "novel lncRNA"         = unname(cat_cols[["Novel lncRNA"]]),
                 "novel non-coding"     = unname(cat_cols[["Novel non-coding"]]),
                 spk_cols)
od_share <- od_x %>%
  mutate(od_bin = cut(overdispersion, c(0, 1.05, 1.25, 1.5, 2, 3, Inf), labels = od_bin_lv)) %>%
  dplyr::count(protocol, category, od_bin, .drop = FALSE) %>%
  group_by(protocol, category) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  filter(is.finite(pct)) %>%
  mutate(spike = category %in% names(spk_cols))
write_csv(od_share, file.path(outroot, "overdispersion_share_by_category.csv"))
p_b2 <- ggplot(od_share, aes(od_bin, pct, colour = category, group = category)) +
  geom_line(aes(linetype = spike), linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = cat_cols_b2, breaks = cat_levels, name = NULL) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
  guides(colour = guide_legend(override.aes = list(
    linetype = ifelse(cat_levels %in% names(spk_cols), "22", "solid")))) +
  facet_wrap(~ protocol, nrow = 1) +
  labs(x = "RTA overdispersion", y = "% of the category's transcripts",
       title = "RTA overdispersion by transcript category",
       subtitle = "Expressed transcripts (>= 10 reads in >= 2 libraries); dashed = spike-ins") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.spacing.x = unit(3, "mm"))

## b3: human transcripts, by isoforms per gene in the whole index
iso_lv <- c("1", "2-3", "4-9", "10-19", ">= 20")
od_iso <- od_x %>%
  filter(!grepl(spike_re, id), !is.na(n_iso)) %>%
  mutate(iso_bin = cut(n_iso, c(0, 1, 3, 9, 19, Inf), labels = iso_lv)) %>%
  group_by(protocol, iso_bin) %>%
  summarise(n = n(), median = median(overdispersion),
            q25 = quantile(overdispersion, 0.25), q75 = quantile(overdispersion, 0.75),
            .groups = "drop")
write_csv(od_iso, file.path(outroot, "overdispersion_by_isoforms.csv"))
## line plot: x = the same overdispersion bins, y = % of the transcripts of each
## isoforms-per-gene bin; one line per bin, light (1 isoform) to dark (>= 20),
## viridis steps without the pale yellow end; one facet per protocol
od_iso_x <- od_x %>%
  filter(!grepl(spike_re, id), !is.na(n_iso)) %>%
  mutate(iso_bin = cut(n_iso, c(0, 1, 3, 9, 19, Inf), labels = iso_lv),
         od_bin  = cut(overdispersion, c(0, 1.05, 1.25, 1.5, 2, 3, Inf), labels = od_bin_lv))
od_iso_share <- od_iso_x %>%
  dplyr::count(protocol, iso_bin, od_bin, .drop = FALSE) %>%
  group_by(protocol, iso_bin) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  filter(is.finite(pct))
write_csv(od_iso_share, file.path(outroot, "overdispersion_share_by_isoforms.csv"))
iso_cols <- setNames(c("#B8DE29", "#55C667", "#1F968B", "#39568C", "#440154"), iso_lv)
p_b3 <- ggplot(od_iso_share, aes(od_bin, pct, colour = iso_bin, group = iso_bin)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = iso_cols, name = "Isoforms per gene") +
  facet_wrap(~ protocol, nrow = 1) +
  labs(x = "RTA overdispersion", y = "% of the group's transcripts",
       title = "Human transcripts: ambiguity grows with isoforms per gene",
       subtitle = "Expressed reference and novel transcripts; isoforms per gene in the oarfish index") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # 35-degree labels would collide at this width
        panel.spacing.x = unit(3, "mm"))

## ====================================== C: human transcripts, Bambu vs oarfish ==
hum <- tx %>% filter(!grepl(spike_re, id))
agree <- hum %>%
  group_by(protocol) %>%
  summarise(n_any        = sum(pmax(bambu_mean, oarfish_mean) >= 1),
            spearman_any = cor(bambu_mean[pmax(bambu_mean, oarfish_mean) >= 1],
                               oarfish_mean[pmax(bambu_mean, oarfish_mean) >= 1], method = "spearman"),
            n_ge10       = sum(pmax(bambu_mean, oarfish_mean) >= 10),
            pct_within_2fold_ge10 = 100 * mean(abs(lr[pmax(bambu_mean, oarfish_mean) >= 10]) <= 1),
            median_abs_lr_ge10    = median(abs(lr[pmax(bambu_mean, oarfish_mean) >= 10])),
            .groups = "drop")
write_csv(agree, file.path(outroot, "human_agreement.csv"))
print(agree, width = Inf)

## c1: x = detection threshold, y = median
## Jaccard of the libraries (Bambu vs oarfish detected sets, human transcripts)
jac_med <- jac %>% group_by(protocol, threshold) %>%
  summarise(jaccard = median(jaccard), .groups = "drop")
p_c1 <- ggplot(jac_med, aes(factor(threshold), jaccard, colour = protocol, group = protocol)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = protocol_cols, breaks = prot_levels, name = NULL) +
  labs(x = "Detected at (>= reads)", y = "Jaccard (Bambu, oarfish)",
       title = "Detection agreement",
       subtitle = "Human transcripts of the curated set; median of the libraries")

c2_lab <- agree %>%
  mutate(lab = sprintf("Spearman = %.2f\n%.0f%% within 2-fold (>= 10 reads)",
                       spearman_any, pct_within_2fold_ge10))
p_c2 <- ggplot(filter(hum, pmax(bambu_mean, oarfish_mean) >= 1),
               aes(log10(bambu_mean + 1), log10(oarfish_mean + 1))) +
  geom_bin_2d(bins = 70) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey30", linewidth = 0.3) +
  geom_text(data = c2_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2,
            lineheight = 1.2, inherit.aes = FALSE) +
  scale_fill_viridis_c(trans = "log10", name = "Transcripts", option = "D") +
  facet_wrap(~ protocol, nrow = 1) +
  coord_equal() +
  labs(x = "Bambu, log10(mean reads + 1)", y = "oarfish, log10(mean reads + 1)",
       title = "Human transcripts: Bambu vs oarfish counts",
       subtitle = "Curated set, mean over libraries; dashed: identity")

od_lv <- c("1-1.05", "1.05-1.25", "1.25-1.5", "1.5-2", "2-3", "> 3")
dis <- hum %>%
  filter(pmax(bambu_mean, oarfish_mean) >= 10, is.finite(overdispersion)) %>%
  mutate(od_bin = cut(overdispersion, c(0, 1.05, 1.25, 1.5, 2, 3, Inf), labels = od_lv)) %>%
  group_by(protocol, od_bin) %>%
  summarise(n = n(), median_abs_lr = median(abs(lr)),
            pct_gt_2fold = 100 * mean(abs(lr) > 1), .groups = "drop")
write_csv(dis, file.path(outroot, "disagreement_by_overdispersion.csv"))
## plain dots with n printed at each; in each bin the highest protocol's n goes
## above its dot, the lowest's below, the middle one's to the right; label
## colours = darker protocol colours
k_lab <- function(x) ifelse(x >= 9950, paste0(round(x / 1000), "k"),   # 9,950+ -> "10k", not "10.0k"
                     ifelse(x >= 1000, sprintf("%.1fk", x / 1000), as.character(x)))   # 1.7k, 227
c3_lab_cols <- vapply(protocol_cols, function(cl) colorRampPalette(c(cl, "black"))(10)[3], "")
c3_lab <- dis %>% group_by(od_bin) %>%
  mutate(rk = rank(-pct_gt_2fold, ties.method = "first"), nb = n()) %>% ungroup() %>%
  mutate(vj = case_when(rk == 1 ~ -0.9, rk == nb ~ 1.9, TRUE ~ 0.5),
         hj = if_else(rk == 1 | rk == nb, 0.5, -0.35))
p_c3 <- ggplot(dis, aes(od_bin, pct_gt_2fold, colour = protocol, group = protocol)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  ## labels repelled from every dot, started above the highest and below the
  ## lowest dot of each bin; short leader lines when a label had to move
  ggrepel::geom_text_repel(data = c3_lab, aes(label = k_lab(n)), size = 7 / .pt,
            family = font_family, colour = c3_lab_cols[as.character(c3_lab$protocol)],
            nudge_y = ifelse(c3_lab$rk == 1, 2.8, ifelse(c3_lab$rk == c3_lab$nb, -2.8, 0)),
            point.padding = 0.7, box.padding = 0.3, min.segment.length = 0.4, segment.size = 0.2,
            segment.colour = "grey50", max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = protocol_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "RTA overdispersion", y = "% with > 2-fold difference\n(Bambu vs oarfish)",
       title = "Quantifiers disagree where reads are ambiguous",
       subtitle = "Human transcripts with >= 10 reads (either quantifier);\nnumbers = transcripts in the bin") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

## ================================================================ save ======
## sizes (mm) = extent on the supplementary figure page (340 x 347 mm, Fig. 1's
## width; 5 mm margins, 4 mm gutters, no empty areas),
## letter in brackets:
##   row 1 (85)  a1 150 [a] | a3 176 [b]
##   row 2 (80)  c2 216 [c] | c1 110 [d]   (c2: square facets)
##   row 3 (80)  b1 150 [e] | b2 176 [f]
##   row 4 (80)  b3 200 [g] | c3 126 [h]
## (a2 = Fig. 1f, drawn by 08_fig1_panels.R)
panels <- list(
  supp_q_a1_spikein_metrics       = list(p = p_a1, size = c(150,  85)),
  supp_q_a3_spikein_quant_metrics = list(p = p_a3, size = c(176,  85)),
  supp_q_c2_human_counts          = list(p = p_c2, size = c(216,  80)),
  supp_q_c1_jaccard               = list(p = p_c1, size = c(110,  80)),
  supp_q_b1_divided_vs_undivided  = list(p = p_b1, size = c(150,  80)),
  supp_q_b2_overdispersion        = list(p = p_b2, size = c(176,  80)),
  supp_q_b3_overdispersion_isoforms = list(p = p_b3, size = c(200,  80)),
  supp_q_c3_disagreement          = list(p = p_c3, size = c(126,  80))
)
for (nm in names(panels)) save_panel(panels[[nm]]$p, nm, panels[[nm]]$size, figdir)
## figure data behind every panel, kept for sharing
source(file.path(repo, "downstream/helpers/source_data.R"))
save_source_data(panels, file.path(outroot, "source_data_ed1.rds"))
## the same panels without title and subtitle (same file names and sizes), in
## figures/quantifier_no_titles/
figdir_nt <- file.path(lb_dir, "spike-ins/figures/quantifier_no_titles")
for (nm in names(panels))
  save_panel(panels[[nm]]$p + labs(title = NULL, subtitle = NULL), nm, panels[[nm]]$size, figdir_nt)
message("Done: ", length(panels), " panels in ", figdir, "; tables in ", outroot)
