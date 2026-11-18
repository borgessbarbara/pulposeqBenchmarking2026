## ===========================================================================
## 03_quantifier_comparison.R - oarfish vs Bambu on the SAME pulposeq runs
##                              (same alignments, annotation and libraries)
##                              -> Supplementary table + figure
##
## Reads what longbench/01, longbench/02 and dong/02 wrote for each quantifier:
##   oarfish  spike-ins/results/{sirv_ercc,sequins}/      (primary, Fig. 1 / 2)
##   bambu    spike-ins/results/{sirv_ercc,sequins}_bambu/ (Rscript <script> bambu)
##   Dong     dong/sequins/results{,_bambu}/
## Rows are quantifier x count mode, so the effect of the quantifier
## (bambu-undivided vs oarfish-undivided) is separated from the effect of the
## RTA division (oarfish-undivided vs oarfish-divided).
##
## Outputs: spike-ins/results/quantifier_comparison/
##   comparison_summary.csv     one row per dataset x protocol x quantifier x count mode
##   sirv_per_isoform_paired.csv  per SIRV: log2(obs/exp) with each quantifier + overdispersion
##   open_questions.csv         isoforms with known problems:
##                              R2_18 (all runs), SIRV612 / SIRV609 (ONT),
##                              PacBio gene-level wrong signs (R1_13, R1_91, R2_41)
##   figures: spike-ins/figures/supp_quantifier_sirv.pdf (+ .png)
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")   # sirv_dir, seq_dir, platforms
source(file.path(shared_helper_dir, "spikein_counts.R"))                    # quant_dir()
plat_labels <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio")
prot_levels <- unname(plat_labels)

dong_seq <- file.path(base, "downstream/dong/sequins/results")   # dong/00_paths.R seq_dir
outroot  <- file.path(spk_dir, "quantifier_comparison")
figdir   <- file.path(lb_dir, "spike-ins/figures")
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir,  recursive = TRUE, showWarnings = FALSE)

quants <- c("oarfish", "bambu")
rd <- function(dir, file) {
  f <- file.path(dir, file)
  if (!file.exists(f)) { message("missing (run that script first): ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE)
}
## the same, with a quantifier column (for tables that do not carry one)
rdq <- function(dir, file, q) { x <- rd(dir, file); if (is.null(x)) NULL else mutate(x, quantifier = q) }

## ---------------------------------------------- LongBench: SIRV / ERCC ----
se <- map_dfr(quants, ~ rd(quant_dir(sirv_dir, .x), "sirv_ercc_summary.csv"))
lc <- map_dfr(quants, ~ rdq(quant_dir(sirv_dir, .x), "length_class_summary.csv", .x))
sirv_tab <- se %>%
  transmute(dataset = "LongBench", quantifier, protocol,
            sirv_detected = sirv_tp, sirv_f1, e0_within_2fold, e0_within_8fold,
            long_quantified, e0_median_overdispersion, sirv_error_od_rho,
            ercc_detected, ercc_r2, ercc_r2_hi, ercc_lod70_amol_ul) %>%
  left_join(lc %>%
              transmute(quantifier, protocol,
                        class = recode(class, "E0 < 1 kb" = "med_e0_lt1kb",
                                       "E0 1-2.5 kb" = "med_e0_1to2.5kb",
                                       "long SIRV (4-12 kb)" = "med_long"),
                        median_log2) %>%
              pivot_wider(names_from = class, values_from = median_log2),
            by = c("quantifier", "protocol"))

## ------------------------------------------------- sequins: both datasets ----
lb_fc <- map_dfr(quants, ~ rd(file.path(quant_dir(seq_dir, .x), "edger"), "fold_change_summary_edger.csv")) %>%
  mutate(dataset = "LongBench", protocol = unname(plat_labels[platform]))
dg_fc <- map_dfr(quants, ~ rd(quant_dir(dong_seq, .x), "summary.csv")) %>%
  mutate(dataset = "Dong", protocol = "ONT cDNA", n_called = called)
fc <- bind_rows(lb_fc, dg_fc)
de_tab <- fc %>%
  transmute(dataset, quantifier, protocol, level, count_mode, tested, truly_de,
            n_called, true_pos, wrong_sign, false_pos, precision, sensitivity,
            sensitivity_in_mix, f1, spearman, r2, r2_no_zero, median_overdispersion) %>%
  pivot_wider(names_from = level,
              values_from = c(tested, truly_de, n_called, true_pos, wrong_sign, false_pos,
                              precision, sensitivity, sensitivity_in_mix, f1, spearman, r2,
                              r2_no_zero, median_overdispersion),
              names_glue = "{level}_{.value}")

lb_id <- map_dfr(quants, ~ rd(quant_dir(seq_dir, .x), "sequin_identification.csv")) %>%
  mutate(dataset = "LongBench")
dg_id <- map_dfr(quants, ~ rd(quant_dir(dong_seq, .x), "sequin_identification.csv")) %>%
  mutate(dataset = "Dong", protocol = "ONT cDNA")
id_tab <- bind_rows(lb_id, dg_id) %>%
  transmute(dataset, quantifier, protocol, sequin_detected = detected, sequin_recall = recall,
            sequin_f1 = f1_validated, sequin_lod70 = lod70_amol_ul,
            sequin_capture_length_rho = capture_length_rho)

comparison <- de_tab %>%
  left_join(id_tab,   by = c("dataset", "quantifier", "protocol")) %>%
  left_join(sirv_tab, by = c("dataset", "quantifier", "protocol")) %>%
  mutate(protocol = factor(protocol, levels = prot_levels),
         dataset  = factor(dataset, levels = c("LongBench", "Dong"))) %>%
  arrange(dataset, protocol, quantifier, count_mode)
write_csv(comparison, file.path(outroot, "comparison_summary.csv"))
print(comparison, n = Inf, width = Inf)

## ------------------------------------------- per-SIRV paired comparison ----
sirv_iso <- map_dfr(quants, function(q) map_dfr(platforms, function(p)
  rdq(file.path(quant_dir(sirv_dir, q), p), "sirv_per_isoform.csv", q)))
paired <- sirv_iso %>%
  dplyr::select(quantifier, protocol, set, id, median, overdispersion) %>%
  pivot_wider(names_from = quantifier, values_from = c(median, overdispersion)) %>%
  dplyr::select(-any_of("overdispersion_bambu")) %>%
  rename(log2_obs_exp_oarfish = median_oarfish, log2_obs_exp_bambu = median_bambu,
         overdispersion = overdispersion_oarfish) %>%
  mutate(protocol = factor(protocol, levels = prot_levels))
write_csv(paired, file.path(outroot, "sirv_per_isoform_paired.csv"))

p_sirv <- ggplot(paired, aes(log2_obs_exp_bambu, log2_obs_exp_oarfish)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_point(aes(colour = set), size = 1) +
  scale_colour_manual(values = c("SIRV E0" = "grey30", "long SIRV" = "#E64B35"), name = NULL) +
  facet_wrap(~ protocol, nrow = 1) +
  coord_equal() +
  labs(x = "Bambu: log2(observed / expected)", y = "oarfish: log2(observed / expected)",
       title = "SIRV abundance, Bambu vs oarfish (same alignments)",
       subtitle = "Each point is one SIRV (median across libraries); dashed: identity") +
  theme_bw(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(file.path(figdir, "supp_quantifier_sirv.pdf"), p_sirv, width = 180, height = 75, units = "mm")
ggsave(file.path(figdir, "supp_quantifier_sirv.png"), p_sirv, width = 180, height = 75, units = "mm",
       dpi = 300)

## ------------------------------------------------------- open questions ----
## R2_18 (both datasets, all protocols), PacBio gene-level wrong signs, SIRV6
lb_de <- map_dfr(quants, ~ rd(file.path(quant_dir(seq_dir, .x), "edger"), "sequin_de_results.csv")) %>%
  mutate(dataset = "LongBench", protocol = unname(plat_labels[platform]))
dg_de <- map_dfr(quants, ~ rd(quant_dir(dong_seq, .x), "sequin_de_results.csv")) %>%
  mutate(dataset = "Dong", protocol = "ONT cDNA")
watch_seq <- c("R2_18", "R2_18_1", "R2_18_2",                       # isoform swap, all runs
               "R1_13", "R1_13_1", "R1_13_2", "R1_91", "R1_91_1", "R1_91_2",
               "R2_41", "R2_41_1", "R2_41_2")                       # PacBio gene-level wrong signs
oq_seq <- bind_rows(lb_de, dg_de) %>%
  filter(id %in% watch_seq) %>%
  transmute(question = if_else(grepl("^R2_18", id), "R2_18 isoform swap",
                               "PacBio gene-level wrong sign (short-isoform capture)"),
            dataset, protocol, quantifier, count_mode, level, id,
            exp_lfc, obs_lfc, padj, outcome, overdispersion)
oq_sirv <- sirv_iso %>%
  filter(id %in% c("SIRV612", "SIRV609", "SIRV311", "SIRV608", "SIRV618")) %>%
  transmute(question = if_else(id %in% c("SIRV612", "SIRV609"),
                               "SIRV6 locus: SIRV612 0 reads in ONT, SIRV609 over-estimated",
                               "Short E0 isoforms missed by PacBio"),
            dataset = "LongBench", protocol, quantifier, count_mode = "undivided",
            level = "transcript", id, log2_obs_exp = median, overdispersion)
open_q <- bind_rows(oq_seq, oq_sirv) %>% arrange(question, id, dataset, protocol, quantifier, count_mode)
write_csv(open_q, file.path(outroot, "open_questions.csv"))
print(open_q, n = Inf, width = Inf)

message("Done: ", outroot)
