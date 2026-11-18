## ===========================================================================
## 02_sirv_ercc.R - SIRV-Set 4 benchmark per protocol (transcript level)  (Fig. 1b, c)
##   oarfish runs; the same script also runs on Bambu's counts of the SAME runs
##   (Rscript 02_sirv_ercc.R bambu) for the quantifier comparison (dong/03)
##
## THREE QUESTIONS (spike-ins are an answer key: amounts are known)
##   Q1 identification : are all spike-ins in the mix recovered, and does the
##                       pipeline avoid inventing models?   (SIRV P / R / F1)
##   Q2 abundance      : are equimolar molecules equally abundant, and how does
##                       that change with length?          (SIRV log2 obs/exp)
##   Q3 dose-response  : does signal rise with input, and what is the lowest
##                       input still detected?              (ERCC)
##
## Counts: UNDIVIDED estimated reads per transcript (nothing is tested here),
##   quantifier "oarfish" (default; downstream/helpers/oarfish_counts.R) or
##   "bambu" (Bambu `counts`), both from pulposeq's curated tables. They are
##   EM-estimated reads (full-length + partial), NOT molecules. With oarfish,
##   each isoform's read-to-transcript ambiguity (RTA) overdispersion from the
##   bootstraps (edgeR::catchOarfish; >= 1, 1 = no ambiguity) is reported
##   next to its abundance error: assignment errors should come with high
##   overdispersion, capture (library) errors with low overdispersion.
##   Needs edgeR >= 4.99.12 (catchOarfish via read_oarfish()).
##
## Expected amounts: Lexogen SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx
##   SIRV E0     69 isoforms, all 60 amol/ul (equimolar)
##   long SIRVs  15 transcripts, 4-12 kb, also 60 amol/ul
##   ERCC        92 transcripts, one mix (SIRV ERCC Mix 1), ~2e6-fold
##               concentration range; every sample received the same mix
##
## SIRV annotation given to pulposeq: Lexogen "correct" (c) annotation =
## exactly the 69 E0 isoforms. The 31 isoforms only in the over-annotated "o"
## set are not in the mix AND not in the GTF used here, so they cannot receive
## reads: the decoy test is reported as not applicable (a stated limitation).
## The "insufficient" (i) annotation (discovery of the 25 missing E0 isoforms)
## was not run either.
##
## Definitions (identical for both quantifiers)
##   * ONE denominator for every SIRV statistic: CPM within the 84 SIRVs
##     (E0 + long), expected share = 1/84 (as in LongBench 03_Bulk_quantification)
##   * two-tier detection, set once below:
##       detected   = >= det_min_reads reads in >= det_min_libs libraries
##                    (counts rounded per library: oarfish / Bambu are fractional)
##       quantified = detected AND |median log2(obs/exp)| <= 3 (within 8-fold,
##                    both directions)
##     F1 therefore measures PRESENCE, not abundance accuracy
##   * ONE value per SIRV isoform everywhere (figure dots, length classes,
##     summary, Fig. 1b label): the median across libraries of the per-library
##     log2(obs/exp)
##   * detection-rule sensitivity grid (detection_grid.csv), incl. a
##     LongBench-like rule (filterByExpr min.count = 5, min.total.count = 15)
##   * length classes (E0 < 1 kb, E0 1-2.5 kb, long SIRV)
##   * ERCC: fit on ERCCs passing the detection rule; hi-concentration fit,
##     Pearson r over ALL ERCCs, residual-vs-length rho, reads per library
##
## Per protocol (ONT cDNA, ONT DRS, PacBio):
##   1. SIRV accuracy (84 SIRVs): log2(observed / expected) per isoform and
##      library; within 2-fold, within 8-fold, lost (median < 1/8 expected),
##      unstable (SD of log2 ratio > 1)
##   2. reads on over-annotation-only SIRV isoforms - only if annotated (else NA)
##   3. novel models (BambuTx*) on spike-in sequences, raw Bambu vs validated
##   4. SIRV identification P / R / F1, raw and validated
##        positives = 69 E0 + 15 long SIRVs;  TP = detected;  FN = not detected
##        FP = novel SIRV models at that stage (+ decoy reads when annotated)
##   5. length bias: E0 + long SIRVs are equimolar -> log2 ratio vs length
##   6. ERCC: dose response (log2 CPM vs log2 input among detected ERCCs),
##      70% detection limit (logistic), detection above fixed concentrations
##
## Inputs : <run>/transcriptome/oarfish/ (catchOarfish) or the curated Bambu
##          table, <run>/transcriptome/{bambu,curated/annotation}/ GTFs,
##          data/SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx
## Outputs: spike-ins/results/sirv_ercc/ (oarfish; read by 08) or
##          spike-ins/results/sirv_ercc_bambu/ (Bambu; read by dong/03):
##   sirv_ercc_summary.csv, length_bias.csv, length_class_summary.csv,
##   detection_grid.csv, ercc_per_transcript.csv, novel_models_on_spikeins.csv,
##   <platform>/sirv_e0_per_isoform.csv (E0 only), <platform>/sirv_per_isoform.csv (all 84,
##   with the RTA overdispersion)
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")   # platforms, run_dir
source(file.path(shared_helper_dir, "oarfish_counts.R"))                    # read_oarfish()
source(file.path(shared_helper_dir, "spikein_counts.R"))                    # read_counts_q()
plat_labels <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio")

## "oarfish" (default, primary) or "bambu": Rscript 02_sirv_ercc.R bambu,
## or set it here when sourcing in RStudio
quantifier <- quantifier_arg()
outroot <- quant_dir(sirv_dir, quantifier)
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
message("Quantifier: ", quantifier, " -> ", outroot)

lex_file  <- sirv_design            # references/SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx

## Novel models are counted at two stages (independent of the quantifier):
##   raw        <run>/transcriptome/bambu/bambu_novel_transcripts.gtf  (Bambu discovery)
##   validated  <run>/transcriptome/curated/annotation/BambuOutput_annotations_validated.gtf
##              (what pulposeq delivers; novel = transcript_id BambuTx*)
gtf_stage <- list(
  raw       = setNames(vapply(platforms, run_file, "", key = "bambu_novel_gtf"), platforms),
  validated = setNames(vapply(platforms, run_file, "", key = "validated_gtf"), platforms))

sirv_annotation <- "correct (c): 69 E0 isoforms"   # recorded in the summary table

## ---- detection rules (decide once, report the grid) ----------------------
det_min_reads   <- 3          # SIRV detected: >= this many reads ...
det_min_libs    <- 2          # ... in >= this many of the libraries
quant_abs_log2  <- 3          # "quantified": |median log2(obs/exp)| <= 3 (within 8-fold, both directions)
ercc_min_reads  <- 1          # ERCC detected: >= this many reads ...
ercc_min_libs   <- 2          # ... in >= this many libraries (single-library hits can be cross-talk)
sirv_min_reads  <- 1          # only for the (not applicable) decoy count
ercc_conc_cuts  <- c(1, 11.71875)                  # amol/ul; detection reported above these
sirv_rule_txt <- sprintf(">=%d reads in >=%d libraries", det_min_reads, det_min_libs)
ercc_rule_txt <- sprintf(">=%d reads in >=%d libraries", ercc_min_reads, ercc_min_libs)

## ----------------------------------------------------- expected amounts ---
## readxl drops empty leading columns, so columns are located relative to
## the ID column (found by its content: SIRV101 / SIRV4001 / ERCC-00002),
## not by fixed position
id_re <- "^(SIRV\\d{3,}|ERCC-\\d+)$"
read_sheet <- function(sheet, conc_off, len_off, flag_off = NULL) {
  d  <- readxl::read_excel(lex_file, sheet = sheet, col_names = FALSE,
                           .name_repair = "minimal")
  ic <- which.max(vapply(d, function(x) sum(grepl(id_re, as.character(x))), numeric(1)))
  out <- tibble(id   = as.character(d[[ic]]),
                conc = suppressWarnings(as.numeric(d[[ic + conc_off]])),
                len  = suppressWarnings(as.numeric(d[[ic + len_off]])))
  if (!is.null(flag_off))                               # annotation flags c / i / o
    out <- mutate(out,
                  in_correct      = !is.na(suppressWarnings(as.numeric(d[[ic + flag_off[1]]]))),
                  in_insufficient = !is.na(suppressWarnings(as.numeric(d[[ic + flag_off[2]]]))),
                  in_over         = !is.na(suppressWarnings(as.numeric(d[[ic + flag_off[3]]]))))
  filter(out, grepl(id_re, id))
}
sirv_all <- read_sheet("SIRV isoform", conc_off = 18, len_off = 21, flag_off = c(2, 3, 4))
e0       <- filter(sirv_all, !is.na(conc))                         # 69 in the mix
over_only <- filter(sirv_all, is.na(conc), in_over, !in_correct)$id # 31, "o" annotation only
long     <- read_sheet("longSIRV", conc_off = 3, len_off = 6)
ercc     <- read_sheet("ERCC",     conc_off = 4, len_off = 7) %>% filter(!is.na(conc))
message(sprintf("Lexogen: %d E0 isoforms, %d over-annotation-only, %d long SIRVs, %d ERCC",
                nrow(e0), length(over_only), nrow(long), nrow(ercc)))
## sanity: counts, E0 + long SIRVs at 60 amol/ul, long SIRVs 4-12 kb, the
## correct annotation is exactly the E0 mix
stopifnot("Lexogen sheet not read as expected" =
            nrow(e0) == 69 && length(over_only) == 31 && nrow(long) == 15 && nrow(ercc) == 92 &&
            all(e0$conc == 60) && all(long$conc == 60) &&
            setequal(e0$id, filter(sirv_all, in_correct)$id) &&
            all(long$len > 3000 & long$len < 13000) &&
            all(e0$len > 100 & e0$len < 5000) && !anyNA(ercc$len))
positives <- c(e0$id, long$id)                                      # 84 SIRV molecules in the mix

## ------------------------------------------------------------- helpers ----
## undivided estimated reads of the SIRV / ERCC rows (curated set); the RTA
## overdispersion (oarfish; NA for Bambu) is kept as an attribute
read_tx <- function(plat) {
  x <- read_counts_q(run_dir[[plat]], "transcript", quantifier, rows = "^(SIRV|ERCC)",
                     cache_dir = oarfish_cache)
  m <- x$reads
  attr(m, "overdispersion") <- x$od
  m
}

## novel Bambu models (BambuTx*) on SIRV / ERCC sequences in one GTF.
## Transcripts are rebuilt from exon records, so it works whether or not the
## GTF has "transcript" lines.
read_novel_spikein <- function(f, plat, stage) {
  if (!file.exists(f))
    stop(sprintf("%s GTF not found for %s: %s\n  edit `gtf_stage` at the top of the script",
                 stage, plat, f))
  g <- read.table(f, sep = "\t", comment.char = "#", quote = "",
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attr"),
                  colClasses = c("character", "character", "character", "integer", "integer",
                                 "character", "character", "character", "character"))
  g %>%
    filter(feature == "exon", grepl("^(SIRV|ERCC)", seqname)) %>%
    mutate(transcript_id = str_match(attr, 'transcript_id "([^"]+)"')[, 2],
           gene_id       = str_match(attr, 'gene_id "([^"]+)"')[, 2]) %>%
    filter(grepl("^BambuTx", transcript_id)) %>%
    group_by(transcript_id, gene_id, seqname, strand) %>%
    ## exonic_len first: once `start` / `end` are summarised they would mask
    ## the per-exon columns
    summarise(exonic_len = sum(end - start + 1L), exons = n(),
              start = min(start), end = max(end), .groups = "drop") %>%
    transmute(protocol = plat_labels[[plat]], stage,
              family = if_else(grepl("^SIRV", seqname), "SIRV", "ERCC"),
              seqname, start, end, strand, exons, exonic_len, gene_id, transcript_id)
}

prf <- function(tp, fn, fp) {
  p  <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  r  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  f1 <- if (is.na(p) || is.na(r) || p + r == 0) NA_real_ else 2 * p * r / (p + r)
  c(precision = p, recall = r, f1 = f1)
}

## rows for `ids`, 0 for anything absent from the matrix (= not detected)
take <- function(m, ids) {
  out <- matrix(0, length(ids), ncol(m), dimnames = list(ids, colnames(m)))
  hit <- intersect(ids, rownames(m))
  out[hit, ] <- m[hit, ]
  out
}
cpm_mat <- function(m) sweep(m, 2, pmax(colSums(m), 1), "/") * 1e6   # not `cpm`: would mask edgeR::cpm

## detected = >= min_reads in >= min_libs libraries (named logical, one per row);
## counts rounded per library first (estimated reads are fractional)
detect <- function(x, min_reads, min_libs) rowSums(round(x) >= min_reads) >= min_libs

## sensitivity of "detected" to the rule: reads x libraries grid + LongBench-like rule
## (edgeR filterByExpr, min.count = 5, each library its own group, min.total.count = 15;
##  raw counts instead of the CPM cut-off, so an approximation)
detect_grid <- function(x, spikein) {
  g <- expand_grid(min_reads = c(1, 2, 3, 5), min_libs = c(1, 2, 4)) %>%
    mutate(rule = sprintf(">=%d reads in >=%d libraries", min_reads, min_libs),
           detected = map2_int(min_reads, min_libs, ~ sum(detect(x, .x, .y))))
  bind_rows(g, tibble(rule = "LongBench-like: >=5 reads in >=1 library, >=15 in total",
                      detected = sum(rowSums(round(x) >= 5) >= 1 & rowSums(round(x)) >= 15))) %>%
    mutate(spikein = spikein, total = nrow(x))
}

## ----------------------------------------------------------------- run ----
len_all <- list(); ercc_all <- list(); summ <- list(); novel_all <- list(); grid_all <- list()

for (plat in platforms) {
  prot   <- plat_labels[[plat]]
  outdir <- file.path(outroot, plat)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  m  <- read_tx(plat)
  od <- attr(m, "overdispersion")                       # NA: Bambu, or not in the curated set

  ## ---- 1. SIRV accuracy: all 84 SIRVs share ONE denominator (expected 1/84) ----
  sirv_m   <- take(m, positives)
  sirv_cpm <- cpm_mat(sirv_m)                          # columns sum to 1e6 over the 84 SIRVs
  sirv_l <- as_tibble(log2((sirv_cpm + 1) / (1e6 / nrow(sirv_m))), rownames = "id") %>%
    pivot_longer(-id, names_to = "sample", values_to = "log2_ratio") %>%
    mutate(protocol = prot,
           set  = if_else(id %in% e0$id, "SIRV E0", "long SIRV"),
           gene = if_else(set == "SIRV E0", sub("^(SIRV\\d)\\d+$", "\\1", id), NA_character_))
  sirv_iso <- sirv_l %>%
    group_by(protocol, set, gene, id) %>%
    summarise(median = median(log2_ratio), sd = sd(log2_ratio), .groups = "drop") %>%
    mutate(overdispersion = unname(od[id]))
  e0_iso <- sirv_iso %>% filter(set == "SIRV E0") %>% select(protocol, gene, id, median, sd, overdispersion)
  write_csv(e0_iso,   file.path(outdir, "sirv_e0_per_isoform.csv"))   # E0 only (Fig. 1b label)
  write_csv(sirv_iso, file.path(outdir, "sirv_per_isoform.csv"))      # all 84

  ## detection (tier 1) and quantification (tier 2)
  det_sirv <- detect(sirv_m, det_min_reads, det_min_libs)
  med_sirv <- setNames(sirv_iso$median, sirv_iso$id)
  quant_ok <- det_sirv & abs(med_sirv[names(det_sirv)]) <= quant_abs_log2

  ## ---- 2. reads on over-annotation-only SIRV isoforms (if annotated) ----
  ## An over-annotation-only ID missing from the curated table was not in the
  ## GTF (or had too few reads to be curated); with the correct annotation used
  ## here this is moot - SIRV104 etc. are not in the Lexogen GTF at all.
  over_in_annot  <- intersect(over_only, rownames(m))
  decoy_testable <- length(over_in_annot) > 0
  if (!decoy_testable)
    message(sprintf("[%s] none of the %d over-annotation-only SIRV isoforms is in the annotation: decoy test not applicable",
                    prot, length(over_only)))
  over_reads    <- if (decoy_testable) sum(m[over_in_annot, ]) else NA_real_
  over_with_rds <- if (decoy_testable) sum(rowSums(round(m[over_in_annot, , drop = FALSE])) >= sirv_min_reads) else NA_integer_
  sirv_tot      <- sum(take(m, c(positives, over_in_annot)))

  ## ---- 3. novel models on spike-in sequences: raw Bambu vs validated ----
  nov <- bind_rows(read_novel_spikein(gtf_stage$raw[[plat]],       plat, "raw"),
                   read_novel_spikein(gtf_stage$validated[[plat]], plat, "validated"))
  novel_all[[plat]] <- nov
  n_nov <- function(fam, st) sum(nov$family == fam & nov$stage == st)
  for (st in c("raw", "validated")) {
    x <- filter(nov, stage == st)
    if (nrow(x))
      message(sprintf("[%s] %s: novel models on spike-ins: %s", prot, st,
                      paste(sprintf("%s %s:%d-%d %s, %d exon(s), %d nt", x$transcript_id, x$seqname,
                                    x$start, x$end, x$strand, x$exons, x$exonic_len), collapse = "; ")))
  }

  ## ---- 4. SIRV identification: precision / recall / F1, per stage ----
  ## TP / FN from the validated counts (annotated SIRVs are not curated away);
  ## FP = novel SIRV models at that stage (+ decoy isoforms with reads, if annotated)
  pos_det  <- det_sirv[positives]
  tp <- sum(pos_det)
  fn <- sum(!pos_det)
  decoy_fp <- if (decoy_testable) over_with_rds else 0L
  fp_raw   <- n_nov("SIRV", "raw")       + decoy_fp
  fp_val   <- n_nov("SIRV", "validated") + decoy_fp
  s_raw    <- prf(tp, fn, fp_raw)
  s_val    <- prf(tp, fn, fp_val)
  if (fn > 0)
    message(sprintf("[%s] SIRV isoforms in the mix not detected (%s): %s", prot, sirv_rule_txt,
                    paste(positives[!pos_det], collapse = ", ")))

  ## detection-rule sensitivity (SIRV: with P/R/F1 after validation; ERCC: counts only)
  g_sirv <- detect_grid(sirv_m, "SIRV (E0 + long)")
  s_grid <- t(vapply(g_sirv$detected, function(d) prf(d, nrow(sirv_m) - d, fp_val), numeric(3)))
  g_sirv <- bind_cols(g_sirv, as_tibble(s_grid))
  g_ercc <- detect_grid(take(m, ercc$id), "ERCC")
  grid_all[[plat]] <- bind_rows(g_sirv, g_ercc) %>% mutate(protocol = prot, .before = 1)

  ## ---- 5. length bias: E0 + long SIRVs, all equimolar ----
  ## same per-isoform value as sirv_*_per_isoform.csv and the summary: median
  ## across libraries of the per-library log2(obs/exp)
  eq      <- bind_rows(mutate(e0, set = "SIRV E0"), mutate(long, set = "long SIRV"))
  eq_rel  <- med_sirv
  len_t   <- tibble(protocol = prot, id = eq$id, set = eq$set, len = eq$len,
                    log2_ratio = unname(eq_rel[eq$id]), overdispersion = unname(od[eq$id]))
  len_all[[plat]] <- len_t
  len_fit <- lm(log2_ratio ~ log2(len), data = len_t)   # descriptive only: PacBio is a step, not a line
  e0_ge1kb <- len_t$set == "SIRV E0" & len_t$len >= 1000

  ## ---- 6. ERCC ----
  er_m    <- take(m, ercc$id)
  er_cpm  <- cpm_mat(er_m)
  det_frac <- rowMeans(round(er_m) >= 1)                  # fraction of libraries with >= 1 read (logistic fit)
  ercc_det     <- detect(er_m, ercc_min_reads, ercc_min_libs)   # rule used for "detected / in mix"
  ercc_det_any <- rowSums(round(er_m) >= 1) >= 1                 # most lenient, for comparison
  er_t <- tibble(protocol = prot, id = ercc$id, conc = ercc$conc, len = ercc$len,
                 mean_count = rowMeans(er_m), detected = det_frac,
                 detected_rule = unname(ercc_det),
                 log2_cpm = log2(rowMeans(er_cpm)))
  ercc_all[[plat]] <- er_t

  ok   <- is.finite(er_t$log2_cpm) & er_t$detected_rule   # fit only ERCCs passing the detection rule
  efit <- if (sum(ok) >= 3) lm(log2_cpm ~ log2(conc), data = er_t[ok, ]) else NULL
  hi   <- ok & er_t$conc >= 11.71875                      # well-detected range: slope without the detection floor
  hfit <- if (sum(hi) >= 3) lm(log2_cpm ~ log2(conc), data = er_t[hi, ]) else NULL
  ## Pearson r over ALL ERCCs and libraries, undetected ones at the pseudocount floor (LongBench style)
  r_all <- suppressWarnings(cor(as.vector(log2(er_cpm + 1)),
                                rep(log2(ercc$conc), times = ncol(er_cpm))))
  ## do residuals follow ERCC length? (PacBio: short ERCCs depleted, like the SIRVs)
  res_rho <- if (is.null(efit)) NA_real_ else
    suppressWarnings(cor(resid(efit), log2(er_t$len[ok]), method = "spearman"))

  k    <- round(det_frac * ncol(m))
  lfit <- suppressWarnings(glm(cbind(k, ncol(m) - k) ~ log2(conc), family = binomial, data = er_t))
  l_b0 <- unname(coef(lfit)[1]); l_b1 <- unname(coef(lfit)[2])
  lod  <- if (is.finite(l_b1) && l_b1 > 0) 2^((qlogis(0.7) - l_b0) / l_b1) else NA_real_  # conc at 70% detection

  ## ERCC detection above fixed concentrations (same rule as `ercc_detected`)
  above <- list()
  for (cc in ercc_conc_cuts) {
    s <- er_t$conc >= cc
    above[[sprintf("ercc_detected_ge_%s", format(cc, trim = TRUE))]] <- sum(ercc_det[s])
    above[[sprintf("ercc_total_ge_%s",    format(cc, trim = TRUE))]] <- sum(s)
  }
  ercc_above <- as_tibble(above)

  ## assignment vs capture: does the abundance error follow the RTA
  ## overdispersion? (detected SIRVs; NA for Bambu, which has no bootstraps)
  od_ok  <- det_sirv[sirv_iso$id] & is.finite(sirv_iso$overdispersion)
  od_rho <- if (sum(od_ok) >= 5)
    suppressWarnings(cor(abs(sirv_iso$median[od_ok]), log2(sirv_iso$overdispersion[od_ok]),
                         method = "spearman")) else NA_real_

  summ[[plat]] <- bind_cols(tibble(
    quantifier = quantifier,
    protocol = prot,
    sirv_annotation      = sirv_annotation,
    sirv_detection_rule  = sirv_rule_txt,
    ## SIRV E0 / long SIRV: detection and accuracy (one 1/84 denominator)
    e0_detected          = sum(det_sirv[e0$id]),
    e0_quantified        = sum(quant_ok[e0$id]),
    e0_median_abs_log2   = median(abs(e0_iso$median)),
    e0_within_2fold      = mean(abs(e0_iso$median) <= 1),
    e0_within_8fold      = mean(abs(e0_iso$median) <= 3),
    e0_lost              = sum(e0_iso$median < -3),
    e0_unstable          = sum(e0_iso$sd > 1),
    e0_median_overdispersion = median(e0_iso$overdispersion, na.rm = TRUE),   # RTA, bootstraps
    sirv_error_od_rho    = od_rho,          # Spearman |log2 obs/exp| vs log2 overdispersion
    long_detected        = sum(det_sirv[long$id]),
    long_quantified      = sum(quant_ok[long$id]),
    ## decoy test (NA when the over-annotation isoforms are not in the GTF)
    decoy_test           = if (decoy_testable) "applicable" else "not applicable (not in annotation)",
    over_only_in_annotation   = length(over_in_annot),
    over_only_with_reads      = over_with_rds,
    over_only_read_pct        = if (decoy_testable) 100 * over_reads / max(sirv_tot, 1) else NA_real_,
    ## assembly: novel models on spike-in sequences, Bambu vs after validation
    novel_sirv_raw       = n_nov("SIRV", "raw"),
    novel_sirv_validated = n_nov("SIRV", "validated"),
    novel_ercc_raw       = n_nov("ERCC", "raw"),
    novel_ercc_validated = n_nov("ERCC", "validated"),
    ## SIRV identification = PRESENCE only (positives = 69 E0 + 15 long SIRVs)
    sirv_tp              = tp,
    sirv_fn              = fn,
    sirv_fp_raw          = fp_raw,
    sirv_fp_validated    = fp_val,
    sirv_recall          = unname(s_val["recall"]),
    sirv_precision_raw   = unname(s_raw["precision"]),
    sirv_f1_raw          = unname(s_raw["f1"]),
    sirv_precision       = unname(s_val["precision"]),     # validated = what pulposeq delivers
    sirv_f1              = unname(s_val["f1"]),
    ## length bias (kept for compatibility; prefer length_class_summary.csv)
    length_slope         = unname(coef(len_fit)[2]),   # log2 ratio per doubling of length (linear: descriptive)
    long_vs_e0_log2      = mean(len_t$log2_ratio[len_t$set == "long SIRV"]) -
                           mean(len_t$log2_ratio[len_t$set == "SIRV E0"]),
    long_vs_e0_ge1kb_log2 = mean(len_t$log2_ratio[len_t$set == "long SIRV"]) -
                           mean(len_t$log2_ratio[e0_ge1kb]),   # fair contrast: E0 >= 1 kb only
    ## ERCC
    ercc_detection_rule  = ercc_rule_txt,
    ercc_detected        = sum(ercc_det),
    ercc_detected_any    = sum(ercc_det_any),            # >= 1 read in any library
    ercc_total           = nrow(ercc),
    ercc_reads_per_library = sum(rowMeans(er_m)),        # depth confounds detection limits
    ercc_slope           = if (is.null(efit)) NA_real_ else unname(coef(efit)[2]),   # detected ERCCs only
    ercc_r2              = if (is.null(efit)) NA_real_ else summary(efit)$r.squared,
    ercc_slope_hi        = if (is.null(hfit)) NA_real_ else unname(coef(hfit)[2]),   # >= 11.7 amol/ul
    ercc_r2_hi           = if (is.null(hfit)) NA_real_ else summary(hfit)$r.squared,
    ercc_pearson_all     = r_all,                        # all ERCCs, undetected at the floor
    ercc_resid_length_rho = res_rho,                     # Spearman: residual vs log2 length
    ercc_lod70_amol_ul   = lod,                          # "70% detection limit" (Anaquin-style logistic, not Anaquin's LOD)
    ercc_lod70_logistic_slope = l_b1), ercc_above)
}

## ------------------------------------------------------------- outputs ----
summ      <- bind_rows(summ)
len_all   <- bind_rows(len_all)
ercc_all  <- bind_rows(ercc_all)
novel_all <- bind_rows(novel_all)
grid_all  <- bind_rows(grid_all)

## length classes: <1 kb E0, 1-2.5 kb E0, long SIRVs (4-12 kb). Values are the plotted
## per-isoform log2 ratios (median across libraries, 1/84 denominator)
len_class <- len_all %>%
  mutate(class = case_when(set == "long SIRV" ~ "long SIRV (4-12 kb)",
                           len < 1000         ~ "E0 < 1 kb",
                           TRUE               ~ "E0 1-2.5 kb")) %>%
  group_by(protocol, class) %>%
  summarise(n = n(), median_log2 = median(log2_ratio),
            within_2x = mean(abs(log2_ratio) <= 1), .groups = "drop")

write_csv(summ,      file.path(outroot, "sirv_ercc_summary.csv"))
write_csv(len_all,   file.path(outroot, "length_bias.csv"))
write_csv(len_class, file.path(outroot, "length_class_summary.csv"))
write_csv(grid_all,  file.path(outroot, "detection_grid.csv"))
write_csv(ercc_all,  file.path(outroot, "ercc_per_transcript.csv"))
write_csv(novel_all, file.path(outroot, "novel_models_on_spikeins.csv"))
print(summ, width = Inf)
print(len_class)

message("Done: ", outroot)
