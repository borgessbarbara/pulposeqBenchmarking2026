## ===========================================================================
## 01_sequin_edger.R - LongBench sequin ground truth with edgeR      (Fig. 1e)
##   oarfish runs; the same script also runs on Bambu's counts of the SAME runs
##   (Rscript 01_sequin_edger.R bambu) for the quantifier comparison (dong/03)
##
## COUNTS (downstream/helpers/oarfish_counts.R + spikein_counts.R)
##   oarfish  transcript level edgeR::catchOarfish(), gene level
##            edgeR::catchOarfishGene() (edgeR >= 4.99.12): the bootstraps give
##            each feature's read-to-transcript ambiguity (RTA) overdispersion.
##            count_mode "divided"   = counts / overdispersion (PRIMARY; not
##                                     rounded; Baldoni et al. 2024 NAR, 2024 NARGAB)
##            count_mode "undivided" = estimated reads (sensitivity check: what
##                                     the division changes)
##   bambu    Bambu `counts`, count_mode "undivided" only (no bootstraps)
##   The DE filter (>= min_count in >= min_samples libraries) uses `filter_counts`:
##   "undivided" (default) = estimated reads, as dong/01_de_edger.R, so
##   the sequins validate exactly the test used for the biology; "divided" = the
##   counts edgeR tests, as filterByExpr in the published pipeline (Baldoni et al.
##   2024 NARGAB). Change it together with dong/01_de_edger.R. Identification always
##   uses the undivided estimated reads, rounded per library.
##
## WHAT SEQUINS TEST (the 4th question, after SIRV / ERCC in 02)
##   Q4 differential abundance: every sequin isoform is spiked at a KNOWN
##      concentration in Mix A and in Mix B, so its true log2FC (B / A) is
##      known. Mix B went into SCLC-A samples, Mix A into LUAD + SCLC-P samples.
##      A per-molecule capture bias (length, RT, size selection) is the same in
##      both mixes, so it cancels in a ratio: sequins test RELATIVE accuracy,
##      SIRV E0 tests ABSOLUTE accuracy.
##
## Test (as dong/02_sequin_edger.R and dong/01_de_edger.R):
##   edgeR v4 glmQLFit(robust = TRUE) -> glmTreat(lfc = 0.5), FDR 0.05;
##   display log2FC = predFC(prior.count = 2); filter >= min_count reads in
##   >= min_samples samples; library sizes AND TMM factors from the ratio-1
##   control sequins (|expected log2FC| < 0.25), so the normalisation never
##   sees the sequins that are meant to change. This also absorbs any change
##   in total sequin input between libraries (LongBench raised the sequin
##   input tenfold for DRS libraries run on secondary flow cells).
##   Contrast: Mix B (SCLC-A) vs Mix A (expected log2FC = log2(B / A)).
##
## Truth and scoring (identical to dong/02_sequin_edger.R)
##   truly DE  = |expected log2FC| > lfc
##   called    = padj < alpha
##   true positive = called, truly DE, correct sign; wrong-sign calls are
##   false discoveries; precision = TP / called = 1 - FDR; F1 from precision
##   and sensitivity. AUC ranks -log10(p) against truly DE.
##
## Identification (no DE tool involved): every sequin isoform in the mix
## (Anaquin, 164) is a positive; FP = novel BambuTx models on chrIS, counted
## in Bambu's GTF (raw) and after pulposeq validation.
##
## Also reported
##   * detected = >= 3 reads in >= 2 libraries, as the SIRVs (02); `detected_any`
##     (>= 1 read summed over all libraries) is kept for comparison
##   * sensitivity among TESTED sequins (conditional on the expression filter) and
##     among ALL truly-DE sequins in the mix (`sensitivity_in_mix`; detection x DE)
##   * robust fold-change agreement: Spearman rho, median |obs - exp|, and R2 /
##     slope (a) without sequins with 0 reads in every library of one mix (their
##     unshrunk logFC is ~ +/-10 and dominates R2) and (b) on the shrunken log2FC
##     (predFC, prior.count = 2)
##   * sequin_fc_outliers.csv: transcript-level sequins with |obs - exp| > 2 or a
##     wrong-sign call, with every sibling isoform (reads per mix) and the
##     gene-level result, to test whether outliers are read-assignment swaps
##     between isoforms; also genes with a gene-level error (a gene total is a
##     capture-weighted sum of its isoforms, so per-isoform bias does NOT cancel
##     at gene level when isoform proportions change between mixes)
##   * sequin dose-response as the ERCCs: reads per library, 70% detection limit
##     (logistic, >= 1 read per library vs that library's mix concentration),
##     detection by concentration band
##   * sequin_per_isoform.csv: expected A / B, reads per mix, detected, exonic
##     length and relative capture efficiency (reads per amol, centred on the
##     median detected isoform); Spearman rho of capture vs length per protocol
##     in sequin_identification.csv
##
## Inputs : <run>/transcriptome/oarfish/ (catchOarfish) or the Bambu tables,
##          <run>/transcriptome/{bambu,curated/annotation}/ GTFs,
##          data/samplesheet_*.csv, Anaquin
## Outputs: spike-ins/results/sequins/ (oarfish; read by 08) or
##          spike-ins/results/sequins_bambu/ (Bambu; read by dong/03):
##   edger/<platform>/<level>/<count_mode>/sequin_de.csv
##   edger/sequin_de_results.csv           every tested sequin, protocol / level / count mode
##   edger/fold_change_summary_edger.csv   one row per protocol x level x count mode
##   edger/sequin_fc_outliers.csv          outliers + sibling isoforms (primary count mode)
##   sequin_identification.csv             identification + dose-response on chrIS
##   sequin_per_isoform.csv                one row per sequin isoform x protocol
##   sequin_detection_by_conc.csv          per-library detection by concentration band
##   novel_models_on_sequins.csv           novel models on chrIS, raw vs validated
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(Anaquin)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
source(file.path(shared_helper_dir, "oarfish_counts.R"))   # read_oarfish()
source(file.path(shared_helper_dir, "spikein_counts.R"))   # read_counts_q()
plat_labels <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio")

## "oarfish" (default, primary) or "bambu": Rscript 01_sequin_edger.R bambu,
## or set it here when sourcing in RStudio
quantifier <- quantifier_arg()
seqroot <- quant_dir(seq_dir, quantifier)
message("Quantifier: ", quantifier, " -> ", seqroot)
outroot <- file.path(seqroot, "edger")
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------- settings ---
## platforms: 00_paths.R
sheets    <- c(ont_cdna = "samplesheet_ont.csv",
               ont_drs  = "samplesheet_ont_drs.csv",
               pb       = "samplesheet_pb.csv")
levels_to_run <- c("gene", "transcript")
alpha       <- 0.05   # LongBench DE settings
lfc         <- 0.5    # LongBench DE settings
min_count   <- 10     # LongBench DE settings
min_samples <- 2      # LongBench DE settings
ctrl_lfc    <- 0.25   # |expected log2FC| below this = ratio-1 control sequin
outlier_lfc <- 2      # |observed - expected log2FC| above this -> sequin_fc_outliers.csv
## sensitivity analysis (Supplementary Table 5, "without R2_18" rows): the same
## test with the R2_18 gene and both its isoforms removed before filtering and
## normalisation (R2_18_2 has no reads in any Mix A library, so R2_18_1 is
## called with the wrong sign in every protocol; see dong/03_quantifier_comparison.R)
sens_drop_re <- "^R2_18(_[0-9]+)?$"
## counts edgeR receives (see header); the first mode is the primary analysis
count_modes  <- if (quantifier == "oarfish")
  list(transcript = c("divided", "undivided"), gene = c("divided", "undivided")) else
  list(transcript = "undivided", gene = "undivided")
primary_mode <- count_modes$transcript[1]
filter_counts <- "undivided"   # DE filter on estimated reads (as dong/01_de_edger.R) or "divided"

## identification (no DE tool involved: reads per sequin + novel models on chrIS)
sequin_seq    <- "chrIS"  # sequin chromosome name in the genome / GTF
det_min_reads <- 3        # detected = >= this many reads ...
det_min_libs  <- 2        # ... in >= this many libraries (same rule as the SIRVs, 02)
det_rule_txt  <- sprintf(">=%d reads in >=%d libraries", det_min_reads, det_min_libs)
lod_p         <- 0.7      # detection limit = concentration detected in 70% of libraries
conc_breaks   <- c(0, 1, 10, 100, 1000, Inf)   # amol/ul bands for sequin_detection_by_conc.csv
gtf_stage <- list(
  raw       = setNames(vapply(platforms, run_file, "", key = "bambu_novel_gtf"), platforms),
  validated = setNames(vapply(platforms, run_file, "", key = "validated_gtf"), platforms))
tool        <- "edgeR"

## ---------------------------------------------- expected concentrations ---
## Anaquin mixture table: concentration per isoform in Mix A and Mix B (amol/ul)
data("RnaQuinIsoformMixture", package = "Anaquin")
iso_mix <- as_tibble(RnaQuinIsoformMixture) %>%
  transmute(id = as.character(Name), A = MixA, B = MixB)
gene_mix <- iso_mix %>%
  mutate(gene = sub("_[0-9]+$", "", id)) %>%
  group_by(gene) %>%
  summarise(A = sum(A), B = sum(B), .groups = "drop") %>%
  transmute(id = gene, A, B)
expected <- list(transcript = iso_mix, gene = gene_mix)
stopifnot("Anaquin mixture: non-positive concentration" =
            all(iso_mix$A > 0 & iso_mix$B > 0))
## sequin length from the Anaquin mixture table when it has one (all 164);
## otherwise from the validated GTFs (below)
anaquin_len <- if ("Length" %in% names(RnaQuinIsoformMixture))
  setNames(as.numeric(RnaQuinIsoformMixture$Length), as.character(RnaQuinIsoformMixture$Name)) else NULL

## -------------------------------------------------------------- inputs ----
## sequin rows of one run: list(reads, divided, od) from the chosen quantifier
## (curated = FALSE: every model before pulposeq's curation)
read_sequins <- function(plat, lv, curated = TRUE)
  read_counts_q(run_dir[[plat]], lv, quantifier, rows = "^R[12]_",
                curated = curated, cache_dir = oarfish_cache)

## rows for `ids`, 0 for anything absent from the matrix (= not in the curated
## set: too few reads)
take <- function(m, ids) {
  out <- matrix(0, length(ids), ncol(m), dimnames = list(ids, colnames(m)))
  hit <- intersect(ids, rownames(m))
  out[hit, ] <- m[hit, ]
  out
}

## detected = >= min_reads in >= min_libs libraries (named logical, one per row);
## estimated reads are rounded per library first
detect <- function(x, min_reads = det_min_reads, min_libs = det_min_libs)
  rowSums(round(x) >= min_reads) >= min_libs

sample_mix <- function(samples, sheet) {
  ss  <- read.csv(file.path(data_dir, sheet))
  grp <- gsub("-", "_", ss$group[match(samples, ss$sample)])
  stopifnot(!anyNA(grp))
  setNames(factor(ifelse(grp == "SCLC_A", "B", "A"), levels = c("A", "B")), samples)
}

auc_fn <- function(score, truth) {
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---------------------------------------------------------- scoring -------
## same definitions as dong/02_sequin_edger.R (do not edit one without the other)
score_calls <- function(t) {
  t %>%
    mutate(truth   = abs(exp_lfc) > lfc,
           called  = !is.na(padj) & padj < alpha,
           outcome = case_when(called & truth & sign(obs_lfc) == sign(exp_lfc) ~ "true positive",
                               called & truth                                  ~ "wrong sign",
                               called & !truth                                 ~ "false positive",
                               !called & truth                                 ~ "false negative",
                               TRUE                                            ~ "true negative"))
}

r2_slope <- function(y, x) {
  if (sum(is.finite(y) & is.finite(x)) < 3)
    return(c(slope = NA_real_, intercept = NA_real_, r2 = NA_real_))
  f <- lm(y ~ x)
  c(slope = unname(coef(f)[2]), intercept = unname(coef(f)[1]), r2 = summary(f)$r.squared)
}

## `truly_de_in_mix`: truly-DE sequins in the whole mix at this level, tested or not
summarise_calls <- function(t, truly_de_in_mix) {
  all_fit    <- r2_slope(t$obs_lfc, t$exp_lfc)
  nz         <- !t$one_mix_zero
  nz_fit     <- r2_slope(t$obs_lfc[nz], t$exp_lfc[nz])
  shr_fit    <- r2_slope(t$disp_lfc, t$exp_lfc)
  t %>%
    summarise(tested      = n(),
              truly_de    = sum(truth),
              n_called    = sum(called),
              true_pos    = sum(outcome == "true positive"),
              false_pos   = sum(outcome == "false positive"),
              wrong_sign  = sum(outcome == "wrong sign"),
              false_neg   = sum(outcome == "false negative"),
              true_neg    = sum(outcome == "true negative"),
              sensitivity = true_pos / max(truly_de, 1),            # among TESTED sequins
              fdr         = if (n_called > 0) (false_pos + wrong_sign) / n_called else NA_real_,
              precision   = if (n_called > 0) true_pos / n_called else NA_real_,
              f1          = if (is.na(precision) || precision + sensitivity == 0) NA_real_ else
                              2 * precision * sensitivity / (precision + sensitivity),
              auc         = auc_fn(-log10(pvalue), truth),
              slope       = unname(all_fit["slope"]),
              intercept   = unname(all_fit["intercept"]),
              r2          = unname(all_fit["r2"]),
              rmse        = sqrt(mean((obs_lfc - exp_lfc)^2)),
              ## robust agreement (outliers below dominate R2 and slope)
              spearman        = cor(obs_lfc, exp_lfc, method = "spearman"),
              median_abs_err  = median(abs(obs_lfc - exp_lfc)),
              n_one_mix_zero  = sum(one_mix_zero),
              slope_no_zero   = unname(nz_fit["slope"]),
              r2_no_zero      = unname(nz_fit["r2"]),
              slope_shrunk    = unname(shr_fit["slope"]),
              r2_shrunk       = unname(shr_fit["r2"]),
              ## detection x DE: TP among every truly-DE sequin in the mix
              truly_de_in_mix    = truly_de_in_mix,
              sensitivity_in_mix = true_pos / max(truly_de_in_mix, 1))
}

## ------------------------------------------------------------- edgeR ------
## cpm() is always called as edgeR::cpm(): a plain cpm() defined by another
## script in the same R session would otherwise mask it
## DGEList with library sizes + TMM factors from the control sequins
norm_dge <- function(m, mix, ctrl, lab) {
  y <- DGEList(m, group = mix)
  if (sum(ctrl) >= 3) {
    yc <- normLibSizes(DGEList(m[ctrl, , drop = FALSE]))
    y$samples$lib.size     <- yc$samples$lib.size
    y$samples$norm.factors <- yc$samples$norm.factors
  } else {
    message(sprintf("[%s] <3 ratio-1 sequins; TMM on all sequins", lab))
    y <- normLibSizes(y)
  }
  y
}

run_edger <- function(m, mix, ctrl, lab) {
  y      <- norm_dge(m, mix, ctrl, lab)
  design <- model.matrix(~ mix)
  cn     <- "mixB"                                         # Mix B vs Mix A
  fit    <- glmQLFit(y, design, robust = TRUE)
  tt     <- topTags(glmTreat(fit, coef = cn, lfc = lfc), n = Inf, sort.by = "none")$table
  shr    <- predFC(estimateDisp(y, design), design, prior.count = 2)[, cn]
  eff    <- y$samples$lib.size * y$samples$norm.factors
  ids    <- rownames(tt)
  a      <- mix == "A"
  tibble(id = ids, obs_lfc = tt$logFC, disp_lfc = unname(shr[ids]),
         pvalue = tt$PValue, padj = tt$FDR,
         mean = rowMeans(edgeR::cpm(y)[ids, , drop = FALSE]) * mean(eff) / 1e6,
         ## raw reads per library in each mix; 0 in every library of one mix
         ## makes the unshrunk logFC ~ +/-10
         mean_A = rowMeans(m[ids, a, drop = FALSE]),
         mean_B = rowMeans(m[ids, !a, drop = FALSE]),
         one_mix_zero = mean_A == 0 | mean_B == 0)
}

## ----------------------------------------------------------------- run ----
res_all <- list(); fc_all <- list(); fc_sens <- list(); mix_reads <- list()

for (plat in platforms) {
  for (lv in levels_to_run) {
    x     <- read_sequins(plat, lv)
    e     <- expected[[lv]]
    r_mix <- take(x$reads, e$id)                       # every sequin in the mix, 0 if absent
    mix   <- sample_mix(colnames(r_mix), sheets[[plat]])
    exp_all  <- setNames(log2(e$B / e$A), e$id)
    ctrl_all <- abs(exp_all) < ctrl_lfc
    det_all  <- detect(r_mix)
    a        <- mix == "A"

    ## estimated reads per mix for every sequin (siblings of outliers, tested or not)
    mix_reads[[length(mix_reads) + 1]] <-
      tibble(platform = plat, level = lv, id = e$id, exp_lfc = unname(exp_all),
             mean_A = rowMeans(r_mix[, a, drop = FALSE]),
             mean_B = rowMeans(r_mix[, !a, drop = FALSE]),
             detected = unname(det_all),
             overdispersion = unname(x$od[e$id]))

    for (mode in count_modes[[lv]]) {
      lab    <- paste(plat_labels[[plat]], lv, mode, "(edgeR)")
      outdir <- file.path(outroot, plat, lv, mode)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

      ## counts edgeR tests, and the counts the filter is applied to (`filter_counts`)
      m_mix <- take(if (mode == "divided") x$divided else x$reads, e$id)
      f_mix <- if (mode == "divided" && filter_counts == "divided") m_mix else r_mix
      keep  <- rowSums(f_mix >= min_count) >= min_samples
      m2    <- m_mix[keep, , drop = FALSE]
      message(sprintf("[%s] %d sequins in the mix, %d detected (%s), %d pass the DE filter, %d ratio-1 controls, %d truly DE",
                      lab, nrow(m_mix), sum(det_all), det_rule_txt, nrow(m2),
                      sum(ctrl_all[keep]), sum(abs(exp_all[keep]) > lfc)))

      ## ---- DE: Mix B vs Mix A ----
      t <- run_edger(m2, mix, ctrl_all[keep], lab) %>%
        mutate(exp_lfc = unname(exp_all[id]), ratio = round(exp_lfc),
               overdispersion = unname(x$od[id])) %>%
        filter(!is.na(pvalue), is.finite(obs_lfc)) %>%
        score_calls() %>%
        mutate(tool = tool, quantifier = quantifier, platform = plat, level = lv,
               count_mode = mode, .before = 1)
      write_csv(t, file.path(outdir, "sequin_de.csv"))
      res_all[[length(res_all) + 1]] <- t

      fc_all[[length(fc_all) + 1]] <- summarise_calls(t, sum(abs(exp_all) > lfc)) %>%
        mutate(tool = tool, quantifier = quantifier, platform = plat, level = lv,
               count_mode = mode,
               detection_rule = det_rule_txt,
               detected       = sum(det_all),                         # >= 3 reads in >= 2 libraries
               detected_any   = sum(rowSums(round(r_mix)) >= 1),      # >= 1 read in total
               median_overdispersion = median(x$od[rownames(m2)], na.rm = TRUE),
               .before = 1)

      ## ---- sensitivity: the same test without R2_18 (gene + isoforms) ----
      drop_s <- grepl(sens_drop_re, e$id)
      keep_s <- keep & !drop_s
      t_s <- run_edger(m_mix[keep_s, , drop = FALSE], mix, ctrl_all[keep_s],
                       paste(lab, "without R2_18")) %>%
        mutate(exp_lfc = unname(exp_all[id]), ratio = round(exp_lfc),
               overdispersion = unname(x$od[id])) %>%
        filter(!is.na(pvalue), is.finite(obs_lfc)) %>%
        score_calls()
      fc_sens[[length(fc_sens) + 1]] <- summarise_calls(t_s, sum(abs(exp_all[!drop_s]) > lfc)) %>%
        mutate(tool = tool, quantifier = quantifier, platform = plat, level = lv,
               count_mode = mode, subset = "without R2_18",
               detection_rule = det_rule_txt,
               detected       = sum(det_all[!drop_s]),
               detected_any   = sum(rowSums(round(r_mix[!drop_s, , drop = FALSE])) >= 1),
               median_overdispersion = median(x$od[rownames(m_mix)[keep_s]], na.rm = TRUE),
               .before = 1)
    }
  }
}

res_all   <- bind_rows(res_all)
fc_all    <- bind_rows(fc_all)
fc_sens   <- bind_rows(fc_sens)
mix_reads <- bind_rows(mix_reads)
write_csv(res_all, file.path(outroot, "sequin_de_results.csv"))
write_csv(fc_all,  file.path(outroot, "fold_change_summary_edger.csv"))
write_csv(fc_sens, file.path(outroot, "fold_change_summary_edger_without_R2_18.csv"))
print(dplyr::select(fc_sens, platform, level, count_mode, tested, truly_de, n_called, wrong_sign,
                    sensitivity, precision, f1, slope, r2), n = Inf, width = Inf)
print(dplyr::select(fc_all, platform, level, count_mode, detected, tested, truly_de, n_called, wrong_sign,
                    sensitivity, sensitivity_in_mix, precision, f1, auc,
                    slope, r2, n_one_mix_zero, r2_no_zero, r2_shrunk, spearman), n = Inf, width = Inf)

## ---------------------------------------------- identification on chrIS ---
## exon records of one GTF on the sequin chromosome
read_gtf_seq <- function(f, seq) {
  if (!file.exists(f))
    stop("GTF not found: ", f, "\n  edit `gtf_stage` at the top of the script")
  g <- read.table(f, sep = "\t", comment.char = "#", quote = "",
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attr"),
                  colClasses = c("character", "character", "character", "integer", "integer",
                                 "character", "character", "character", "character"))
  g %>%
    filter(feature == "exon", seqname == seq) %>%
    mutate(transcript_id = str_match(attr, 'transcript_id "([^"]+)"')[, 2],
           gene_id       = str_match(attr, 'gene_id "([^"]+)"')[, 2])
}

## one row per transcript (exon_len first: summarised start / end would mask it)
tx_spans <- function(ex) {
  ex %>%
    group_by(transcript_id, gene_id, strand) %>%
    summarise(exonic_len = sum(end - start + 1L), exons = n(),
              start = min(start), end = max(end), .groups = "drop")
}

prf <- function(tp, fn, fp) {
  p  <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  r  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  f1 <- if (is.na(p) || is.na(r) || p + r == 0) NA_real_ else 2 * p * r / (p + r)
  c(precision = p, recall = r, f1 = f1)
}

## header kept when no novel model exists anywhere
novel_template <- tibble(platform = character(), protocol = character(), stage = character(),
                         transcript_id = character(), gene_id = character(), strand = character(),
                         exonic_len = integer(), exons = integer(), start = integer(), end = integer(),
                         same_strand_overlap = character(), opposite_strand_overlap = character(),
                         candidate_mirror = logical())

## sequin lengths: Anaquin's table when it has them, else the union of all
## protocols' validated GTFs (a validated GTF only keeps sequins with reads)
seq_len_tbl <- if (!is.null(anaquin_len))
  tibble(transcript_id = names(anaquin_len), exonic_len = unname(anaquin_len)) else
  map_dfr(platforms, function(pl)
    tx_spans(filter(read_gtf_seq(gtf_stage$validated[[pl]], sequin_seq),
                    !grepl("^BambuTx", transcript_id)))) %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  dplyr::select(transcript_id, exonic_len)
n_nolen <- sum(!iso_mix$id %in% seq_len_tbl$transcript_id)
if (n_nolen > 0)
  message(sprintf("%d sequin isoforms have no length: %s", n_nolen,
                  paste(setdiff(iso_mix$id, seq_len_tbl$transcript_id), collapse = ", ")))

ident_all <- list(); novel_seq_all <- list(); per_iso_all <- list(); by_conc_all <- list()
for (plat in platforms) {
  prot <- plat_labels[[plat]]

  ## positives = ALL sequin isoforms in the mix (Anaquin, 164); one missing
  ## from a count table had no reads there
  xv    <- read_sequins(plat, "transcript")              # curated set, undivided reads
  m_val <- xv$reads
  sm    <- round(take(m_val, iso_mix$id))                # rows in iso_mix order, reads rounded
  det   <- detect(sm)
  tp    <- sum(det); fn <- sum(!det)
  if (fn > 0)
    message(sprintf("[%s] sequin isoforms not detected (%s; %d of %d): %s", prot, det_rule_txt,
                    fn, length(det), paste(head(names(det)[!det], 30), collapse = ", ")))

  ## does curation remove annotated sequins? the same quantifier before curation
  ## (oarfish: every model in the index; Bambu: its raw table)
  m_rawt  <- tryCatch(read_sequins(plat, "transcript", curated = FALSE)$reads,
                      error = function(e) { message("[", prot, "] ", conditionMessage(e)); NULL })
  det_raw <- if (is.null(m_rawt)) NA_integer_ else sum(detect(take(m_rawt, iso_mix$id)))
  if (!is.null(m_rawt) && det_raw > tp)
    message(sprintf("[%s] %d sequin isoforms are detected before curation but not after",
                    prot, det_raw - tp))

  ## ---- dose-response: per library, >= 1 read vs that library's mix concentration ----
  mix  <- sample_mix(colnames(sm), sheets[[plat]])
  cmat <- vapply(as.character(mix), function(mx) if (mx == "A") iso_mix$A else iso_mix$B,
                 numeric(nrow(iso_mix)))                 # sequins x libraries, same layout as sm
  hit  <- as.vector(sm >= 1)
  lc   <- log2(as.vector(cmat))
  lfit <- suppressWarnings(glm(hit ~ lc, family = binomial))
  l_b0 <- unname(coef(lfit)[1]); l_b1 <- unname(coef(lfit)[2])
  lod  <- if (is.finite(l_b1) && l_b1 > 0) 2^((qlogis(lod_p) - l_b0) / l_b1) else NA_real_
  by_conc_all[[plat]] <- tibble(conc = as.vector(cmat), hit = hit) %>%
    mutate(band = cut(conc, conc_breaks, right = FALSE, dig.lab = 5)) %>%
    group_by(band) %>%
    summarise(pairs = n(), frac_libraries_with_reads = mean(hit), .groups = "drop") %>%
    mutate(platform = plat, protocol = prot, .before = 1)

  ## annotated sequins on chrIS (from the validated GTF) - also checks the name
  val_ex <- read_gtf_seq(gtf_stage$validated[[plat]], sequin_seq)
  ann    <- tx_spans(filter(val_ex, !grepl("^BambuTx", transcript_id)))
  if (!nrow(ann))
    stop(sprintf("[%s] no annotated transcripts on '%s' in the validated GTF - set `sequin_seq`",
                 prot, sequin_seq))

  ## capture efficiency = reads per amol/ul, mean of the two mixes, log2,
  ## centred on the median DETECTED isoform (0 = typical; -5 = 32x fewer reads
  ## than a typical sequin at the same input). Undetected isoforms: NA.
  rA  <- rowMeans(sm[, mix == "A", drop = FALSE])
  rB  <- rowMeans(sm[, mix == "B", drop = FALSE])
  cap <- log2((rA / iso_mix$A + rB / iso_mix$B) / 2)
  cap[!det | !is.finite(cap)] <- NA
  cap <- cap - median(cap, na.rm = TRUE)
  len <- seq_len_tbl$exonic_len[match(iso_mix$id, seq_len_tbl$transcript_id)]
  cap_rho <- suppressWarnings(cor(cap, log2(len), method = "spearman", use = "complete.obs"))

  per_iso_all[[plat]] <- tibble(platform = plat, protocol = prot, id = iso_mix$id,
                                gene = sub("_[0-9]+$", "", iso_mix$id),
                                len = len,
                                conc_A = iso_mix$A, conc_B = iso_mix$B,
                                exp_lfc = log2(iso_mix$B / iso_mix$A),
                                reads_A = rA, reads_B = rB,
                                libs_with_reads = rowSums(sm >= 1),
                                detected = unname(det),
                                detected_any = rowSums(sm) >= 1,
                                capture_log2 = cap,
                                overdispersion = unname(xv$od[iso_mix$id]))

  row <- tibble(quantifier = quantifier, platform = plat, protocol = prot,
                annotated = nrow(iso_mix),                     # 164 (Anaquin)
                detection_rule = det_rule_txt,
                in_curated_table = length(intersect(rownames(m_val), iso_mix$id)),
                in_uncurated_table = if (is.null(m_rawt)) NA_integer_ else
                                   length(intersect(rownames(m_rawt), iso_mix$id)),
                detected_before_curation = det_raw,            # same rule and quantifier
                detected = tp, not_detected = fn,
                detected_any = sum(rowSums(sm) >= 1),  # >= 1 read in total
                sequin_reads_per_library = mean(colSums(sm)),  # depth confounds the detection limit
                lod70_amol_ul = lod,                           # >= 1 read in 70% of libraries
                lod70_logistic_slope = l_b1,
                capture_length_rho = cap_rho)                  # Spearman, detected isoforms
  for (st in c("raw", "validated")) {
    ex  <- if (st == "validated") val_ex else read_gtf_seq(gtf_stage$raw[[plat]], sequin_seq)
    nov <- tx_spans(filter(ex, grepl("^BambuTx", transcript_id)))
    if (nrow(nov)) {
      ## overlapping annotated sequins, same vs opposite strand
      ov <- function(i, same) {
        hit <- ann$start <= nov$end[i] & ann$end >= nov$start[i] &
               (if (same) ann$strand == nov$strand[i] else ann$strand != nov$strand[i])
        paste(ann$transcript_id[hit], collapse = ";")
      }
      nov <- nov %>%
        mutate(same_strand_overlap     = map_chr(seq_len(n()), ov, same = TRUE),
               opposite_strand_overlap = map_chr(seq_len(n()), ov, same = FALSE),
               candidate_mirror        = opposite_strand_overlap != "" & same_strand_overlap == "")
      message(sprintf("[%s] %s: %d novel model(s) on %s (%d candidate mirror(s))", prot, st,
                      nrow(nov), sequin_seq, sum(nov$candidate_mirror)))
      novel_seq_all[[length(novel_seq_all) + 1]] <-
        mutate(nov, platform = plat, protocol = prot, stage = st, .before = 1)
    }
    s <- prf(tp, fn, nrow(nov))
    row[[paste0("novel_", st)]]     <- nrow(nov)
    row[[paste0("precision_", st)]] <- unname(s["precision"])
    row[[paste0("f1_", st)]]        <- unname(s["f1"])
    row[["recall"]]                 <- unname(s["recall"])
  }
  ident_all[[plat]] <- row
}
ident_all     <- bind_rows(ident_all)
novel_seq_all <- bind_rows(novel_template, novel_seq_all)
per_iso_all   <- bind_rows(per_iso_all)
by_conc_all   <- bind_rows(by_conc_all)
write_csv(ident_all,     file.path(seqroot, "sequin_identification.csv"))
write_csv(novel_seq_all, file.path(seqroot, "novel_models_on_sequins.csv"))   # header only if none
write_csv(per_iso_all,   file.path(seqroot, "sequin_per_isoform.csv"))
write_csv(by_conc_all,   file.path(seqroot, "sequin_detection_by_conc.csv"))
print(ident_all, width = Inf)
print(by_conc_all, n = Inf)

## ---- fold-change outliers and their sibling isoforms (transcript + gene level) ----
## Two ways a fold change can go wrong although a per-molecule bias cancels
## in a ratio:
##   isoform level - reads assigned to the wrong sibling isoform, and the
##     error depends on the siblings' relative abundance, which differs
##     between Mix A and Mix B (R2_18 in every protocol)
##   gene level    - the gene total is a capture-weighted sum of its isoforms;
##     when siblings with very different capture (capture_log2, len) change in
##     opposite directions, the gene follows the well-captured one (PacBio:
##     R1_13, R1_91, R2_41)
tx_res <- res_all %>% filter(level == "transcript", count_mode == primary_mode) %>%
  mutate(gene = sub("_[0-9]+$", "", id))
out_genes <- bind_rows(
  tx_res %>% filter(abs(obs_lfc - exp_lfc) > outlier_lfc | outcome == "wrong sign"),
  res_all %>% filter(level == "gene", count_mode == primary_mode,
                     abs(obs_lfc - exp_lfc) > outlier_lfc | outcome %in% c("wrong sign", "false positive")) %>%
    mutate(gene = id)) %>%
  distinct(platform, gene)
gene_res <- res_all %>% filter(level == "gene", count_mode == primary_mode) %>%
  transmute(platform, gene = id, gene_exp_lfc = exp_lfc, gene_obs_lfc = obs_lfc,
            gene_outcome = outcome)
outliers <- mix_reads %>%
  filter(level == "transcript") %>%
  mutate(gene = sub("_[0-9]+$", "", id)) %>%
  semi_join(out_genes, by = c("platform", "gene")) %>%
  left_join(dplyr::select(tx_res, platform, id, obs_lfc, padj, outcome), by = c("platform", "id")) %>%
  left_join(gene_res, by = c("platform", "gene")) %>%
  left_join(dplyr::select(per_iso_all, platform, id, len, capture_log2), by = c("platform", "id")) %>%
  mutate(protocol = unname(plat_labels[platform]), tested = !is.na(outcome)) %>%
  dplyr::select(gene, id, protocol, len, capture_log2, overdispersion, exp_lfc, obs_lfc, padj, outcome, tested,
                mean_A, mean_B, gene_exp_lfc, gene_obs_lfc, gene_outcome) %>%
  arrange(gene, id, protocol)
write_csv(outliers, file.path(outroot, "sequin_fc_outliers.csv"))
if (nrow(outliers)) {
  message("Fold-change outliers (|obs - exp| > ", outlier_lfc, " or wrong sign) and their siblings:")
  print(outliers, n = Inf, width = Inf)
}

message("Done: ", outroot)
