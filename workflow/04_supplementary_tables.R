## ===========================================================================
## 04_supplementary_tables.R - builds the Supplementary Tables workbook from the
## current results (oarfish-integrated runs). Rerunnable: every table is read
## from the files the analysis scripts write; nothing is typed in by hand
## except titles, descriptions and the SAMtools and MultiQC versions of ST3
## (added only while the pipeline's versions file lacks them).
## Descriptions refer to other tables by their id below, as {STn}; after the drop
## and renumbering, {STn} becomes "Supplementary Table <its new number>".
##
## ST1  library metrics (NanoComp: raw, restranded, aligned) + Restrander stats
## ST2  computational resources per run and per step (01_computational_resources.R)
## ST3  software versions (pipeline versions file + downstream R packages)
## ST4  spike-in identification (LongBench, Dong)
## ST5  SIRV abundance per isoform (LongBench)
## ST6  ERCC dose response per transcript (LongBench)
## ST7  sequin differential abundance, with "without R2_18" rows (LongBench, Dong)
## ST8  quantifier comparison, Bambu vs oarfish (LongBench, Dong)
## ST9  read-level support (LongBench, Dong; Dong also by GffCompare class)
## ST10 novel transcripts by biotype and class code (LongBench, Dong)
## ST11 unique intron chain sharing and transcript-end differences (LongBench)
## ST12 classification agreement of shared intron chains (LongBench)
## ST13 DE library sizes and TMM factors (Dong)
## ST14 DE results, gene and transcript level (Dong)
## ST15 cis pairs: every DE-tested pair within 1 Mb with rho and DE calls (Dong)
##
## Output: <project>/pulposeq_supplementary_tables.xlsx (rewritten on every run)
## Needs: tidyverse, openxlsx, yaml-free (the versions file is parsed here)
## ===========================================================================
suppressPackageStartupMessages({ library(tidyverse); library(openxlsx) })
project <- "C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026"
ds      <- file.path(project, "workflow/downstream")
lb      <- file.path(ds, "longbench"); dg <- file.path(ds, "dong")
out     <- file.path(project, "pulposeq_supplementary_tables.xlsx")
options(openxlsx.datetimeFormat = "yyyy-mm-dd hh:mm:ss")  # dates would otherwise show as ####
runs    <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio", dong = "ONT cDNA")
dset    <- c(ont_cdna = "LongBench", ont_drs = "LongBench", pb = "LongBench", dong = "Dong")
rd      <- function(...) read_csv(file.path(...), show_col_types = FALSE, guess_max = 1e5)
rt      <- function(...) read_tsv(file.path(...), show_col_types = FALSE, guess_max = 1e5)
tag     <- function(d, dataset) bind_cols(tibble(Dataset = rep(dataset, nrow(d))), d)
## the workbook uses the labels of the figures and the text; the result files
## keep their own and are relabelled here only
prot_lab <- function(x) {            # UIC scripts: ONTd / ONTc / PB -> ONT DRS / ONT cDNA / PacBio
  x <- sub("_only$", " only", x)
  x <- sub("^shared_all$", "all 3 protocols", x)
  x <- gsub("\\bONTd\\b", "ONT DRS", x, perl = TRUE)
  x <- gsub("\\bONTc\\b", "ONT cDNA", x, perl = TRUE)
  x <- gsub("\\bPB\\b", "PacBio", x, perl = TRUE)
  gsub("+", " + ", x, fixed = TRUE)
}
uic_tab <- function(d) {             # a table of the UIC scripts: "UIC" (not "UJC") and protocol names
  names(d) <- prot_lab(sub("^UJC ", "UIC ", names(d)))
  mutate(d, across(any_of(c("Protocol", "Protocols", "Category", "Group", "Sharing category")), prot_lab))
}
rs_lab <- function(d) {              # read support: biotypes and support classes as in Fig. 1g and 2a
  bt <- c("predicted protein-coding" = "novel protein-coding", "novel lncRNA candidate" = "novel lncRNA",
          "predicted non-coding" = "novel non-coding", "other biotypes" = "annotated other")
  sp <- c("full-length + unique" = "full-length + unique reads", "unique only" = "unique reads",
          "full-length only" = "full-length reads", "neither" = "partial & ambiguous reads")
  mutate(d, across(any_of("biotype"), ~ ifelse(.x %in% names(bt), unname(bt[.x]), .x)),
            across(any_of("support"), ~ ifelse(.x %in% names(sp), unname(sp[.x]), .x)))
}

ST <- list()                       # id -> list(title, desc, blocks = list(name = df))
## the workbook keeps summaries only; per-feature tables are in the benchmarking
## repository with the other results. Remove an id or a block letter here to
## bring it back. Tables are renumbered in order after the drop.
drop_tables <- c("ST5", "ST6", "ST14", "ST15")      # SIRV / ERCC per isoform, DE results, cis pairs
drop_blocks <- list(ST7 = c("c", "d"),               # per-sequin results
                    ST8 = c("e", "f", "g", "h"),     # Jaccard, overdispersion, disagreement, Dong per isoform
                    ST12 = "b")                      # itemised disagreements
add <- function(id, title, desc, ...) {
  if (id %in% drop_tables) return(invisible())
  b <- list(...)
  if (!is.null(drop_blocks[[id]])) b <- b[!sub(" \\|.*$", "", names(b)) %in% drop_blocks[[id]]]
  ST[[id]] <<- list(title = title, desc = desc, blocks = b)
}

## ---- ST1: NanoComp NanoStats (general summary) per library and stage -------
## parsed once by 03_sequencing_metrics.R (which also draws Extended Data Fig. 2)
sm   <- file.path(ds, "sequencing-metrics/results")
prot <- function(run) sub("^(LongBench|Dong) ", "", run)
st1  <- rd(sm, "sequencing_metrics_per_library.csv") %>%
  transmute(Dataset = dataset, Protocol = prot(run), Stage = stage, library, sample, across(-c(dataset, run, stage, library, sample)))
rs1  <- rd(sm, "restrander_orientation.csv") %>%
  mutate(Dataset = ifelse(grepl("^Dong", run), "Dong", "LongBench"), Protocol = prot(run), .before = 1) %>% select(-run)
## Restrander's total_reads = plus + minus + undetermined: every read of the raw
## FASTQ. NanoComp's raw count in (a) is a few reads lower: NanoComp leaves some
## records out, all of them undetermined reads (plus + minus equals NanoComp's
## count of the oriented reads). Checked here; the description quotes the gap.
nc_n <- st1 %>% select(Dataset, Protocol, sample, Stage, n = `Number of reads`) %>%
  pivot_wider(names_from = Stage, values_from = n)
chk1 <- rs1 %>% inner_join(nc_n, by = c("Dataset", "Protocol", "sample")) %>%
  mutate(gap = total_reads - `Raw reads`, gap_pct = 100 * gap / total_reads)
if (nrow(chk1) < nrow(rs1)) message("  ! ST1: ", nrow(rs1) - nrow(chk1), " Restrander libraries not matched to NanoComp counts")
if (any(chk1$plus + chk1$minus + chk1$undetermined != chk1$total_reads, na.rm = TRUE)) warning("ST1: total_reads differs from plus + minus + undetermined")
if (any(chk1$plus + chk1$minus != chk1$`Restranded reads`, na.rm = TRUE)) warning("ST1: plus + minus differs from NanoComp's count of oriented reads")
if (any(chk1$gap < 0, na.rm = TRUE)) warning("ST1: NanoComp counts more raw reads than Restrander for some libraries")
gap_txt <- chk1 %>% group_by(Dataset) %>%
  summarise(m = max(gap_pct, na.rm = TRUE),
            t = sprintf("%s %s-%s reads, %s of total_reads", first(Dataset),
                        formatC(min(gap, na.rm = TRUE), format = "d", big.mark = ","),
                        formatC(max(gap, na.rm = TRUE), format = "d", big.mark = ","),
                        if (m < 0.01) "<0.01%" else sprintf("up to %.2f%%", m)), .groups = "drop") %>%
  pull(t) %>% paste(collapse = "; ")
rs1 <- rs1 %>% mutate(oriented_fraction = round(oriented_fraction, 3))
add("ST1", "Sequencing, orientation and alignment metrics per library (LongBench and Dong).",
    paste("Unit: sequencing libraries. (a) NanoComp NanoStats general summary per library at each stage the pipeline reports:",
          "raw reads as supplied, reads oriented by Restrander (ONT cDNA only; undetermined reads are not in this file) and minimap2 alignments.",
          "Number of reads counts reads at the first two stages and, for the alignments, BAM records of mapped reads (primary and supplementary alignments; secondary alignments and unmapped reads are not counted).",
          "library is the file name at each stage; sample identifies the same library across stages.",
          "Fraction of bases aligned = Total bases aligned / Total bases, recomputed from these totals (NanoComp prints it to one decimal).",
          "LongBench: 8 cell lines x 3 protocols, one library each; Dong: H1975 and HCC827, 3 cultures each, ONT PCR-cDNA.",
          "(b) Restrander orientation per ONT cDNA library. total_reads = plus + minus + undetermined, every read of the raw FASTQ:",
          "plus and minus are the reads oriented to the + or - strand (together, the oriented reads of (a)); undetermined (?) are the reads Restrander could not orient, written to a separate file and not analysed;",
          "oriented_fraction = (plus + minus) / total_reads, checked against restrand_min_frac.",
          paste0("total_reads exceeds the raw read count in (a) (", gap_txt, "): NanoComp leaves these records out, and all of them are undetermined reads, since plus + minus equals NanoComp's count of oriented reads."),
          "TSO-TSO and RTP-RTP count artefactual reads (the same primer at both ends) within these classes. Extended Data Fig. 2."),
    `a | NanoComp summaries` = st1, `b | Restrander orientation (ONT cDNA)` = rs1)

## ---- ST2: computational resources ------------------------------------------
res <- file.path(ds, "resources/results")
st2a <- rd(res, "resources_per_run.csv")
names(st2a)[names(st2a) == "critical_path_h"] <- "execution_time_h"   # column named as in the text and ED1
resumed <- st2a$run[st2a$cached_tasks > 0]                              # runs resumed with cached tasks
add("ST2", "Computational resources of pulposeq per run and per pipeline step (LongBench ONT cDNA, ONT DRS, PacBio; Dong ONT cDNA).",
    paste("From the newest Nextflow execution trace of each run (by the date and time in its name). Queue time is excluded: task start = submit + duration - realtime.",
          "Execution time (execution_time_h) = idealized critical-path (uncontended) time: each step takes its slowest task's realtime (all its tasks in parallel) and the run the longest chain of steps through its Nextflow DAG (critical_path, the chain itself); this removes HPC queue and partition limits.",
          "For transparency the time observed on the cluster is also given (elapsed_h: time with at least one task executing, queue excluded; span_h: first start to last end; idle_h = span_h - elapsed_h, time with no task executing;",
          "elapsed_h_with_queue: first submission to last end, queue included). Per step: slowest_task_h, realtime_h (union of its task intervals) and realtime_sum_h.",
          "CPU time = realtime x %CPU; peak memory = largest task peak RSS; CPU cores = the busiest task's %CPU / 100; data read / written = summed rchar / wchar (intermediate files included).",
          "Per run: final_tasks = completed + cached tasks; failed_attempts = failed task attempts, with their CPU time in cpu_h_failed (empty when the trace records no realtime or %CPU for any of them; failed_without_metrics counts the attempts without these values);",
          "queue_h_summed = summed queue time of the tasks (submit to start); peak_step and peak_cores_step = the steps of the largest task by peak RSS and by %CPU.",
          if (length(resumed)) paste0(paste(resumed, collapse = " and "), if (length(resumed) > 1) " were" else " was",
                                      " resumed: cached tasks are timed where they first ran, so span_h, idle_h and elapsed_h_with_queue include the interval before the resume, which elapsed_h and execution_time_h do not.") else "",
          "Per step: completed and cached tasks (cached ones with the metrics of the run that computed them), summed over samples; first_start and last_end in the cluster's clock time, as written in the trace (Extended Data Fig. 1)."),
    `a | Per run` = st2a,
    `b | Per pipeline step` = rd(res, "resources_per_step.csv"))

## ---- ST3: software versions -------------------------------------------------
yml <- function(f) {
  l <- readLines(f, warn = FALSE); proc <- NA_character_; out <- list()
  for (x in l) {
    if (grepl("^[^ ].*:\\s*$", x)) proc <- sub(":\\s*$", "", x)
    ## any indent: the Workflow block (Nextflow, pulposeq) is indented by 4
    else if (grepl("^ +[^ ].*:", x)) out[[length(out) + 1]] <- tibble(Process = proc,
      Software = trimws(sub(":.*$", "", x)), Version = trimws(sub("^[^:]*:", "", x)))
  }
  ## NanoComp's version string ends in a stray ")" ("1.25.6)") in the versions
  ## file: drop a closing bracket that has no opening one
  bind_rows(out) %>% mutate(Version = ifelse(grepl("\\)$", Version) & !grepl("\\(", Version),
                                             sub("\\)$", "", Version), Version))
}
v_pipe <- map_dfr(names(runs), ~ yml(file.path(project, .x, "pipeline_info", "pulposeq_software_mqc_versions.yml"))) %>%
  distinct() %>% arrange(Process, Software)
## oarfish is in the versions file and edgeR in the downstream rows. SAMtools
## (sorting and indexing in MINIMAP2_ALIGN) and MultiQC are not recorded there:
## added by hand, and only while the file lacks them, so they never appear twice
v_hand <- tibble(Process = c("MINIMAP2_ALIGN", "MULTIQC"), Software = c("samtools", "multiqc"),
                 Version = c("1.23.1", "1.35")) %>%
  filter(!tolower(Software) %in% tolower(v_pipe$Software))
v_pipe <- bind_rows(v_pipe, v_hand) %>% arrange(Process, Software)
pk <- c("edgeR", "limma", "tidyverse", "ggplot2", "patchwork", "ggh4x", "ggrepel", "circlize", "ComplexHeatmap",
        "plotgardener", "rtracklayer", "GenomicRanges", "GenomicAlignments", "Rsamtools", "SummarizedExperiment",
        "isoformic", "svglite", "showtext", "magick", "openxlsx", "data.table")
v_ds <- tibble(Process = "downstream (R)", Software = pk,
               Version = vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), "")) %>%
  filter(!is.na(Version))
v_ds <- bind_rows(tibble(Process = "downstream (R)", Software = "R", Version = paste(R.version$major, R.version$minor, sep = ".")),
                  v_ds)  # R first
## pipeline processes named on their first row only (blank below); downstream
## rows keep the name on every row
v_pipe <- v_pipe %>% mutate(Process = ifelse(Process == lag(Process, default = ""), NA_character_, Process))
v_all <- bind_rows(v_pipe, v_ds)
add("ST3", "Software versions of the pulposeq pipeline and of the downstream analyses.",
    paste(paste0("Pipeline rows: software of each Nextflow process as recorded in the consolidated versions file every pulposeq run writes (union of the four runs)",
                 if (nrow(v_hand)) paste0("; ", paste(sprintf("%s %s (%s)", c(samtools = "SAMtools", multiqc = "MultiQC")[v_hand$Software],
                                                               v_hand$Version, v_hand$Process), collapse = " and "),
                                          ", which the file does not record, were added by hand.") else "."),
          "oarfish v0.10.3 runs through the nf-core oarfish/readmode module; edgeR 4.99.12 (catchOarfish, catchOarfishGene) was used downstream.",
          "Downstream rows: R packages of the benchmarking and case-study analyses (Bioconductor 3.24, development branch), as installed when this workbook was built."),
    `Software versions` = v_all)

## ---- ST4: spike-in identification -------------------------------------------
id_lb <- rd(lb, "figure-panels/results/figure_tables/spikein_identification_by_set.csv")
## E0 and long SIRVs together (84 isoforms), with the same formulas as the
## per-set rows
id_sirv <- id_lb %>% filter(set %in% c("SIRV E0", "long SIRV")) %>%
  group_by(protocol) %>%
  summarise(across(c(detected, expected, fp_raw, fp_val), sum), .groups = "drop") %>%
  mutate(set = "SIRV (E0 + long)", recall = detected / expected,
         precision = detected / (detected + fp_val), f1 = 2 * precision * recall / (precision + recall),
         precision_raw = detected / (detected + fp_raw), f1_raw = 2 * precision_raw * recall / (precision_raw + recall))
st4a <- bind_rows(id_lb, id_sirv[names(id_lb)]) %>%
  mutate(set = sub("^Sequins$", "sequin", set)) %>%
  arrange(match(set, c("SIRV E0", "long SIRV", "SIRV (E0 + long)", "ERCC", "sequin")),
          match(protocol, c("ONT cDNA", "ONT DRS", "PacBio")))
## Dong: identification columns named as in (a); the differential abundance
## columns of the same row get _de (sensitivity is the recall of DE)
st4b <- rd(dg, "figure-panels/results/figure_tables/sequin_identification_dong.csv") %>%
  mutate(set = sub("^Sequins$", "sequin", set)) %>%
  rename(any_of(c(precision = "precision_id", f1 = "f1_id", recall_de = "sensitivity", r2_de = "r2")))
if (!all(c("recall_de", "precision_de", "f1_de", "r2_de") %in% names(st4b)))
  message("  ! ST4: Dong differential abundance columns not found to relabel")
add("ST4", "Spike-in identification per protocol: detected against expected isoforms, false positives, precision, recall and F1 (LongBench and Dong).",
    paste("Unit: spike-in set x protocol. Detected = >=3 estimated reads (oarfish, rounded per library) in >=2 libraries (ERCCs: >=1 read in >=2 libraries).",
          "False positives = novel models assembled on spike-in sequences, in Bambu's output and after pulposeq validation.",
          "Fig. 1b and the text report the two SIRV sets separately (F1 = 0.95-1.00); the combined row SIRV (E0 + long, 84 isoforms) is given for completeness.",
          "Sequins: the 164 isoforms of the Anaquin mixture table (R2_18 and R2_66 included; Supplementary Note 1).",
          "Recall = detected / expected (the sensitivity of identification).",
          "In (b), recall_de, precision_de, f1_de and r2_de are sequin isoform differential abundance (recall_de = sensitivity, as in {ST7}; r2_de = R2 of observed against expected log2FC),",
          "with undivided counts for both quantifiers; the sequin differential abundance F1 quoted in the text uses oarfish counts divided by RTA overdispersion ({ST7})."),
    `a | LongBench` = st4a,
    `b | Dong (sequins; Bambu and oarfish counts)` = st4b)

## ---- ST5, ST6: SIRV and ERCC -------------------------------------------------
sirv <- map_dfr(c("ont_cdna", "ont_drs", "pb"), ~ rd(lb, "spike-ins/results/sirv_ercc", .x, "sirv_per_isoform.csv") %>%
                  mutate(Protocol = runs[[.x]], .before = 1))
add("ST5", "SIRV abundance accuracy and length bias per isoform (LongBench).",
    paste("Unit: SIRV isoforms (69 E0 and 15 long SIRVs, all 60 amol/ul) x protocol. Per library, log2 ratio = log2((CPM + 1) / expected CPM), CPM within these 84 transcripts (expected share 1/84);",
          "summarized per isoform as the median and SD across libraries. Within two-fold: |median| <= 1; lost: median < -3; unstable: SD > 1. Counts: oarfish, undivided.",
          "(b) Length classes; (c) regression of the log2 ratio on log2 length."),
    `a | Per isoform` = sirv,
    `b | Length classes` = rd(lb, "spike-ins/results/sirv_ercc/length_class_summary.csv"),
    `c | Length bias` = rd(lb, "spike-ins/results/sirv_ercc/length_bias.csv"))
add("ST6", "ERCC dose response per transcript (LongBench).",
    paste("Unit: ERCC transcripts (92, SIRV ERCC Mix 1) x protocol. Detected = >=1 read in >=2 libraries. Dose response: linear fit of log2 CPM on log2 input among detected ERCCs;",
          "limit of detection (LOD70) = input at which 70% of ERCCs are detected (logistic fit, as Anaquin). (b) Detection by concentration."),
    `a | Per transcript` = rd(lb, "spike-ins/results/sirv_ercc/ercc_per_transcript.csv"),
    `b | Detection grid` = rd(lb, "spike-ins/results/sirv_ercc/detection_grid.csv"))

## ---- ST7: sequin differential abundance -------------------------------------
lbsq <- file.path(lb, "spike-ins/results/sequins/edger")
## the result files' quantifier column already says which count table was used
## ("bambu" written as in the other tables)
dgs  <- function(dir) bind_rows(rd(dg, dir, "summary.csv") %>% mutate(sequins = "all 164"),
                                rd(dg, dir, "summary_without_R2_18.csv") %>% mutate(sequins = "without R2_18")) %>%
  mutate(quantifier = recode(quantifier, bambu = "Bambu"))
st7a <- bind_rows(rd(lbsq, "fold_change_summary_edger.csv") %>% mutate(sequins = "all 164"),
                  rd(lbsq, "fold_change_summary_edger_without_R2_18.csv") %>% mutate(sequins = "without R2_18")) %>%
  mutate(platform = unname(runs[platform])) %>% rename(protocol = platform)
## Dong's rmse of the tested log2FC is the `rmse` of (a); that of the shrunken
## log2FC is named like slope_shrunk / r2_shrunk
st7b <- bind_rows(dgs("sequins/results"), dgs("sequins/results_bambu")) %>%
  rename(any_of(c(rmse = "rmse_test_lfc", rmse_shrunk = "rmse_display_lfc")))
## a `subset` column of the result files that labels the same two groups as
## `sequins` is a duplicate: dropped
drop_subset <- function(d) {
  if ("subset" %in% names(d) && nrow(distinct(d, sequins, subset)) == n_distinct(d$sequins) &&
      n_distinct(d$subset) == n_distinct(d$sequins)) select(d, -subset) else d
}
## (a) and (b) with the same column names and order: a count named without the
## "n_" prefix in one block and with it in the other takes the prefixed name
## (called -> n_called); (b) then follows the column order of (a)
n_names <- function(d, other) {
  for (x in setdiff(names(d), names(other)))
    if (paste0("n_", x) %in% setdiff(names(other), names(d))) names(d)[names(d) == x] <- paste0("n_", x)
  d
}
## accuracy of the calls, as in the sequin DE rows of {ST8} (a), for both datasets
add_acc <- function(d) mutate(d, accuracy = (true_pos + true_neg) / tested, .after = f1)
st7a <- add_acc(st7a); st7b <- add_acc(st7b)
st7a <- drop_subset(st7a); st7b <- drop_subset(st7b)
st7a <- n_names(st7a, st7b); st7b <- n_names(st7b, st7a)
st7b <- select(st7b, any_of(names(st7a)), everything())
only7 <- c(setdiff(names(st7a), names(st7b)), setdiff(names(st7b), names(st7a)))
if (length(only7)) message("  ST7: columns in only one of (a) and (b): ", paste(only7, collapse = ", "))
add("ST7", "Sequin differential abundance against the known mixture ratios (LongBench and Dong), with and without R2_18.",
    paste("edgeR v4 glmQLFit (robust) -> glmTreat (H0: |log2FC| <= 0.5), FDR < 0.05; library sizes and TMM factors from the ratio-1 control sequins (|expected log2FC| < 0.25).",
          "Tested: sequins with >=10 estimated reads (before division) in >=2 (LongBench) or >=3 (Dong) libraries. Truly differential: |expected log2FC| > 0.5; wrong-sign calls are false discoveries.",
          "Counts: oarfish divided by RTA overdispersion (count_mode divided) and undivided; Bambu (undivided). LongBench: Mix B (n = 3) vs Mix A (n = 5); Dong: Mix A (H1975) vs Mix B (HCC827).",
          "'without R2_18' rows: the same tests without R2_18 (its two isoforms at transcript level, the R2_18 gene at gene level). Anaquin lists R2_18_2 at the same concentration as R2_18_1, but no read covered its unique first exon in any LongBench Mix A library, and its coverage there was at most 2 in the Dong Mix A libraries (Supplementary Note 1).",
          "Columns: detected = >=3 reads in >=2 libraries (estimated reads rounded per library, as {ST4}); detected_any = >=1 read summed over all libraries (empty in the Dong 'without R2_18' rows);",
          "median_overdispersion = median RTA overdispersion of the tested sequins (empty for Bambu); tested = sequins passing the filter; truly_de = tested sequins that are truly differential; n_called = FDR < 0.05;",
          "true_pos = called, truly differential and of the right sign; false_pos = called, not truly differential; wrong_sign = called, truly differential, wrong sign; false_neg = truly differential, not called; true_neg = neither;",
          "sensitivity (= recall) = true_pos / truly_de; fdr = (false_pos + wrong_sign) / n_called; precision = true_pos / n_called; f1 from precision and sensitivity; accuracy = (true_pos + true_neg) / tested; auc = area under the ROC curve of -log10(p) for truly differential sequins;",
          "slope, intercept (LongBench only) and r2 = linear regression of observed (unshrunk) on expected log2FC; rmse = root mean square of observed - expected; spearman = Spearman correlation of observed and expected; median_abs_err = median |observed - expected|;",
          "n_one_mix_zero = tested sequins with no reads in every library of one mix (unshrunk log2FC near +/-10), left out of slope_no_zero and r2_no_zero; slope_shrunk and r2_shrunk = the same fit on the shrunken log2FC (predFC, prior count 2);",
          "truly_de_in_mix = truly differential sequins in the whole mix, tested or not, and sensitivity_in_mix = true_pos / truly_de_in_mix. In (b), rmse_shrunk = the rmse of the shrunken log2FC.",
          "Per-sequin results are in the benchmarking repository."),
    `a | LongBench summary` = st7a,
    `b | Dong summary` = st7b,
    `c | LongBench per sequin` = rd(lbsq, "sequin_de_results.csv"),
    `d | Dong per sequin` = bind_rows(rd(dg, "sequins/results/sequin_de_results.csv"), rd(dg, "sequins/results_bambu/sequin_de_results.csv")))

## ---- ST8: quantifier comparison ----------------------------------------------
lq <- file.path(lb, "spike-ins/results/quantifier_comparison"); dq <- file.path(dg, "figure-panels/results/quantifier_comparison")
## ERCC reads per library and LOD70 on the linear scale, as in Fig. 1d and the
## text (the result file holds log10, the scale Extended Data Fig. 4b plots)
q_spk <- rd(lq, "q_metrics_spikeins.csv") %>%
  mutate(lg = grepl("(log10)", metric, fixed = TRUE),
         value = ifelse(lg, 10^value, value),
         metric = trimws(gsub("[[:space:]]+", " ", sub("(log10)", "", metric, fixed = TRUE)))) %>%
  select(-lg) %>%
  ## the SIRV F1 here is that of E0 and long SIRVs together (84 isoforms): named
  ## so, as it is not the per-set F1 of the text (0.95-1.00)
  mutate(metric = ifelse(grepl("^SIRV F1", metric) & !grepl("E0", metric, fixed = TRUE),
                         sub("^SIRV F1", "SIRV F1 (E0 + long)", metric), metric))
if (!any(grepl("SIRV F1 (E0 + long)", q_spk$metric, fixed = TRUE))) message("  ! ST8: no SIRV F1 metric found to relabel")
add("ST8", "Bambu and oarfish quantification of the same runs (LongBench and Dong).",
    paste("Same alignments, annotation and libraries. Quantifier comparisons use undivided counts for both quantifiers, so only the quantifier differs (Bambu has no RTA division);",
          "the count_mode column of (a) marks them, and the oarfish counts divided by RTA overdispersion (the primary analysis, Fig. 1e and the text) are added for sequin DE.",
          "(a, b) LongBench spike-in metrics by set and per spike-in. In (a), r2 is the dose-response R2 (log2 CPM on log2 input) in ERCC rows and the R2 of observed against expected log2FC in sequin DE rows.",
          "In (b), SIRV F1 (E0 + long) is the F1 of the two SIRV sets together (84 isoforms; the text and Fig. 1b report the sets separately), and ERCC reads per library and the limit of detection (LOD70, amol/ul) are on the linear scale (Extended Data Fig. 4b plots them as log10).",
          "In (a), recall in sequin DE rows is the sensitivity, and accuracy = (true positives + true negatives) / tested sequins (sequin DE rows only: identification has no true negatives).",
          "(c) LongBench and (d) Dong human transcripts, mean reads over libraries: Spearman correlation among transcripts with a mean of >=1 read with either quantifier (n_any, spearman_any), and, among those with >=10 (n_ge10),",
          "the share within two-fold and the median |log2((oarfish + 1) / (Bambu + 1))| (pct_within_2fold_ge10, median_abs_lr_ge10).",
          "Extended Data Fig. 4; Fig. 2c-g. Jaccard indices, overdispersion shares, disagreement by overdispersion and Dong sequins per isoform are in the benchmarking repository."),
    `a | LongBench metrics by set` = rd(lq, "q_metrics_by_set.csv"),
    `b | LongBench spike-ins` = q_spk,
    `c | LongBench human agreement` = rd(lq, "human_agreement.csv"),
    `d | Dong human agreement` = rd(dq, "human_agreement.csv"),
    `e | Jaccard by threshold` = bind_rows(rd(lq, "jaccard_per_sample.csv") %>% mutate(Dataset = "LongBench", .before = 1),
                                           rd(dq, "jaccard_per_sample.csv") %>% mutate(Dataset = "Dong", .before = 1)),
    `f | Overdispersion by category` = bind_rows(rd(lq, "overdispersion_share_by_category.csv") %>% mutate(Dataset = "LongBench", .before = 1),
                                                 rd(dq, "overdispersion_share_by_category.csv") %>% mutate(Dataset = "Dong", .before = 1)),
    `g | Disagreement by overdispersion` = bind_rows(rd(lq, "disagreement_by_overdispersion.csv") %>% mutate(Dataset = "LongBench", .before = 1),
                                                     rd(dq, "disagreement_by_overdispersion.csv") %>% mutate(Dataset = "Dong", .before = 1)),
    `h | Dong sequins per isoform` = rd(dq, "sequin_quantifiers.csv"))

## ---- ST9: read support ---------------------------------------------------------
add("ST9", "Read-level support of the validated transcriptome (LongBench and Dong).",
    paste("Unit: transcripts of the validated assembly with Bambu read classes (transcripts kept on oarfish counts alone have none).",
          "Full-length = >=1 full-length read in any sample; unique = >=1 uniquely assigned read in any sample (not necessarily the same read). Reference transcripts by reference biotype, novel ones by pulposeq biotype.",
          "Support classes as in Fig. 1g and 2a: full-length + unique reads; unique reads (no full-length read); full-length reads (no uniquely assigned read); partial & ambiguous reads (neither).",
          "(c) Dong, by GffCompare class code of the novel isoforms (reference = annotated isoforms), per biotype."),
    `a | LongBench` = rs_lab(rd(lb, "read-support/results/read_support_all_runs_long.csv")),
    `b | Dong` = rs_lab(rd(dg, "read-support/results/read_support_all_runs_long.csv")),
    `c | Dong by class code` = rs_lab(rd(dg, "figure-panels/results/figure_tables/read_support_by_class_dong.csv")))

## ---- ST10: novel composition ---------------------------------------------------
## the Protocol column of the result files reads "Dong ONT cDNA" for the case
## study: protocol names as in the other tables, with the dataset in its own column
with_dataset <- function(d, dataset) if (any(tolower(names(d)) == "dataset")) d else tag(d, dataset)
st10 <- bind_rows(with_dataset(rd(lb, "figure-panels/results/figure_tables/supp_novel_biotype_classcode_longbench.csv"), "LongBench"),
                  with_dataset(rd(dg, "figure-panels/results/figure_tables/supp_novel_biotype_classcode_dong.csv"), "Dong")) %>%
  mutate(across(any_of(c("Protocol", "protocol")), ~ sub("^(LongBench|Dong) ", "", .x)))
add("ST10", "Novel transcripts by pulposeq biotype and GffCompare class code (LongBench and Dong).",
    paste("Unit: validated novel transcript models. Class codes grouped as junction match (j), contains reference (k), retained intron (m, n), sense overlap (sense intronic i, o, y),",
          "antisense (antisense intronic i, x) and intergenic (u). Percentages of all novel transcripts of the run and within each biotype."),
    `Novel transcripts` = st10)

## ---- ST11, ST12: unique intron chains (07_supp_tables_ujc.R) ------------------
uj <- file.path(lb, "unique_intron_chains/results/supp_tables")
## the length-stratified sharing of Fig. 1j (05_ujc_length.R); bins of fewer
## than ten transcripts are omitted, as in the figure
uic_len <- rd(lb, "unique_intron_chains/results/length/corroboration_by_length.csv") %>%
  filter(n >= 10) %>%
  transmute(`Source protocol` = from, `Also found by` = found_in, `Spliced length in the source protocol` = bin,
            `Multi-exon novel transcripts` = n, `% with UIC also found` = pct_found)
add("ST11", "Sharing of novel intron chains between protocols and transcript-end differences within shared chains (LongBench).",
    paste("(a) Novel models per protocol and (b) their composition, before cross-protocol matching. (c) Unit: unique intron chain (UIC) entities, one row per distinct chain;",
          "mono-exon models have no UIC and are never matched; transcript ends are not part of the match. (d) Unit: protocol pairs (an entity of three protocols contributes three pairs);",
          "|d5| and |d3| = distances between the two protocols' 5' and 3' ends.",
          "(e) Unit: multi-exon novel transcripts of the source protocol; percentage whose UIC was also found by each other protocol, by spliced length in the source protocol (bins of fewer than ten transcripts omitted; Fig. 1j)."),
    `a | Novel models per protocol` = uic_tab(rt(uj, "ST2_discovery_per_protocol.tsv")),
    `b | Composition` = uic_tab(rt(uj, "ST3_composition_of_novel_calls.tsv")),
    `c | UIC sharing` = uic_tab(rt(uj, "ST4_ujc_sharing.tsv")),
    `d | Transcript-end differences` = uic_tab(rt(uj, "ST6_transcript_end_differences.tsv")),
    `e | UIC sharing by length` = uic_len)
add("ST12", "Agreement between protocols on the classification of shared intron chains (LongBench).",
    paste("UICs reconstructed by >=2 protocols. Unit: UIC entities; concordant when every protocol assigned the same class code, biotype or reference host gene; Wilson 95% intervals.",
          "Every disagreement, itemized, is in the benchmarking repository."),
    `a | Concordance` = uic_tab(rt(uj, "ST5_classification_concordance.tsv")),
    `b | Disagreements` = uic_tab(rt(uj, "ST7_itemised_disagreements.tsv")))

## ---- ST13-ST15: case-study DE and cis pairs -----------------------------------
add("ST13", "Library sizes and TMM normalization factors of the differential expression analysis (Dong).",
    paste("Unit: samples, per level. lib.size = reads in the filtered count table (oarfish counts divided by RTA overdispersion); norm.factors = TMM (edgeR); eff.lib.size = lib.size x norm.factors;",
          "rel.scale = mean(eff.lib.size) / eff.lib.size."),
    `Normalization` = rd(dg, "isoformic/results/tables/00_normalisation_factors.csv"))
de <- bind_rows(rd(dg, "diff-exp/results/gene/de_TKI_resistant_vs_TKI_sensitive.csv") %>% mutate(level = "gene", .before = 1),
                rd(dg, "diff-exp/results/transcript/de_TKI_resistant_vs_TKI_sensitive.csv") %>% mutate(level = "transcript", .before = 1))
add("ST14", "Differential expression at gene and transcript level: H1975 vs HCC827 (Dong).",
    paste("Unit: genes or transcripts (annotated and novel) passing the filter (>=10 estimated reads before division in >=3 samples; rRNA, Mt_rRNA, processed pseudogenes and chrY excluded).",
          "Counts: oarfish divided by RTA overdispersion; TMM ({ST13}); edgeR glmQLFit (robust) -> glmTreat (H0: |log2FC| <= 0.5); call at FDR < 0.05. Positive log2FC = higher in H1975."),
    `DE results` = de)
add("ST15", "Novel transcripts and every tested transcript within 1 Mb: cis co-expression (Dong).",
    paste("Unit: pairs of a novel transcript and a neighbouring tested transcript (annotated or novel) on the same chromosome within 1 Mb (distance_bp = 0: overlapping).",
          "rho = Spearman correlation across the 6 samples of logCPM; DE calls and padj from {ST14} (transcript level); same_gene / same_strand compare the two transcripts.",
          "Pairs used in the text: both DE, reference partner of another gene (same_gene = FALSE), |rho| > 0.5 (478 novel transcripts, 233 lncRNA candidates). Figure pairs: Fig. 2q,r (15, any biotype) and Extended Data Fig. 7 (10 novel lncRNAs, partners of other genes)."),
    `All pairs within 1 Mb` = rd(dg, "correlation-cis-1Mb/results/cis_pairs_all.csv"))

## renumber after the drop: sheet names follow the table numbers, and {STn} in a
## description becomes the new number of STn ("the benchmarking repository" when
## STn was dropped, as dropped tables are kept there)
new_no <- setNames(seq_along(ST), names(ST))
xref <- function(s) {
  for (k in unique(regmatches(s, gregexpr("\\{ST[0-9]+\\}", s))[[1]])) {
    id <- gsub("[{}]", "", k)
    s <- gsub(k, if (id %in% names(new_no)) paste("Supplementary Table", new_no[[id]]) else "the benchmarking repository", s, fixed = TRUE)
  }
  gsub(" {2,}", " ", s)
}
ST <- lapply(ST, function(x) { x$desc <- xref(x$desc); x })
names(ST) <- paste0("ST", seq_along(ST))

## ------------------------------------------------------------------- write ---
wb <- createWorkbook()
st_title <- createStyle(fontSize = 12, textDecoration = "bold")
st_desc  <- createStyle(wrapText = TRUE, valign = "top")
st_blk   <- createStyle(textDecoration = "bold", fontColour = "#1F4E79")
st_hdr   <- createStyle(textDecoration = "bold", border = "bottom", fgFill = "#F2F2F2")
addWorksheet(wb, "Descriptions")
idx <- tibble(Sheet = names(ST), Table = paste("Supplementary Table", seq_along(ST)),
              Title = vapply(ST, `[[`, "", "title"), Description = vapply(ST, `[[`, "", "desc"))
writeData(wb, "Descriptions", "pulposeq - Supplementary Tables", startRow = 1); addStyle(wb, "Descriptions", st_title, 1, 1)
writeData(wb, "Descriptions", "One row per table; click a sheet name to open it. Built by 04_supplementary_tables.R from the oarfish-integrated runs. Per-feature results (SIRV and ERCC per isoform, per-sequin DE, quantifier details, itemized intron-chain disagreements, DE results and cis pairs) are in the benchmarking repository.", startRow = 2)
writeData(wb, "Descriptions", idx, startRow = 4, headerStyle = st_hdr)
for (i in seq_len(nrow(idx))) writeFormula(wb, "Descriptions", startRow = 4 + i, startCol = 1,
                                           x = makeHyperlinkString(sheet = idx$Sheet[i], row = 1, col = 1, text = idx$Sheet[i]))
setColWidths(wb, "Descriptions", 1:4, c(8, 22, 60, 120)); addStyle(wb, "Descriptions", st_desc, 5:(4 + nrow(idx)), 3:4, gridExpand = TRUE)
for (i in seq_along(ST)) {
  id <- names(ST)[i]; s <- ST[[i]]
  addWorksheet(wb, id)
  writeData(wb, id, sprintf("Supplementary Table %d | %s", i, s$title), startRow = 1); addStyle(wb, id, st_title, 1, 1)
  writeData(wb, id, s$desc, startRow = 2); addStyle(wb, id, st_desc, 2, 1)
  mergeCells(wb, id, cols = 1:12, rows = 2); setRowHeights(wb, id, 2, 15 * ceiling(nchar(s$desc) / 160))
  r <- 4; w <- integer(0)
  for (b in names(s$blocks)) {
    d <- as.data.frame(s$blocks[[b]])
    if (length(s$blocks) > 1) { writeData(wb, id, b, startRow = r); addStyle(wb, id, st_blk, r, 1); r <- r + 1 }
    writeData(wb, id, d, startRow = r, headerStyle = st_hdr, keepNA = FALSE)
    r <- r + nrow(d) + 3
    ## widths from the table itself (header and values), not from the title and
    ## description in A1:A2
    for (j in seq_along(d)) {
      v <- d[[j]]
      v <- if (inherits(v, "POSIXt")) format(v, "%Y-%m-%d %H:%M:%S") else if (is.numeric(v)) format(signif(v, 7)) else as.character(v)
      w[j] <- max(c(w[j], nchar(names(d)[j]), nchar(v[!is.na(v)])), na.rm = TRUE)
    }
  }
  setColWidths(wb, id, seq_along(w), pmin(pmax(w + 2, 8), 45))
  message(sprintf("  %s: %s (%s rows)", id, s$title, paste(vapply(s$blocks, nrow, 0), collapse = " + ")))
}
if (file.exists(out)) message("  overwriting ", basename(out))

## openxlsx 4.2.9 declares xl/drawings/drawingN.xml in [Content_Types].xml
## and in each sheet's rels without writing the files (strict readers such as
## openpyxl refuse the workbook). Drop every reference to a part not in the file.
repair_xlsx <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)   # zip::zip works inside root
  tmp <- tempfile("xlsx_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = tmp)
  parts <- list.files(tmp, recursive = TRUE, all.files = TRUE)
  rd_txt <- function(f) readChar(f, file.size(f), useBytes = TRUE)
  ct_f <- file.path(tmp, "[Content_Types].xml"); ct <- rd_txt(ct_f)
  ov   <- regmatches(ct, gregexpr('<Override PartName="/[^"]*"[^>]*/>', ct))[[1]]
  miss <- ov[!sub('^<Override PartName="/([^"]*)".*$', "\\1", ov) %in% parts]
  for (m in miss) ct <- sub(m, "", ct, fixed = TRUE)
  writeChar(ct, ct_f, eos = NULL, useBytes = TRUE)
  n_rel <- 0L
  for (rf in list.files(file.path(tmp, "xl/worksheets/_rels"), full.names = TRUE)) {
    r    <- rd_txt(rf)
    rel  <- regmatches(r, gregexpr("<Relationship [^>]*/>", r))[[1]]
    tgt  <- sub('^.*Target="\\.\\./([^"]*)".*$', "xl/\\1", rel)
    drop <- rel[grepl('Target="\\.\\./', rel) & !tgt %in% parts]
    for (d in drop) r <- sub(d, "", r, fixed = TRUE)
    n_rel <- n_rel + length(drop)
    writeChar(r, rf, eos = NULL, useBytes = TRUE)
  }
  files <- c("[Content_Types].xml", setdiff(parts, "[Content_Types].xml"))
  unlink(path)
  zip::zip(path, files, root = tmp, mode = "mirror")
  message(sprintf("  workbook repaired: %d content-type and %d sheet-relationship entries for missing parts removed",
                  length(miss), n_rel))
  invisible(path)
}
saveWorkbook(wb, out, overwrite = TRUE)
repair_xlsx(out)
message("Done: ", out)
