## ===========================================================================
## 06_read_support.R - read-level support per transcript biotype, LongBench   (Fig. 1g)
##
##   Universe     every transcript of the curated assembly
##                (BambuOutput_counts_transcript_validated.txt). In these
##                oarfish runs the curated set is chosen on oarfish counts
##                (>= min_total_counts summed across samples), so transcripts
##                kept on oarfish counts alone have no Bambu reads and fall in
##                "neither". Full-length and unique support are Bambu's.
##   Full-length  >= 1 full-length read in any sample
##                (BambuOutput_fullLengthCounts_transcript_validated.txt)
##   Unique       >= 1 uniquely assigned read in any sample
##                (BambuOutput_uniqueCounts_transcript_validated.txt)
##                A transcript missing from a table has no such read.
##   Classes      full-length + unique | unique only | full-length only |
##                neither (supported only by partial, multi-assigned reads)
##   Biotypes     annotated: reference metadata, using the column that
##                reproduces the report's Bambu_lncRNA_PC_summary.csv counts
##                (protein-coding, lncRNA, other biotypes);
##                novel: pulposeq_transcript_biotype of the validated novel
##                metadata (predicted protein-coding, novel lncRNA candidate,
##                predicted non-coding)
##   Percentages  full precision in the tables, one decimal in pct_1dp (the
##                transcriptome report rounds to whole numbers)
##
## Checks per run (message; stop only if a table cannot be read):
##   transcripts / full-length / uniquely mapped totals (annotated + novel)
##   against the report's Bambu_assembly_summary.csv, when that file exists.
##
## Outputs: downstream/longbench/read-support/results/
##   read_support_<run>.csv            one row per biotype x support class
##   read_support_all_runs_long.csv    the runs stacked (input of the panel)
##   read_support_all_runs_wide.csv    one row per run x biotype (n and %)
##   read_support_tables.xlsx          one sheet per run + "All runs"
##                                     (only if openxlsx is installed)
## Needs: tidyverse; optional openxlsx
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

## ---------------------------------------------------------------- paths ---
source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")   # base, rs_dir
out_dir <- rs_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## one row per run: folder under `base` and the label used in the tables.
## Every run folder (run_dir in 00_paths.R) is a pulposeq output directory
## with transcriptome/curated/ and transcriptome_report/.
runs <- tribble(
  ~run,       ~label,
  "ont_cdna", "LongBench ONT cDNA",
  "ont_drs",  "LongBench ONT DRS",
  "pb",       "LongBench PacBio"
)

## relative to each run folder (edit here if a run is laid out differently)
paths <- list(
  counts    = "transcriptome/curated/counts/BambuOutput_counts_transcript_validated.txt",
  full      = "transcriptome/curated/counts/BambuOutput_fullLengthCounts_transcript_validated.txt",
  unique    = "transcriptome/curated/counts/BambuOutput_uniqueCounts_transcript_validated.txt",
  ann_meta  = "transcriptome/curated/metadata/annotated_transcriptome_metadata.csv",
  nov_dir   = "transcriptome/curated/metadata",    # *validated_metadata*.csv inside
  rep_pc    = "transcriptome_report/Bambu_lncRNA_PC_summary.csv",
  rep_asm   = "transcriptome_report/Bambu_assembly_summary.csv"
)

## ------------------------------------------------------------- settings ---
sup_lv <- c("full-length + unique", "unique only", "full-length only", "neither")
bt_lv  <- c("annotated protein-coding", "annotated lncRNA", "other biotypes",
            "predicted protein-coding", "novel lncRNA candidate", "predicted non-coding")
## summary rows added under the six biotypes
sum_lv <- c("all annotated", "all novel", "all transcripts")

## ------------------------------------------------------------ helpers -----
bare_tx <- function(x) sub("\\.[0-9]+$", "", sub("_[A-Za-z=]+$", "", x))

pick_col <- function(df, cands, what, file) {
  hit <- intersect(cands, names(df))
  if (!length(hit)) stop("no ", what, " column in ", file, " (looked for ",
                         paste(cands, collapse = ", "), "); columns: ",
                         paste(names(df), collapse = ", "))
  hit
}

need <- function(f) {
  if (!file.exists(f)) stop("file not found: ", f)
  f
}

## TXNAME of transcripts with >= 1 read in any sample of a Bambu count table
## (columns 1-2 are TXNAME and GENEID, the rest one per sample)
ids_with_reads <- function(f) {
  d <- read_tsv(need(f), show_col_types = FALSE, progress = FALSE)
  d$TXNAME[rowSums(as.matrix(d[, -(1:2)]), na.rm = TRUE) > 0]
}

## counts per biotype x support class, + the three summary rows
tabulate_support <- function(sup) {
  by_bt <- sup %>% dplyr::count(biotype, support, name = "n")
  add   <- bind_rows(
    sup %>% filter(status == "annotated") %>% dplyr::count(support, name = "n") %>% mutate(biotype = sum_lv[1]),
    sup %>% filter(status == "novel")     %>% dplyr::count(support, name = "n") %>% mutate(biotype = sum_lv[2]),
    sup %>%                                  dplyr::count(support, name = "n") %>% mutate(biotype = sum_lv[3]))
  bind_rows(by_bt, add) %>%
    mutate(biotype = factor(biotype, levels = c(bt_lv, sum_lv)),
           support = factor(support, levels = sup_lv)) %>%
    complete(biotype, support, fill = list(n = 0L)) %>%        # every class shown, also 0
    group_by(biotype) %>%
    mutate(total = sum(n), pct = if_else(total > 0, 100 * n / total, NA_real_)) %>%
    ungroup() %>%
    filter(total > 0) %>%                                      # drop biotypes absent from a run
    mutate(pct_1dp = sprintf("%.1f%%", pct)) %>%
    arrange(biotype, support)
}

## ------------------------------------------------------------ one run -----
read_support_run <- function(run, label) {
  rd <- run_dir[[run]]
  f  <- function(k) file.path(rd, paths[[k]])
  message("== ", label, " (", rd, ")")

  cnt    <- read_tsv(need(f("counts")), show_col_types = FALSE, progress = FALSE,
                     col_select = c(TXNAME))
  fl_ids <- ids_with_reads(f("full"))
  un_ids <- ids_with_reads(f("unique"))
  is_ann <- !grepl("^BambuTx", cnt$TXNAME)

  ## annotated biotype: candidate columns tried in order; keep the one that
  ## reproduces the report's annotated protein-coding / lncRNA counts
  ann_meta <- read_csv(need(f("ann_meta")), show_col_types = FALSE, progress = FALSE)
  ann_id   <- pick_col(ann_meta, c("transcript_id_version", "transcript_id", "TXNAME"),
                       "transcript id", f("ann_meta"))[1]
  ann_cols <- pick_col(ann_meta, c("gene_type", "gene_biotype", "transcript_type", "transcript_biotype"),
                       "biotype", f("ann_meta"))
  ann_key  <- match(bare_tx(cnt$TXNAME[is_ann]), bare_tx(ann_meta[[ann_id]]))
  target   <- NULL
  if (file.exists(f("rep_pc"))) {
    rp <- read_csv(f("rep_pc"), show_col_types = FALSE) %>%
      filter(Bambu == "transcripts", Status == "Annotated")
    target <- setNames(rp$Count, rp$Biotype)[c("Protein-coding", "lncRNA")]
  }
  ann_bt <- NULL; ann_src <- ann_cols[1]
  for (cc in ann_cols) {
    b   <- ann_meta[[cc]][ann_key]
    got <- c(sum(b %in% "protein_coding"), sum(b %in% "lncRNA"))
    if (is.null(ann_bt)) ann_bt <- b                           # fallback: first candidate
    if (!is.null(target) && isTRUE(all(got == target))) { ann_bt <- b; ann_src <- cc; break }
  }
  if (is.null(target)) {
    message("  ! no Bambu_lncRNA_PC_summary.csv: annotated biotype from '", ann_src, "' (not checked)")
  } else if (!isTRUE(all(c(sum(ann_bt %in% "protein_coding"), sum(ann_bt %in% "lncRNA")) == target))) {
    message("  ! annotated protein-coding / lncRNA counts differ from the report with every column - ",
            "using '", ann_src, "'; check the biotype column")
  } else {
    message("  annotated biotype from '", ann_src, "' (matches the report)")
  }
  n_unmatched <- sum(is.na(ann_key))
  if (n_unmatched) message("  ! ", n_unmatched, " annotated transcripts not in the reference metadata (-> other biotypes)")

  ## novel biotype from the pulposeq metadata of the validated novel transcripts
  nov_f <- list.files(file.path(rd, paths$nov_dir), pattern = "validated_metadata.*\\.csv$",
                      full.names = TRUE)[1]
  if (is.na(nov_f)) stop("no ", file.path(rd, paths$nov_dir), "/*validated_metadata*.csv")
  nov_meta <- read_csv(nov_f, show_col_types = FALSE, progress = FALSE)
  nov_id   <- pick_col(nov_meta, c("qry_id", "transcript_id", "TXNAME", "id", "transcript_name"),
                       "transcript id", nov_f)[1]
  nov_bc   <- pick_col(nov_meta, c("pulposeq_transcript_biotype", "transcript_biotype", "biotype"),
                       "biotype", nov_f)[1]
  nov_key  <- match(bare_tx(cnt$TXNAME[!is_ann]), bare_tx(nov_meta[[nov_id]]))
  message(sprintf("  novel: %d / %d found in %s", sum(!is.na(nov_key)), sum(!is_ann), basename(nov_f)))
  nov_bt   <- nov_meta[[nov_bc]][nov_key]

  ## per transcript: status, biotype, support class
  sup <- tibble(id = cnt$TXNAME, status = if_else(is_ann, "annotated", "novel"),
                biotype = NA_character_, fl = cnt$TXNAME %in% fl_ids, un = cnt$TXNAME %in% un_ids) %>%
    mutate(support = case_when(fl & un ~ sup_lv[1], un ~ sup_lv[2], fl ~ sup_lv[3], TRUE ~ sup_lv[4]))
  sup$biotype[is_ann]  <- case_when(ann_bt %in% "protein_coding" ~ bt_lv[1],
                                    ann_bt %in% "lncRNA"         ~ bt_lv[2],
                                    TRUE                         ~ bt_lv[3])
  sup$biotype[!is_ann] <- case_when(grepl("protein", nov_bt) ~ bt_lv[4],
                                    grepl("lncRNA",  nov_bt) ~ bt_lv[5],
                                    TRUE                     ~ bt_lv[6])

  ## check against the report's assembly summary (annotated / novel totals)
  if (file.exists(f("rep_asm"))) {
    asm <- read_csv(f("rep_asm"), show_col_types = FALSE) %>%
      filter(Status %in% c("Annotated", "Novel"))
    mine <- tribble(
      ~Bambu,           ~Status,      ~mine,
      "transcripts",    "Annotated",  sum(is_ann),
      "transcripts",    "Novel",      sum(!is_ann),
      "fullLength",     "Annotated",  sum(sup$fl & is_ann),
      "fullLength",     "Novel",      sum(sup$fl & !is_ann),
      "uniquelyMapped", "Annotated",  sum(sup$un & is_ann),
      "uniquelyMapped", "Novel",      sum(sup$un & !is_ann))
    chk <- inner_join(mine, asm, by = c("Bambu", "Status"))
    bad <- filter(chk, mine != Count)
    if (nrow(bad)) {
      message("  ! differs from Bambu_assembly_summary.csv:\n",
              paste(sprintf("    %s %s: %d here vs %d in the report", bad$Bambu, bad$Status, bad$mine, bad$Count),
                    collapse = "\n"))
    } else {
      message(sprintf("  matches the report: %s transcripts, %s full-length, %s uniquely mapped",
                      format(sum(chk$mine[chk$Bambu == "transcripts"]), big.mark = ","),
                      format(sum(chk$mine[chk$Bambu == "fullLength"]), big.mark = ","),
                      format(sum(chk$mine[chk$Bambu == "uniquelyMapped"]), big.mark = ",")))
    }
  } else {
    message("  ! no Bambu_assembly_summary.csv: totals not checked")
  }

  tab <- tabulate_support(sup) %>% mutate(run = label, .before = 1)
  write_csv(tab, file.path(out_dir, paste0("read_support_", run, ".csv")))
  tab
}

## ------------------------------------------------------------- all runs ---
long <- map2_dfr(runs$run, runs$label, read_support_run) %>%
  mutate(run = factor(run, levels = runs$label))
write_csv(long, file.path(out_dir, "read_support_all_runs_long.csv"))

## wide: one row per run x biotype, n and % for each support class, plus the
## share with any full-length / any unique read
wide <- long %>%
  dplyr::select(run, biotype, total, support, n, pct) %>%
  pivot_wider(names_from = support, values_from = c(n, pct), names_glue = "{support} ({.value})") %>%
  mutate(`any full-length (%)` = `full-length + unique (pct)` + `full-length only (pct)`,
         `any unique (%)`      = `full-length + unique (pct)` + `unique only (pct)`) %>%
  rename_with(~ sub(" \\(pct\\)$", " (%)", .x)) %>%
  rename(Run = run, Biotype = biotype, Transcripts = total) %>%
  dplyr::select(Run, Biotype, Transcripts,
                all_of(as.vector(rbind(paste0(sup_lv, " (n)"), paste0(sup_lv, " (%)")))),
                `any full-length (%)`, `any unique (%)`) %>%
  arrange(Run, Biotype)
write_csv(wide, file.path(out_dir, "read_support_all_runs_wide.csv"))

## ------------------------------------------------------ Excel (optional) ---
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb  <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#EDEDED",
                               border = "bottom", wrapText = TRUE, valign = "top")
  pct <- openxlsx::createStyle(numFmt = "0.0")
  int <- openxlsx::createStyle(numFmt = "#,##0")
  sumrow <- openxlsx::createStyle(textDecoration = "bold", border = "top")
  write_sheet <- function(name, d) {
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, d, headerStyle = hdr)
    pc <- grep("\\(%\\)$", names(d)); nc <- grep("\\(n\\)$|^Transcripts$", names(d))
    rows <- 2:(nrow(d) + 1)
    if (length(pc)) openxlsx::addStyle(wb, name, pct, rows = rows, cols = pc, gridExpand = TRUE, stack = TRUE)
    if (length(nc)) openxlsx::addStyle(wb, name, int, rows = rows, cols = nc, gridExpand = TRUE, stack = TRUE)
    s <- which(as.character(d$Biotype) %in% sum_lv) + 1
    if (length(s)) openxlsx::addStyle(wb, name, sumrow, rows = s, cols = seq_along(d), gridExpand = TRUE, stack = TRUE)
    openxlsx::setColWidths(wb, name, cols = seq_along(d), widths = c(20, 24, rep(12, ncol(d) - 2)))
    openxlsx::freezePane(wb, name, firstRow = TRUE)
  }
  for (i in seq_len(nrow(runs))) {
    write_sheet(substr(runs$label[i], 1, 31), filter(wide, Run == runs$label[i]))
  }
  write_sheet("All runs", wide)
  openxlsx::saveWorkbook(wb, file.path(out_dir, "read_support_tables.xlsx"), overwrite = TRUE)
  message("Excel: read_support_tables.xlsx")
} else {
  message("openxlsx not installed: CSV only (install.packages(\"openxlsx\") for the .xlsx)")
}

## console summary: any full-length / any unique per run, all transcripts
wide %>%
  filter(Biotype %in% c("all annotated", "all novel")) %>%
  transmute(Run, Biotype, Transcripts,
            `full-length + unique` = sprintf("%.1f%%", `full-length + unique (%)`),
            `neither`              = sprintf("%.1f%%", `neither (%)`)) %>%
  print(n = Inf)
message("Done: ", out_dir)
