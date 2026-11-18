#!/usr/bin/env Rscript
# =============================================================================
# run_all.R - runs every downstream script in order:
#   downstream/longbench/01-09   Fig. 1, Extended Data Fig. 4 a-h (LongBench)
#   downstream/dong/01-11        Fig. 2, Extended Data Figs. 4 i-l, 5-7 (Dong)
#   01-04 in this folder         both datasets: Extended Data Figs. 1-3,
#                                Supplementary Tables
#
# Each script runs in its OWN fresh R process (Rscript), with the script's folder
# as working directory, so nothing one script defines can mask a function in the
# next one. Everything a script prints is shown here as it runs; nothing is
# saved to a log file. A summary table is printed at the end.
#
# Use: Rscript run_all.R from this folder, after editing the settings in section 0.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Settings -- edit here
# -----------------------------------------------------------------------------
BASE <- "C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow"

RUN_FOLDERS <- c("longbench", "dong", "main") # which folders, in this order, to run
SKIP <- character(0)   # scripts not to run, written <folder>/<script>
ONLY <- character(0)
STOP_ON_ERROR <- FALSE # FALSE: a failed script only skips the scripts that need it
# TRUE : stop at the first failure

## where each folder's scripts are, relative to BASE
FOLDERS <- c(longbench = "downstream/longbench", dong = "downstream/dong", main = ".")

# Scripts in run order, each with the scripts it needs (written <folder>/<script>).
# "<script>@<arg>" runs the script with that argument: "@bambu" repeats a
# spike-in script on Bambu's counts of the same runs (quantifier comparison);
# "@pass2" reruns dong/08 once dong/10 has saved its panel sizes (the argument is ignored).
PIPELINE <- list(
  longbench = list(
    "01_sequin_edger.R" = character(0),
    "02_sirv_ercc.R" = character(0),
    "01_sequin_edger.R@bambu" = character(0),
    "02_sirv_ercc.R@bambu" = character(0),
    "03_quantifier_plots.R" = c(
      "longbench/01_sequin_edger.R", "longbench/02_sirv_ercc.R",
      "longbench/01_sequin_edger.R@bambu", "longbench/02_sirv_ercc.R@bambu"
    ),
    "04_ujc_report.R" = character(0),
    "05_ujc_length.R" = "longbench/04_ujc_report.R",
    "06_read_support.R" = character(0),
    "07_supp_tables_ujc.R" = "longbench/04_ujc_report.R",
    "08_fig1_panels.R" = c(
      "longbench/01_sequin_edger.R", "longbench/02_sirv_ercc.R",
      "longbench/03_quantifier_plots.R", "longbench/04_ujc_report.R",
      "longbench/05_ujc_length.R", "longbench/06_read_support.R"
    ),
    "09_assemble_figures.R" = c("longbench/08_fig1_panels.R", "longbench/03_quantifier_plots.R")
  ),
  dong = list(
    "01_de_edger.R" = character(0),
    "02_sequin_edger.R" = character(0),
    "02_sequin_edger.R@bambu" = character(0),
    ## oarfish vs Bambu on the spike-ins of both datasets
    "03_quantifier_comparison.R" = c(
      "dong/02_sequin_edger.R", "dong/02_sequin_edger.R@bambu",
      "longbench/01_sequin_edger.R", "longbench/02_sirv_ercc.R",
      "longbench/01_sequin_edger.R@bambu", "longbench/02_sirv_ercc.R@bambu"
    ),
    "04_cis_circos.R" = "dong/01_de_edger.R",
    "05_cis_lncRNA_top10.R" = c("dong/01_de_edger.R", "dong/04_cis_circos.R"),
    "06_isoformic.R" = c("dong/01_de_edger.R", "dong/04_cis_circos.R"),
    "07_read_support.R" = character(0),
    ## Fig. 2 rows 1-2 share plot boxes across 08 and 10 (helpers/fig2_rows.R):
    ## 08, 10, then 08 again with both sizes known; 11 assembles Fig. 2 and the Dong EDs
    "08_fig2_panels.R" = c(
      "dong/01_de_edger.R", "dong/02_sequin_edger.R",
      "dong/02_sequin_edger.R@bambu", "dong/07_read_support.R"
    ),
    "09_circos_coverage_biotypes.R" = character(0),
    "10_quantifier_panels.R" = c("dong/03_quantifier_comparison.R", "dong/08_fig2_panels.R"),
    "08_fig2_panels.R@pass2" = c("dong/10_quantifier_panels.R", "dong/06_isoformic.R"),
    "11_assemble_figures.R" = c(
      "dong/04_cis_circos.R", "dong/05_cis_lncRNA_top10.R",
      "dong/06_isoformic.R", "dong/10_quantifier_panels.R",
      "dong/08_fig2_panels.R@pass2", "longbench/09_assemble_figures.R"
    )
  ),
  main = list(
    "01_computational_resources.R" = character(0),
    "02_transcript_structure.R" = character(0),
    "03_sequencing_metrics.R" = character(0),  # Extended Data Fig. 2, tables for ST1
    "04_supplementary_tables.R" = c(
      "main/01_computational_resources.R", "main/03_sequencing_metrics.R",
      "longbench/01_sequin_edger.R", "longbench/02_sirv_ercc.R",
      "longbench/03_quantifier_plots.R",
      "longbench/06_read_support.R", "longbench/07_supp_tables_ujc.R",
      "longbench/08_fig1_panels.R",
      "dong/02_sequin_edger.R", "dong/02_sequin_edger.R@bambu",
      "dong/03_quantifier_comparison.R",
      "dong/06_isoformic.R", "dong/07_read_support.R",
      "dong/08_fig2_panels.R@pass2", "dong/10_quantifier_panels.R"
    )
  )
)

# -----------------------------------------------------------------------------
# 1. Set-up
# -----------------------------------------------------------------------------
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("Rscript not found at ", rscript)
if (!dir.exists(BASE)) stop("BASE folder not found: ", BASE)

unknown <- setdiff(RUN_FOLDERS, names(PIPELINE))
if (length(unknown)) stop("Unknown folder in RUN_FOLDERS: ", paste(unknown, collapse = ", "))
steps <- do.call(rbind, lapply(RUN_FOLDERS, function(fd) {
  data.frame(
    folder = fd, script = names(PIPELINE[[fd]]), id = paste0(fd, "/", names(PIPELINE[[fd]])),
    stringsAsFactors = FALSE
  )
}))
needs_of <- function(id) {
  i <- match(id, steps$id)
  PIPELINE[[steps$folder[i]]][[steps$script[i]]]
}
all_ids <- unlist(lapply(names(PIPELINE), function(fd) paste0(fd, "/", names(PIPELINE[[fd]]))))
bad <- setdiff(c(SKIP, ONLY, unlist(PIPELINE)), all_ids)
if (length(bad)) {
  stop(
    "Not in the pipeline: ", paste(bad, collapse = ", "),
    " (write them as <folder>/<script>.R, folder = longbench, dong or main)"
  )
}

fmt_time <- function(s) {
  sprintf(
    "%02d:%02d:%02d", as.integer(s) %/% 3600,
    (as.integer(s) %% 3600) %/% 60, as.integer(s) %% 60
  )
}
status <- setNames(rep(NA_character_, nrow(steps)), steps$id)
elapsed <- setNames(rep(NA_real_, nrow(steps)), steps$id)
rule <- strrep("=", 78)

cat(sprintf(
  "run_all.R  |  R %s  |  %d scripts  |  no log files are written\n",
  getRversion(), nrow(steps)
))

# -----------------------------------------------------------------------------
# 2. Run
# -----------------------------------------------------------------------------
t_all <- Sys.time()
for (i in seq_len(nrow(steps))) {
  id <- steps$id[i]
  part <- strsplit(steps$script[i], "@", fixed = TRUE)[[1]] # "<script>@<arg>"
  scr <- normalizePath(file.path(BASE, FOLDERS[[steps$folder[i]]], part[1]),
    winslash = "/", mustWork = FALSE
  )
  args <- part[-1]

  ## not selected / skipped on purpose / missing / a needed script did not succeed
  if (length(ONLY) && !id %in% ONLY) {
    status[id] <- "not selected"
    next
  }
  cat(sprintf(
    "\n%s\n[%2d/%d] %s   (started %s)\n%s\n", rule, i, nrow(steps), id,
    format(Sys.time(), "%H:%M:%S"), rule
  ))
  if (id %in% SKIP) {
    status[id] <- "SKIPPED (setting)"
    cat("SKIPPED (setting)\n")
    next
  }
  if (!file.exists(scr)) {
    status[id] <- "FAILED (file missing)"
    cat("FAILED: file missing:", scr, "\n")
  } else {
    failed_need <- intersect(needs_of(id), names(status)[grepl("^(FAILED|SKIPPED \\(needs)", status)])
    if (length(failed_need)) {
      status[id] <- "SKIPPED (needs a failed script)"
      cat("SKIPPED: needs", paste(failed_need, collapse = ", "), "\n")
      next
    }
    t0 <- Sys.time()
    owd <- setwd(dirname(scr))
    ## output goes straight to this console (stdout and stderr), nothing is saved
    code <- tryCatch(
      system2(rscript, c("--no-save", "--no-restore", shQuote(scr), args), stdout = "", stderr = ""),
      error = function(e) {
        message(conditionMessage(e))
        1L
      },
      finally = setwd(owd)
    )
    elapsed[id] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    status[id] <- if (identical(as.integer(code), 0L)) "OK" else "FAILED"
    cat(sprintf("--> %s %s  %s\n", id, status[id], fmt_time(elapsed[id])))
  }
  if (STOP_ON_ERROR && grepl("^FAILED", status[id])) {
    cat("\nSTOP_ON_ERROR = TRUE: stopping here.\n")
    break
  }
}

# -----------------------------------------------------------------------------
# 3. Summary (printed only)
# -----------------------------------------------------------------------------
res <- data.frame(
  script = names(status), status = ifelse(is.na(status), "not run", status),
  time = ifelse(is.na(elapsed), "", fmt_time(elapsed)), stringsAsFactors = FALSE
)
cat(sprintf("\n%s\nSummary\n%s\n", rule, rule))
cat(sprintf("  %-45s %-32s %s\n", res$script, res$status, res$time), sep = "")
cat(sprintf(
  "\nDone in %s.  OK: %d  |  failed: %d  |  skipped: %d\n",
  fmt_time(as.numeric(difftime(Sys.time(), t_all, units = "secs"))),
  sum(res$status == "OK"), sum(grepl("^FAILED", res$status)),
  sum(grepl("^SKIPPED", res$status))
))
if (any(grepl("^FAILED", res$status))) {
  cat("Failed: ", paste(res$script[grepl("^FAILED", res$status)], collapse = ", "),
    "  -> scroll up to its block for the error\n",
    sep = ""
  )
}
