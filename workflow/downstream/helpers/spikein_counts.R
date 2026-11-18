## ===========================================================================
## spikein_counts.R - spike-in rows of one pulposeq run, from either quantifier,
##                    so the spike-in scripts can compare oarfish with Bambu on
##                    the SAME alignments and annotation
## Sourced after oarfish_counts.R (needs read_oarfish()).
##
##   quantifier = "oarfish" : read_oarfish() -> undivided counts, counts divided
##                            by the read-to-transcript ambiguity (RTA)
##                            overdispersion, and that overdispersion
##   quantifier = "bambu"   : Bambu's `counts` table; no bootstraps, so
##                            divided = undivided and overdispersion = NA
## Both are EM-estimated reads per transcript (or gene), not molecules.
##
## curated = TRUE : the curated set pulposeq delivers
## curated = FALSE: every model the quantifier saw (oarfish: the whole index;
##                  Bambu: its own raw table) - to show what curation removes
## ===========================================================================

quantifier_arg <- function(default = "oarfish") {
  a <- commandArgs(trailingOnly = TRUE)
  q <- if (length(a)) a[1] else default
  match.arg(q, c("oarfish", "bambu"))
}

## results folder for a quantifier: oarfish = the primary folder the figure
## scripts read; bambu = "<folder>_bambu" (comparison only)
quant_dir <- function(dir, quantifier) if (quantifier == "oarfish") dir else paste0(dir, "_bambu")

bambu_table <- function(run_dir, level, curated) {
  if (curated)
    file.path(run_dir, "transcriptome", "curated", "counts",
              sprintf("BambuOutput_counts_%s_validated.txt", level))
  else
    file.path(run_dir, "transcriptome", "bambu", sprintf("BambuOutput_counts_%s.txt", level))
}

## list(reads, divided, od) for the rows matching `rows` (regex on the IDs)
read_counts_q <- function(run_dir, level = c("transcript", "gene"),
                          quantifier = c("oarfish", "bambu"), rows,
                          curated = TRUE, cache_dir = NULL) {
  level <- match.arg(level); quantifier <- match.arg(quantifier)
  if (quantifier == "oarfish") {
    x <- read_oarfish(run_dir, level, curated = curated, cache_dir = cache_dir)
    s <- grepl(rows, rownames(x$counts))
    return(list(reads   = x$counts[s, , drop = FALSE],
                divided = x$divided[s, , drop = FALSE],
                od      = x$overdispersion[s]))
  }
  f <- bambu_table(run_dir, level, curated)
  if (!file.exists(f)) stop("Missing Bambu table: ", f, call. = FALSE)
  x <- utils::read.table(f, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  x$GENEID <- NULL                                    # transcript table: TXNAME, GENEID, samples
  m <- as.matrix(x[grepl(rows, rownames(x)), , drop = FALSE])
  colnames(m) <- sub("[.]bam$", "", colnames(m))
  list(reads = m, divided = m, od = setNames(rep(NA_real_, nrow(m)), rownames(m)))
}
