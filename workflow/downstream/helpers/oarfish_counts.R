## ===========================================================================
## oarfish_counts.R - oarfish counts with read-to-transcript ambiguity (RTA)
##                    overdispersion, for every downstream script that models counts
##
## Method (edgeR reference manual, catchSalmon / catchOarfish page; Baldoni et al.
## 2024a NAR 52:e13, 2024b NARGAB lqae151, 2025 NAR gkaf1305):
##   - edgeR::catchOarfish() reads each sample's <sample>.quant, .meta_info.json and
##     .infreps.pq (the bootstrap replicates), and estimates one RTA overdispersion
##     per transcript (>= 1; 1 = no variance inflation), moderated towards a prior.
##   - Divided counts (counts / overdispersion) "can be input into standard
##     differential expression pipelines designed for gene-level counts"; the edgeR
##     quasi-likelihood pipeline performs well with them. They are not rounded:
##     they "behave much like negative binomial counts".
##   - Gene level: edgeR::catchOarfishGene() sums each sample's bootstraps per gene
##     and estimates a gene-level overdispersion the same way (edgeR >= 4.99.12,
##     i.e. Bioconductor devel / edgeR 5.0.0).
##
## Both are run on EVERY model oarfish quantified (so the overdispersion prior uses
## the whole index), then restricted to pulposeq's curated set. Gene counts are
## the sum over all of a gene's models, as in pulposeq's own oarfish gene tables.
##
## Needs: edgeR >= 4.99.12, readr, nanoparquet, jsonlite.
## Install edgeR from Bioconductor devel until 5.0.0 is released:
##   BiocManager::install(version = "devel"); BiocManager::install("edgeR")
## ===========================================================================

## pulposeq output files used here, relative to a run's --outdir
oarfish_files <- function(run_dir) list(
  oarfish  = file.path(run_dir, "transcriptome", "oarfish"),
  bambu_tx = file.path(run_dir, "transcriptome", "bambu", "BambuOutput_counts_transcript.txt"),
  cur_tx   = file.path(run_dir, "transcriptome", "curated", "counts",
                       "oarfish_counts_transcript_validated.txt"),
  cur_gene = file.path(run_dir, "transcriptome", "curated", "counts",
                       "oarfish_counts_gene_validated.txt")
)

check_oarfish_deps <- function(level = "transcript") {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR is not installed.", call. = FALSE)
  need <- if (level == "gene") "catchOarfishGene" else "catchOarfish"
  if (!exists(need, envir = asNamespace("edgeR"), inherits = FALSE))
    stop(sprintf(paste0("edgeR %s has no %s(). Install edgeR >= 4.99.12 (Bioconductor devel,",
                        " to be released as 5.0.0):\n  BiocManager::install(version = \"devel\");",
                        " BiocManager::install(\"edgeR\")"),
                 as.character(utils::packageVersion("edgeR")), need), call. = FALSE)
  for (p in c("readr", "nanoparquet", "jsonlite"))
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is needed to read oarfish output: install.packages(\"", p, "\")",
           call. = FALSE)
  invisible(TRUE)
}

## first column of a pulposeq table (TXNAME or GENEID)
first_col <- function(path) {
  if (!file.exists(path)) stop("Missing pulposeq output: ", path, call. = FALSE)
  as.character(readr::read_tsv(path, col_select = 1, show_col_types = FALSE, progress = FALSE)[[1]])
}

## ---------------------------------------------------------------------------
## read_oarfish(run_dir, level)
##   run_dir   a pulposeq --outdir that ran with --run_oarfish true
##   level     "transcript" or "gene"
##   curated   TRUE: keep the curated set only (pulposeq's curated oarfish tables)
##   cache_dir optional folder for an .rds copy, reused while the oarfish files
##             are unchanged (reading the bootstraps is the slow part)
## Returns a list:
##   counts          undivided counts (curated rows x samples)
##   divided         counts / overdispersion (not rounded)
##   overdispersion  named vector, one value per row of counts
##   prior           the overdispersion prior catchOarfish/catchOarfishGene used
##   gene            transcript level only: GENEID of each transcript
## ---------------------------------------------------------------------------
read_oarfish <- function(run_dir, level = c("transcript", "gene"), curated = TRUE,
                         cache_dir = NULL, verbose = FALSE) {
  level <- match.arg(level)
  f <- oarfish_files(run_dir)
  if (!dir.exists(f$oarfish))
    stop("No oarfish output in ", f$oarfish, " - was the run made with --run_oarfish true?",
         call. = FALSE)
  check_oarfish_deps(level)

  ## optional cache, invalidated when any oarfish output is newer
  cache <- NULL
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    ## short name (Windows path-length limit): the run's last two folders,
    ## e.g. pulposeq_oarfish_ont_cdna
    run_tag <- gsub("[^A-Za-z0-9]+", "_",
                    paste(basename(dirname(run_dir)), basename(run_dir), sep = "_"))
    cache <- file.path(cache_dir, sprintf("oarfish_%s_%s_%s.rds", level,
                                          if (curated) "curated" else "all", run_tag))
    src <- list.files(f$oarfish, "\\.(quant|infreps\\.pq|meta_info\\.json)$", full.names = TRUE)
    if (file.exists(cache) && all(file.mtime(src) <= file.mtime(cache)))
      return(readRDS(cache))
  }

  if (level == "transcript") {
    x <- edgeR::catchOarfish(parent.dir = f$oarfish, DGEList = FALSE, divide = FALSE,
                             verbose = verbose)
    od <- setNames(x$annotation$Overdispersion, rownames(x$counts))
    tx2gene <- utils::read.delim(f$bambu_tx, check.names = FALSE, stringsAsFactors = FALSE)[, 1:2]
    gene <- setNames(tx2gene[[2]], tx2gene[[1]])
  } else {
    ## tx2gene from Bambu's own table (every model; IDs exactly as pulposeq writes
    ## them, so remove.version = FALSE keeps them matching the curated tables)
    tx2gene <- utils::read.delim(f$bambu_tx, check.names = FALSE, stringsAsFactors = FALSE)[, 1:2]
    x <- edgeR::catchOarfishGene(parent.dir = f$oarfish, tx2gene = tx2gene,
                                 remove.version = FALSE, DGEList = FALSE, divide = FALSE,
                                 verbose = verbose)
    od <- setNames(x$annotation$Overdispersion, rownames(x$counts))
  }
  counts <- x$counts
  colnames(counts) <- basename(colnames(counts))

  if (curated) {
    keep <- first_col(if (level == "transcript") f$cur_tx else f$cur_gene)
    miss <- setdiff(keep, rownames(counts))
    if (length(miss))
      stop(sprintf("%d curated %ss are not in the oarfish output, e.g. %s", length(miss), level,
                   paste(head(miss, 3), collapse = ", ")), call. = FALSE)
    counts <- counts[keep, , drop = FALSE]
    od     <- od[keep]
  }

  out <- list(counts = counts, divided = counts / od, overdispersion = od,
              prior = x$overdispersion.prior, level = level,
              gene = if (level == "transcript") unname(gene[rownames(counts)]) else NULL)
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

## pulposeq names Bambu's sample columns "<sample>.bam"; oarfish uses "<sample>"
strip_bam <- function(x) sub("[.]bam$", "", x)
