## ---------------------------------------------------------------------------
## Readable labels for the Bambu validated count matrices.
##
## Adapted to the pulposeq output schema as of main @ c0362cb (2026-09-14).
## Both metadata tables were renamed upstream; the old column names are listed
## against the new ones in `PULPOSEQ_SCHEMA` below, and require_cols() fails
## loudly (naming the offending file) rather than silently producing all-NA
## labels.
##
##   gene matrices       rownames = GENEID  (ENSG*.v or BambuGene*)
##   transcript matrices rownames = TXNAME  (ENST*.v or BambuTx*) + GENEID column
##
## Two metadata tables per dataset, both under the pipeline outdir:
##
##   ref_transcripts/annotated_transcriptome_metadata.csv
##     transcript_id_version -> transcript_name
##     gene_id_version       -> gene_name
##     plus gene_biotype / transcript_biotype
##
##   novel_transcripts/novel_transcripts_validated_metadata.csv
##     the validated novel set: one row per BambuTx, with class_code, its
##     human-readable `classification` ("multi-exon SJ match (j)"), the query
##     gene it belongs to (qry_gene_id / qry_gene_name / qry_gene_biotype),
##     the reference gene it was compared against (compared_gene_id /
##     compared_gene_name / compared_gene_biotype) and a
##     pulposeq_transcript_biotype of novel_lncRNA / novel_protein_coding /
##     novel_non_coding.
##
##     qry_gene_name and qry_gene_biotype are blank exactly for the BambuGene
##     rows (a novel gene has no reference name or biotype); for ENSG rows
##     qry_gene_id is versioned and equals compared_gene_id.
##
## Every dataset is quantified against the same reference release, so the
## versioned ID is the key throughout: no version is ever stripped.
##
## Naming rules
##   annotated feature           -> its symbol (SNHG5-225, MALAT1)
##   symbol shared by >1 feature -> SYMBOL_GENEID, not SYMBOL_1
##   gene with only novel        -> symbol from the novel table, since its ENSG
##     transcripts                  is absent from the annotated one
##   novel transcript            -> never renamed, only tagged with its class
##                                  code: BambuTx1_j
##   novel gene                  -> left exactly as it is: BambuGene9345
##
## NOTE: both metadata files are comma separated with NO row-name column, so
## read them with read.csv() and WITHOUT row.names = 1 (row.names = 1 shifts
## every column and errors on duplicate row names).
## ---------------------------------------------------------------------------

## Current column names, and the pre-2026-09 names they replaced. Kept as data
## so the error message can tell you which rename you are looking at instead of
## just reporting a missing column.
PULPOSEQ_SCHEMA <- list(
  ref = c(transcript_id_version = "ensembl_transcript_id_version",
          transcript_name       = "external_transcript_name",
          gene_id_version       = "ensembl_gene_id_version",
          gene_name             = "external_gene_name",
          gene_biotype          = "gene_biotype",
          transcript_biotype    = "transcript_biotype"),
  nov = c(qry_id                      = "qry_id",
          qry_gene_id                 = "qry_gene_id",
          qry_gene_name               = "gene_name",
          qry_gene_biotype            = NA_character_,
          compared_gene_id            = "ref_gene_id",
          compared_gene_name          = "ref_gene_name",
          compared_gene_biotype       = "ref_gene_biotype",
          class_code                  = "class_code",
          classification              = NA_character_,
          pulposeq_transcript_biotype = "transcript_biotype")
)

require_cols <- function(x, which = c("ref", "nov"), what = which) {
  which <- match.arg(which)
  need  <- names(PULPOSEQ_SCHEMA[[which]])
  miss  <- setdiff(need, colnames(x))
  if (!length(miss)) return(invisible(x))

  was <- PULPOSEQ_SCHEMA[[which]][miss]
  hint <- ifelse(is.na(was) | !nzchar(was), "(new column)",
                 paste0("(was `", was, "`)"))
  stop(sprintf(
    paste0("`%s` is not a current pulposeq %s metadata table.\n",
           "Missing column(s):\n%s\n",
           "Present: %s\n",
           "Expected the schema written by pulposeq main @ c0362cb or later."),
    what, which,
    paste0("  - ", miss, " ", hint, collapse = "\n"),
    paste(colnames(x), collapse = ", ")), call. = FALSE)
}

is_novel <- function(x) grepl("^Bambu", as.character(x))

## "" and the literal string "NA" both mean "no value" in these tables
blank <- function(x) {
  x <- as.character(x)
  is.na(x) | !nzchar(x) | x == "NA"
}

## first non-blank of a and b, element by element
coalesce_chr <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  a[blank(a)] <- b[blank(a)]
  a
}

## named vector: names = versioned Ensembl IDs, values = symbol ---------------
build_name_map <- function(ref, level = c("transcript", "gene")) {
  level <- match.arg(level)
  require_cols(ref, "ref")

  if (level == "transcript") {
    key <- ref$transcript_id_version
    val <- ref$transcript_name
  } else {
    key <- ref$gene_id_version
    val <- ref$gene_name
  }
  key <- as.character(key)
  val <- as.character(val)

  ok  <- !blank(key) & !blank(val)
  key <- key[ok]; val <- val[ok]

  keep <- !duplicated(key)
  stats::setNames(val[keep], key[keep])
}

## Reference genes whose only validated feature is a novel transcript never
## reach the annotated table, but the novel one names the gene they came from.
## Ensembl keys only: novel genes are never renamed.
##
## qry_gene_name is blank for the BambuGene rows, so it contributes nothing
## there; compared_gene_name covers the reference gene a novel transcript of a
## novel gene was compared against.
build_novel_gene_map <- function(nov) {
  require_cols(nov, "nov")

  key <- c(as.character(nov$qry_gene_id), as.character(nov$compared_gene_id))
  val <- c(coalesce_chr(nov$qry_gene_name, nov$compared_gene_name),
           as.character(nov$compared_gene_name))

  ok  <- !blank(key) & !blank(val) & grepl("^ENS", key)
  key <- key[ok]; val <- val[ok]

  keep <- !duplicated(key)
  stats::setNames(val[keep], key[keep])
}

## class code of each novel transcript. Only transcripts carry one: a novel
## gene can hold several novel transcripts with different codes, so BambuGene
## IDs are left exactly as they are.
build_class_map <- function(nov) {
  require_cols(nov, "nov")
  cc   <- as.character(nov$class_code)
  ok   <- !blank(cc)
  key  <- as.character(nov$qry_id)[ok]
  val  <- cc[ok]
  keep <- !duplicated(key)
  stats::setNames(val[keep], key[keep])
}

## human-readable class, e.g. "multi-exon SJ match (j)". Too long for a plot
## label, so it is carried as a column in the result tables rather than being
## pasted onto the feature name.
build_classification_map <- function(nov) {
  require_cols(nov, "nov")
  cl   <- as.character(nov$classification)
  ok   <- !blank(cl)
  key  <- as.character(nov$qry_id)[ok]
  val  <- cl[ok]
  keep <- !duplicated(key)
  stats::setNames(val[keep], key[keep])
}

## exact lookup on the versioned ID, otherwise the ID unchanged
id2name <- function(ids, map) {
  ids  <- as.character(ids)
  out  <- unname(map[ids])
  miss <- is.na(out)
  out[miss] <- ids[miss]
  out
}

## Label a vector of feature IDs (no uniquifying) -----------------------------
## annotated -> symbol; novel -> ID_classcode; anything unknown -> ID
feature_label <- function(ids, ref, nov = NULL, level = c("transcript", "gene")) {
  level <- match.arg(level)
  ids   <- as.character(ids)
  out   <- ids
  novel <- is_novel(ids)

  if (any(!novel))
    out[!novel] <- id2name(ids[!novel], build_name_map(ref, level = level))

  if (!is.null(nov)) {
    ## genes carrying only novel transcripts are missing from the annotated table
    if (level == "gene") {
      still <- !novel & out == ids
      if (any(still))
        out[still] <- id2name(ids[still], build_novel_gene_map(nov))
    }
    ## novel transcripts carry their class code; novel genes stay untouched
    if (level == "transcript" && any(novel)) {
      sfx <- unname(build_class_map(nov)[ids[novel]])
      out[novel] <- ifelse(blank(sfx), ids[novel], paste0(ids[novel], "_", sfx))
    }
  }
  out
}

## gene-level label, used for the GENEID column of the transcript matrices
gene_label <- function(ids, ref, nov = NULL) feature_label(ids, ref, nov, "gene")

## Make labels unique ---------------------------------------------------------
## A symbol carried by more than one feature becomes SYMBOL_GENEID. make.unique
## is only a backstop for the rare case where the gene ID does not separate them
## either (two transcripts of one gene sharing a name).
unique_labels <- function(labels, gene_ids) {
  dup <- labels %in% labels[duplicated(labels)]
  if (any(dup)) {
    gid <- as.character(gene_ids)
    gid[blank(gid)] <- labels[blank(gid)]
    labels[dup] <- paste0(labels[dup], "_", gid[dup])
  }
  make.unique(labels, sep = "_")
}

## Which features are lncRNAs -------------------------------------------------
## named vector: feature ID -> "lncRNA" | "novel lncRNA" | "other".
## Reference genes seen only through a novel transcript are typed from the
## novel table's compared_gene_biotype.
build_biotype_map <- function(ref, nov = NULL, level = c("transcript", "gene")) {
  level <- match.arg(level)
  require_cols(ref, "ref")

  if (level == "transcript") {
    key <- as.character(ref$transcript_id_version)
    bt  <- as.character(ref$transcript_biotype)
  } else {
    key <- as.character(ref$gene_id_version)
    bt  <- as.character(ref$gene_biotype)
  }
  ok   <- !blank(key)
  key  <- key[ok]
  val  <- ifelse(bt[ok] == "lncRNA", "lncRNA", "other")
  keep <- !duplicated(key)
  out  <- stats::setNames(val[keep], key[keep])

  if (is.null(nov)) return(out)
  require_cols(nov, "nov")

  nbt <- as.character(nov$pulposeq_transcript_biotype)
  if (level == "transcript") {
    nkey <- as.character(nov$qry_id)
    nval <- ifelse(nbt == "novel_lncRNA", "novel lncRNA", "other")
  } else {
    ## a novel gene counts as a novel lncRNA if any of its novel transcripts is
    qkey <- as.character(nov$qry_gene_id)
    qok  <- is_novel(qkey)
    lnc  <- tapply((nbt == "novel_lncRNA")[qok], qkey[qok], any)
    nkey <- names(lnc)
    nval <- ifelse(unname(lnc), "novel lncRNA", "other")
    ## reference genes reached only through a novel transcript keep their own
    ## reference biotype: the gene is annotated, only its transcript is novel.
    ## qry_gene_biotype is the biotype of the ENSG rows' own gene, and
    ## compared_gene_biotype that of the gene each novel transcript was
    ## compared against - both are reference biotypes, so both are used.
    rkey <- c(qkey[!qok], as.character(nov$compared_gene_id))
    rbt  <- c(as.character(nov$qry_gene_biotype)[!qok],
              as.character(nov$compared_gene_biotype))
    rval <- ifelse(rbt == "lncRNA", "lncRNA", "other")
    rok  <- !blank(rkey) & grepl("^ENS", rkey)
    nkey <- c(nkey, rkey[rok]); nval <- c(nval, rval[rok])
  }
  nok  <- !blank(nkey)
  nkey <- nkey[nok]; nval <- nval[nok]
  keep <- !duplicated(nkey) & !(nkey %in% names(out))
  c(out, stats::setNames(nval[keep], nkey[keep]))
}

## Main entry point -----------------------------------------------------------
## x        : count matrix / data.frame with feature IDs as rownames
## ref      : annotated_transcriptome_metadata.csv
## nov      : novel_transcripts_validated_metadata.csv (NULL to skip novel handling)
## level    : "transcript" for *_tcounts, "gene" for *_gcounts
## gene_col : GENEID column of the transcript matrices; relabelled in place, so
##            a novel isoform of an annotated gene reports that gene's symbol.
##            Set to NULL to leave it alone.
swap_ids <- function(x,
                     ref,
                     nov      = NULL,
                     level    = c("transcript", "gene"),
                     gene_col = "GENEID",
                     verbose  = TRUE) {

  level <- match.arg(level)

  ids <- rownames(x)
  if (is.null(ids)) stop("`x` has no rownames to translate.")

  ## gene ID of each row, used only to break symbol collisions
  gene_ids <- if (level == "gene") ids
              else if (!is.null(gene_col) && gene_col %in% colnames(x)) x[[gene_col]]
              else ids

  lab <- feature_label(ids, ref, nov, level = level)
  new <- unique_labels(lab, gene_ids)
  rownames(x) <- new

  if (!is.null(gene_col) && nzchar(gene_col) && gene_col %in% colnames(x))
    x[[gene_col]] <- gene_label(x[[gene_col]], ref, nov)

  if (verbose) {
    novel <- is_novel(ids)
    message(sprintf(
      "%-10s: %d/%d annotated named, %d kept as ID, %d novel (%d tagged), %d disambiguated",
      level, sum(!novel & lab != ids), sum(!novel), sum(!novel & lab == ids),
      sum(novel), sum(novel & lab != ids), sum(new != lab)))
  }
  x
}
