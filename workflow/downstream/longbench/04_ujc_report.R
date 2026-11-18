# =============================================================================
# 04_ujc_report.R - cross-protocol comparison of novel transcripts by unique
#                   intron chain (tables for Fig. 1i, k, 05 and 07)
#
# WHAT THIS SCRIPT DOES AND DOES NOT CLAIM
#   It groups novel transcript models by Unique Junction Chain (UJC; Nanni et
#   al. 2024) and DESCRIBES how the protocols agree. It does NOT build a
#   consolidated ID space for quantification, does not pick a representative
#   model, and does not assert that a UJC is one molecule. Counts are never
#   touched: each protocol's matrix stays on its own transcript IDs.
#
# MONO-EXON TRANSCRIPTS ARE NEVER MATCHED
#   A UJC is defined by junctions. A mono-exon model has none, so it cannot
#   have a UJC and is always reported as protocol-exclusive. This is stricter
#   than overlap-based matching on purpose - mono-exon models replicated across
#   protocols at ~1% in this dataset, i.e. essentially never.
#
# ID SCHEME (column UJC_ID)
#   UJC<n>    intron chain seen in >= 2 protocols   (numbered in genomic order)
#   ONTd<n>   multi- or mono-exon model seen only in ONTd
#   ONTc<n>   ... only in ONTc
#   PB<n>     ... only in PB
#   (the exclusive prefixes are taken from --labels, so they follow your names)
#
# OUTPUTS (--outdir)
#   ujc_table.tsv              one row per novel entity: UJC_ID, category, and
#                              per protocol transcript_id / gene_id / biotype /
#                              class_code / spliced_length / TSS / TES; plus
#                              d5_range, d3_range, spliced_range and the
#                              concordance flags for shared entities
#   ujc_pair_deltas.tsv        one row per (shared entity x protocol pair):
#                              signed d5 / d3, delta spliced, agreement flags
#   summary_sharing.tsv        n and % per sharing category (+ TOTAL rows)
#   summary_disagreements.tsv  attribute x category: compared / discordant / %
#   summary_deltas.tsv         mean, median, IQR, p90, max of |d5|, |d3| and
#                              |delta spliced|, by category, by protocol pair,
#                              and by each kind of disagreement
#   sessionInfo.txt
#
# SIGN CONVENTION (pair deltas; protocol_2 relative to protocol_1)
#   d5 > 0  -> protocol_2 reaches further UPSTREAM   (longer 5' end)
#   d3 > 0  -> protocol_2 reaches further DOWNSTREAM (longer 3' end)
#     '+':  d5 = start_1 - start_2      d3 = end_2   - end_1
#     '-':  d5 = end_2   - end_1        d3 = start_1 - start_2
#   Per-entity *_range columns are max-min across the protocols present, so
#   they are order-independent and always >= 0.
#
# Inputs : <run>/transcriptome/curated/annotation/novel_transcripts_validated.gtf (ONT DRS, ONT cDNA, PacBio)
# Outputs: unique_intron_chains/results/
# =============================================================================

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")

## ---- CONFIG ---------------------------------------------------------------------
opt <- list(
  ## validated novel GTFs, same order as `labels` (labels are also the
  ## exclusive-ID prefixes: ONTd1, ONTc1, PB1 ...)
  gtfs   = paste(file.path(run_dir[c("ont_drs", "ont_cdna", "pb")], pp$novel_gtf), collapse = ","),
  labels = "ONTd,ONTc,PB",
  attrs  = "gene_id,pulposeq_transcript_biotype,class_code",   # GTF attributes carried through
  novel_gene_regex = "^BambuGene",   # gene_id values matching this are NOVEL loci (never compared)
  outdir = ujc_dir
)

splitc <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
gtfs <- splitc(opt$gtfs); labels <- splitc(opt$labels); keys <- splitc(opt$attrs)
if (length(gtfs) != length(labels))
  stop("--gtfs and --labels must be the same length (", length(gtfs), " vs ", length(labels), ")")
if (anyDuplicated(labels)) stop("--labels must be unique")
miss <- gtfs[!file.exists(gtfs)]
if (length(miss)) stop("missing GTF(s):\n  ", paste(miss, collapse = "\n  "))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

BIOTYPE <- "pulposeq_transcript_biotype"
CLASS   <- "class_code"
GENE    <- "gene_id"

have_dt <- requireNamespace("data.table", quietly = TRUE)
write_tsv <- function(df, f) {
  p <- file.path(opt$outdir, f)
  if (have_dt) data.table::fwrite(df, p, sep = "\t", quote = FALSE, na = "NA")
  else write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("wrote ", f, "  (", nrow(df), " rows)")
}

# ---- 1. GTF parsing ----------------------------------------------------------
# read.delim defaults to quote = "\"", which mangles field 9 -> quote = "" is
# mandatory. Attributes are parsed by splitting on ';' and matching the key
# EXACTLY; a greedy '.*transcript_id "' regex would capture cmp_ref_transcript_id.
parse_one <- function(a) {
  tk <- trimws(strsplit(a, ";", fixed = TRUE)[[1]])
  tk <- tk[nzchar(tk) & grepl("^[^[:space:]]+[[:space:]]", tk)]
  if (!length(tk)) return(setNames(character(0), character(0)))
  k <- sub("^([^[:space:]]+)[[:space:]].*$", "\\1", tk)
  v <- trimws(sub("^[^[:space:]]+[[:space:]]+", "", tk))
  v <- sub('^"', "", v); v <- sub('"$', "", v)
  setNames(v, k)
}
pull <- function(av, k) {
  L <- lapply(av, parse_one)
  vapply(L, function(x) unname(x[k])[1], NA_character_)
}
chain_of <- function(s, e) {                      # s, e = exon starts/ends, sorted
  if (length(s) < 2L) return(NA_character_)
  paste0(utils::head(e, -1) + 1, "-", utils::tail(s, -1) - 1, collapse = ",")
}

read_novel_gtf <- function(path, label) {
  g <- read.delim(path, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE, quote = "")
  if (ncol(g) < 9) stop("[", label, "] fewer than 9 columns in ", path)
  names(g)[1:9] <- c("chr","src","type","start","end","score","strand","frame","attr")
  g$start <- as.numeric(g$start); g$end <- as.numeric(g$end)
  ex <- g[g$type == "exon", , drop = FALSE]
  tr <- g[g$type == "transcript", , drop = FALSE]
  if (!nrow(ex)) stop("[", label, "] no exon records in ", path)
  ex$tid <- pull(ex$attr, "transcript_id")
  if (anyNA(ex$tid)) stop("[", label, "] ", sum(is.na(ex$tid)), " exons have no transcript_id")

  ex <- ex[order(ex$tid, ex$start), , drop = FALSE]   # chain_of needs sorted exons
  sp <- split(seq_len(nrow(ex)), ex$tid)
  tids <- names(sp)
  out <- data.frame(
    protocol       = label,
    transcript_id  = tids,
    chr            = vapply(sp, function(i) ex$chr[i[1]], ""),
    strand         = vapply(sp, function(i) ex$strand[i[1]], ""),
    start          = vapply(sp, function(i) min(ex$start[i]), 0),
    end            = vapply(sp, function(i) max(ex$end[i]), 0),
    n_exons        = vapply(sp, length, 0L),
    spliced_length = vapply(sp, function(i) sum(ex$end[i] - ex$start[i] + 1), 0),
    junctions      = vapply(sp, function(i) chain_of(ex$start[i], ex$end[i]), NA_character_),
    stringsAsFactors = FALSE)
  rownames(out) <- NULL

  # attributes: prefer transcript records, fall back to the first exon
  if (nrow(tr)) {
    ttid <- pull(tr$attr, "transcript_id")
    src_attr <- tr$attr[match(out$transcript_id, ttid)]
  } else {
    src_attr <- ex$attr[vapply(sp, function(i) i[1], 0L)][match(out$transcript_id, tids)]
  }
  for (k in keys) out[[k]] <- pull(src_attr, k)

  message(sprintf("[%s] %d novel transcripts (%d multi-exon, %d mono-exon)",
                  label, nrow(out), sum(out$n_exons > 1), sum(out$n_exons == 1)))
  for (k in keys)
    message(sprintf("       %-32s %5d/%d present", k, sum(!is.na(out[[k]])), nrow(out)))
  out
}

tx <- do.call(rbind, Map(read_novel_gtf, gtfs, labels)); rownames(tx) <- NULL
tx$protocol <- factor(tx$protocol, levels = labels)
tx$structure <- ifelse(tx$n_exons > 1, "multi-exon", "mono-exon")
# transcript-oriented ends
tx$tss <- ifelse(tx$strand == "+", tx$start, tx$end)
tx$tes <- ifelse(tx$strand == "+", tx$end,   tx$start)
# host gene = a reference locus; BambuGene* is protocol-local and never compared
tx$host_gene <- if (GENE %in% keys)
  ifelse(grepl(opt$novel_gene_regex, tx[[GENE]]), NA_character_, tx[[GENE]]) else NA_character_

n_in <- table(tx$protocol)
cat(sprintf("\ninput: %d novel transcripts across %d protocols (%s)\n",
            nrow(tx), length(labels),
            paste(sprintf("%s=%d", labels, as.integer(n_in[labels])), collapse = ", ")))

# ---- 2. UJC grouping ---------------------------------------------------------
# multi-exon: key = chr:strand:intron chain. mono-exon: no key -> exclusive.
tx$ujc_key <- ifelse(is.na(tx$junctions), NA_character_,
                     paste(tx$chr, tx$strand, tx$junctions, sep = ":"))
tx$group <- ifelse(is.na(tx$ujc_key),
                   paste0("SOLO:", tx$protocol, ":", tx$transcript_id),
                   tx$ujc_key)

gi <- split(seq_len(nrow(tx)), tx$group)
n_prot  <- vapply(gi, function(i) length(unique(as.character(tx$protocol[i]))), 0L)
n_model <- vapply(gi, length, 0L)
dupes <- sum(n_model > n_prot)
if (dupes)
  cat(sprintf("  NOTE: %d UJCs contain >1 model from the same protocol; the longest is kept per protocol\n", dupes))

g_chr   <- vapply(gi, function(i) tx$chr[i[1]], "")
g_start <- vapply(gi, function(i) min(tx$start[i]), 0)
g_lab   <- vapply(gi, function(i) as.character(tx$protocol[i[1]]), "")

shared <- n_prot >= 2L
ord <- order(factor(g_chr, levels = sort(unique(g_chr))), g_start)
UJC_ID <- character(length(gi)); names(UJC_ID) <- names(gi)
k_shared <- 0L; k_excl <- setNames(integer(length(labels)), labels)
for (j in ord) {
  if (shared[j]) { k_shared <- k_shared + 1L; UJC_ID[j] <- sprintf("UJC%d", k_shared) }
  else { L <- g_lab[j]; k_excl[[L]] <- k_excl[[L]] + 1L
         UJC_ID[j] <- sprintf("%s%d", L, k_excl[[L]]) }
}
tx$UJC_ID <- UJC_ID[tx$group]

# ---- 3. entity (wide) table --------------------------------------------------
cat_of <- function(present) {
  p <- labels[labels %in% present]
  if (length(p) == length(labels)) "shared_all"
  else if (length(p) == 1L) paste0(p, "_only")
  else paste(p, collapse = "+")
}
ent_ids <- UJC_ID[ord]
ent <- data.frame(UJC_ID = ent_ids, stringsAsFactors = FALSE)
gi_ord <- gi[ord]
ent$category    <- vapply(gi_ord, function(i) cat_of(unique(as.character(tx$protocol[i]))), "")
ent$n_protocols <- as.integer(n_prot[ord])
ent$structure   <- vapply(gi_ord, function(i) tx$structure[i[1]], "")
ent$chr         <- vapply(gi_ord, function(i) tx$chr[i[1]], "")
ent$strand      <- vapply(gi_ord, function(i) tx$strand[i[1]], "")
ent$start       <- vapply(gi_ord, function(i) min(tx$start[i]), 0)
ent$end         <- vapply(gi_ord, function(i) max(tx$end[i]), 0)
ent$n_exons     <- vapply(gi_ord, function(i) tx$n_exons[i[1]], 0L)

per_cols <- c("transcript_id", keys, "host_gene", "spliced_length", "tss", "tes")
for (L in labels) {
  # index of this protocol's model in each group (longest if duplicated)
  idx <- vapply(gi_ord, function(i) {
    k <- i[as.character(tx$protocol[i]) == L]
    if (!length(k)) NA_integer_ else k[which.max(tx$spliced_length[k])]
  }, NA_integer_)
  for (cc in per_cols) ent[[paste0(cc, "_", L)]] <- tx[[cc]][idx]
  ent[[paste0("present_", L)]] <- !is.na(idx)
}

rng <- function(m) vapply(seq_len(nrow(m)), function(i) {
  v <- m[i, ]; v <- v[!is.na(v)]
  if (length(v) < 2L) NA_real_ else max(v) - min(v) }, NA_real_)
ent$d5_range      <- rng(as.matrix(ent[, paste0("tss_", labels), drop = FALSE]))
ent$d3_range      <- rng(as.matrix(ent[, paste0("tes_", labels), drop = FALSE]))
ent$spliced_range <- rng(as.matrix(ent[, paste0("spliced_length_", labels), drop = FALSE]))

ndist <- function(m) vapply(seq_len(nrow(m)), function(i) {
  v <- m[i, ]; length(unique(v[!is.na(v)])) }, 0L)
conc_attrs <- intersect(c(BIOTYPE, CLASS, "host_gene"), names(tx))
for (k in conc_attrs) {
  m <- as.matrix(ent[, paste0(k, "_", labels), drop = FALSE])
  nd <- ndist(m)
  ent[[paste0(k, "_n_distinct")]] <- nd
  ent[[paste0(k, "_concordant")]] <- ifelse(ent$n_protocols >= 2 & nd >= 1, nd == 1L, NA)
}
write_tsv(ent, "ujc_table.tsv")

# ---- 4. pair deltas ----------------------------------------------------------
prs <- if (length(labels) >= 2) utils::combn(labels, 2, simplify = FALSE) else list()
dl <- list()
for (pr in prs) {
  a <- pr[1]; b <- pr[2]
  ok <- ent[[paste0("present_", a)]] & ent[[paste0("present_", b)]]
  if (!any(ok)) next
  sa <- ent[[paste0("tss_", a)]][ok]; ea <- ent[[paste0("tes_", a)]][ok]
  sb <- ent[[paste0("tss_", b)]][ok]; eb <- ent[[paste0("tes_", b)]][ok]
  plus <- ent$strand[ok] == "+"
  d5 <- ifelse(plus, sa - sb, sb - sa)     # tss is start on '+', end on '-'
  d3 <- ifelse(plus, eb - ea, ea - eb)     # tes is end on '+',   start on '-'
  d <- data.frame(UJC_ID = ent$UJC_ID[ok], category = ent$category[ok],
                  structure = ent$structure[ok],
                  protocol_1 = a, protocol_2 = b, strand = ent$strand[ok],
                  d5 = d5, d3 = d3, abs_d5 = abs(d5), abs_d3 = abs(d3),
                  spliced_1 = ent[[paste0("spliced_length_", a)]][ok],
                  spliced_2 = ent[[paste0("spliced_length_", b)]][ok],
                  stringsAsFactors = FALSE)
  d$delta_spliced <- d$spliced_2 - d$spliced_1
  d$abs_delta_spliced <- abs(d$delta_spliced)
  for (k in conc_attrs) {
    v1 <- ent[[paste0(k, "_", a)]][ok]; v2 <- ent[[paste0(k, "_", b)]][ok]
    d[[paste0(k, "_agree")]] <- ifelse(is.na(v1) | is.na(v2), NA, v1 == v2)
  }
  dl[[length(dl) + 1L]] <- d
}
pd <- if (length(dl)) do.call(rbind, dl) else NULL
if (!is.null(pd)) {
  pd$protocol_pair <- paste(pd$protocol_1, pd$protocol_2, sep = " vs ")
  write_tsv(pd, "ujc_pair_deltas.tsv")
  # spliced_2 - spliced_1 must equal d5 + d3 exactly: with an identical intron
  # chain the only length difference is in the two terminal exons. A failure
  # means the grouping is wrong, not that the data are odd.
  bad <- sum(abs(pd$delta_spliced - (pd$d5 + pd$d3)) > 0.5, na.rm = TRUE)
  cat(sprintf("\nidentity check  spliced_2 - spliced_1 == d5 + d3 : %d/%d exact\n",
              nrow(pd) - bad, nrow(pd)))
  if (bad) cat("  *** failures mean the intron chains are NOT identical - the UJC grouping is broken\n")
}

# ---- 5. summaries ------------------------------------------------------------
cat_levels <- c("shared_all",
                unlist(lapply(prs, function(p) paste(p, collapse = "+"))),
                paste0(labels, "_only"))
cat_levels <- cat_levels[cat_levels %in% unique(ent$category)]
ent$category <- factor(ent$category, levels = cat_levels)

tot <- nrow(ent); n_sh <- sum(ent$n_protocols >= 2)
share <- do.call(rbind, lapply(cat_levels, function(cl) {
  s <- ent[ent$category == cl, , drop = FALSE]
  data.frame(category = cl, n = nrow(s),
             pct_of_all_entities = round(100 * nrow(s) / tot, 2),
             pct_of_shared = if (all(s$n_protocols >= 2)) round(100 * nrow(s) / n_sh, 2) else NA_real_,
             n_multiexon = sum(s$structure == "multi-exon"),
             n_monoexon  = sum(s$structure == "mono-exon"),
             stringsAsFactors = FALSE)
}))
share <- rbind(share,
  data.frame(category = "TOTAL_shared(>=2 protocols)", n = n_sh,
             pct_of_all_entities = round(100 * n_sh / tot, 2), pct_of_shared = 100,
             n_multiexon = sum(ent$n_protocols >= 2 & ent$structure == "multi-exon"),
             n_monoexon  = sum(ent$n_protocols >= 2 & ent$structure == "mono-exon"),
             stringsAsFactors = FALSE),
  data.frame(category = "TOTAL_entities", n = tot, pct_of_all_entities = 100,
             pct_of_shared = NA_real_,
             n_multiexon = sum(ent$structure == "multi-exon"),
             n_monoexon  = sum(ent$structure == "mono-exon"),
             stringsAsFactors = FALSE),
  local({
    # explicit levels so the column exists even when a protocol has none
    st <- table(tx$protocol,
                factor(tx$structure, levels = c("multi-exon", "mono-exon")))
    data.frame(category = paste0("INPUT_", labels), n = as.integer(n_in[labels]),
               pct_of_all_entities = NA_real_, pct_of_shared = NA_real_,
               n_multiexon = as.integer(st[labels, "multi-exon"]),
               n_monoexon  = as.integer(st[labels, "mono-exon"]),
               stringsAsFactors = FALSE)
  }))
write_tsv(share, "summary_sharing.tsv")
cat("\n=== sharing ===\n"); print(share, row.names = FALSE)

dis <- NULL
if (!is.null(pd)) {
  rows <- list()
  for (k in conc_attrs) {
    ag <- paste0(k, "_concordant")
    for (cl in c("ALL_shared", cat_levels[grepl("shared_all|\\+", cat_levels)])) {
      s <- if (cl == "ALL_shared") ent[ent$n_protocols >= 2, ] else ent[ent$category == cl, ]
      v <- s[[ag]]; v <- v[!is.na(v)]
      if (!length(v)) next
      rows[[length(rows) + 1L]] <- data.frame(
        attribute = k, category = cl, n_compared = length(v),
        n_concordant = sum(v), n_discordant = sum(!v),
        pct_discordant = round(100 * mean(!v), 2), stringsAsFactors = FALSE)
    }
  }
  dis <- do.call(rbind, rows)
  write_tsv(dis, "summary_disagreements.tsv")
  cat("\n=== disagreements ===\n"); print(dis, row.names = FALSE)
}

if (!is.null(pd)) {
  qstat <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NULL)
    data.frame(n = length(x), mean = round(mean(x), 1),
               median = stats::median(x), sd = round(stats::sd(x), 1),
               q25 = unname(stats::quantile(x, .25)), q75 = unname(stats::quantile(x, .75)),
               p90 = unname(stats::quantile(x, .90)), max = max(x))
  }
  metrics <- c("abs_d5", "abs_d3", "abs_delta_spliced")
  blocks <- list(list("overall", "all_shared_pairs", rep(TRUE, nrow(pd))))
  for (cl in cat_levels[cat_levels %in% pd$category])
    blocks[[length(blocks) + 1L]] <- list("category", cl, pd$category == cl)
  for (p in unique(pd$protocol_pair))
    blocks[[length(blocks) + 1L]] <- list("protocol_pair", p, pd$protocol_pair == p)
  for (k in conc_attrs) {
    ag <- pd[[paste0(k, "_agree")]]
    blocks[[length(blocks) + 1L]] <- list(paste0(k, "_agreement"), "agree",
                                          !is.na(ag) & ag)
    blocks[[length(blocks) + 1L]] <- list(paste0(k, "_agreement"), "DISAGREE",
                                          !is.na(ag) & !ag)
  }
  drows <- list()
  for (b in blocks) for (m in metrics) {
    q <- qstat(pd[[m]][b[[3]]])
    if (is.null(q)) next
    drows[[length(drows) + 1L]] <- data.frame(grouping = b[[1]], group = b[[2]],
                                              metric = m, q, stringsAsFactors = FALSE)
  }
  dsum <- do.call(rbind, drows)
  write_tsv(dsum, "summary_deltas.tsv")
  cat("\n=== end deltas (bp) ===\n"); print(dsum, row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))
cat("\nDone ->", normalizePath(opt$outdir), "\n")
