#!/usr/bin/env Rscript
# =============================================================================
# build_spikein_gtfs.R
#
# Builds two GENCODE-format spike-in annotations from their original sources:
#
#   spikein_longbench.gtf   Lexogen SIRV-Set 4 in full (SIRV + long SIRV + ERCC,
#                           Lexogen coordinates) + RNA sequins v2.4
#                           196 genes / 336 transcripts / 1639 exons
#
#   spikein_encode.gtf      Lexogen SIRV-Set 4 SIRVs + long SIRVs (ENCODE's
#                           sirv4.gtf) + Thermo Fisher ERCC92
#                           120 genes / 176 transcripts / 464 exons
#
# The ERCC blocks DIFFER between the two on purpose: Lexogen and Thermo use
# different end definitions (ERCC-00002 is 1045 nt vs 1061 nt). Each file
# matches the FASTA its dataset was aligned to. Do not harmonise them.
#
# SIRV gene grouping is taken from ENCODE's sirv4.gtf in both files, because
# the Lexogen multi-fasta GTF puts sense and antisense transcripts in the SAME
# gene_id for 6 of 7 loci. Coordinates are identical between the two sources
# (verified: 0 differences across all 84 transcripts).
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- "C:/Users/mbarb/OneDrive/paper/refs/spikeins_processed"

read_gtf <- function(path) {
  x <- read.delim(path, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE, quote = "")
  names(x) <- c("seqid","source","feature","start","end",
                "score","strand","frame","attr")
  x
}

get_attr <- function(attr, key) {
  pat <- sprintf('%s "[^"]*"', key)
  out <- rep(NA_character_, length(attr))
  hit <- regexpr(pat, attr) > 0
  out[hit] <- sub(sprintf('^%s "(.*)"$', key), "\\1",
                  regmatches(attr, regexpr(pat, attr)))
  out
}

classify <- function(seqid) {
  ifelse(grepl("^ERCC-", seqid), "ERCC",
  ifelse(seqid == "chrIS",       "sequin",
  ifelse(grepl("^SIRV[0-9]{4,5}$", seqid), "long_SIRV", "SIRV_isoform")))
}
gene_class <- function(tt) ifelse(tt %in% c("SIRV_isoform","long_SIRV"), "SIRV", tt)

# ---- SIRV gene grouping + ENCODE IDs, from sirv4.gtf -------------------------
s4 <- read_gtf(file.path("C:/Users/mbarb/OneDrive/paper/refs/sirv4_rush_original/SIRV4_ENCFF873AAI.gtf/sirv4.gtf"))
s4t <- s4[s4$feature == "transcript", ]
s4t$encode  <- get_attr(s4t$attr, "transcript_id")
s4t$lexogen <- sub("^.*-(SIRV[0-9]+)$", "\\1", s4t$encode)
SIRV_GENE   <- setNames(get_attr(s4t$attr, "gene_id"), s4t$lexogen)
SIRV_ENCODE <- setNames(s4t$encode, s4t$lexogen)

# ---- collect exon records into a common shape -------------------------------
#  tx  gene  seqid  source  strand  start  end  alt_key  alt_val
exon_table <- function(gtf, tx_fun, gene_fun, alt_fun) {
  e <- gtf[gtf$feature == "exon", ]
  tx <- tx_fun(e)
  out <- data.frame(tx = tx, gene = gene_fun(e, tx), seqid = e$seqid, source = e$source,
             strand = e$strand, start = e$start, end = e$end,
             alt_key = alt_fun(e, tx)$key, alt_val = alt_fun(e, tx)$val,
             stringsAsFactors = FALSE)
  # Guard: an unresolved ID (failed regex, failed lookup) must not reach the GTF.
  if (anyNA(out$tx) || anyNA(out$gene) || any(!grepl("[A-Za-z0-9]", out$tx)))
    stop("exon_table: unresolved transcript_id or gene_id - check the regex backreferences (\\\\1, not /1) and the SIRV_GENE lookup")
  out
}

write_gencode <- function(tab, path, description, provider) {
  tab$tt <- classify(tab$seqid)
  tab$gt <- gene_class(tab$tt)
  tab <- tab[order(tab$gene, tab$tx,
                   ifelse(tab$strand == "-", -tab$start, tab$start)), ]
  tab$exon_number <- ave(seq_len(nrow(tab)), tab$tx, FUN = seq_along)

  a_gene <- function(g, gt) sprintf('gene_id "%s"; gene_type "%s"; gene_name "%s"; level 3;', g, gt, g)
  a_core <- function(r) sprintf('gene_id "%s"; transcript_id "%s"; gene_type "%s"; gene_name "%s"; transcript_type "%s"; transcript_name "%s";',
                                r$gene, r$tx, r$gt, r$gene, r$tt, r$tx)
  a_alt  <- function(r) ifelse(is.na(r$alt_val), "",
                               sprintf(' %s "%s";', r$alt_key, r$alt_val))

  # gene records
  gs <- tapply(tab$start, tab$gene, min); ge <- tapply(tab$end, tab$gene, max)
  gmeta <- tab[!duplicated(tab$gene), c("gene","seqid","source","strand","gt")]
  gene_rows <- data.frame(
    seqid = gmeta$seqid, source = gmeta$source, feature = "gene",
    start = gs[gmeta$gene], end = ge[gmeta$gene], score = ".",
    strand = gmeta$strand, frame = ".",
    attr = a_gene(gmeta$gene, gmeta$gt), stringsAsFactors = FALSE)

  # transcript records
  tmeta <- tab[!duplicated(tab$tx), ]
  ts <- tapply(tab$start, tab$tx, min); te <- tapply(tab$end, tab$tx, max)
  tx_rows <- data.frame(
    seqid = tmeta$seqid, source = tmeta$source, feature = "transcript",
    start = ts[tmeta$tx], end = te[tmeta$tx], score = ".",
    strand = tmeta$strand, frame = ".",
    attr = paste0(a_core(tmeta), " level 3; tag \"basic\";", a_alt(tmeta)),
    stringsAsFactors = FALSE)

  # exon records
  ex_rows <- data.frame(
    seqid = tab$seqid, source = tab$source, feature = "exon",
    start = tab$start, end = tab$end, score = ".",
    strand = tab$strand, frame = ".",
    attr = paste0(a_core(tab),
                  sprintf(' exon_number %d; exon_id "%s.E%d";', tab$exon_number, tab$tx, tab$exon_number),
                  " level 3; tag \"basic\";", a_alt(tab)),
    stringsAsFactors = FALSE)

  all <- rbind(gene_rows, tx_rows, ex_rows)
  rk  <- match(all$feature, c("gene","transcript","exon"))
  all <- all[order(all$seqid, all$start, rk, all$end), ]

  stopifnot(!anyNA(tab$tx), !anyNA(tab$gene),
            !anyDuplicated(tmeta$tx),
            setequal(get_attr(gene_rows$attr, "gene_id"),
                     get_attr(rbind(tx_rows, ex_rows)$attr, "gene_id")))

  writeLines(c(sprintf("##description: %s", description),
               sprintf("##provider: %s", provider),
               "##format: gtf",
               sprintf("##date: %s", Sys.Date())), path)
  write.table(all, path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE, append = TRUE)
  message(sprintf("%s: %d genes / %d transcripts / %d exons",
                  basename(path), nrow(gene_rows), nrow(tx_rows), nrow(ex_rows)))
}

# ============================== LongBench ====================================
lex <- read_gtf(file.path("C:/Users/mbarb/OneDrive/paper/refs/sirv4_original/SIRV_Set4_Norm_Sequences_20210507/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf"))
lex_tab <- exon_table(
  lex,
  tx_fun   = function(e) ifelse(grepl("^ERCC-", e$seqid), e$seqid,
                                get_attr(e$attr, "transcript_id")),
  gene_fun = function(e, tx) ifelse(grepl("^ERCC-", e$seqid), e$seqid,
                                    unname(SIRV_GENE[tx])),
  alt_fun  = function(e, tx) list(
     key = ifelse(grepl("^ERCC-", e$seqid), "genbank_id", "encode_transcript_id"),
     val = ifelse(grepl("^ERCC-", e$seqid), get_attr(e$attr, "transcript_id"),
                  unname(SIRV_ENCODE[tx]))))

seq <- read_gtf(file.path("C:/Users/mbarb/OneDrive/paper/refs/sequins_original/rnasequin_annotation_2.4.gtf"))
seq_tab <- exon_table(
  seq,
  tx_fun   = function(e) get_attr(e$attr, "transcript_id"),
  gene_fun = function(e, tx) get_attr(e$attr, "gene_id"),
  alt_fun  = function(e, tx) list(key = NA_character_, val = NA_character_))

write_gencode(rbind(lex_tab, seq_tab),
  file.path(out_dir, "spikein_longbench.gtf"),
  "Lexogen SIRV-Set 4 (SIRV + long SIRV + ERCC) and RNA sequins v2.4 - LongBench",
  "SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf + rnasequin_annotation_2.4.gtf")

# ================================ ENCODE =====================================
enc_tab <- exon_table(
  s4,
  tx_fun   = function(e) sub("^.*-(SIRV[0-9]+)$", "\\1", get_attr(e$attr, "transcript_id")),
  gene_fun = function(e, tx) get_attr(e$attr, "gene_id"),
  alt_fun  = function(e, tx) list(key = "encode_transcript_id",
                                  val = get_attr(e$attr, "transcript_id")))

ercc <- read_gtf(file.path("C:/Users/mbarb/OneDrive/paper/refs/ambion ercc spike ins/ERCC92.gtf"))
ercc_tab <- exon_table(
  ercc,
  tx_fun   = function(e) e$seqid,
  gene_fun = function(e, tx) e$seqid,
  alt_fun  = function(e, tx) list(key = "genbank_id",
                                  val = get_attr(e$attr, "transcript_id")))

write_gencode(rbind(enc_tab, ercc_tab),
  file.path(out_dir, "spikein_encode.gtf"),
  "Lexogen SIRV-Set 4 (SIRV + long SIRV) and Thermo Fisher ERCC92 - ENCODE RUSH AD",
  "ENCODE sirv4.gtf + ERCC92.gtf")
