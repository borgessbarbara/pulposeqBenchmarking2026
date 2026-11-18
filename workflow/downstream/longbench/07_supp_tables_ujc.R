# =============================================================================
# 07_supp_tables_ujc.R - builds the UJC supplementary tables (drafted as ST2-ST7)
# NOTE: the ST numbers in the sheet names and captions are the draft numbers;
# the final workbook (pulposeq_supplementary_tables.xlsx) renumbers them.
# Emits one .tsv per table, a single .xlsx with a sheet
# per table when openxlsx or writexl is available, and a draft caption for each.
# UNITS DIFFER BETWEEN TABLES AND THIS IS THE EASIEST THING TO GET WRONG:
#   ST2-ST5 count UJC ENTITIES (one row per distinct intron chain)
#   ST6-ST7 count PROTOCOL PAIRS (a 3-protocol entity contributes three)
# Every caption states its unit. Do not mix them in the text.
# Inputs : unique_intron_chains/results/ujc_table.tsv, unique_intron_chains/results/ujc_pair_deltas.tsv
# Outputs: unique_intron_chains/results/supp_tables/ (one .tsv per table, captions .txt, .xlsx)
# =============================================================================
source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")

opt <- list(ujc    = file.path(ujc_dir, "ujc_table.tsv"),
            pairs  = file.path(ujc_dir, "ujc_pair_deltas.tsv"),
            outdir = supp_dir)
for (f in c(opt$ujc, opt$pairs)) if (!file.exists(f)) stop("no such file: ", f)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

have_dt <- requireNamespace("data.table", quietly = TRUE)
# NB the braces are load-bearing: with a brace-less function body, R closes the
# definition at the end of the `if` and the `else` on the next line is a syntax
# error at top level.
rd <- function(p) {
  if (have_dt) as.data.frame(data.table::fread(p, sep = "\t", data.table = FALSE))
  else read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
}
ent <- rd(opt$ujc); pd <- rd(opt$pairs)

# protocol labels, in the order the analysis used them
labels <- sub("^present_", "", grep("^present_", names(ent), value = TRUE))
message("protocols: ", paste(labels, collapse = ", "))
BIOTYPE <- "pulposeq_transcript_biotype"; CLASS <- "class_code"

# ---- formatting helpers ------------------------------------------------------
cm  <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
np  <- function(n, d) ifelse(is.na(d) | d == 0, cm(n),
                             sprintf("%s (%.1f%%)", cm(n), 100 * n / d))
wilson <- function(x, n, conf = 0.95) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2); ph <- x / n; dn <- 1 + z^2 / n
  ctr <- (ph + z^2 / (2 * n)) / dn
  hw  <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / dn
  c(max(0, ctr - hw), min(1, ctr + hw))
}
pct_ci <- function(x, n) {
  ci <- wilson(x, n)
  sprintf("%.2f%% (%.2f-%.2f)", 100 * x / n, 100 * ci[1], 100 * ci[2])
}
tabs <- list(); caps <- character()
add <- function(id, df, cap) {
  tabs[[id]] <<- df; caps[[id]] <<- cap
  write.table(df, file.path(opt$outdir, paste0(id, ".tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "")
  message("built ", id, " (", nrow(df), " rows)")
}

# ============================== ST2 ===========================================
st2 <- do.call(rbind, lapply(labels, function(L) {
  pres <- ent[[paste0("present_", L)]]
  s <- ent[pres, , drop = FALSE]
  sl <- s[[paste0("spliced_length_", L)]]
  hg <- s[[paste0("host_gene_", L)]]
  multi <- sum(s$structure == "multi-exon"); mono <- sum(s$structure == "mono-exon")
  data.frame(
    Protocol = L,
    `Novel transcript models` = cm(nrow(s)),
    `Multi-exon` = np(multi, nrow(s)),
    `Mono-exon`  = np(mono,  nrow(s)),
    `Median spliced length (nt)` = cm(round(stats::median(sl, na.rm = TRUE))),
    `Spliced length IQR (nt)` = sprintf("%s-%s",
        cm(round(stats::quantile(sl, .25, na.rm = TRUE))),
        cm(round(stats::quantile(sl, .75, na.rm = TRUE)))),
    `Assigned to a reference locus` = np(sum(!is.na(hg)), nrow(s)),
    `Assigned to a novel locus` = np(sum(is.na(hg)), nrow(s)),
    check.names = FALSE, stringsAsFactors = FALSE)
}))
add("ST2_discovery_per_protocol", st2, paste0(
  "Supplementary Table 2 | Novel transcript models called by each long-read protocol. ",
  "Unit: novel transcript models, before any cross-protocol matching. Percentages are of that ",
  "protocol's own novel models. 'Reference locus' means Bambu assigned the model to an annotated ",
  "gene; 'novel locus' means it received a new BambuGene identifier. Bambu NDR = 0.1, stranded; ",
  "models were filtered to >=1 count in >=1 sample and to the pulposeq class-code set before this table."))

# ============================== ST3 ===========================================
comp_block <- function(attr, nice) {
  vals <- sort(unique(unlist(lapply(labels, function(L) {
    v <- ent[[paste0(attr, "_", L)]][ent[[paste0("present_", L)]]]; v[!is.na(v)] }))))
  tot <- vapply(labels, function(L) sum(ent[[paste0("present_", L)]]), 0L)
  out <- data.frame(Attribute = nice, Value = vals, check.names = FALSE,
                    stringsAsFactors = FALSE)
  for (i in seq_along(labels)) {
    L <- labels[i]
    v <- ent[[paste0(attr, "_", L)]][ent[[paste0("present_", L)]]]
    out[[L]] <- vapply(vals, function(x) np(sum(v == x, na.rm = TRUE), tot[i]), "")
  }
  ord <- order(-vapply(seq_along(vals), function(j) sum(vapply(labels, function(L) {
    v <- ent[[paste0(attr, "_", L)]][ent[[paste0("present_", L)]]]
    sum(v == vals[j], na.rm = TRUE) }, 0L)), 0))
  out[ord, , drop = FALSE]
}
st3 <- rbind(comp_block(CLASS, "GffCompare class code"),
             comp_block(BIOTYPE, "pulposeq transcript biotype"))
add("ST3_composition_of_novel_calls", st3, paste0(
  "Supplementary Table 3 | Composition of each protocol's novel transcript models. ",
  "Unit: novel transcript models. Percentages are of that protocol's own novel models, so each ",
  "block sums to 100% per column. Class codes are GffCompare's; biotypes are pulposeq's ",
  "assignment. Rows are ordered by total frequency across protocols."))

# ============================== ST4 ===========================================
cats <- unique(ent$category)
ord_cat <- c("shared_all",
             unlist(lapply(utils::combn(labels, 2, simplify = FALSE), paste, collapse = "+")),
             paste0(labels, "_only"))
cats <- ord_cat[ord_cat %in% cats]
tot <- nrow(ent); nsh <- sum(ent$n_protocols >= 2)
memb <- function(cl) {
  if (cl == "shared_all") return(paste(labels, collapse = ", "))
  if (grepl("_only$", cl)) return(sub("_only$", "", cl))
  gsub("\\+", ", ", cl)
}
st4 <- do.call(rbind, lapply(cats, function(cl) {
  s <- ent[ent$category == cl, , drop = FALSE]
  data.frame(Category = cl, Protocols = memb(cl),
             `UJC entities` = cm(nrow(s)),
             `% of all entities` = sprintf("%.2f%%", 100 * nrow(s) / tot),
             `% of shared entities` = if (all(s$n_protocols >= 2))
               sprintf("%.2f%%", 100 * nrow(s) / nsh) else "-",
             `Multi-exon` = cm(sum(s$structure == "multi-exon")),
             `Mono-exon`  = cm(sum(s$structure == "mono-exon")),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
st4 <- rbind(st4,
  data.frame(Category = "Shared (>=2 protocols)", Protocols = "-",
             `UJC entities` = cm(nsh), `% of all entities` = sprintf("%.2f%%", 100 * nsh / tot),
             `% of shared entities` = "100.00%",
             `Multi-exon` = cm(sum(ent$n_protocols >= 2 & ent$structure == "multi-exon")),
             `Mono-exon`  = cm(sum(ent$n_protocols >= 2 & ent$structure == "mono-exon")),
             check.names = FALSE, stringsAsFactors = FALSE),
  data.frame(Category = "All entities", Protocols = "-", `UJC entities` = cm(tot),
             `% of all entities` = "100.00%", `% of shared entities` = "-",
             `Multi-exon` = cm(sum(ent$structure == "multi-exon")),
             `Mono-exon`  = cm(sum(ent$structure == "mono-exon")),
             check.names = FALSE, stringsAsFactors = FALSE))
add("ST4_ujc_sharing", st4, paste0(
  "Supplementary Table 4 | Novel transcript models grouped by unique junction chain (UJC). ",
  "Unit: UJC entities - one row per distinct intron chain, not per model. A UJC is defined by ",
  "junctions alone, so MONO-EXON MODELS CARRY NO UJC AND ARE NEVER MATCHED; each is reported as ",
  "protocol-exclusive. Transcript ends are not part of the match (see ST6)."))

# ============================== ST5 ===========================================
attr_nice <- c(setNames("GffCompare class code", CLASS),
               setNames("pulposeq transcript biotype", BIOTYPE),
               host_gene = "Reference host gene")
rows <- list()
for (k in names(attr_nice)) {
  cc <- paste0(k, "_concordant"); if (!cc %in% names(ent)) next
  grp <- c("All shared" = NA, setNames(cats[!grepl("_only$", cats)],
                                       cats[!grepl("_only$", cats)]))
  for (nm in names(grp)) {
    s <- if (is.na(grp[[nm]])) ent[ent$n_protocols >= 2, ] else ent[ent$category == grp[[nm]], ]
    v <- s[[cc]]; v <- v[!is.na(v)]
    if (!length(v)) next
    rows[[length(rows) + 1L]] <- data.frame(
      Attribute = attr_nice[[k]], `Sharing category` = nm,
      `UJC entities compared` = cm(length(v)),
      Concordant = cm(sum(v)), Discordant = cm(sum(!v)),
      `% discordant (95% CI)` = pct_ci(sum(!v), length(v)),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}
st5 <- do.call(rbind, rows)
add("ST5_classification_concordance", st5, paste0(
  "Supplementary Table 5 | Agreement between protocols on the classification of the same unique ",
  "junction chain. Unit: UJC entities reconstructed by >=2 protocols. An entity is concordant when ",
  "every protocol that reconstructed it assigned the same value. Entities for which no protocol ",
  "assigned a reference host gene are excluded from that row. Intervals are Wilson 95% - several ",
  "rates rest on fewer than ten events and the interval, not the point estimate, is the result."))

# ============================== ST6 ===========================================
pd$protocol_pair <- paste(pd$protocol_1, pd$protocol_2, sep = " vs ")
mets <- c(abs_d5 = "5' end |d5|", abs_d3 = "3' end |d3|",
          abs_delta_spliced = "spliced length |delta|")
mets <- mets[names(mets) %in% names(pd)]
blocks <- c(list(list("All shared pairs", rep(TRUE, nrow(pd)))),
            lapply(cats[!grepl("_only$", cats)], function(cl) list(cl, pd$category == cl)),
            lapply(unique(pd$protocol_pair), function(p) list(p, pd$protocol_pair == p)))
rows <- list()
for (b in blocks) for (m in names(mets)) {
  x <- pd[[m]][b[[2]]]; x <- x[!is.na(x)]
  if (!length(x)) next
  rows[[length(rows) + 1L]] <- data.frame(
    Group = b[[1]], Metric = mets[[m]], `Protocol pairs` = cm(length(x)),
    `Identical (0 bp)` = sprintf("%.1f%%", 100 * mean(x == 0)),
    `<=10 bp` = sprintf("%.1f%%", 100 * mean(x <= 10)),
    `<=50 bp` = sprintf("%.1f%%", 100 * mean(x <= 50)),
    `<=500 bp` = sprintf("%.1f%%", 100 * mean(x <= 500)),
    Median = cm(stats::median(x)),
    `90th pct` = cm(round(unname(stats::quantile(x, .90)))),
    Maximum = cm(max(x)), check.names = FALSE, stringsAsFactors = FALSE)
}
st6 <- do.call(rbind, rows)
add("ST6_transcript_end_differences", st6, paste0(
  "Supplementary Table 6 | Transcript-end differences between protocols reconstructing an identical ",
  "intron chain. Unit: PROTOCOL PAIRS - an entity found by all three contributes three rows, so the ",
  "totals here exceed the entity counts in ST4-ST5. Because the intron chain is identical by ",
  "construction, every difference is terminal-exon extent, and spliced_2 - spliced_1 = d5 + d3 exactly. ",
  "|delta spliced| is therefore largely redundant with |d3| and is given for completeness."))

# ============================== ST7 ===========================================
# ujc_pair_deltas.tsv carries only the *_agree flags, not the attribute VALUES,
# so the values are looked up per row from ujc_table.tsv using that row's own
# protocol_1 / protocol_2. (Row-wise, because the pair of protocols differs
# between rows - a plain column subset would take the wrong protocol's value.)
getval <- function(k, prot, idx) {
  out <- rep(NA_character_, length(prot))
  for (L in labels) {
    col <- paste0(k, "_", L)
    if (!col %in% names(ent)) next
    sel <- which(prot == L)
    if (length(sel)) out[sel] <- as.character(ent[[col]][idx[sel]])
  }
  out
}
rows <- list()
for (k in names(attr_nice)) {
  ag <- paste0(k, "_agree"); if (!ag %in% names(pd)) next
  s <- pd[!is.na(pd[[ag]]) & !pd[[ag]], , drop = FALSE]
  if (!nrow(s)) next
  idx <- match(s$UJC_ID, ent$UJC_ID)
  if (anyNA(idx)) stop("ST7: ", sum(is.na(idx)), " UJC_IDs in the pair file are absent from ",
                       opt$ujc, " - are the two files from the same run?")
  val1 <- getval(k, s$protocol_1, idx)
  val2 <- getval(k, s$protocol_2, idx)
  h1 <- getval("host_gene", s$protocol_1, idx)
  h2 <- getval("host_gene", s$protocol_2, idx)
  host_switch <- (!is.na(h1) | !is.na(h2)) &
                 (is.na(h1) != is.na(h2) | (!is.na(h1) & !is.na(h2) & h1 != h2))
  cross500 <- !is.na(s$spliced_1) & !is.na(s$spliced_2) &
              ((s$spliced_1 < 500) != (s$spliced_2 < 500))
  mech <- if (k == BIOTYPE)
    ifelse(host_switch, "host gene switch",
      ifelse(cross500, "crosses the 500 nt threshold", "coding potential")) else
    ifelse(abs(s$d3) > abs(s$d5), "3' end", "5' end")
  rows[[length(rows) + 1L]] <- data.frame(
    UJC_ID = s$UJC_ID, Attribute = attr_nice[[k]], `Sharing category` = s$category,
    `Protocol 1` = s$protocol_1, `Value 1` = val1,
    `Protocol 2` = s$protocol_2, `Value 2` = val2,
    `Host gene 1` = h1, `Host gene 2` = h2,
    Strand = s$strand, `d5 (bp)` = s$d5, `d3 (bp)` = s$d3,
    `Spliced 1 (nt)` = s$spliced_1, `Spliced 2 (nt)` = s$spliced_2,
    `Delta spliced (nt)` = s$delta_spliced, `Likely driver` = mech,
    check.names = FALSE, stringsAsFactors = FALSE)
}
st7 <- do.call(rbind, rows)
st7 <- st7[order(st7$Attribute, -abs(st7$`Delta spliced (nt)`)), , drop = FALSE]
add("ST7_itemised_disagreements", st7, paste0(
  "Supplementary Table 7 | Every classification disagreement between protocols, itemised. ",
  "Unit: PROTOCOL PAIRS. d5 > 0 means protocol 2 reaches further upstream; d3 > 0 means it reaches ",
  "further downstream. 'Likely driver' for biotype distinguishes a change of reference host gene, a ",
  "crossing of the 500 nt lncRNA minimum, and a coding-potential flip at constant host; for the other ",
  "attributes it names the end with the larger absolute difference. These are inferences from the ",
  "coordinates, not from re-running the classifier."))

# ---- write ------------------------------------------------------------------
writeLines(unlist(lapply(names(caps), function(i)
  paste0(caps[[i]], "\n"))), file.path(opt$outdir, "supplementary_table_captions.txt"))
xlsx <- file.path(opt$outdir, "supplementary_tables_ST2-ST7.xlsx")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(textDecoration = "bold", border = "bottom",
                               valign = "center", wrapText = TRUE)
  for (id in names(tabs)) {
    sh <- substr(sub("_.*$", "", id), 1, 31)
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, tabs[[id]], headerStyle = hdr)
    openxlsx::freezePane(wb, sh, firstRow = TRUE)
    openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(tabs[[id]])), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE); message("wrote ", basename(xlsx))
} else if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(setNames(tabs, sub("_.*$", "", names(tabs))), xlsx)
  message("wrote ", basename(xlsx), " (install openxlsx for formatted headers)")
} else message("neither openxlsx nor writexl installed - .tsv files only")
cat("\nDone ->", normalizePath(opt$outdir), "\n")
