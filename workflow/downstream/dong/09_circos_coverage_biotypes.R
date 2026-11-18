#!/usr/bin/env Rscript
# =============================================================================
# 09_circos_coverage_biotypes.R - genome-wide circos of the Dong assembly:
#   RNA coverage per cell line, annotated transcripts (density) and pulposeq /
#   Bambu novel transcripts (exact positions), one ring per strand

# Layout (outside -> inside):
#   chromosome names
#   ideogram (hg38 cytobands if downloadable, otherwise plain chromosomes)
#   rings C1-C2  RNA coverage, mean of replicates per cell line (bigWig)
#                (HCC827 outside, H1975 inside; label "RNA cov." in the gap)
#   ring + strand  annotated stacked bars per bin (protein-coding / lncRNA / other)
#                  + novel ticks in the outer band (one tick = one transcript)
#   ring - strand  same, same scale
#   centre       legend: "Transcript types" (ref / novel colour pairs) and
#                "RNA coverage" (HCC827 / H1975)
#
# Inputs (read only): <run>/transcriptome/curated/metadata/annotated_transcriptome_metadata.csv,
#   <run>/transcriptome/curated/metadata/novel_transcripts_validated_metadata.csv,
#   <run>/coverage/*.bw (one per sample)
# Outputs: figures/circos_coverage/circos_coverage_biotypes.{pdf,png,svg}
#          results/circos_coverage/circos_coverage_biotypes_{annotated_density,coverage}_bins.tsv
# Needs  : data.table, circlize (CRAN); rtracklayer (Bioconductor) for the
#          coverage rings; svglite for editable SVG text; internet for the hg38
#          cytobands (otherwise a plain ideogram is drawn)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
})

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")   # data_dir, helper_dir, res_dir, fig_dir

# -----------------------------------------------------------------------------
# 0. Configuration -- edit here
# -----------------------------------------------------------------------------
ANNOT_CSV  <- run_path("ann_meta")
NOVEL_CSV  <- run_path("novel_meta")
BW_DIR     <- run_path("coverage_dir")                    # H1975_1.bw ... HCC827_3.bw
OUT_FIG    <- file.path(fig_dir, "circos_coverage")
OUT_TAB    <- file.path(res_dir, "circos_coverage")
OUT_NAME   <- "circos_coverage_biotypes"
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# font: Source Sans 3, registered by the project theme (stops if it can't be)
THEME_FILE <- file.path(helper_dir, "plot_theme.R")
source(THEME_FILE)
FONT <- font_family

CHROMS       <- paste0("chr", c(1:22, "X"))  # chrY dropped: female samples, chrY
                                              # transcripts are PAR copies (also on chrX)
BIN_SIZE     <- 1e6       # density / coverage bin (1 Mb = convention in the papers cited)
POSITION     <- "tss"     # "tss" (strand-aware 5' end) or "mid" (transcript midpoint)
ANNOT_GROUP_BY <- "transcript_biotype"  # or "gene_biotype" (see notes at bottom)
DROP_PSEUDOGENES <- TRUE  # TRUE: *pseudogene* transcripts removed entirely
                          # (FALSE: they are counted in "Other")
DENSITY_CAP_Q  <- 0.995   # cap each ring at this quantile so 1-2 extreme bins
                          # (e.g. chr6 HLA / chr19 ZNF clusters) don't flatten it
USE_CYTOBAND <- TRUE      # try to download hg38 cytobands from UCSC (needs internet)
FAI          <- NULL      # optional: path to genome .fai to take chromosome lengths from
GAP_START    <- 30        # degrees of empty space at 12 o'clock for the ring labels

# --- coverage rings (bigWig) ----------------------------------------------------
COVERAGE     <- TRUE      # FALSE = no coverage rings (rtracklayer not needed)
BW_PATTERN   <- "\\.bw$"
BW_GROUP_RE  <- "_[0-9]+$"            # replicate suffix removed to get the cell line
COV_GROUPS   <- c("HCC827", "H1975")  # ring order, outside -> inside
COV_NORMALISE <- TRUE     # scale each sample to the same total binned signal
                          # (harmless if the bigWigs are already CPM-normalised)
COV_COLS     <- c(HCC827 = "#4DBBD5", H1975 = "#E64B35")

# Colours: the same hue links an annotated class to its novel counterpart
# (light = annotated density, saturated = novel positions); the novel colours
# are the cis-acting circos colours
## the biotype palette of every figure (as make_palette in 06_isoformic.R)
ANNOT_COLS <- c(
  "Protein-coding" = "#A6D854",
  "lncRNA"         = "#E78AC3",
  "Other"          = "#CFCFCF"
)
NOVEL_COLS <- c(
  "Novel protein-coding" = "#238B45",
  "Novel lncRNA"         = "#AE017E",
  "Novel non-coding"     = "#6A51A3"
)
KEY_LABELS <- c("ref/novel protein-coding", "ref/novel lncRNA", "ref other / novel non-coding")
INK       <- "#0b0b0b"
RING_BG   <- "#f3f2ef"
STRAND_BG <- c("+" = "#f6f5f2", "-" = "#ebe9e4")
BAR_FRAC  <- 0.68         # annotated bars use the inner 68% of a strand ring ...
TICK_FROM <- 0.72         # ... novel ticks run from 72% to the outer edge
NOVEL_LWD <- 1.5

# --- page geometry (as the final figure) -----------------------------------------
# radius 1 = outer edge of the chromosome-name track; track heights are
# fractions of that radius. The canvas leaves the space the figure has around
# the circle (more room below than above).
FIG_IN    <- c(8, 7.49)                 # width, height (in)
CANVAS_X  <- c(-1.136, 1.136)
CANVAS_Y  <- c(-1.141, 0.986)
H_CHR     <- 0.048                      # chromosome names
H_IDEO    <- 0.031                      # ideogram
H_COV     <- 0.045                      # each coverage ring
H_STRAND  <- 0.18                       # each strand ring
# text (pt): chromosome names, ring labels in the gap, centre legend
PT_CHR    <- 13.4
PT_RING   <- 9.3
PT_KEY    <- 10.5
KEY_FIT   <- 0.94    # centre legend must fit inside this fraction of the hole radius

# GRCh38 primary-assembly lengths (fallback when no .fai / cytoband available); chrY omitted
HG38_LEN <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
  chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
  chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
  chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
  chr21 = 46709983, chr22 = 50818468, chrX = 156040895
)

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

norm_strand <- function(s) {
  s <- as.character(s)
  fcase(s %in% c("+", "1"), "+", s %in% c("-", "-1"), "-", default = NA_character_)
}

anchor_pos <- function(start, end, strand, mode = POSITION) {
  if (mode == "tss") fifelse(strand == "-", end, start)
  else as.numeric(round((start + end) / 2))
}

group_annotated <- function(bt) {
  fcase(
    bt == "protein_coding",      "Protein-coding",
    bt == "lncRNA",              "lncRNA",
    grepl("pseudogene", bt),     "Pseudogene",   # processed, unprocessed, transcribed_*,
                                                 # unitary, rRNA_, IG_*/TR_* pseudogenes
    default =                    "Other"
  )
}

group_novel <- function(bt) {
  fcase(
    bt == "novel_protein_coding", "Novel protein-coding",
    bt == "novel_lncRNA",         "Novel lncRNA",
    bt == "novel_non_coding",     "Novel non-coding",
    default = NA_character_
  )
}

check_cols <- function(dt, cols, file) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(sprintf("%s is missing column(s): %s", file, paste(miss, collapse = ", ")))
}

cap_of <- function(v) {                     # quantile cap over non-zero values
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(1)
  quantile(v, DENSITY_CAP_Q, names = FALSE)
}

# -----------------------------------------------------------------------------
# 2. Chromosome lengths + bin grid
# -----------------------------------------------------------------------------
if (!is.null(FAI)) {
  fai     <- fread(FAI, header = FALSE, select = 1:2, col.names = c("chr", "len"))
  chr_len <- setNames(fai$len, fai$chr)[CHROMS]
} else {
  chr_len <- HG38_LEN[CHROMS]
}
if (anyNA(chr_len)) stop("No length for: ", paste(CHROMS[is.na(chr_len)], collapse = ", "))
chr_dt <- data.table(chr = factor(CHROMS, levels = CHROMS), start = 0, end = as.numeric(chr_len))

bins <- chr_dt[, .(bin = 0:(ceiling(end / BIN_SIZE) - 1), len = end), by = chr]
bins[, `:=`(start = bin * BIN_SIZE, end = pmin((bin + 1) * BIN_SIZE, len))]

cytoband <- NULL
if (USE_CYTOBAND) {
  cytoband <- tryCatch(
    read.cytoband(species = "hg38", chromosome.index = CHROMS)$df,
    error = function(e) { message("Cytoband download failed -> plain ideogram (", conditionMessage(e), ")"); NULL }
  )
}

# -----------------------------------------------------------------------------
# 3. Read + tidy annotated transcripts
# -----------------------------------------------------------------------------
ann <- fread(ANNOT_CSV)
check_cols(ann, c("chromosome_name", "strand", "transcript_start", "transcript_end", ANNOT_GROUP_BY), ANNOT_CSV)

n_ann_all <- nrow(ann)
ann <- ann[chromosome_name %in% CHROMS]
ann[, `:=`(
  chr    = factor(chromosome_name, levels = CHROMS),
  strand = norm_strand(strand),
  group  = group_annotated(get(ANNOT_GROUP_BY))
)]
n_pseudo <- ann[group == "Pseudogene", .N]
if (DROP_PSEUDOGENES) ann <- ann[group != "Pseudogene"] else ann[group == "Pseudogene", group := "Other"]
ann <- ann[!is.na(strand)]
ann[, pos := anchor_pos(transcript_start, transcript_end, strand)]
message(sprintf("Annotated: %s transcripts; %s plotted (%s-%s, pseudogenes %s: %s)",
                fmt_n(n_ann_all), fmt_n(nrow(ann)), CHROMS[1], CHROMS[length(CHROMS)],
                if (DROP_PSEUDOGENES) "removed" else "in Other", fmt_n(n_pseudo)))
print(ann[, .N, by = .(group, strand)][order(group, strand)])

# -----------------------------------------------------------------------------
# 4. Read + tidy novel transcripts
# -----------------------------------------------------------------------------
nov <- fread(NOVEL_CSV)
check_cols(nov, c("seqnames", "start", "end", "strand", "pulposeq_transcript_biotype"), NOVEL_CSV)

n_nov_all <- nrow(nov)
nov <- nov[seqnames %in% CHROMS]
nov[, `:=`(
  chr    = factor(seqnames, levels = CHROMS),
  strand = norm_strand(strand),
  group  = group_novel(pulposeq_transcript_biotype)
)]
if (anyNA(nov$group)) {
  warning("Unrecognised pulposeq_transcript_biotype values dropped: ",
          paste(unique(nov[is.na(group), pulposeq_transcript_biotype]), collapse = ", "))
  nov <- nov[!is.na(group)]
}
nov <- nov[!is.na(strand)]
nov[, pos := anchor_pos(start, end, strand)]
message(sprintf("Novel: %s transcripts, %s plotted (dropped %s on chrY/scaffolds/unstranded)",
                fmt_n(n_nov_all), fmt_n(nrow(nov)), fmt_n(n_nov_all - nrow(nov))))
print(nov[, .N, by = .(group, strand)][order(group, strand)])

# -----------------------------------------------------------------------------
# 5. Bin annotated transcripts per group x strand (zero-filled)
# -----------------------------------------------------------------------------
grid <- bins[, .(group = rep(names(ANNOT_COLS), each = 2),
                 strand = rep(c("+", "-"), times = length(ANNOT_COLS))),
             by = .(chr, bin, start, end)]

ann[, bin := as.integer(pos %/% BIN_SIZE)]
counts <- ann[, .(n = .N), by = .(chr, group, strand, bin)]

dens <- counts[grid, on = .(chr, group, strand, bin)]
dens[is.na(n), n := 0L]
# stacked bars: protein-coding at the base, then lncRNA, then other. One cap on
# the bin TOTAL, shared by both strands; a capped bin keeps its composition.
dens[, group := factor(group, levels = names(ANNOT_COLS))]
dens[, total := sum(n), by = .(chr, bin, strand)]
stack_cap <- max(1, ceiling(cap_of(unique(dens[, .(chr, bin, strand, total)])$total)))
dens[, k := fifelse(total > stack_cap, stack_cap / total, 1)]
setorder(dens, chr, bin, strand, group)
dens[, ytop := cumsum(n * k), by = .(chr, bin, strand)]
dens[, ybot := ytop - n * k]
dens[, group := as.character(group)]

# -----------------------------------------------------------------------------
# 6. Coverage per bin from bigWigs (mean over ALL bases of the bin)
# -----------------------------------------------------------------------------
cov <- NULL
if (COVERAGE) {
  # bare file names: the bigWig reader (UCSC library) takes the text before the first
  # ":" as a URL protocol, so Windows paths like C:\... fail ("Unrecognized protocol C");
  # the files are opened from inside BW_DIR instead
  bw_files <- list.files(BW_DIR, pattern = BW_PATTERN, full.names = FALSE)
  if (!length(bw_files)) {
    message("No bigWig files in ", BW_DIR, " -> coverage rings skipped")
  } else if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("rtracklayer not installed -> coverage rings skipped ",
            "(BiocManager::install(\"rtracklayer\"))")
  } else {
    # attached (not ::) so summary() dispatches to rtracklayer's BigWigFile method;
    # all data.table work above is already done
    suppressPackageStartupMessages(library(rtracklayer))
    bin_signal <- function(f) {
      bwf <- BigWigFile(f)
      sl  <- GenomeInfoDb::seqlengths(bwf)
      use_chr <- any(grepl("^chr", names(sl)))            # "chr1" vs "1" naming
      sq  <- if (use_chr) as.character(bins$chr) else sub("^chr", "", as.character(bins$chr))
      ok  <- sq %in% names(sl)
      gr  <- GenomicRanges::GRanges(sq[ok], IRanges::IRanges(bins$start[ok] + 1,
                                                             pmin(bins$end[ok], sl[sq[ok]])))
      summ <- function(type) {
        s <- summary(bwf, which = gr, size = 1L, type = type, defaultValue = 0)
        v <- if (methods::is(s, "GRangesList")) unlist(s, use.names = FALSE)$score
             else vapply(as.list(s), function(g) g$score[1], numeric(1))
        v[is.na(v)] <- 0
        v
      }
      out <- rep(0, nrow(bins))
      # bigWig "mean" is over covered bases only -> multiply by the covered fraction
      out[ok] <- summ("mean") * summ("coverage")
      out
    }
    message("Reading ", length(bw_files), " bigWig files from ", normalizePath(BW_DIR), " ...")
    owd <- setwd(BW_DIR)                              # restored even if reading fails
    sig <- tryCatch(vapply(bw_files, bin_signal, numeric(nrow(bins))),
                    finally = setwd(owd))
    colnames(sig) <- sub(BW_PATTERN, "", basename(bw_files))
    if (COV_NORMALISE) {
      tot <- colSums(sig * (bins$end - bins$start))
      sig <- sweep(sig, 2, tot / mean(tot), "/")
    }
    smp_group <- sub(BW_GROUP_RE, "", colnames(sig))
    miss <- setdiff(COV_GROUPS, smp_group)
    if (length(miss)) stop("No bigWig for cell line(s): ", paste(miss, collapse = ", "),
                           " (found: ", paste(unique(smp_group), collapse = ", "), ")")
    cov <- rbindlist(lapply(COV_GROUPS, function(g) {
      bins[, .(chr, start, end, group = g,
               value = log10(1 + rowMeans(sig[, smp_group == g, drop = FALSE])))]
    }))
    message("Coverage replicates per cell line: ",
            paste(sprintf("%s = %d", COV_GROUPS, as.integer(table(smp_group)[COV_GROUPS])), collapse = ", "))
    cov_cap <- cap_of(cov$value)                      # shared scale for both cell lines
    cov[, y := pmin(value, cov_cap)]
  }
}

# -----------------------------------------------------------------------------
# 7. Plot
# -----------------------------------------------------------------------------
cex_of <- function(pt) pt / 12          # base device font = 12 pt

gap_x <- function() {
  # x coordinate (in the chr1 cell) of the centre of the empty gap at 12 o'clock
  deg_per_x <- (CELL_META$cell.start.degree - CELL_META$cell.end.degree) / CELL_META$xrange
  CELL_META$cell.xlim[1] - (GAP_START / 2) / deg_per_x
}
ring_label <- function(txt, y = CELL_META$ycenter) {
  if (CELL_META$sector.index != CHROMS[1]) return(invisible())
  circos.text(gap_x(), y, txt, facing = "downward", adj = c(0.5, 0.5),
              cex = cex_of(PT_RING), col = INK, font = 1)
}

## legend in the empty centre of the circle (canvas centre = 0, 0), left-aligned
## block: "Transcript types" + one row per biotype (ref swatch, novel swatch,
## label), a blank line, "RNA coverage" + one row (HCC827 / H1975 swatches).
## Sized from the cap height; shrunk if it would leave the hole (r_max).
centre_key <- function(cex, r_max = Inf) {
  rows <- c(list(list(head = "Transcript types")),
            lapply(seq_along(KEY_LABELS), function(i)
              list(cols = c(ANNOT_COLS[[i]], NOVEL_COLS[[i]]), text = KEY_LABELS[i])))
  if (!is.null(cov))
    rows <- c(rows, list(list(gap = TRUE), list(head = "RNA coverage"),
                         list(cols = unname(COV_COLS[COV_GROUPS]),
                              text = paste(COV_GROUPS, collapse = " / "))))
  ## user units per inch in x and y differ (the device is not square)
  usr <- par("usr"); pin <- par("pin")
  x_per_y <- (diff(usr[1:2]) / pin[1]) / (diff(usr[3:4]) / pin[2])
  geom <- function(cex) {
    ch <- strheight("H", cex = cex)                  # cap height (y units)
    g <- list(ch = ch, pitch = 2.65 * ch, gap = 0.65 * 2.65 * ch,
              sw_w = 1.77 * ch * x_per_y, sw_h = 1.42 * ch,
              sw_sep = 0.68 * ch * x_per_y, txt_sep = 0.97 * ch * x_per_y)
    g$h <- sum(vapply(rows, function(r) if (isTRUE(r$gap)) g$gap else g$pitch, numeric(1)))
    g$w <- max(vapply(rows, function(r) {
      if (!is.null(r$head)) strwidth(r$head, cex = cex)
      else if (!is.null(r$text)) 2 * g$sw_w + g$sw_sep + g$txt_sep + strwidth(r$text, cex = cex)
      else 0
    }, numeric(1)))
    g
  }
  g <- geom(cex)
  half_diag <- function(g) sqrt((g$w / 2 / x_per_y)^2 + (g$h / 2)^2)   # in y units
  if (is.finite(r_max) && half_diag(g) > r_max) {
    cex <- cex * r_max / half_diag(g); g <- geom(cex)
    message(sprintf("centre legend shrunk to %.1f pt to fit the hole", cex * 12))
  }
  x0 <- -g$w / 2; y <- g$h / 2
  for (r in rows) {
    if (isTRUE(r$gap)) { y <- y - g$gap; next }
    yc <- y - g$pitch / 2
    if (!is.null(r$head)) {
      text(x0, yc, r$head, adj = c(0, 0.5), cex = cex, col = INK)
    } else {
      for (j in 1:2) {
        xs <- x0 + (j - 1) * (g$sw_w + g$sw_sep)
        rect(xs, yc - g$sw_h / 2, xs + g$sw_w, yc + g$sw_h / 2, col = r$cols[j], border = NA)
      }
      text(x0 + 2 * g$sw_w + g$sw_sep + g$txt_sep, yc, r$text, adj = c(0, 0.5), cex = cex, col = INK)
    }
    y <- y - g$pitch
  }
}

draw_circos <- function(fam = "sans") {
  par(family = fam, mar = c(0, 0, 0, 0))
  circos.clear()
  circos.par(
    start.degree   = 90 - GAP_START / 2,
    gap.after      = c(rep(1, length(CHROMS) - 1), GAP_START),
    cell.padding   = c(0, 0, 0, 0),
    track.margin   = c(0.003, 0.003),
    canvas.xlim    = CANVAS_X,
    canvas.ylim    = CANVAS_Y,
    points.overflow.warning = FALSE
  )
  circos.genomicInitialize(as.data.frame(chr_dt), plotType = NULL)

  # --- chromosome names -----------------------------------------------------
  circos.track(ylim = c(0, 1), track.height = H_CHR, bg.border = NA,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, 0.3, sub("^chr", "", CELL_META$sector.index),
                             cex = cex_of(PT_CHR), col = INK, facing = "bending.inside",
                             niceFacing = TRUE)
               })

  # --- ideogram -------------------------------------------------------------
  if (!is.null(cytoband)) {
    circos.genomicIdeogram(cytoband, track.height = H_IDEO)
  } else {
    circos.track(ylim = c(0, 1), track.height = H_IDEO, bg.col = "#d9d8d2", bg.border = INK)
  }

  # --- coverage rings (one per cell line, shared log scale) -----------------
  if (!is.null(cov)) {
    for (g in COV_GROUPS) {
      d <- cov[group == g, .(chr, start, end, y)]
      circos.genomicTrack(
        as.data.frame(d), ylim = c(0, cov_cap), track.height = H_COV,
        bg.col = RING_BG, bg.border = NA,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(region, value, ytop = value$y, ybottom = rep(0, nrow(region)),
                             col = COV_COLS[[g]], border = NA)
          if (g == COV_GROUPS[1]) ring_label("RNA cov.")           # one label for both rings
        })
    }
  }

  # --- one ring per strand: annotated + novel together ----------------------
  for (s in c("+", "-")) {
    d  <- as.data.frame(dens[strand == s & n > 0, .(chr, start, end, ybot, ytop, group)])
    nd <- as.data.frame(nov[strand == s, .(chr, pos, group)])
    circos.track(ylim = c(0, 1), track.height = H_STRAND, bg.col = STRAND_BG[[s]],
                 bg.border = NA, track.margin = c(0.01, 0.01), panel.fun = function(x, y) {
      dd <- d[d$chr == CELL_META$sector.index, , drop = FALSE]
      if (nrow(dd))                                   # annotated: light stacked bars
        circos.rect(dd$start, BAR_FRAC * dd$ybot / stack_cap, dd$end,
                    BAR_FRAC * dd$ytop / stack_cap, col = ANNOT_COLS[dd$group], border = NA)
      nn <- nd[nd$chr == CELL_META$sector.index, , drop = FALSE]
      if (nrow(nn))                                   # novel: ticks, outer band
        circos.segments(nn$pos, TICK_FROM, nn$pos, 1, col = NOVEL_COLS[nn$group],
                        lwd = NOVEL_LWD)
      ring_label(if (s == "+") "+  strand" else "-  strand")
    })
  }

  ## radius of the empty hole = inner edge of the last (- strand) ring
  r_hole <- get.cell.meta.data("cell.bottom.radius", sector.index = CHROMS[1],
                               track.index = get.current.track.index())
  circos.clear()
  centre_key(cex = cex_of(PT_KEY), r_max = KEY_FIT * r_hole)
  invisible()
}

# -----------------------------------------------------------------------------
# 8. Write outputs
# -----------------------------------------------------------------------------
# SVG: svglite (showtext off) keeps text editable; base svg() outlines it
open_svg <- function(file, w, h) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    fam <- setdiff(unique(c("sans", FONT)), "Source Sans 3")
    svglite::svglite(file, width = w, height = h, fix_text_size = FALSE,
                     system_fonts = setNames(rep(list("Source Sans 3"), length(fam)), fam))
  } else stop("svglite is needed for SVG with editable text (install.packages(\"svglite\")); base svg() would turn the text into glyph outlines")
}
cairo_ok <- capabilities("cairo")
# one device per file; showtext needs 72 dpi for vector output, 300 for PNG and
# is switched off for SVG (editable text); if the font fails, default sans
render <- function(file, w, h) {
  ext <- tolower(tools::file_ext(file))
  st  <- "showtext" %in% loadedNamespaces()
  if (st) {
    showtext::showtext_opts(dpi = if (ext == "png") 300 else 72)
    if (ext == "svg") showtext::showtext_auto(FALSE)
  }
  go <- function(fam) {
    if (ext == "pdf") {
      if (cairo_ok) cairo_pdf(file, width = w, height = h)
      else pdf(file, width = w, height = h, useDingbats = FALSE)
    } else if (ext == "png") {
      png(file, width = w, height = h, units = "in", res = 300,
          type = if (cairo_ok) "cairo" else getOption("bitmapType"))
    } else open_svg(file, w, h)
    on.exit(dev.off())
    tryCatch({ draw_circos(fam); TRUE }, error = function(e) {
      message("  ", basename(file), " with font '", fam, "' failed: ", conditionMessage(e)); FALSE })
  }
  if (!go(FONT)) stop("could not draw ", basename(file), " in ", FONT)
  if (st) {
    showtext::showtext_opts(dpi = 300)
    if (ext == "svg") showtext::showtext_auto(TRUE)
  }
  invisible(file)
}

for (ext in c("pdf", "png", "svg"))
  render(file.path(OUT_FIG, paste0(OUT_NAME, ".", ext)), FIG_IN[1], FIG_IN[2])

fwrite(dens[, .(chr, start, end, group, strand, n)],
       file.path(OUT_TAB, paste0(OUT_NAME, "_annotated_density_bins.tsv")), sep = "\t")
if (!is.null(cov))
  fwrite(cov[, .(chr, start, end, group, log10_1p_signal = value)],
         file.path(OUT_TAB, paste0(OUT_NAME, "_coverage_bins.tsv")), sep = "\t")
message("Done: ", file.path(OUT_FIG, OUT_NAME), ".pdf / .png / .svg; bin tables in ", OUT_TAB)

# -----------------------------------------------------------------------------
# Notes
# - "Other" with ANNOT_GROUP_BY = "transcript_biotype" absorbs nonsense_mediated_decay,
#   retained_intron and protein_coding_CDS_not_defined (~70k transcripts, mostly
#   from protein-coding genes). Set ANNOT_GROUP_BY = "gene_biotype" to classify by
#   the parent gene instead (then those go to "Protein-coding").
# - chrY (PAR copies; samples are female), chrIS (sequins), chrM and
#   unplaced/unlocalised scaffolds are not drawn.
# - Coverage: the bigWigs are unstranded (one file per sample), so the coverage
#   rings are not split by strand.
# =============================================================================
