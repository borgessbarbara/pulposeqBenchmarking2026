## ---------------------------------------------------------------------------
## 04_cis_circos.R - cis co-expression of novel transcripts (Dong)   (Fig. 2q, r)
##   Fig. 2r = paper_cis_circos.{pdf,svg,png}, Fig. 2q = paper_cis_heatmap.{pdf,svg,png},
##   both in the layout of the printed Figure 2 (A3 page); circos_cis_de_pairs.csv is the
##   cis-pair table read by 06_isoformic.R
##   Inputs : diff-exp/results/transcript/{logcpm.rds, de_*.csv} (01_de_edger.R; in this
##            logCPM from the RTA-divided oarfish counts),
##            <run>/transcriptome/curated/annotation/bambu_annotated_transcriptome.gtf,
##            <run>/transcriptome/curated/annotation/novel_transcripts_validated.gtf
##   Outputs: correlation-cis-1Mb/results/
## Journal panels, no bold text: Fig. 2r = the circle (6 pt) with a centred
##     title + two-column key below it (8 pt); Fig. 2q = the cis-pair heatmap
##     (pairs with an annotated protein-coding / lncRNA partner, numbered as on
##     the circle; names 8 pt, rho / locus 7 pt, colour bar top right).
##     Positions were measured on the printed A3 page. The variant with the key on
##     the right of the circle: write_v5_keyright <- TRUE (paper_cis_circos_keyright.*)
## Notes:
##   - fixed paths under dong/ (no working-directory search, no download)
##   - annotated coordinates from the Bambu annotated transcriptome GTF
##     (instead of GENCODE v50), read in chunks (the file is ~750 MB)
##   - chrY removed: female samples, so chrY features are PAR copies
##   - no "Novel unclassified" category (every novel transcript is classified)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(circlize)
  library(RColorBrewer)
})

## ------------------------------------------------------------- settings ---
source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")   # data_dir, de_dir, cis_dir
gtf_ref   <- run_path("ref_tx_gtf")
gtf_novel <- run_path("novel_gtf")
drop_chr  <- "chrY"      # female samples: chrY features are PAR copies

window_bp  <- 1e6        # max gap between novel transcript and neighbour
rho_min    <- 0.5        # keep pairs with |rho| > rho_min
exclude_same_gene <- FALSE
partners_from <- "all"   # "all" or "annotated"

## cis heatmap
only_de_novel <- FALSE
max_novel     <- 40
max_partners  <- 8

## circos
de_filter      <- "both"   # "both": novel AND neighbour DE; "novel": novel DE only
## pairs drawn (as in the Methods and 05_cis_lncRNA_top10.R): reference partner
## of another gene; ranked by |rho| (2 decimals), then by the novel transcript's
## DE FDR. With 3 + 3 samples |rho| saturates (many DE novel transcripts have a
## DE neighbour at |rho| = 1), so the FDR breaks the ties
drop_same_gene   <- TRUE
ref_partner_only <- TRUE     # partner = annotated (reference) transcript
one_per_novel  <- TRUE
max_pairs      <- 15
## optional variant (set by a wrapper such as 05_cis_lncRNA_top10.R before it
## sources this script): list(outdir, max_pairs, filter = function(sel) ..., hm_keep)
cis_variant <- get0("cis_variant", ifnotfound = list())
if (!is.null(cis_variant$max_pairs)) max_pairs <- cis_variant$max_pairs
lfc_cap        <- 8
lab_cex        <- 0.5      # size of the "pair N . distance" labels in the inner circle
## strand-track lines, tuned for print: novel lines thick and opaque on a
## lighter annotated background
lwd_novel_de   <- 1.8  # novel, DE
lwd_novel_ns   <- 1.0  # novel, not DE
lwd_annot      <- 0.3  # annotated
alpha_novel_ns <- 0.6  # novel, not DE
col_annot_de   <- "grey65"  # annotated, DE
col_annot_ns   <- "grey88"  # annotated, not DE
## layout / chrome (visual only)
gap_top    <- 14  # degrees at 12 o'clock for axis + strand names
strand_h   <- 0.14  # height of each strand track
track_bg   <- "#F4F4F2"  # strand-track background
ideo_cols  <- c("grey45", "grey75")  # alternating chromosomes
ink        <- "grey15"                # all text
zoom_gap   <- 3                       # degrees between zoomed pair sectors (inner circle)
zoom_min_w <- 0.02                    # min drawn width of a transcript body (fraction of window)
zoom_names <- FALSE                   # TRUE adds a 2nd label line: "novel / neighbour" names

grp_ref <- "TKI_sensitive"; grp_num <- "TKI_resistant"
ct_name <- paste0(grp_num, "_vs_", grp_ref)
line_lab <- c(TKI_sensitive = "HCC827", TKI_resistant = "H1975")

cat_levels <- c("Annotated protein-coding", "Annotated lncRNA", "Annotated other",
                "Novel protein-coding", "Novel lncRNA", "Novel non-coding")
## the biotype palette of every figure (as make_palette in 06_isoformic.R);
## "other" in light grey
cat_cols <- setNames(c("#A6D854", "#E78AC3", "#C9C9C9",
                       "#238B45", "#AE017E", "#6A51A3"), cat_levels)
nov_cols <- cat_cols[grepl("^Novel", names(cat_cols))]
de_cols  <- c(Up = "darkred", Down = "darkblue", NS = "grey85")
hm_anchor <- c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2", "#F2CBB7", "#EE8468", "#B40426")
hm_cols  <- colorRampPalette(hm_anchor, space = "Lab")(100)
rho_col  <- colorRamp2(seq(-1, 1, length.out = 7),
                       c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2",
                         "#F2CBB7", "#EE8468", "#B40426"))

## hg38 chromosome lengths
chr_len <- c(chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
             chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
             chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
             chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
             chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
             chr21 = 46709983, chr22 = 50818468, chrX = 156040895)   # chrY dropped

## -------------------------------------------------------------- paths -----
source(file.path(helper_dir, "id_helpers.R"))
source(file.path(helper_dir, "plot_theme.R"))
## ---- editable SVG text ------------------------------------------------------
## plot_theme.R switches showtext on, and showtext draws every letter as a glyph
## outline on EVERY device, SVG included. These helpers switch showtext off while
## an SVG is written, so svglite writes real <text> elements in the figure font:
## text boxes you can edit in Inkscape / Illustrator / Affinity / Vectornator.
##   svg_font        font name written into the SVG. It must be installed on the
##                   computer that edits the file (Google Fonts "Source Sans 3",
##                   static TTFs, "Install for all users"), or the editor swaps it.
##   fix_text_size   FALSE: no textLength attribute, so edited text is not
##                   stretched or squeezed to the old width.
## PNG/PDF are unchanged (showtext stays on for them).
if (!exists("svg_editable_helpers", inherits = TRUE)) {
  svg_editable_helpers <- TRUE
  svg_font <- get0("svg_font", ifnotfound = "Source Sans 3")
  ## R family names used in the plots -> the real font name
  svg_aliases <- function() {
    fams <- unique(c("sans", get0("font_family", ifnotfound = NULL), get0("theme_font", ifnotfound = NULL)))
    fams <- setdiff(fams, c("serif", "mono", "symbol", svg_font))
    stats::setNames(rep(list(svg_font), length(fams)), fams)
  }
  ## is the font known to systemfonts (the text measurer svglite uses)? checked once
  svg_font_ok <- local({
    ok <- NULL
    function() {
      if (!is.null(ok)) return(ok)
      have <- tryCatch(svg_font %in% systemfonts::system_fonts()$family, error = function(e) FALSE)
      if (!have && "require_font" %in% getNamespaceExports("systemfonts")) {
        ## systemfonts >= 1.2: fetch from Google Fonts for this R session (metrics only)
        have <- isTRUE(tryCatch({ systemfonts::require_font(svg_font); TRUE },
                                error = function(e) FALSE))
      }
      if (!have) message("  ! '", svg_font, "' is not installed: SVG text is still editable, ",
                         "but measured with a fallback font. Install it and restart R.")
      ok <<- have
      ok
    }
  })
  svg_device <- function(filename, width, height, bg = "white", ...) {
    svg_font_ok()
    svglite::svglite(filename, width = width, height = height, bg = bg,
                     system_fonts = svg_aliases(), fix_text_size = FALSE)
  }
  showtext_pause  <- function() {
    on <- "showtext" %in% loadedNamespaces()
    if (on) showtext::showtext_auto(FALSE)
    on
  }
  showtext_resume <- function(was_on) if (isTRUE(was_on)) showtext::showtext_auto(TRUE)
  ## ggsave() for SVG with editable text
  ggsave_svg <- function(filename, plot, width, height, ...) {
    st <- showtext_pause(); on.exit(showtext_resume(st))
    ggplot2::ggsave(filename, plot, width = width, height = height, device = svg_device, ...)
  }
  ## grid / base graphics: st <- svg_open(file, w, h); draw ...; svg_close(st)
  svg_open  <- function(filename, width, height) {
    st <- showtext_pause(); svg_device(filename, width, height); invisible(st)
  }
  svg_close <- function(st) { invisible(grDevices::dev.off()); showtext_resume(st) }
}
## ---------------------------------------------------------------------------
outdir <- cis_variant$outdir %||% cis_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------ functions ---
## read in chunks and keep only `keep` IDs, because
## the Bambu annotated GTF is ~750 MB
read_gtf_tx <- function(file, keep = NULL, feat = "transcript") {
  message("Reading ", basename(file), " ...")
  cb <- function(g, pos) {
    g <- g[g$X3 == feat, , drop = FALSE]
    if (!nrow(g)) return(NULL)
    attr <- function(key) str_match(g$X9, paste0(key, ' "([^"]+)"'))[, 2]
    out <- tibble(id = attr("transcript_id"), gene_id = attr("gene_id"),
                  chr = g$X1, start = g$X4, end = g$X5, strand = g$X7) %>%
      filter(!is.na(id))
    if (!is.null(keep)) out <- filter(out, id %in% keep)
    out
  }
  g <- read_tsv_chunked(file, DataFrameCallback$new(cb), chunk_size = 250000,
                        comment = "#", quote = "", col_names = FALSE, progress = FALSE,
                        col_types = cols_only(X1 = "c", X3 = "c", X4 = "d", X5 = "d",
                                              X7 = "c", X9 = "c"))
  if ((is.null(g) || !nrow(g)) && feat == "transcript")    # GTF with exon lines only
    return(read_gtf_tx(file, keep, feat = "exon"))
  g %>%
    group_by(id) %>%
    summarise(gene_id = first(gene_id), chr = first(chr), start = min(start),
              end = max(end), strand = first(strand), .groups = "drop")
}

std_rows  <- function(m) { m <- m - rowMeans(m); m / sqrt(rowSums(m^2)) }
rank_rows <- function(m) { r <- t(apply(m, 1, rank)); dimnames(r) <- dimnames(m); r }

## Cairo device; showtext at 72 dpi for PDF (else text ~4x too big); font fallback
cairo_ok <- capabilities("cairo") && isTRUE(tryCatch({
  f <- tempfile(fileext = ".png"); png(f, type = "cairo"); dev.off(); unlink(f); TRUE
}, warning = function(w) { try(dev.off(), silent = TRUE); FALSE },
   error = function(e) FALSE))
save_plot <- function(file, w, h, draw, res = 300) {
  is_pdf <- grepl("[.]pdf$", file)
  is_svg <- grepl("[.]svg$", file)
  if ("showtext" %in% loadedNamespaces()) {
    showtext::showtext_opts(dpi = if (is_pdf || is_svg) 72 else res)
    on.exit(showtext::showtext_opts(dpi = 300), add = TRUE)
  }
  ## draw into a temporary file, copy at the end: a destination that is open in a
  ## viewer or being synced by OneDrive can't break the drawing itself
  if (is_svg) { st_svg <- showtext_pause(); on.exit(showtext_resume(st_svg), add = TRUE) }
  tmp <- tempfile(fileext = if (is_pdf) ".pdf" else if (is_svg) ".svg" else ".png")
  open <- function() {
    if (is_svg) {
      ## svglite + showtext off = editable text (Inkscape / Illustrator); base svg() outlines it
      if (requireNamespace("svglite", quietly = TRUE)) svg_device(tmp, width = w, height = h)
      else stop("svglite is needed for SVG with editable text (install.packages(\"svglite\")); base svg() would turn the text into glyph outlines")
    } else if (is_pdf) {
      if (cairo_ok) cairo_pdf(tmp, width = w, height = h) else pdf(tmp, width = w, height = h)
    } else {
      type <- if (cairo_ok) "cairo" else if (capabilities("aqua")) "quartz" else getOption("bitmapType")
      png(tmp, width = w, height = h, units = "in", res = res, type = type)
    }
  }
  try_draw <- function(fam) {
    ## spacer device: circlize::circos.nested() does a dry run on pdf(NULL) and then
    ## returns to "the previous device", but jumps to the IDE plot pane when that pane
    ## sits right before ours (RStudio/Positron). A null device in between keeps it on ours.
    pdf(NULL); spacer <- dev.cur()
    open(); dev_out <- dev.cur()
    on.exit({
      if (dev_out %in% dev.list()) dev.off(dev_out)
      if (spacer  %in% dev.list()) dev.off(spacer)
    })
    tryCatch({ draw(fam); TRUE }, error = function(e) {
      message("  ", basename(file), " with font '", fam, "' failed: ", conditionMessage(e)); FALSE })
  }
  ok <- try_draw(font_family)                      # Source Sans 3 only, no Arial fallback
  if (ok) {
    if (!isTRUE(file.copy(tmp, file, overwrite = TRUE)))
      message("  could not write ", file, " (open in a viewer / OneDrive syncing?). ",
              "The figure is at: ", tmp)
  }
  invisible(file)
}

## ======================================================== CIS 1 Mb ========
rds <- readRDS(file.path(de_dir, "transcript", "logcpm.rds"))
ids <- setNames(as.character(rds$id), rownames(rds$norm))
xs  <- log2(rds$norm + 1)
xs  <- xs[apply(xs, 1, sd) > 0, , drop = FALSE]
Sa  <- std_rows(rank_rows(xs))

de <- read_csv(file.path(de_dir, "transcript",
                         paste0("de_", ct_name, ".csv")), show_col_types = FALSE)
info <- tibble(feature = rownames(xs), id = unname(ids[rownames(xs)])) %>%
  left_join(dplyr::select(de, id, category, call, padj, log2FC), by = "id") %>%
  mutate(novel = is_novel(id), call = replace_na(call, "NS"),
         category = factor(category, levels = cat_levels))

coords <- bind_rows(read_gtf_tx(gtf_novel, keep = info$id[info$novel]),
                    read_gtf_tx(gtf_ref,   keep = info$id[!info$novel])) %>%
  distinct(id, .keep_all = TRUE)
info <- left_join(info, coords, by = "id")
miss <- filter(info, is.na(chr))
write_csv(dplyr::select(miss, feature, id, category), file.path(outdir, "missing_coordinates.csv"))
message(sprintf("%d / %d tested transcripts with coordinates (%d novel missing, %d annotated missing)",
                sum(!is.na(info$chr)), nrow(info), sum(miss$novel), sum(!miss$novel)))
message(sum(info$chr %in% drop_chr, na.rm = TRUE), " tested transcripts on ",
        paste(drop_chr, collapse = ", "), " removed")
info <- filter(info, !is.na(chr), !chr %in% drop_chr)

nov_tab <- filter(info, novel)
by_chr  <- split(info, info$chr)
pairs <- map_dfr(seq_len(nrow(nov_tab)), function(i) {
  n  <- nov_tab[i, ]
  cc <- by_chr[[n$chr]]
  gap <- pmax(0, pmax(cc$start, n$start) - pmin(cc$end, n$end))
  hit <- gap <= window_bp & cc$feature != n$feature
  if (!any(hit)) return(NULL)
  tibble(novel = n$feature, partner = cc$feature[hit], distance_bp = gap[hit])
})
if (!nrow(pairs)) stop("No neighbour within ", window_bp, " bp of any novel transcript.")

pairs <- pairs %>%
  mutate(rho = rowSums(Sa[novel, , drop = FALSE] * Sa[partner, , drop = FALSE])) %>%
  left_join(dplyr::select(info, novel = feature, novel_chr = chr, novel_start = start,
                          novel_end = end, novel_strand = strand, novel_gene = gene_id,
                          novel_category = category, novel_call = call, novel_padj = padj),
            by = "novel") %>%
  left_join(dplyr::select(info, partner = feature, partner_start = start, partner_end = end,
                          partner_strand = strand, partner_gene = gene_id,
                          partner_novel = novel, partner_category = category,
                          partner_call = call), by = "partner") %>%
  mutate(same_gene = novel_gene == partner_gene,
         same_strand = novel_strand == partner_strand) %>%
  arrange(novel_chr, novel_start, desc(abs(rho)))
write_csv(pairs, file.path(outdir, "cis_pairs_all.csv"))

cor_pairs <- pairs %>%
  filter(abs(rho) > rho_min, !(exclude_same_gene & same_gene),
         partners_from == "all" | !partner_novel)
write_csv(cor_pairs, file.path(outdir, "cis_pairs_correlated.csv"))
message(sprintf("%d pairs within %.1f Mb; %d with |rho| > %.2f, involving %d novel transcripts",
                nrow(pairs), window_bp / 1e6, nrow(cor_pairs), rho_min,
                n_distinct(cor_pairs$novel)))

## ---- cis heatmap: novel (rows) x correlated neighbours (cols), genome order
rows <- cor_pairs %>%
  { if (only_de_novel) filter(., novel_call %in% c("Up", "Down")) else . } %>%
  group_by(novel, novel_call, novel_padj) %>%
  summarise(best = max(abs(rho)), .groups = "drop") %>%
  arrange(!novel_call %in% c("Up", "Down"), novel_padj, desc(best)) %>%
  slice_head(n = max_novel) %>%
  pull(novel)
hp <- cor_pairs %>%
  filter(novel %in% rows) %>%
  group_by(novel) %>%
  slice_max(abs(rho), n = max_partners, with_ties = FALSE) %>%
  ungroup()

if (length(rows) >= 1 && nrow(hp) >= 2) {
  chr_ord <- function(x) {
    k <- sub("^chr", "", x); num <- suppressWarnings(as.numeric(k))
    ifelse(!is.na(num), num, match(k, c("X", "Y", "M", "MT"), nomatch = 99) + 100)
  }
  pos <- dplyr::select(info, feature, chr, start)
  ord <- function(f) { p <- pos[match(f, pos$feature), ]; f[order(chr_ord(p$chr), p$chr, p$start)] }
  r_ord <- ord(unique(hp$novel)); c_ord <- ord(unique(hp$partner))
  m <- matrix(NA_real_, length(r_ord), length(c_ord), dimnames = list(r_ord, c_ord))
  m[cbind(match(hp$novel, r_ord), match(hp$partner, c_ord))] <- hp$rho

  chr_of <- function(f) info$chr[match(f, info$feature)]
  gaps   <- function(f) which(head(chr_of(f), -1) != tail(chr_of(f), -1))
  ann <- function(f) data.frame(
    category = droplevels(info$category[match(f, info$feature)]),
    DE = factor(info$call[match(f, info$feature)], levels = names(de_cols)),
    row.names = f)
  ar <- ann(r_ord); ac <- ann(c_ord)
  acol <- list(category = cat_cols[union(levels(ar$category), levels(ac$category))],
               DE = de_cols)
  lab <- function(f) paste0(f, "  (", chr_of(f), ")")
  args <- list(mat = m, color = hm_cols, breaks = seq(-1, 1, length.out = 101),
               na_col = "white", cluster_rows = FALSE, cluster_cols = FALSE,
               gaps_row = gaps(r_ord), gaps_col = gaps(c_ord),
               labels_row = lab(r_ord), labels_col = lab(c_ord),
               annotation_row = ar, annotation_col = ac, annotation_colors = acol,
               border_color = "grey90", fontsize = 8, fontsize_row = 7,
               fontsize_col = 6.5, angle_col = 90,
               main = paste0("Novel transcripts vs correlated neighbours (\u2264 ",
                             window_bp / 1e6, " Mb, |\u03c1| > ", rho_min, ")\n",
                             "Spearman, n = ", ncol(xs),
                             " samples; blank = not a local correlated pair"))
  for (ext in c("png", "pdf"))
    save_plot(file.path(outdir, paste0("heatmap_cis_1mb.", ext)),
              4.5 + 0.14 * ncol(m), 3.5 + 0.15 * nrow(m),
              function(fam) {
                ## pheatmap measures text on its own pdf(NULL) and, after closing it,
                ## draws on whatever device R picks next: save_plot's spacer, not our file
                ## (empty PDF, no PNG). Build the layout silently, then draw it on our device.
                dev_draw <- dev.cur()
                g <- do.call(pheatmap, c(args, list(fontfamily = fam, filename = NA, silent = TRUE)))
                dev.set(dev_draw)
                grid::grid.newpage(); grid::grid.draw(g$gtable)
              })
  write_csv(as_tibble(m, rownames = "novel"), file.path(outdir, "heatmap_cis_1mb.csv"))
}

## ========================================================== CIRCOS ========
## every novel transcript with coordinates, for the strand tracks
tx_all <- info %>%
  filter(chr %in% names(chr_len), !is.na(log2FC), strand %in% c("+", "-")) %>%
  mutate(cat = as.character(category),
         de  = call %in% c("Up", "Down"),
         x   = (start + end) / 2,
         y   = pmax(pmin(log2FC, lfc_cap), -lfc_cap),
         col = case_when(novel &  de ~ unname(nov_cols[cat]),
                         novel & !de ~ adjustcolor(unname(nov_cols[cat]), alpha.f = alpha_novel_ns),
                         !novel & de ~ col_annot_de,
                         TRUE        ~ col_annot_ns),
         lwd   = case_when(novel & de ~ lwd_novel_de, novel ~ lwd_novel_ns, TRUE ~ lwd_annot),
         layer = 2 * novel + de) %>%
  filter(!is.na(col)) %>%
  arrange(layer)
## pairs to show
is_de <- function(x) x %in% c("Up", "Down")
sel <- cor_pairs %>%
  filter(novel_chr %in% names(chr_len), is_de(novel_call),
         de_filter != "both" | is_de(partner_call),
         !(drop_same_gene & same_gene),
         !(ref_partner_only & partner_novel)) %>%
  arrange(desc(round(abs(rho), 2)), novel_padj)   # ties in |rho|: smallest FDR first
if (!is.null(cis_variant$filter)) sel <- cis_variant$filter(sel)   # variant: extra pair filter
if (one_per_novel) sel <- distinct(sel, novel, .keep_all = TRUE)
## drop reciprocal duplicates (A-B and B-A both kept when each is the other's best
## partner); keeps the copy with the larger |rho|, before max_pairs is applied
sel <- sel %>%
  mutate(pair_key = paste(pmin(novel, partner), pmax(novel, partner))) %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  dplyr::select(-pair_key)
if (!nrow(sel)) stop("No correlated pair passes the DE filter - try de_filter <- \"novel\"",
                     call. = FALSE)
sel <- sel %>%
  slice_head(n = max_pairs) %>%
  arrange(match(novel_chr, names(chr_len)), novel_start) %>%
  mutate(mid_n = (novel_start + novel_end) / 2,
         pair_lab = paste0(novel, " \u2013 ", partner))
write_csv(sel, file.path(outdir, "circos_cis_de_pairs.csv"))
message(nrow(sel), " pairs in the circos")

## pair heatmap: z-score, novel row then neighbour row, blank row between pairs
grp  <- setNames(as.character(rds$group), colnames(rds$norm))
smp  <- c(names(grp)[grp == grp_ref], names(grp)[grp == grp_num])
z    <- t(scale(t(log2(rds$norm[unique(c(sel$novel, sel$partner)), smp, drop = FALSE] + 1))))
hrows <- list(); labs <- character(0)
for (i in seq_len(nrow(sel))) {
  hrows <- c(hrows, list(z[sel$novel[i], ], z[sel$partner[i], ], rep(NA, length(smp))))
  labs  <- c(labs, sel$novel[i], sel$partner[i], "")
}
hm <- do.call(rbind, head(hrows, -1)); labs <- head(labs, -1)
z_lim <- 2

## ---------------------------------------------------------------- figure ---
## Text sizes are in cex units of a 12 pt device font (cex 1 = 12 pt).
## FS_FULL  = the combined screen figure (14 x 8.5 in).
## FS_POSTER = poster panels at their real print size (cm below); sizes in pt / 12
## so labels print at 11-16 pt, in line with the poster's 15-25 pt text.
FS_FULL <- list(chr = 0.6, tick = 0.38, strand = 0.5, zoom = lab_cex, sign = 0.65,
                key = 0.62, note = 0.58, hm = 0.58, hm_num = 0.75, hm_small = 0.5,
                hm_group = 0.65, hm_title = 0.75, chip = 0.58, cbar = 0.55, cbar_note = 0.5,
                gap_top = gap_top, zoom_dist = TRUE, notes = TRUE, canvas = 1)
FS_POSTER <- list(chr = 12, tick = 10, strand = 12, zoom = 12, sign = 14,
                  key = 12, note = 12, hm = 11, hm_num = 13, hm_small = 10,
                  hm_group = 13, hm_title = 14, chip = 12, cbar = 11, cbar_note = 11)
FS_POSTER <- c(lapply(FS_POSTER, function(pt) pt / 12),
               list(gap_top = 22, zoom_dist = FALSE, notes = FALSE,
                    canvas = 1.12,    # room for chromosome labels inside a tight panel
                    centre_key = TRUE, ckey = 10 / 12))  # DE key inside the circle
FS_FULL$centre_key <- FALSE
## FS_PAPER = journal panel (Fig. 2) at its print size, PAPER_CIRCOS_MM: every
## text 6 pt as the other panels (5 pt axis ticks), lines thinner for the
## smaller circle, key on the right so the centre stays blank
FS_PAPER <- c(lapply(list(chr = 6, tick = 5, strand = 6, zoom = 6, sign = 6, key = 6, note = 6,
                          hm = 6, hm_num = 6, hm_small = 6, hm_group = 6, hm_title = 7,
                          chip = 6, cbar = 6, cbar_note = 6), function(pt) pt / 12),
              list(gap_top = 22,                 # as the poster: strand names + log2FC axis
                   zoom_dist = TRUE,  # "pair N . distance"
                   notes = FALSE, canvas = 1.08,
                   strand_h = 0.09,                      # thinner strand tracks -> larger inner
                                                         # circle, so the pair labels fit at 6 pt
                   centre_key = FALSE, key_side = "right", circ_frac = 0.78, key_frac = 0.12,
                   lwd_k = 0.7))
## full-width Fig. 2 row: 171.5 mm (174 mm box - 2.5 mm letter) x 133 mm. The
## circle (~123 mm) takes the left 78 %, the key the right 22 % (~38 mm).
PAPER_CIRCOS_MM <- c(171.5, 133)
POSTER_CIRCOS_CM  <- c(19, 19)  # square (key inside the circle)
LINK_H_RATIO <- 0.8              # depth of the pair links (circlize h.ratio)
KEY_FIT      <- 0.94             # centre key must fit inside this fraction of the free hole
KEY_MIN_PT   <- 7                # message if fitting pushes the key below this size
POSTER_HEATMAP_CM <- c(14, 24)   # width, height of the poster heatmap panel
cm2in <- function(x) x / 2.54

## ---- Fig. 2r / 2q as printed (Figure 2, A3 page) -----------------------------
## Every position below was measured on the printed page (mm). Offsets for 2r
## are from the circle centre; 2q positions are from the panel's top-left corner.
##
## Fig. 2r: the journal circle at exactly the same size. circlize draws the canvas
## (+/- canvas) with plot()'s 4 % padding per side, so +/- canvas spans
## 2 * canvas * 1.08 units; the journal circle is 133 mm for canvas 1.08, i.e. 57.0 mm per unit.
## The circle sits in a W x W square at the top of the panel; title and key below.
CIRCOS_FIG_MM    <- c(123, 163.5)                         # width, height (mm)
MM_PER_UNIT <- PAPER_CIRCOS_MM[2] / (2 * 1.08 * 1.08)  # journal circle scale
FS_CIRCOS_FIG <- modifyList(FS_PAPER, list(
  key_side = "paper",
  canvas   = CIRCOS_FIG_MM[1] / (2 * 1.08 * MM_PER_UNIT),  # same mm per unit as the journal circle
  key_frac = 1 - CIRCOS_FIG_MM[1] / CIRCOS_FIG_MM[2],          # bottom share of the page = key
  key      = 8 / 12,                                 # title, subtitle and key: 8 pt
  dev_mm   = CIRCOS_FIG_MM))
CIRCOS_FIG_KEY <- list(
  title_dx = -0.4, title_dy = 67.1, sub_dy = 70.8,   # title / subtitle centres
  row1_dy  = 78.4, row_h = 3.45,                      # 7 key rows
  x_lab    = -35.6,                                   # left column labels (left edge)
  x_de     = -16.65, x_ns = -8.1, sw_len = 6.05,      # DE / not-DE swatch centres, length
  de_col   = c("grey20", "grey20"), ns_col = c("grey60", "grey65"),   # annotated, novel
  de_lwd   = c(0.6, 2.8),         ns_lwd = c(0.6, 1.6),
  link_len = 3.7, link_txt = 5.4, link_lty_neg = "64", # link swatch, text start, dashes
  x_box    = 5.3, box_w = 3.4, box_h = 2.1, x_box_txt = 10.4)   # right column (colours)

## Fig. 2q: pairs whose partner is annotated protein-coding or lncRNA (as in
## 06_isoformic.R); sample columns HCC827_1-3 | H1975_1-3; z-score of
## log2(normalised count + 1) per row, colour limits = +/- the largest |z| shown
HM_FIG_KEEP <- c("Annotated protein-coding", "Annotated lncRNA")
if (!is.null(cis_variant$hm_keep)) HM_FIG_KEEP <- cis_variant$hm_keep
HM_FIG <- list(
  width   = 82.7,                       # mm; the height follows the number of pairs
  y_tile0 = 3.66,                       # top of the first tile row
  row_dy  = 3.52, pair_dy = 9.146, tile_h = 3.2,     # novel -> partner row, pair -> pair
  x_tile1 = 32.82, x_tile2 = 41.96, tile_dx = 2.94, tile_w = 2.78,   # HCC827 | H1975 tiles
  x_cat   = 29.94, cat_w = 2.14,        # biotype square
  x_name  = 28.45,                      # right edge of the transcript names
  x_num   = 3.15,                       # right edge of the pair numbers
  x_side  = 52.95,                      # left edge of rho / locus
  y_head  = 1.54,                       # centre of HCC827 / H1975
  bar_x   = 72.25, bar_w = 2.73, bar_top = 13.75, bar_h = 18.9,   # colour bar
  x_bar_title = 71.7, y_bar_title = c(6.15, 9.55), x_bar_lab = 76.5,
  pt_name = 8, pt_side = 7, pt_head = 8, pt_bar = 8)

## part = "both" (combined figure), "circos" or "heatmap" (poster panels).
## Pairs are numbered 1..n in genome order: the same number marks the zoomed
## pair on the circle and its rows in the heatmap.
draw <- function(fam, part = "both", fs = FS_FULL) {
  par(family = fam, mar = c(0, 0, 0, 0), xpd = NA)
  pid  <- seq_len(nrow(sel))
  dist_txt <- ifelse(sel$distance_bp == 0, "overlap", sprintf("%.0f kb", sel$distance_bp / 1e3))
  chr1 <- names(chr_len)[1]
  chrN <- names(chr_len)[length(chr_len)]
  new_panel <- FALSE                       # TRUE once something is on the page

  ## ---- circos -------------------------------------------------------------
  if (part %in% c("both", "circos")) {
    ## Outer circle = genome (ideogram + strand tracks). Inner circle = one zoomed
    ## sector per pair (circlize::circos.nested), because 1 Mb is ~0.1 degree on
    ## the genome circle and links there would be invisible.
    key_right <- identical(fs$key_side, "right")
    key_below <- identical(fs$key_side, "below")
    key_paper <- identical(fs$key_side, "paper")  # Fig. 2r as printed
    lk    <- if (is.null(fs$lwd_k)) 1 else fs$lwd_k        # line-width factor
    key_h <- if (part == "both") 0.12 else if (isTRUE(fs$centre_key) || key_right || key_below || key_paper) 0 else 0.16
    r_in  <- NA                    # inner radius of the zoom transcripts (set in f_inner)
    par(fig = if (part == "both") c(0, 0.58, key_h, 1) else if (key_right) c(0, fs$circ_frac, 0, 1)
              else if (key_below || key_paper) c(0, 1, fs$key_frac, 1) else c(0, 1, key_h, 1))
    on.exit(circos.clear(), add = TRUE)
    circos.clear()

    ## zoom window per pair: novel + neighbour span, padded 5 % (min 2 kb)
    win_lo <- pmin(sel$novel_start, sel$partner_start)
    win_hi <- pmax(sel$novel_end, sel$partner_end)
    win_pd <- pmax(0.05 * (win_hi - win_lo), 2000)
    win_s  <- pmax(0, win_lo - win_pd)
    win_e  <- win_hi + win_pd
    zoom_sec <- paste0("p", pid)
    corr <- data.frame(chr1 = sel$novel_chr, start1 = win_s, end1 = win_e,
                       chr2 = zoom_sec, start2 = win_s, end2 = win_e)

    f_outer <- function() {
      circos.par(start.degree = 90, gap.after = c(rep(1, length(chr_len) - 1), fs$gap_top),
                 cell.padding = c(0, 0, 0, 0), track.margin = c(0.005, 0.005),
                 canvas.xlim = c(-1, 1) * fs$canvas, canvas.ylim = c(-1, 1) * fs$canvas)
      circos.genomicInitialize(data.frame(chr = names(chr_len), start = 0, end = chr_len),
                               plotType = NULL)

      ## ideogram: thin, alternating mid-greys, names outside
      circos.track(ylim = c(0, 1), track.height = 0.03, bg.border = NA,
                   panel.fun = function(x, y) {
        i <- CELL_META$sector.numeric.index
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, border = NA,
                    col = ideo_cols[2 - i %% 2])
        ## "downward" = always upright and horizontal (no flipped label at 9 o'clock)
        circos.text(CELL_META$xcenter, 1 + mm_y(3 + 5 * (fs$chr - 0.6)), sub("^chr", "", CELL_META$sector.index),
                    cex = fs$chr, col = ink, facing = "downward", adj = c(0.5, 0.5))
      })

      ## strand tracks: +/-cap/2 gridlines, a log2FC axis (on chrX, in the gap)
      ## and the strand name (beside chr1)
      for (s in c("+", "-")) {
        d <- filter(tx_all, strand == s)
        circos.track(ylim = c(-lfc_cap, lfc_cap), track.height = if (is.null(fs$strand_h)) strand_h else fs$strand_h,
                     bg.border = NA,
                     bg.col = track_bg, panel.fun = function(x, y) {
          for (h in c(-1, 1) * lfc_cap / 2)
            circos.lines(CELL_META$xlim, c(h, h), col = "white", lwd = 0.8)
          dd <- d[d$chr == CELL_META$sector.index, ]
          if (nrow(dd)) circos.segments(dd$x, 0, dd$x, dd$y, col = dd$col, lwd = dd$lwd * lk)
          circos.lines(CELL_META$xlim, c(0, 0), col = "grey55", lwd = 0.4)
          if (CELL_META$sector.index == chrN)
            circos.yaxis(side = "right", at = c(-lfc_cap, 0, lfc_cap) / 2,
                         labels = c(-lfc_cap / 2, 0, paste0("+", lfc_cap / 2)),
                         labels.cex = fs$tick, col = "grey45", labels.col = "grey30", lwd = 0.4)
          if (CELL_META$sector.index == chr1) {
            deg_per_x <- (CELL_META$cell.start.degree - CELL_META$cell.end.degree) / CELL_META$xrange
            circos.text(CELL_META$xlim[1] - 0.8 / deg_per_x, 0, paste(s, "strand"),
                        facing = "downward", adj = c(1, 0.5), cex = fs$strand, col = ink)
          }
        })
      }
      ## no circos.clear() here: circos.nested() handles it
    }

    f_inner <- function() {
      ## no start.degree / canvas here: circos.nested() sets them
      circos.par(gap.after = zoom_gap, cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0.01))
      circos.genomicInitialize(data.frame(chr = factor(zoom_sec, levels = zoom_sec),
                                          start = win_s, end = win_e),
                               plotType = NULL, sector.width = rep(1, length(zoom_sec)))

      ## sector label along the arc: "pair N . distance" (poster: "pair N"; the
      ## distance is in the heatmap) and, optionally, the two names
      circos.track(ylim = c(0, 1), track.height = if (zoom_names) 0.11 else 0.07,
                   bg.border = NA, panel.fun = function(x, y) {
        i <- CELL_META$sector.numeric.index
        circos.text(CELL_META$xcenter, if (zoom_names) 0.68 else 0.35,
                    if (fs$zoom_dist) sprintf("pair %d \u00b7 %s", i, dist_txt[i]) else sprintf("pair %d", i),
                    cex = fs$zoom, font = 1, col = ink,
                    facing = "bending.inside", niceFacing = TRUE)
        if (zoom_names)
          circos.text(CELL_META$xcenter, 0.2, paste(sel$novel[i], "/", sel$partner[i]),
                      cex = fs$zoom * 0.8, col = "grey35",
                      facing = "bending.inside", niceFacing = TRUE)
      })

      ## transcript bodies: outer lane = novel, inner lane = neighbour, colour = category
      circos.track(ylim = c(0, 1), track.height = 0.16, bg.col = track_bg, bg.border = NA,
                   panel.fun = function(x, y) {
        i    <- CELL_META$sector.numeric.index
        wmin <- zoom_min_w * CELL_META$xrange          # keep short transcripts visible
        body <- function(s, e) { m <- (s + e) / 2; h <- max(e - s, wmin) / 2; c(m - h, m + h) }
        ## transcript as a bar + its strand sign (+ / -) on the roomier side of the bar
        tx_bar <- function(s, e, strand, y, col) {
          b <- body(s, e)
          circos.rect(b[1], y - 0.17, b[2], y + 0.17, col = col, border = NA)
          right <- mean(b) < CELL_META$xcenter
          circos.text(if (right) b[2] + mm_x(1.2 * fs$sign / 0.65) else b[1] - mm_x(1.2 * fs$sign / 0.65), y,
                      if (identical(strand, "-")) "\u2212" else "+",
                      facing = "downward", adj = c(0.5, 0.5), cex = fs$sign, font = 1, col = ink)
        }
        tx_bar(sel$novel_start[i], sel$novel_end[i], sel$novel_strand[i], 0.73,
               cat_cols[as.character(sel$novel_category[i])])
        tx_bar(sel$partner_start[i], sel$partner_end[i], sel$partner_strand[i], 0.27,
               cat_cols[as.character(sel$partner_category[i])])
      })
      ## inner edge of the transcript track = where the links start
      r_in <<- get.cell.meta.data("cell.bottom.radius", sector.index = zoom_sec[1],
                                  track.index = get.current.track.index())

      ## cis links (innermost): solid = positive rho, dashed = negative rho
      for (i in pid) {
        m1 <- (sel$novel_start[i] + sel$novel_end[i]) / 2
        m2 <- (sel$partner_start[i] + sel$partner_end[i]) / 2
        if (m1 != m2)
          circos.link(zoom_sec[i], m1, zoom_sec[i], m2, h.ratio = LINK_H_RATIO,
                      col = ink, lwd = 1.3 * lk, lty = if (sel$rho[i] < 0) 2 else 1)
      }
    }

    circos.nested(f_outer, f_inner, corr, connection_height = mm_h(6),
                  connection_col = adjustcolor("grey55", alpha.f = 0.12),
                  connection_border = adjustcolor("grey35", alpha.f = 0.5),
                  connection_lwd = 0.3)

    k_de  <- c(col_annot_de, unname(nov_cols))
    k_ns  <- c(col_annot_ns, adjustcolor(unname(nov_cols), alpha.f = alpha_novel_ns))
    k_lwd_de <- 2 * c(lwd_annot, rep(lwd_novel_de, length(nov_cols)))  # x2 so thin lines show
    k_lwd_ns <- 2 * c(lwd_annot, rep(lwd_novel_ns, length(nov_cols)))

    ## ---- poster: DE / not-DE key in the empty centre of the zoom circle --------
    ## (coordinates here = the inner circle's canvas, centre = 0, 0). The links dip
    ## at most rou * (1 - h.ratio * cos(beta / 2)) inwards (circlize getQuadraticPoints),
    ## so everything inside rou * h.ratio * cos(beta / 2) is free.
    if (isTRUE(fs$centre_key)) {
      beta   <- (360 - length(zoom_sec) * zoom_gap) / length(zoom_sec)   # sector width, deg
      r_free <- KEY_FIT * r_in * LINK_H_RATIO * cos(beta / 2 * pi / 180)
      rows <- c("annotated", "novel", "protein-coding", "lncRNA", "non-coding")
      idx  <- c(1, NA, 2, 3, 4)                  # which k_* colour each row uses
      box <- function(cex) {
        sw <- max(strwidth(c("not DE", "MMM"), cex = cex)); gap <- strwidth("m", cex = cex) * 0.8
        lh <- strheight("Mg", cex = cex) * 1.6
        w  <- max(strwidth(rows, cex = cex)) + 2 * (gap + sw)
        list(sw = sw, gap = gap, lh = lh, w = w, h = (length(rows) + 1) * lh,
             r = sqrt((w / 2)^2 + (((length(rows) + 1) * lh) / 2)^2))
      }
      cex <- fs$ckey
      if (box(cex)$r > r_free) {
        cex <- cex * r_free / box(cex)$r
        if (cex * 12 < KEY_MIN_PT)
          message(sprintf("cis centre key shrunk to %.1f pt to fit (lower LINK_H_RATIO to enlarge the hole)", cex * 12))
      }
      b  <- box(cex)
      x0 <- -b$w / 2                               # left edge of the labels
      xs <- x0 + max(strwidth(rows, cex = cex)) + b$gap + b$sw / 2   # DE swatch centre
      xn <- xs + b$sw + b$gap                      # not-DE swatch centre
      ys <- b$h / 2 - b$lh * (seq_len(length(rows) + 1) - 0.5)       # header + rows
      text(c(xs, xn), ys[1], c("DE", "not DE"), cex = cex, font = 1, col = ink)
      for (r in seq_along(rows)) {
        y <- ys[r + 1]
        is_sub <- r > 2                            # novel biotypes, indented
        text(x0 + if (is_sub) strwidth("  ", cex = cex) else 0, y, rows[r], adj = c(0, 0.5),
             cex = cex, col = ink, font = if (rows[r] == "novel") 3 else 1)
        k <- idx[r]
        if (!is.na(k)) {
          segments(xs - b$sw / 2, y, xs + b$sw / 2, y, col = k_de[k], lwd = k_lwd_de[k], lend = 1)
          segments(xn - b$sw / 2, y, xn + b$sw / 2, y, col = k_ns[k], lwd = k_lwd_ns[k], lend = 1)
        }
      }
    }

    ## ---- journal: key on the right (centre blank) ------------------------------
    ## block 1 = strand tracks (DE / not DE); block 2 = inner circle: the bars are
    ## coloured by the six transcript categories and the links by the sign of rho
    if (key_right) {
      par(fig = c(fs$circ_frac, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
      plot.new(); plot.window(c(0, 1), c(0, 1))
      cex <- fs$key
      lh  <- strheight("Mg", cex = cex) * 1.75           # line height
      sw  <- 0.12                                        # swatch width
      k_lab <- c("annotated", sub("^Novel ", "novel ", names(nov_cols)))
      x0  <- 0.02
      ## swatch columns placed after the longest label
      xs  <- max(0.60, x0 + max(strwidth(k_lab, cex = cex)) + 0.03 + sw / 2)
      xn  <- min(xs + sw + 0.05, 1 - sw / 2)
      y   <- 0.5 + lh * 9                                # block centred on the circle
      text(x0, y, "Strand tracks", adj = c(0, 0.5), cex = cex, col = ink); y <- y - lh
      text(x0, y, "line = log2FC H1975 / HCC827", adj = c(0, 0.5), cex = cex, col = "grey35"); y <- y - lh
      text(c(xs, xn), y, c("DE", "not DE"), cex = cex, col = ink); y <- y - lh
      for (r in seq_along(k_lab)) {
        text(x0, y, k_lab[r], adj = c(0, 0.5), cex = cex, col = ink)
        segments(xs - sw / 2, y, xs + sw / 2, y, col = k_de[r], lwd = k_lwd_de[r] * lk, lend = 1)
        segments(xn - sw / 2, y, xn + sw / 2, y, col = k_ns[r], lwd = k_lwd_ns[r] * lk, lend = 1)
        y <- y - lh
      }
      y <- y - 0.6 * lh
      text(x0, y, "Cis-acting pairs (inner circle)", adj = c(0, 0.5), cex = cex, col = ink); y <- y - lh

      c_lab <- sub("^Novel ", "novel ", sub("^Annotated ", "annotated ", names(cat_cols)))
      for (r in seq_along(cat_cols)) {
        rect(x0, y - 0.32 * lh, x0 + 0.07, y + 0.32 * lh, col = cat_cols[r], border = NA)
        text(x0 + 0.10, y, c_lab[r], adj = c(0, 0.5), cex = cex, col = ink); y <- y - lh
      }
      for (neg in c(FALSE, TRUE)) {
        segments(x0, y, x0 + 0.07, y, col = ink, lwd = 1.3 * lk, lty = if (neg) 2 else 1)
        text(x0 + 0.10, y, if (neg) "link: \u03c1 < 0" else "link: \u03c1 > 0", adj = c(0, 0.5),
             cex = cex, col = ink)
        y <- y - lh
      }
      text(x0, y, "+ / \u2212 = transcript strand", adj = c(0, 0.5), cex = cex, col = "grey35")
    }

    ## ---- journal, full width: key below the circle, three columns ------------
    if (key_below) {
      par(fig = c(0, 1, 0, fs$key_frac), new = TRUE, mar = c(0, 0, 0, 0))
      plot.new(); plot.window(c(0, 1), c(0, 1))
      cex <- fs$key
      lh  <- strheight("Mg", cex = cex) * 1.6
      top <- 1 - 0.6 * lh
      ## column 1: strand tracks
      x0 <- 0.03; xs <- 0.23; xn <- 0.30; sw <- 0.045
      y  <- top
      text(x0, y, "Strand tracks: line = log2FC H1975 / HCC827", adj = c(0, 0.5), cex = cex, col = ink)
      y <- y - lh
      text(c(xs, xn), y, c("DE", "not DE"), cex = cex, col = ink); y <- y - lh
      k_lab <- c("annotated", sub("^Novel ", "novel ", names(nov_cols)))
      for (r in seq_along(k_lab)) {
        text(x0, y, k_lab[r], adj = c(0, 0.5), cex = cex, col = ink)
        segments(xs - sw / 2, y, xs + sw / 2, y, col = k_de[r], lwd = k_lwd_de[r] * lk, lend = 1)
        segments(xn - sw / 2, y, xn + sw / 2, y, col = k_ns[r], lwd = k_lwd_ns[r] * lk, lend = 1)
        y <- y - lh
      }
      ## column 2: inner-circle bars, the six transcript categories (2 x 3)
      x1 <- 0.40; y <- top
      text(x1, y, "Cis-acting pairs (inner circle)", adj = c(0, 0.5), cex = cex, col = ink)
      c_lab <- sub("^Novel ", "novel ", sub("^Annotated ", "annotated ", names(cat_cols)))
      bw <- 0.018
      for (r in seq_along(cat_cols)) {
        cx <- x1 + if (r > 3) 0.21 else 0
        cy <- top - lh * (1 + (r - 1) %% 3)
        rect(cx, cy - 0.3 * lh, cx + bw, cy + 0.3 * lh, col = cat_cols[r], border = NA)
        text(cx + bw + 0.008, cy, c_lab[r], adj = c(0, 0.5), cex = cex, col = ink)
      }
      ## column 3: links and strand signs
      x2 <- 0.83; y <- top - lh
      for (neg in c(FALSE, TRUE)) {
        segments(x2, y, x2 + 0.03, y, col = ink, lwd = 1.3 * lk, lty = if (neg) 2 else 1)
        text(x2 + 0.04, y, if (neg) "link: \u03c1 < 0" else "link: \u03c1 > 0", adj = c(0, 0.5),
             cex = cex, col = ink)
        y <- y - lh
      }
      text(x2, y, "+ / \u2212 = transcript strand", adj = c(0, 0.5), cex = cex, col = ink)
    }

    ## ---- Fig. 2r as printed: title + two-column key under the circle ------------
    ## coordinates in mm (x from the left, y from the top of the panel); the
    ## circle centre is the centre of the W x W square at the top
    if (key_paper) {
      W  <- fs$dev_mm[1]; H <- fs$dev_mm[2]; K <- CIRCOS_FIG_KEY
      par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
      plot.new(); plot.window(c(0, W), c(H, 0), xaxs = "i", yaxs = "i")
      cx  <- W / 2; cy <- W / 2
      cex <- fs$key
      text(cx + K$title_dx, cy + K$title_dy, cis_variant$title %||% "Cis-acting regulation pairs", cex = cex, col = ink)
      text(cx + K$title_dx, cy + K$sub_dy, "line = log2FC H1975 / HCC827", cex = cex, col = "grey35")
      ry <- cy + K$row1_dy + K$row_h * (0:6)            # row centres
      xl <- cx + K$x_lab
      ## left column: strand-track lines (rows 1-3), cis links (5-6), strand signs (7)
      text(cx + c(K$x_de, K$x_ns), ry[1], c("DE", "not DE"), cex = cex, col = ink)
      h  <- K$sw_len / 2
      for (r in 1:2) {
        y <- ry[r + 1]
        text(xl, y, c("annotated", "novel")[r], adj = c(0, 0.5), cex = cex, col = ink)
        segments(cx + K$x_de - h, y, cx + K$x_de + h, y, col = K$de_col[r], lwd = K$de_lwd[r], lend = 1)
        segments(cx + K$x_ns - h, y, cx + K$x_ns + h, y, col = K$ns_col[r], lwd = K$ns_lwd[r], lend = 1)
      }
      for (k in 1:2) {
        y <- ry[4 + k]
        segments(xl, y, xl + K$link_len, y, col = ink, lwd = 1.3 * lk, lend = 1,
                 lty = if (k == 1) "solid" else K$link_lty_neg)
        text(xl + K$link_txt, y, c("link: ρ > 0", "link: ρ < 0")[k], adj = c(0, 0.5),
             cex = cex, col = ink)
      }
      text(xl, ry[7], "+ / − = transcript strand", adj = c(0, 0.5), cex = cex, col = "grey35")
      ## right column: colours of the transcript bars in the inner circle
      c_lab <- sub("^Novel ", "novel ", sub("^Annotated ", "", names(cat_cols)))
      for (r in seq_along(cat_cols)) {
        rect(cx + K$x_box, ry[r] - K$box_h / 2, cx + K$x_box + K$box_w, ry[r] + K$box_h / 2,
             col = cat_cols[r], border = NA)
        text(cx + K$x_box_txt, ry[r], c_lab[r], adj = c(0, 0.5), cex = cex, col = ink)
      }
    }

    ## ---- key under the circle: the DE / not-DE encoding ----------------------
    if (!isTRUE(fs$centre_key) && !key_right && !key_below && !key_paper) {
    par(fig = if (part == "both") c(0, 0.58, 0, key_h) else c(0, 1, 0, key_h),
        new = TRUE, mar = c(0, 0, 0, 0))
    plot.new(); plot.window(c(0, 1), c(0, 1))
    if (fs$notes) {  # combined figure
      k_lab <- c("annotated", sub("^Novel ", "novel ", names(nov_cols)))
      kx <- seq(0.34, 0.88, length.out = length(k_lab)); kw <- 0.035
      y_lab <- 0.92; y_de <- 0.72; y_ns <- 0.54
      text(0.04, 0.92, "Strand tracks", adj = 0, cex = fs$key, font = 1, col = ink)
    } else {                                 # poster: 2-line labels, no captions
      k_lab <- c("annotated", sub("^Novel ", "novel\n", names(nov_cols)))
      kx <- seq(0.3, 0.9, length.out = length(k_lab)); kw <- 0.05
      y_lab <- 0.8; y_de <- 0.5; y_ns <- 0.3
    }
    text(kx, y_lab, k_lab, cex = fs$key, col = ink)
    text(if (fs$notes) 0.27 else kx[1] - 0.1, c(y_de, y_ns), c("DE", "not DE"),
         adj = 1, cex = fs$key, col = ink)
    segments(kx - kw, y_de, kx + kw, y_de, col = k_de, lwd = k_lwd_de)
    segments(kx - kw, y_ns, kx + kw, y_ns, col = k_ns, lwd = k_lwd_ns)
    if (fs$notes) {
      text(0.5, 0.36, paste0("One line per tested transcript; length = log2FC H1975 / HCC827 ",
                             "(outward = higher in H1975, capped at \u00b1", lfc_cap, ")."),
           cex = fs$note, col = "grey25")
      text(0.5, 0.21, sprintf("Inner circle: pair %s %sneighbour \u2264 %g Mb (|\u03c1| > %.1f), zoomed; outer lane = novel, inner lane = neighbour.",
                              if (one_per_novel) "= DE novel transcript + its most correlated" else "= DE novel transcript + a correlated",
                              if (de_filter == "both") "DE " else "",
                              window_bp / 1e6, rho_min),
           cex = fs$note, col = "grey25")
      text(0.5, 0.06, "Bar = transcript, +/\u2212 = its strand; link = cis pair (solid: \u03c1 > 0, dashed: \u03c1 < 0); pair numbers as in the heatmap.",
           cex = fs$note, col = "grey25")
    } else {
      text(0.5, 0.07, "line length = log2FC H1975 / HCC827 (outward = higher in H1975)",
           cex = fs$note, col = "grey25")
    }
    }                                            # end: key under the circle
    new_panel <- TRUE
  }

  ## ---- pair heatmap ------------------------------------------------------------
  if (part %in% c("both", "heatmap")) {
    rho_txt   <- sprintf("\u03c1 = %+.2f", sel$rho)
    locus_txt <- paste0(sel$novel_chr, " \u00b7 ", dist_txt)
    row_lab   <- c(rbind(sel$novel, sel$partner))
    if (!new_panel) { plot.new(); new_panel <- TRUE }   # strwidth() needs an open page
    ## margins (in lines) sized to the labels, so larger poster text never clips
    csi   <- par("csi")
    mar_l <- (max(strwidth(row_lab, units = "inches", cex = fs$hm)) + 0.18 + 0.12 +
              strwidth("15", units = "inches", cex = fs$hm_num, font = 1) + 0.1) / csi
    mar_r <- (max(strwidth(rho_txt, units = "inches", cex = fs$hm),
                  strwidth(locus_txt, units = "inches", cex = fs$hm_small)) + 0.2) / csi
    mar_t <- if (part == "both") 4 else 3 + 2.2 * fs$hm_title
    hm_fig <- if (part == "both") c(0.6, 0.98, 0.15, 0.93) else c(0, 1, 0.2, 1)
    par(fig = hm_fig, new = new_panel, mar = c(1, mar_l, mar_t, mar_r))
    nr <- nrow(hm); nc <- ncol(hm)
    image(seq_len(nc), seq_len(nr), t(pmax(pmin(hm, z_lim), -z_lim)[nr:1, , drop = FALSE]),
          col = hm_cols, zlim = c(-z_lim, z_lim), axes = FALSE, xlab = "", ylab = "")
    n_ref <- sum(grp[smp] == grp_ref)
    segments(seq_len(nc - 1) + 0.5, 0.5, seq_len(nc - 1) + 0.5, nr + 0.5, col = "white", lwd = 0.6)
    segments(n_ref + 0.5, 0.5, n_ref + 0.5, nr + 0.5, col = "white", lwd = 2.5)

    u_in    <- function(inch) inch * diff(par("usr")[1:2]) / par("pin")[1]  # inches -> x units
    y_nov   <- nr - 3 * (pid - 1)             # novel row; its partner is the row below
    row_y   <- c(rbind(y_nov, y_nov - 1))
    row_cat <- c(rbind(as.character(sel$novel_category), as.character(sel$partner_category)))
    ## category chip (same colours as the circle), name, pair number
    rect(0.5 - u_in(0.13), row_y - 0.42, 0.5 - u_in(0.05), row_y + 0.42,
         col = cat_cols[row_cat], border = NA)
    text(0.5 - u_in(0.18), row_y, row_lab, adj = c(1, 0.5), cex = fs$hm, col = ink)
    x_num <- 0.5 - u_in(0.18) - max(strwidth(row_lab, cex = fs$hm)) - u_in(0.12)
    text(x_num, y_nov - 0.5, pid, adj = c(1, 0.5), cex = fs$hm_num, font = 1, col = ink)
    ## rho as neutral text (the red/blue scale belongs to the z-score only)
    text(nc + 0.5 + u_in(0.1), y_nov, rho_txt, adj = 0, cex = fs$hm, col = ink)
    ## locus of the pair: chromosome (find its wedge) + distance
    text(nc + 0.5 + u_in(0.1), y_nov - 1, locus_txt, adj = 0, cex = fs$hm_small, col = "grey40")

    rect(seq_len(nc) - 0.5, nr + 0.7, seq_len(nc) + 0.5, nr + 1.3,
         col = group_cols[grp[smp]], border = "white")
    for (g in c(grp_ref, grp_num))
      text(mean(which(grp[smp] == g)), nr + 1.6, line_lab[[g]], adj = c(0.5, 0),
           cex = fs$hm_group, font = 1, col = ink)
    mtext("Expression of each pair", side = 3,
          line = if (part == "both") 2.4 else mar_t - 1.2 - fs$hm_title,
          cex = fs$hm_title, font = 1, adj = 0)

    ## chip legend: the six transcript categories (same colours as the circle)
    par(fig = if (part == "both") c(0.62, 0.98, 0.075, 0.145) else c(0, 1, 0.08, 0.19),
        new = TRUE, mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", bty = "n", ncol = if (part == "both") 3 else 2, cex = fs$chip,
           x.intersp = 0.6, text.col = ink,
           legend = sub("^Novel ", "novel ", sub("^Annotated ", "annotated ", names(cat_cols))),
           fill = cat_cols, border = NA)

    ## colour bar: z-score only
    par(fig = if (part == "both") c(0.74, 0.9, 0.02, 0.055) else c(0.2, 0.8, 0, 0.07),
        new = TRUE, mar = c(if (part == "both") 1.2 else 1.8, 0, 0, 0))
    image(seq(-z_lim, z_lim, length.out = 100), 1, matrix(seq(-z_lim, z_lim, length.out = 100)),
          col = hm_cols, axes = FALSE, xlab = "", ylab = "")
    axis(1, at = c(-z_lim, 0, z_lim), cex.axis = fs$cbar, tick = FALSE, line = -0.8)
    mtext(sprintf("z-score of log2 normalized count  |  \u03c1 = Spearman, n = %d samples",
                  ncol(xs)), side = 3, cex = fs$cbar_note, line = 0.15)
  }
}

cairo_ok <- capabilities("cairo") && isTRUE(tryCatch({
  f <- tempfile(fileext = ".png"); png(f, type = "cairo"); dev.off(); unlink(f); TRUE
}, warning = function(w) { try(dev.off(), silent = TRUE); FALSE },
   error = function(e) FALSE))

## ---- Fig. 2q as printed: cis-pair heatmap ---------------------------------------
## rows = the circos pairs whose partner is annotated protein-coding / lncRNA
## (HM_FIG_KEEP), numbered as on the circle; columns = HCC827_1-3 | H1975_1-3
m_idx <- which(as.character(sel$partner_category) %in% HM_FIG_KEEP & !sel$partner_novel)
message(sprintf("Fig. 2n: %d of %d circos pairs (left out: %s)", length(m_idx), nrow(sel),
                paste(setdiff(seq_len(nrow(sel)), m_idx), collapse = ", ")))
smp_m <- c(sort(names(grp)[grp == grp_ref]), sort(names(grp)[grp == grp_num]))
z_m <- log2(rds$norm[unique(c(sel$novel[m_idx], sel$partner[m_idx])), smp_m, drop = FALSE] + 1)
z_m <- t(scale(t(z_m)))
z_m[!is.finite(z_m)] <- 0                       # a row without variance = 0
z_m_lim  <- max(abs(z_m))                       # colour limits = +/- the largest |z|
dist_m   <- ifelse(sel$distance_bp == 0, "overlapping",
                   paste(prettyNum(round(sel$distance_bp / 1e3), big.mark = ","), "kb"))
hm_fig_height <- function(n) with(HM_FIG, y_tile0 + (n - 1) * pair_dy + row_dy + tile_h + 0.56)

draw_hm_fig <- function(fam) {
  G <- HM_FIG; n <- length(m_idx); H <- hm_fig_height(n)
  par(family = fam, mar = c(0, 0, 0, 0), xpd = NA)
  plot.new(); plot.window(c(0, G$width), c(H, 0), xaxs = "i", yaxs = "i")   # mm, y down
  ramp <- grDevices::colorRamp(hm_anchor, space = "Lab")
  zcol <- function(v) grDevices::rgb(ramp(pmin(1, pmax(0, (v + z_m_lim) / (2 * z_m_lim)))),
                                     maxColorValue = 255)
  n_ref  <- sum(grp[smp_m] == grp_ref)
  x_left <- c(G$x_tile1 + G$tile_dx * (seq_len(n_ref) - 1),
              G$x_tile2 + G$tile_dx * (seq_len(length(smp_m) - n_ref) - 1))
  ## line names, centred over their tiles
  for (g in 1:2) {
    k <- if (g == 1) seq_len(n_ref) else setdiff(seq_along(smp_m), seq_len(n_ref))
    text(mean(range(x_left[k]) + c(0, G$tile_w)), G$y_head, line_lab[[c(grp_ref, grp_num)[g]]],
         cex = G$pt_head / 12, col = "black")
  }
  for (j in seq_len(n)) {
    i   <- m_idx[j]
    top <- G$y_tile0 + (j - 1) * G$pair_dy
    yc  <- top + c(0, G$row_dy) + G$tile_h / 2          # novel row, partner row
    f   <- c(sel$novel[i], sel$partner[i])
    cat <- as.character(c(sel$novel_category[i], sel$partner_category[i]))
    for (r in 1:2) {
      yt <- top + (r - 1) * G$row_dy
      rect(x_left, yt, x_left + G$tile_w, yt + G$tile_h, col = zcol(z_m[f[r], smp_m]), border = NA)
      rect(G$x_cat, yt, G$x_cat + G$cat_w, yt + G$tile_h, col = cat_cols[cat[r]], border = NA)
      text(G$x_name, yc[r], f[r], adj = c(1, 0.5), cex = G$pt_name / 12, col = "black")
    }
    text(G$x_num, mean(yc), i, adj = c(1, 0.5), cex = G$pt_name / 12, col = "black")
    text(G$x_side, yc[1], sprintf("ρ = %+.2f", sel$rho[i]), adj = c(0, 0.5),
         cex = G$pt_side / 12, col = "black")
    text(G$x_side, yc[2], sprintf("%s · %s", sel$novel_chr[i], dist_m[i]), adj = c(0, 0.5),
         cex = G$pt_side / 12, col = "grey40")
  }
  ## colour bar (top right): -lim at the bottom, +lim at the top
  ymap <- function(v) G$bar_top + G$bar_h * (z_m_lim - v) / (2 * z_m_lim)
  nb <- 100; v0 <- seq(-z_m_lim, z_m_lim, length.out = nb + 1)
  bc <- zcol((head(v0, -1) + tail(v0, -1)) / 2)
  rect(G$bar_x, ymap(head(v0, -1)), G$bar_x + G$bar_w, ymap(tail(v0, -1)),
       col = bc, border = bc, lwd = 0.25)                 # border = fill: no seams between slices
  tk <- seq(-floor(2 * z_m_lim) / 2, floor(2 * z_m_lim) / 2, by = 0.5)
  segments(G$bar_x, ymap(tk), G$bar_x + 0.2 * G$bar_w, ymap(tk), col = "white", lwd = 0.5)
  segments(G$bar_x + 0.8 * G$bar_w, ymap(tk), G$bar_x + G$bar_w, ymap(tk), col = "white", lwd = 0.5)
  text(G$x_bar_lab, ymap(tk), formatC(tk, format = "f", digits = 1), adj = c(0, 0.5),
       cex = G$pt_bar / 12, col = "black")
  text(G$x_bar_title, G$y_bar_title, c("z-score", "(logCPM)"), adj = c(0, 0.5),
       cex = G$pt_bar / 12, col = "black")
}

## combined figure (screen / paper)
for (ext in c("png", "pdf", "svg"))
  save_plot(file.path(outdir, paste0("circos_cis_de.", ext)), 14, 8.5, draw)

## Fig. 2r and 2q at their print size, as in the printed Figure 2 (place at 100 %)
for (ext in c("png", "pdf", "svg")) {
  save_plot(file.path(outdir, paste0("paper_cis_circos.", ext)),
            CIRCOS_FIG_MM[1] / 25.4, CIRCOS_FIG_MM[2] / 25.4,
            function(fam) draw(fam, part = "circos", fs = FS_CIRCOS_FIG), res = 600)
  if (length(m_idx))
    save_plot(file.path(outdir, paste0("paper_cis_heatmap.", ext)),
              HM_FIG$width / 25.4, hm_fig_height(length(m_idx)) / 25.4, draw_hm_fig, res = 600)
}
## journal panel with the key on the right of the circle; TRUE to write it as well
write_v5_keyright <- FALSE
if (write_v5_keyright) for (ext in c("png", "pdf", "svg"))
  save_plot(file.path(outdir, paste0("paper_cis_circos_keyright.", ext)),
            PAPER_CIRCOS_MM[1] / 25.4, PAPER_CIRCOS_MM[2] / 25.4,
            function(fam) draw(fam, part = "circos", fs = FS_PAPER), res = 600)

## poster panels at their real print size (place at 100 % on the poster)
for (ext in c("png", "pdf", "svg")) {
  save_plot(file.path(outdir, paste0("poster_cis_circos.", ext)),
            cm2in(POSTER_CIRCOS_CM[1]), cm2in(POSTER_CIRCOS_CM[2]),
            function(fam) draw(fam, part = "circos", fs = FS_POSTER))
  save_plot(file.path(outdir, paste0("poster_cis_heatmap.", ext)),
            cm2in(POSTER_HEATMAP_CM[1]), cm2in(POSTER_HEATMAP_CM[2]),
            function(fam) draw(fam, part = "heatmap", fs = FS_POSTER))
}
