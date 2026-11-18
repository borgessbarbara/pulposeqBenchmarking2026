## ===========================================================================
## fig2_rows.R - Fig. 2 rows 1 and 2: aligned plot boxes across 08 and 10
## (rows 1-2 aligned; row 2 = equal square boxes sized to fit the 340 mm
## page; all text in these rows 8 pt, in-plot notes 7 pt)
##
##   row 1  a read support | b novel transcripts by class code | c sequin DE Bambu | oarfish |
##          d sequin length ratio
##   row 2  e human counts | f overdispersion by category | g overdispersion by
##          isoforms per gene | h sample correlation
##   one square size s for every box of rows 1-2 (row 2 fills the page);
##   row 1 boxes s tall, a's bars take the width that is left
## Within a row every file gets the row's largest top / bottom room, so the
## boxes share their top and bottom edges and every file of a row is equally tall.
##
## How: each script measures the room (mm) its panels need around the box
## (axis text, titles, legend; independent of the box size) and stores it in
## results/fig2_rooms.rds; the layout below is computed from all of them.
## 08 (a, b, h) and 10 (c-g) each write their own entries. Run 08, 10, then 08
## once more the first time (or whenever a legend / label width changes);
## 11 checks that the saved files match the layout and says which to rerun.
## ===========================================================================

fig2_W <- 340; fig2_mar <- 5; fig2_gut <- 4; fig2_slack <- 0.5
## every scatter / line / tile box of rows 1-2 is one square size s, set by
## row 2's four squares filling the page; row 1: a's bars fig2_a_w (legend on
## the right), b's bars fig2_b_w, c's two facets and d share the rest (w1
## each), all s tall
fig2_rows <- list(
  r1 = c(a = "fig2a_read_support", b = "fig2b_novel_biotype", c = "ed_q_b_sequin_de",
         d = "ed_q_a_sequin_length"),
  r2 = c(e = "ed_q_c_human_counts", f = "ed_q_f_overdispersion",
         g = "ed_q_g_overdispersion_isoforms", h = "fig2e_sample_correlation"))
fig2_b_w    <- 22  # row 1, b: bar area width
fig2_a_w    <- 36  # row 1, a: bar area width (legend on the right)
fig2_c_gap  <- 3                          # c: gap between the two facets (10's b_gap)
fig2_h_scale <- 0.85  # row 2, h (sample-correlation heatmap): square box this share of s
fig2_rooms_file <- file.path(res_dir, "fig2_rooms.rds")

## 8 pt everywhere in rows 1-2 (10's theme has 6.4 pt ticks / legend text)
fig2_text <- ggplot2::theme(
  axis.text    = ggplot2::element_text(size = 8),  axis.title   = ggplot2::element_text(size = 8),
  legend.text  = ggplot2::element_text(size = 8),  legend.title = ggplot2::element_text(size = 8),
  strip.text   = ggplot2::element_text(size = 8))

## the plot with its box(es) fixed; c has two facets of w each
fig2_fix <- function(p, w, h) p + ggh4x::force_panelsizes(rows = grid::unit(h, "mm"), cols = grid::unit(w, "mm"))

## room (mm) around the box(es): left / right / top / bottom, without plot margin
fig2_room <- function(p, box = 30) {
  g  <- ggplot2::ggplotGrob(fig2_fix(p, box, box) + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0)))
  pl <- g$layout[grepl("^panel", g$layout$name), ]
  w  <- grid::convertWidth(g$widths, "mm", valueOnly = TRUE)
  h  <- grid::convertHeight(g$heights, "mm", valueOnly = TRUE)
  c(L = sum(w[seq_len(min(pl$l) - 1)]), R = sum(w[-seq_len(max(pl$r))]),
    T = sum(h[seq_len(min(pl$t) - 1)]), B = sum(h[-seq_len(max(pl$b))]))
}

## measure on the device that writes the SVG (svglite, showtext paused):
## pdf(NULL) has other font metrics, so the boxes would drift from their row
fig2_measure_open <- function() {
  st <- if (exists("showtext_pause")) showtext_pause() else FALSE
  if (exists("svg_device")) svg_device(tempfile(fileext = ".svg"), 7, 7) else grDevices::pdf(NULL)
  assign(".fig2_st", st, envir = globalenv())
}
fig2_measure_close <- function() {
  grDevices::dev.off()
  if (exists("showtext_resume")) showtext_resume(get0(".fig2_st", envir = globalenv()))
}

## measure the given plots (named by file name) and store them with the others
fig2_update_rooms <- function(plots) {
  r <- if (file.exists(fig2_rooms_file)) readRDS(fig2_rooms_file) else list()
  fig2_measure_open(); on.exit(fig2_measure_close())
  for (nm in names(plots)) r[[nm]] <- fig2_room(plots[[nm]])
  dir.create(dirname(fig2_rooms_file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(r, fig2_rooms_file)
  r
}

## the layout from the stored rooms. A panel not measured yet gets a rough
## guess (and a message): rerun its script, then this one.
fig2_layout <- function(rooms = readRDS(fig2_rooms_file)) {
  guess <- c(L = 12, R = 35, T = 2, B = 10)
  for (nm in unlist(fig2_rows)) if (is.null(rooms[[nm]])) {
    message("  ! fig2_rows: no room stored yet for ", nm, " (rough guess used; rerun its script, then this one)")
    rooms[[nm]] <- guess
  }
  usable <- fig2_W - 2 * fig2_mar
  lr <- function(nm) unname(rooms[[nm]]["L"] + rooms[[nm]]["R"] + 2 * fig2_slack)
  r1 <- fig2_rows$r1; r2 <- fig2_rows$r2
  ## row 2: equal squares filling the row (h at fig2_h_scale of them); s is also row 1's box height
  box2 <- (usable - fig2_gut * (length(r2) - 1) - sum(sapply(r2, lr))) / (length(r2) - 1 + fig2_h_scale)
  ## row 1 (narrower c and d so a's legend fits beside its bars):
  ## a's bars fig2_a_w, b's fig2_b_w; c's two facets and d share what is left,
  ## each w1 wide, all s tall (row 1 and row 2 boxes keep one height)
  a_w  <- fig2_a_w
  w1   <- (usable - fig2_gut * (length(r1) - 1) - sum(sapply(r1, lr)) - fig2_b_w - a_w - fig2_c_gap) / 3
  c_in <- 2 * w1 + fig2_c_gap
  box <- list()                                   # per file: col width, box height, inner width
  box[[r1[["a"]]]] <- c(col = a_w, h = box2, inner = a_w)
  box[[r1[["b"]]]] <- c(col = fig2_b_w, h = box2, inner = fig2_b_w)
  box[[r1[["c"]]]] <- c(col = w1, h = box2, inner = c_in)
  box[[r1[["d"]]]] <- c(col = w1, h = box2, inner = w1)
  if (w1 < 25) message(sprintf("  ! fig2_rows: row 1, c and d boxes only %.1f mm wide", w1))
  for (nm in r2) box[[nm]] <- c(col = box2, h = box2, inner = box2)
  hs <- fig2_h_scale * box2                       # h: smaller square, top-aligned with the row
  box[[r2[["h"]]]] <- c(col = hs, h = hs, inner = hs)
  tb <- lapply(fig2_rows, function(rw) c(T = max(sapply(rw, function(nm) rooms[[nm]]["T"])) + fig2_slack,
                                          B = max(sapply(rw, function(nm) rooms[[nm]]["B"])) + fig2_slack))
  size <- list()
  for (i in seq_along(fig2_rows)) for (nm in fig2_rows[[i]])
    size[[nm]] <- c(w = unname(rooms[[nm]]["L"] + box[[nm]]["inner"] + rooms[[nm]]["R"] + 2 * fig2_slack),
                    h = unname(tb[[i]]["T"] + box2 + tb[[i]]["B"]))     # every file of a row equally tall
  if (box2 < 25) message(sprintf("  ! fig2_rows: boxes only %.1f mm", box2))
  if (a_w < 25)  message(sprintf("  ! fig2_rows: row 1, a's bars only %.1f mm wide", a_w))
  list(rooms = rooms, box = box, tb = tb, size = size, box2 = box2, a_w = a_w, w1 = w1)
}

## save one row-1/2 panel at its layout size, box fixed, top / bottom padded to the row's room
fig2_save <- function(p, nm, lay, dir) {
  i  <- which(sapply(fig2_rows, function(rw) nm %in% rw))
  r  <- lay$rooms[[nm]]; bx <- lay$box[[nm]]; tb <- lay$tb[[i]]
  p  <- fig2_fix(p, bx[["col"]], bx[["h"]]) +
    ggplot2::theme(plot.margin = ggplot2::margin(tb[["T"]] - r[["T"]], fig2_slack,
                                                  tb[["B"]] - r[["B"]] + lay$box2 - bx[["h"]], fig2_slack, "mm"))   # a smaller box (h): extra room below
  save_panel(p, nm, lay$size[[nm]], dir)
  message(sprintf("    %s: box %.1f x %.1f mm", nm, bx[["inner"]], bx[["h"]]))
}

## row 3: volcano | these four profiles (06, section 7e);
## below: DTE heatmap | EGFR locus | LY6D locus, then circos | cis-pair heatmap
fig2_row3_profiles <- c("fig2h_profile_EGFR", "fig2j_profile_UBE2I", "fig2k_profile_RAB11B-AS1",
                        "fig2m_profile_LY6D")
fig2_loci <- c("EGFR", "LY6D")
