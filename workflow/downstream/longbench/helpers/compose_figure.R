## ===========================================================================
## compose_figure.R - assemble single-panel files into one figure page
## (sourced by 09_assemble_figures.R)
##
##   compose_figure(lay, page, panel_dir, out_stem, dpi = 300)
##     lay        data frame, one row per slot: letter, file (panel file name
##                without extension; NA = empty slot, e.g. the workflow drawn
##                in the vector editor), x, y, w, h (mm from the page's top-left)
##                optional column ly: y of the panel letter (default y)
##     page       c(w = , h = ) in mm
##     panel_dir  folder with <file>.svg and <file>.png (save_panel() output)
##     out_stem   output path without extension: writes <out_stem>.svg and .png
##
## SVG: each panel's own SVG nested at its slot (ids prefixed with the panel
## letter so they stay unique), text stays editable (svglite <text>, font
## "Source Sans 3"); panel letters in their own group, 12 pt plain, in the
## gutter above each slot's top-left corner. PNG: the panels' own PNGs
## composited at their slots (a preview; the SVG is the figure).
## Warns when a panel file's size differs from its slot (re-render the panel).
## ===========================================================================

compose_figure <- function(lay, page, panel_dir, out_stem, dpi = 300) {
  stopifnot(all(c("letter", "file", "x", "y", "w", "h") %in% names(lay)))
  pt2mm <- 25.4 / 72
  letter_mm <- 12 * pt2mm                                     # 12 pt letters
  ly <- if (!is.null(lay$ly)) lay$ly else lay$y              # optional: letter height apart from the slot (Fig. 2 row 3)
  hdr <- sprintf(paste0('<?xml version="1.0" encoding="UTF-8"?>\n',
                        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" ',
                        'width="%gmm" height="%gmm" viewBox="0 0 %g %g">\n',
                        '<rect width="100%%" height="100%%" fill="white"/>\n'),
                 page[["w"]], page[["h"]], page[["w"]], page[["h"]])
  letters_svg <- paste0('<g id="panel_letters" font-family="Source Sans 3" font-size="', round(letter_mm, 3),
                        '" fill="black">\n',
                        paste(sprintf('  <text x="%g" y="%g">%s</text>', lay$x, ly - 0.7, lay$letter),
                              collapse = "\n"), '\n</g>\n')

  nest <- function(r) {
    f <- file.path(panel_dir, paste0(r$file, ".svg"))
    if (!file.exists(f)) stop("panel not found: ", f)
    s <- paste(readLines(f, warn = FALSE), collapse = "\n")
    s <- sub("^<\\?xml[^>]*>\\s*", "", s)
    root <- regmatches(s, regexpr("^<svg[^>]*>", s))
    vb   <- regmatches(root, regexpr("viewBox=['\"][^'\"]+['\"]", root))
    dims <- as.numeric(regmatches(root, gregexpr("(?<=(width|height)=['\"])[0-9.]+(?=pt)", root, perl = TRUE))[[1]])
    if (length(dims) == 2 && any(abs(dims * pt2mm - c(r$w, r$h)) > 0.5))
      warning(sprintf("%s is %.1f x %.1f mm but its slot is %g x %g mm: re-render the panel",
                      r$file, dims[1] * pt2mm, dims[2] * pt2mm, r$w, r$h), call. = FALSE)
    inner <- sub("^<svg[^>]*>", "", s)
    inner <- sub("</svg>\\s*$", "", inner)
    pf <- paste0(r$letter, "_")
    inner <- gsub("id=(['\"])", paste0("id=\\1", pf), inner)
    inner <- gsub("url\\(#", paste0("url(#", pf), inner)
    inner <- gsub("href=(['\"])#", paste0("href=\\1#", pf), inner)
    sprintf('<svg id="panel_%s" x="%g" y="%g" width="%g" height="%g" %s preserveAspectRatio="none">%s</svg>',
            r$letter, r$x, r$y, r$w, r$h, vb, inner)
  }
  filled <- lay[!is.na(lay$file), , drop = FALSE]
  body <- paste(vapply(seq_len(nrow(filled)), function(i) nest(filled[i, ]), ""), collapse = "\n")
  writeLines(paste0(hdr, body, "\n", letters_svg, "</svg>\n"), paste0(out_stem, ".svg"))

  ## PNG preview
  px <- function(mm) round(mm / 25.4 * dpi)
  canvas <- magick::image_blank(px(page[["w"]]), px(page[["h"]]), "white")
  for (i in seq_len(nrow(filled))) {
    r  <- filled[i, ]
    im <- magick::image_resize(magick::image_read(file.path(panel_dir, paste0(r$file, ".png"))),
                               sprintf("%dx%d!", px(r$w), px(r$h)))
    canvas <- magick::image_composite(canvas, im, offset = sprintf("+%d+%d", px(r$x), px(r$y)))
  }
  fs <- round(12 / 72 * dpi)
  for (i in seq_len(nrow(lay)))
    canvas <- magick::image_annotate(canvas, lay$letter[i], size = fs, font = "Source Sans 3",
                                     location = sprintf("+%d+%d", px(lay$x[i]), px(ly[i] - 0.7) - round(fs * 0.8)))
  magick::image_write(canvas, paste0(out_stem, ".png"))
  message(sprintf("  assembled %s.{svg,png} (%g x %g mm, %d panels)",
                  basename(out_stem), page[["w"]], page[["h"]], nrow(filled)))
  invisible(lay)
}

## slot table from rows: each row = list(h = height, cells = list(...)); a cell
## is list(letter, file, w) or, for slots stacked in one column,
## list(w = , stack = list(list(letter, file, h), ...)). Rows and cells are
## separated by `gutter`; `margin` around the page.
layout_rows <- function(rows, margin = 5, gutter = 4) {
  out <- list(); y <- margin
  for (rw in rows) {
    x <- margin
    for (cl in rw$cells) {
      if (!is.null(cl$stack)) {
        yy <- y
        for (s in cl$stack) {
          out[[length(out) + 1]] <- data.frame(letter = s[[1]], file = s[[2]], x = x, y = yy, w = cl$w, h = s[[3]])
          yy <- yy + s[[3]] + gutter
        }
        stopifnot("stacked slots must fill the row height" = abs(yy - gutter - (y + rw$h)) < 1e-9)
        x <- x + cl$w + gutter
      } else {
        out[[length(out) + 1]] <- data.frame(letter = cl[[1]], file = cl[[2]], x = x, y = y, w = cl[[3]], h = rw$h)
        x <- x + cl[[3]] + gutter
      }
    }
    y <- y + rw$h + gutter
  }
  lay <- do.call(rbind, out)
  lay$file[lay$file %in% c("", "NA")] <- NA
  lay
}

## every row must end at the right margin and the last row at the bottom margin
check_page <- function(lay, page, margin = 5) {
  right <- tapply(lay$x + lay$w, lay$y, max)
  if (any(abs(right - (page[["w"]] - margin)) > 1e-9))
    stop("a row does not end at the right margin: ", paste(names(right), right, collapse = "; "))
  if (abs(max(lay$y + lay$h) - (page[["h"]] - margin)) > 1e-9)
    stop(sprintf("the last row ends at %g mm, the bottom margin at %g mm",
                 max(lay$y + lay$h), page[["h"]] - margin))
  invisible(TRUE)
}
