## ===========================================================================
## figure_helpers.R - saving panels at print size (sourced by the panel scripts)
##
##   ggsave_svg()  SVG with editable text (showtext paused while svglite writes)
##   save_panel()  one panel -> <name>.pdf, <name>.svg (editable), <name>.png (600 dpi)
##   keep_level()  keep one facet level in a plot and in every layer's own data
##
## Editable SVG text: plot_theme.R switches showtext on, and showtext draws
## every letter as an outline on EVERY device, SVG included. svg_device()
## switches showtext off while an SVG is written, so svglite writes real
## <text> elements in the figure font ("Source Sans 3" must be installed on the
## computer that edits the SVG). fix_text_size = FALSE: edited text is not
## stretched to the old width. PDF / PNG keep showtext.
## ===========================================================================

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

## one panel: PDF + editable SVG + PNG; size in mm
save_panel <- function(p, name, size_mm, dir, dpi = 600) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(dir, name); w <- size_mm[1] / 25.4; h <- size_mm[2] / 25.4
  ggplot2::ggsave(paste0(out, ".pdf"), p, width = w, height = h)
  ggsave_svg(paste0(out, ".svg"), p, width = w, height = h, bg = "white")
  showtext::showtext_opts(dpi = dpi)             # showtext must know the PNG dpi
  ggplot2::ggsave(paste0(out, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  showtext::showtext_opts(dpi = 300)             # back to the plot_theme.R default
  message(sprintf("  saved %s (%g x %g mm)", name, size_mm[1], size_mm[2]))
}

## keep one facet level (e.g. isoform) in a plot and in every layer's own data
keep_level <- function(p, lev = "isoform", var = "level") {
  p$data <- dplyr::filter(p$data, .data[[var]] == lev)
  for (k in seq_along(p$layers)) {
    d <- p$layers[[k]]$data
    if (is.data.frame(d) && var %in% names(d)) p$layers[[k]]$data <- dplyr::filter(d, .data[[var]] == lev)
  }
  p
}

## read-level support bars (Fig. 2a): horizontal stacked
## bars, one per biotype; `rs` = rows of read_support_*_long.csv (one run or
## several, with a `protocol` column for facets)
sup_lv   <- c("full-length + unique", "unique only", "full-length only", "neither")
sup_lab  <- c("full-length + unique" = "full-length + unique reads", "unique only" = "unique reads",
              "full-length only" = "full-length reads", "neither" = "partial & ambiguous reads")
sup_cols <- c("full-length + unique" = "#4459CC", "unique only" = "#2EEA9D",
              "full-length only" = "#F0C93A", "neither" = "#C22A05")
sup_txt  <- c("full-length + unique" = "white", "unique only" = "black",
              "full-length only" = "black", "neither" = "white")
bt_disp  <- c("annotated protein-coding" = "protein-coding", "annotated lncRNA" = "lncRNA",
              "other biotypes" = "other", "predicted protein-coding" = "novel protein-coding",
              "novel lncRNA candidate" = "novel lncRNA", "predicted non-coding" = "novel non-coding")

read_support_plot <- function(rs, label_fs, min_lab = 5, n_sep = " ") {
  has_prot <- "protocol" %in% names(rs)
  rs <- rs %>%
    dplyr::filter(biotype %in% names(bt_disp)) %>%           # no summary rows
    dplyr::mutate(bt      = factor(bt_disp[biotype], levels = bt_disp),
                  support = factor(support, levels = sup_lv),
                  ylab    = sprintf("%s%s(%s)", bt, n_sep, scales::comma(total))) %>%
    dplyr::arrange(dplyr::across(dplyr::any_of("protocol")), bt)
  ## one y position per (protocol x) biotype, labelled "biotype (n)"; first biotype on top
  rs$key <- if (has_prot) paste(rs$protocol, rs$bt) else as.character(rs$bt)
  rs$key <- factor(rs$key, levels = rev(unique(rs$key)))
  y_labs <- stats::setNames(rs$ylab, as.character(rs$key))
  ggplot(rs, aes(pct, key, fill = support)) +
    geom_col(width = 0.78, position = position_stack(reverse = TRUE)) +
    geom_text(data = dplyr::filter(rs, pct >= min_lab),
              aes(label = sprintf("%.1f%%", pct), colour = support),
              position = position_stack(vjust = 0.5, reverse = TRUE), size = label_fs / .pt,
              show.legend = FALSE) +
    scale_fill_manual(values = sup_cols, labels = sup_lab, name = NULL, drop = FALSE) +
    scale_colour_manual(values = sup_txt, guide = "none") +
    scale_y_discrete(labels = y_labs) +
    scale_x_continuous(breaks = c(0, 50, 100), labels = c("0", "50%", "100%"), expand = c(0, 0)) +
    labs(x = "Share of isoforms", y = NULL, title = "Read-level support across biotypes") +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank())
}
