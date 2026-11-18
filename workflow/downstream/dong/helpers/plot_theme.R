## ===========================================================================
## plot_theme.R - shared look for every plot (dong/de copy)
## ===========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(showtext)
})

# ---------------------------------------------------------------- theme -----
## Every figure is drawn in Source Sans 3. The font is registered for showtext
## (PDF / PNG) from, in this order:
##   1. the static TTF/OTF files installed on this computer (also needed for the
##      editable SVG text; Windows: "Install for all users" or per-user install)
##   2. Google Fonts (needs internet)
## If neither works the script STOPS instead of silently drawing in Arial.
font_family <- "Source Sans 3"
register_figure_font <- function(fam = font_family) {
  if (fam %in% sysfonts::font_families()) return(invisible("already registered"))
  ## 1. installed fonts (Windows per-user fonts are not in the default search path)
  user_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "Microsoft", "Windows", "Fonts")
  if (nzchar(Sys.getenv("LOCALAPPDATA")) && dir.exists(user_dir))
    try(sysfonts::font_paths(user_dir), silent = TRUE)
  ff <- tryCatch(sysfonts::font_files(), error = function(e) NULL)
  if (!is.null(ff) && nrow(ff)) {
    ff <- ff[ff$family %in% fam & !grepl("\\[|VF|Variable", ff$file), , drop = FALSE]  # static files only
    pick <- function(face) {
      f <- ff[ff$face == face, , drop = FALSE]
      if (nrow(f)) file.path(f$path[1], f$file[1]) else NULL
    }
    reg <- pick("Regular")
    if (!is.null(reg)) {
      ok <- tryCatch({
        sysfonts::font_add(fam, regular = reg, bold = pick("Bold"), italic = pick("Italic"),
                           bolditalic = pick("Bold Italic"))
        TRUE
      }, error = function(e) FALSE)
      if (ok) return(invisible(paste("installed font", reg)))
    }
  }
  ## 2. Google Fonts
  ok <- tryCatch({ sysfonts::font_add_google(fam, fam); TRUE }, error = function(e) FALSE)
  if (ok) return(invisible("Google Fonts"))
  stop("Font '", fam, "' not found: install it (Google Fonts 'Source Sans 3', static TTF ",
       "files) or connect to the internet, then restart R.", call. = FALSE)
}
message("Figure font: ", font_family, " (", register_figure_font(), ")")
showtext_auto()
showtext_opts(dpi = 300)

base_theme <- theme(
  strip.background = element_rect(fill = "#f0f0f0", color = "black"),
  plot.background  = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
)
theme_set(theme_bw(base_family = font_family, base_size = 10) + base_theme)
update_geom_defaults("text",  list(family = font_family, size = 8 / .pt))
update_geom_defaults("label", list(family = font_family, size = 8 / .pt))
## ggrepel labels (volcanos, PCA) do not inherit the "text" defaults
update_geom_defaults(GeomTextRepel,  list(family = font_family))
update_geom_defaults(GeomLabelRepel, list(family = font_family))
## base R plots (plotDispEsts) and pheatmap get font_family passed explicitly

# --------------------------------------------------------------- colours ----
## dong dataset: H1975 (TKI-resistant) vs HCC827 (TKI-sensitive)
group_cols    <- c(TKI_sensitive = "#4DBBD5", TKI_resistant = "#E64B35")
protocol_cols <- c("ONT cDNA" = "#27408B", "ONT DRS" = "#008B45", "PacBio" = "#D02090")
## NPG palette (ggsci pal_npg "nrc", Nature Reviews Cancer-inspired), for other categories
npg_cols <- c(red = "#E64B35", lightblue = "#4DBBD5", green = "#00A087", navy = "#3C5488",
              salmon = "#F39B7F", greyblue = "#8491B4", mint = "#91D1C2", darkred = "#DC0000",
              brown = "#7E6148", beige = "#B09C85")
plat_labels   <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio")
