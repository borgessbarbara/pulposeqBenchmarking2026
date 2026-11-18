## ===========================================================================
## plot_theme.R - shared look for every plot; sourced by the analysis scripts
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
  if (fam %in% sysfonts::font_families()) {
    return(invisible("already registered"))
  }
  ## 1. installed fonts (Windows per-user fonts are not in the default search path)
  user_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "Microsoft", "Windows", "Fonts")
  if (nzchar(Sys.getenv("LOCALAPPDATA")) && dir.exists(user_dir)) {
    try(sysfonts::font_paths(user_dir), silent = TRUE)
  }
  ff <- tryCatch(sysfonts::font_files(), error = function(e) NULL)
  if (!is.null(ff) && nrow(ff)) {
    ff <- ff[ff$family %in% fam & !grepl("\\[|VF|Variable", ff$file), , drop = FALSE] # static files only
    pick <- function(face) {
      f <- ff[ff$face == face, , drop = FALSE]
      if (nrow(f)) file.path(f$path[1], f$file[1]) else NULL
    }
    reg <- pick("Regular")
    if (!is.null(reg)) {
      ok <- tryCatch(
        {
          sysfonts::font_add(fam,
            regular = reg, bold = pick("Bold"), italic = pick("Italic"),
            bolditalic = pick("Bold Italic")
          )
          TRUE
        },
        error = function(e) FALSE
      )
      if (ok) {
        return(invisible(paste("installed font", reg)))
      }
    }
  }
  ## 2. Google Fonts
  ok <- tryCatch(
    {
      sysfonts::font_add_google(fam, fam)
      TRUE
    },
    error = function(e) FALSE
  )
  if (ok) {
    return(invisible("Google Fonts"))
  }
  stop("Font '", fam, "' not found: install it (Google Fonts 'Source Sans 3', static TTF ",
    "files) or connect to the internet, then restart R.",
    call. = FALSE
  )
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
update_geom_defaults("text", list(family = font_family, size = 8 / .pt))
update_geom_defaults("label", list(family = font_family, size = 8 / .pt))
## ggrepel labels (volcanos, PCA) do not inherit the "text" defaults
update_geom_defaults(GeomTextRepel, list(family = font_family))
update_geom_defaults(GeomLabelRepel, list(family = font_family))
## base R plots (plotDispEsts) and pheatmap get font_family passed explicitly

# --------------------------------------------------------------- colours ----
group_cols <- c(LUAD = "#5858ca", SCLC_A = "#d95f02", SCLC_P = "#1a8e2c")
protocol_cols <- c("ONT cDNA" = "#1f8d75", "ONT DRS" = "#4a59e7", "PacBio" = "#d93b9f")
## NPG palette (ggsci pal_npg "nrc", Nature Reviews Cancer-inspired), for other categories
npg_cols <- c(
  red = "#E64B35", lightblue = "#4DBBD5", green = "#00A087", navy = "#3C5488",
  salmon = "#F39B7F", greyblue = "#8491B4", mint = "#91D1C2", darkred = "#DC0000",
  brown = "#7E6148", beige = "#B09C85"
)
plat_labels <- c(ont_cdna = "ONT cDNA", ont_drs = "ONT DRS", pb = "PacBio")

# ------------------------------------------------ shared DE plot styling ----
## (same choices as the dong DE analysis)
## DE heatmaps (row z-scores): Moreland "cool-warm" blue / light grey / red,
## interpolated in Lab space so the transition is perceptually even
hm_cols <- colorRampPalette(c(
  "#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2",
  "#F2CBB7", "#EE8468", "#B40426"
), space = "Lab")(100)
## sample correlation heatmaps: sequential (all values positive), so they
## never look like a DE heatmap
cor_cols <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(50)

## annotated vs novel lncRNA: isoformic pink / magenta (see cat_cols below)
class_cols <- c(lncRNA = "#E78AC3", novel_lncRNA = "#AE017E")

## feature categories (DE dot plot, DE composition); protein-coding / other in
## greys so the lncRNA and novel classes stand out; lncRNA colours = class_cols
cat_levels <- c(
  "Annotated protein-coding", "Annotated lncRNA", "Annotated other",
  "Novel protein-coding", "Novel lncRNA", "Novel non-coding"
)
## the biotype palette of every figure (as make_palette in dong/06_isoformic.R);
## "other" in light grey
cat_cols <- setNames(c(
  "#A6D854", "#E78AC3", "#C9C9C9",
  "#238B45", "#AE017E", "#6A51A3"
), cat_levels)

## annotated x protein-coding / lncRNA / other;
## novel x protein-coding / lncRNA / non-coding (pulposeq classes; every
## novel feature has one, so anything else stops the run)
feature_category <- function(biotype, novel) {
  type <- dplyr::case_when(
    biotype %in% c("protein_coding", "novel_protein_coding") ~ "protein-coding",
    biotype %in% c("lncRNA", "novel_lncRNA") ~ "lncRNA",
    biotype == "novel_non_coding" ~ "non-coding",
    TRUE ~ "other"
  )
  out <- factor(paste(ifelse(novel, "Novel", "Annotated"), type), levels = cat_levels)
  if (anyNA(out)) stop(sum(is.na(out)), " novel feature(s) without a pulposeq biotype")
  out
}

## ggrepel labels: white box, thin black border, leader arrow to the point.
## data = NULL inherits the plot data; extra arguments go to geom_label_repel
repel_labels <- function(data, mapping, size = 2.5, box.padding = 0.5,
                         label.padding = unit(0.12, "lines"), ...) {
  geom_label_repel(
    data = data, mapping = mapping, size = size, colour = "black",
    fill = "white", label.size = 0.2, label.r = unit(0.05, "lines"),
    label.padding = label.padding, box.padding = box.padding,
    point.padding = 0.2, min.segment.length = 0, max.overlaps = Inf,
    segment.colour = "black", segment.size = 0.3,
    arrow = arrow(length = unit(0.012, "npc"), type = "open"),
    show.legend = FALSE, seed = 1, ...
  )
}
