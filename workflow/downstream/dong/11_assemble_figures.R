## ===========================================================================
## 11_assemble_figures.R - assemble Extended Data Fig. 4, Fig. 2 and Extended Data Fig. 7
##
##   Extended Data Fig. 4:
##            a-h LongBench quantifier comparison (layout saved by
##            longbench/09_assemble_figures.R, panels from
##            longbench/03_quantifier_plots.R), then i-l the Dong
##            quantifier panels of 10 in one row filling the page width
##            (i sequin counts, j divided vs undivided, k Jaccard, l disagreement)
##   Fig. 2 (middle): a-r; rows 1-2 from
##            helpers/fig2_rows.R, row 3 aligned on 06's fig2_profile_rooms.rds
##   Extended Data Fig. 7 (bottom): the novel lncRNAs in cis pairs (05)
## 5 mm margins, 4 mm gutters; panel letters 12 pt, plain, in the gutter above
## each panel's top-left corner. Slot sizes here must match the sizes the panel
## scripts save (a mismatch is reported).
##
## Needs (run first): 10_quantifier_panels.R and longbench/09_assemble_figures.R
## (ED4); for Fig. 2 also 04, 06 and 08
## Output (in <project>/figures/, next to the vector-editor files, which are
## never touched): "Extended data figure 4".{svg,png}, fig2_assembled.{svg,png},
## "Extended data figure 7".{svg,png}
## Needs: magick
## ===========================================================================

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")
source(file.path(repo, "downstream/longbench/helpers/compose_figure.R"))   # compose_figure, layout_rows, check_page
out_dir <- file.path(project, "figures")
dir.create(out_dir, showWarnings = FALSE)

## ------------------------------------------------ Extended Data Fig. 4 ----
## LongBench a-h as longbench/09 laid them out (340 x 347 mm), then the Dong row i-l:
## 10 saves the four with one top / bottom room and boxes sized so the row
## fills 340 - 2 x 5 mm; packed from the left with 4 mm gutters.
lb  <- readRDS(file.path(repo, "downstream/longbench/spike-ins/results/quantifier_comparison/ed1_longbench_layout.rds"))
ps  <- readRDS(file.path(res_dir, "quantifier_comparison", "ed_panel_sizes_mm.rds"))
if (length(ps$rows) != 1) stop("10 saved ", length(ps$rows), " ED rows; ED4 wants the Dong panels in one row: rerun 10")
gut <- 4; mar <- 5
ed1_stage <- file.path(fig_dir, "ed1_assembly_panels")
unlink(ed1_stage, recursive = TRUE); dir.create(ed1_stage, showWarnings = FALSE)
lb_panel_dir <- file.path(repo, "downstream/longbench/spike-ins/figures/quantifier_no_titles")
for (f in lb$lay$file[!is.na(lb$lay$file)]) for (ext in c("svg", "png")) {
  src <- file.path(lb_panel_dir, paste0(f, ".", ext))
  if (!file.exists(src)) stop("LongBench ED4 panel not found: ", src, " (rerun longbench/03_quantifier_plots.R)")
  if (!isTRUE(file.copy(src, ed1_stage, overwrite = TRUE)))
    stop("could not copy LongBench ED4 panel: ", src)
}
rw <- ps$rows[[1]]
for (f in rw) for (ext in c("svg", "png")) {
  src <- file.path(fig_dir, "ed_quantifier_panels_no_titles", paste0(f, ".", ext))
  if (!file.exists(src)) stop("Dong ED4 panel not found: ", src, " (rerun dong/10_quantifier_panels.R)")
  if (!isTRUE(file.copy(src, ed1_stage, overwrite = TRUE)))
    stop("could not copy Dong ED4 panel: ", src)
}
w  <- unname(sapply(ps$sizes[rw], `[`, 1)); h <- unname(ps$sizes[[rw[1]]][2])
y5 <- lb$page[["h"]] - mar + gut                          # under the LongBench rows
x  <- mar + c(0, cumsum(w + gut))[seq_along(rw)]
dong_row <- data.frame(letter = letters[nrow(lb$lay) + seq_along(rw)], file = rw, x = x, y = y5, w = w, h = h)
if (max(dong_row$x + dong_row$w) > lb$page[["w"]] - mar + 0.3)
  stop(sprintf("ED4 Dong row is %.1f mm past the right margin: rerun 10",max(dong_row$x + dong_row$w) - (lb$page[["w"]] - mar)))
ed1 <- rbind(lb$lay, dong_row)
ed1_page <- c(w = lb$page[["w"]], h = y5 + h + mar)
print(ed1, row.names = FALSE, digits = 4)
compose_figure(ed1, ed1_page, ed1_stage, file.path(out_dir, "Extended data figure 4"))
message("Done: Extended Data Fig. 4 in ", out_dir)

## ============================================================== Fig. 2 =====
## 340 mm wide, 5 mm margins, 4 mm gutters; height follows the content.
##   row 1  a read support | b novel transcripts by class code | c sequin DE Bambu | oarfish |
##          d sequin length ratio
##   row 2  e human counts | f overdispersion by category | g overdispersion by
##          isoforms per gene | h sample correlation
##          rows 1-2: one square box size, aligned boxes, 8 pt (helpers/fig2_rows.R;
##          08 and 10 save them at these sizes, checked below)
##   row 3  i volcano | j EGFR | k UBE2I | l RAB11B-AS1 | m LY6D profiles (06, 7e):
##          plot boxes share top and bottom edges (fig2_profile_rooms.rds; 08 saves
##          the volcano with the profiles' box height in the width they leave)
##   then   n DTE heatmap (left) | o EGFR locus | p LY6D locus;
##          under the loci: q cis-pair heatmap (under o) | r circos (right margin)
## MUC1 is not in Fig. 2 (06 still writes its profile).
## Every panel at its own saved size (no scaling here, so the font sizes stay as drawn).
source(file.path(helper_dir, "plot_theme.R"))        # font_family
source(file.path(helper_dir, "figure_helpers.R"))    # save_panel
source(file.path(helper_dir, "fig2_rows.R"))         # fig2_layout (rows 1-2), fig2_row3_profiles, fig2_loci
stage <- file.path(fig_dir, "fig2_assembly_panels")
unlink(stage, recursive = TRUE); dir.create(stage, showWarnings = FALSE)
pn  <- file.path(fig_dir, "fig2_panels_no_titles")
pq  <- file.path(fig_dir, "fig2_quantifier_panels")
isp <- file.path(iso_paper, "fig2_panels")
loc <- function(gene, ext) file.path(iso_paper, paste0("locus_tracks_", toupper(ext)), paste0(gene, ".", ext))
src <- c(
  a = file.path(pn, "fig2a_read_support"),
  b = file.path(pn, "fig2b_novel_biotype"),
  c = file.path(pq, "ed_q_b_sequin_de"),
  d = file.path(pq, "ed_q_a_sequin_length"),
  e = file.path(pq, "ed_q_c_human_counts"),
  f = file.path(pq, "ed_q_f_overdispersion"),
  g = file.path(pq, "ed_q_g_overdispersion_isoforms"),
  h = file.path(pn, "fig2e_sample_correlation"),
  i = file.path(pn, "fig2f_volcano"),
  setNames(file.path(isp, fig2_row3_profiles), c("j", "k", "l", "m")),
  n = file.path(pn, "fig2g_isoform_heatmap"),
  o = fig2_loci[1], p = fig2_loci[2],                         # loci: locus_tracks_{SVG,PNG}/<gene>.* (06)
  q = file.path(cis_dir, "paper_cis_heatmap"),  # q heatmap, r circos
  r = file.path(cis_dir, "paper_cis_circos"))
for (k in names(src)) for (ext in c("svg", "png")) {
  f <- if (k %in% c("o", "p")) loc(src[[k]], ext) else paste0(src[[k]], ".", ext)
  if (!file.exists(f)) stop("Fig. 2", k, ": panel not found (run its script first; profiles and loci: 04): ", f)
  file.copy(f, file.path(stage, paste0("fig2", k, ".", ext)), overwrite = TRUE)
}
## size (mm) of a panel file, from its SVG root (svglite: pt)
svg_mm <- function(f) {
  l <- readLines(f, n = 5, warn = FALSE)
  root <- regmatches(l, regexpr("<svg[^>]*>", l))[1]
  as.numeric(regmatches(root, gregexpr("(?<=(width|height)=['\"])[0-9.]+(?=pt)", root, perl = TRUE))[[1]]) * 25.4 / 72
}
ids <- names(src)
sz  <- lapply(setNames(ids, ids), function(k) svg_mm(file.path(stage, paste0("fig2", k, ".svg"))))
mar <- fig2_mar; gut <- fig2_gut; W <- fig2_W

## rows 1-2: the files must have the sizes of the shared layout (else a script ran on old rooms)
lay <- fig2_layout()
for (k in letters[1:8]) {
  nm <- basename(src[[k]]); want <- lay$size[[nm]]
  if (any(abs(sz[[k]] - want) > 0.3))
    stop(sprintf("Fig. 2%s (%s) is %.1f x %.1f mm, the rows 1-2 layout wants %.1f x %.1f: rerun %s, then this script",
                 k, nm, sz[[k]][1], sz[[k]][2], want[1], want[2], if (grepl("^ed_q", nm)) "10" else "08"))
}
## row 3: profiles as 06 measured them, the volcano as 08 sized it from those
pr_file <- file.path(isp, "fig2_profile_rooms.rds")
if (!file.exists(pr_file)) stop("no ", pr_file, ": run 04 (then 06) first")
pr <- readRDS(pr_file)[fig2_row3_profiles]
if (any(sapply(pr, is.null))) stop("fig2_profile_rooms.rds lacks ",
                                   paste(fig2_row3_profiles[sapply(pr, is.null)], collapse = ", "), ": rerun 04")
for (q in seq_along(pr)) {
  k <- c("j", "k", "l", "m")[q]
  if (any(abs(sz[[k]] - pr[[q]][c("W", "H")]) > 0.3))
    stop(sprintf("Fig. 2%s (%s) is %.1f x %.1f mm, fig2_profile_rooms.rds says %.1f x %.1f: rerun 04",
                 k, fig2_row3_profiles[q], sz[[k]][1], sz[[k]][2], pr[[q]][["W"]], pr[[q]][["H"]]))
}
Tp <- max(sapply(pr, `[[`, "T")); Bp <- max(sapply(pr, `[[`, "B")); bh <- pr[[1]][["box_h"]]
## 06's shared legend of the profiles, after panel m (no letter)
lg <- readRDS(pr_file)[["fig2_profile_legend"]]
if (!is.null(lg)) for (ext in c("svg", "png")) {
  f <- file.path(isp, paste0("fig2_profile_legend.", ext))
  if (!file.exists(f)) stop("Fig. 2 row-3 legend not found: ", f, " (rerun 06)")
  file.copy(f, file.path(stage, paste0("fig2legend.", ext)), overwrite = TRUE)
}
W3 <- W - 2 * mar - gut * length(pr) - sum(sapply(pr, `[[`, "W")) - (if (is.null(lg)) 0 else lg[["W"]] + gut)
if (any(abs(sz$i - c(W3, Tp + bh + Bp)) > 0.3))
  stop(sprintf("Fig. 2i (volcano) is %.1f x %.1f mm, row 3 wants %.1f x %.1f: rerun 06 (after 04), then this script",
               sz$i[1], sz$i[2], W3, Tp + bh + Bp))
## lower block: the two loci beside the DTE heatmap must fit the page
if (sz$n[1] + sz$o[1] + sz$p[1] + 2 * gut > W - 2 * mar + 0.3)
  stop(sprintf("Fig. 2 n + o + p are %.1f mm wide, the page has %.1f: rerun 04 (paper_locus_tracks_mm = 119.5)",
               sz$n[1] + sz$o[1] + sz$p[1] + 2 * gut, W - 2 * mar))

pl <- function(letter, x, y, ly = y) data.frame(letter = letter, file = paste0("fig2", letter), x = x, y = y,
                                                w = sz[[letter]][1], h = sz[[letter]][2], ly = ly)
row <- function(ks, y) {                                     # left-packed row
  w <- vapply(ks, function(k) sz[[k]][1], 0)
  x <- mar + c(0, cumsum(w + gut))[seq_along(ks)]
  do.call(rbind, lapply(seq_along(ks), function(q) pl(ks[q], x[q], y)))
}
y2 <- mar + sz$a[2] + gut
y3 <- y2 + sz$e[2] + gut
## row 3: volcano box top at Tp below y3; each profile file shifted so its box top is there too;
## the letters stay on one line (y3)
r3 <- pl("i", mar, y3); x <- mar + sz$i[1] + gut
for (q in seq_along(pr)) {
  k  <- c("j", "k", "l", "m")[q]
  r3 <- rbind(r3, pl(k, x, y3 + Tp - pr[[q]][["T"]], ly = y3)); x <- x + sz[[k]][1] + gut
}
if (!is.null(lg)) {                                          # legend top on the profiles' box top
  lw <- svg_mm(file.path(stage, "fig2legend.svg"))
  r3 <- rbind(r3, data.frame(letter = "", file = "fig2legend", x = x, y = y3 + Tp, w = lw[1], h = lw[2], ly = y3))
}
y4 <- max(r3$y + r3$h) + gut
xo <- mar + sz$n[1] + gut; xp <- xo + sz$o[1] + gut
y5 <- y4 + max(sz$o[2], sz$p[2]) + gut                       # under the taller locus
fig2 <- rbind(row(c("a", "b", "c", "d"), mar), row(c("e", "f", "g", "h"), y2), r3,
              pl("n", mar, y4), pl("o", xo, y4), pl("p", xp, y4),
              pl("q", xo, y5), pl("r", W - mar - sz$r[1], y5))
fig2 <- fig2[order(fig2$letter), ]
fig2_page <- c(w = W, h = max(fig2$y + fig2$h) + mar)
## checks: inside the margins; no overlaps
if (any(fig2$x + fig2$w > W - mar + 0.3))
  stop("Fig. 2 panels past the right margin: ",
       paste(sprintf("%s by %.1f mm", fig2$letter, (fig2$x + fig2$w) - (W - mar))[fig2$x + fig2$w > W - mar + 0.3], collapse = ", "))
pairs <- combn(nrow(fig2), 2)
ov <- apply(pairs, 2, function(ij) {
  a <- fig2[ij[1], ]; b <- fig2[ij[2], ]
  a$x < b$x + b$w && b$x < a$x + a$w && a$y < b$y + b$h && b$y < a$y + a$h
})
if (any(ov)) stop("overlapping slots in the Fig. 2 layout: ",
                  paste(apply(pairs[, ov, drop = FALSE], 2, function(ij) paste(fig2$letter[ij], collapse = "/")), collapse = ", "))
print(fig2, row.names = FALSE, digits = 4)
message(sprintf("Fig. 2 page: %.0f x %.0f mm; rows 1-2 boxes %.1f mm, a's bars %.1f mm; row 3 boxes %.1f mm tall, volcano %.1f mm wide",
                fig2_page[["w"]], fig2_page[["h"]], lay$box2, lay$a_w, bh, W3))
compose_figure(fig2, fig2_page, stage, file.path(out_dir, "fig2_assembled"))
message("Done: Fig. 2 in ", out_dir)

## ------------------ Extended Data figure 7: the novel lncRNAs in cis pairs ------
## 05's circos (the 10 differentially
## expressed novel lncRNAs with the strongest correlation to a differentially
## expressed reference transcript of another gene) beside its heatmap.
## The files are copied, nothing is redrawn.
ed5_dir <- file.path(fig_dir, "ed4_assembly_panels")
unlink(file.path(fig_dir, c("ed5_assembly_panels", "ed6_assembly_panels")), recursive = TRUE)
dir.create(ed5_dir, showWarnings = FALSE)
lnc_dir <- file.path(dirname(cis_dir), "results_lncRNA_top10")
for (nm in c("paper_cis_circos", "paper_cis_heatmap")) for (ext in c("svg", "png"))
  file.copy(file.path(lnc_dir, paste0(nm, ".", ext)), file.path(ed5_dir, paste0(nm, ".", ext)), overwrite = TRUE)
s  <- function(letter, file, x, y, w, h) data.frame(letter = letter, file = file, x = x, y = y, w = w, h = h)
c5 <- svg_mm(file.path(ed5_dir, "paper_cis_circos.svg")); hm5 <- svg_mm(file.path(ed5_dir, "paper_cis_heatmap.svg"))
ed5 <- rbind(s("a", "paper_cis_circos", 5, 5, c5[1], c5[2]),
             s("b", "paper_cis_heatmap", 5 + c5[1] + 4, 5, hm5[1], hm5[2]))
compose_figure(ed5, c(w = 5 + c5[1] + 4 + hm5[1] + 5, h = 5 + max(c5[2], hm5[2]) + 5), ed5_dir,
               file.path(out_dir, "Extended data figure 7"))
