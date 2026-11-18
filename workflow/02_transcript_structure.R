## ===========================================================================
## 02_transcript_structure.R - transcript structure of the final transcriptome by
## biotype, across the LongBench protocols (ONT cDNA, ONT DRS, PacBio) and the
## Dong case study (ONT cDNA). Extended Data Fig. 3.
##
##   a  transcript length   } multi-exonic (solid) and mono-exonic (pale) side by side
##   c  exon length         } within each biotype: a wide translucent violin, a narrow
##                            opaque one and a white notched box; groups of < 10 as points
##   b  median number of exons per biotype, one line per dataset
##   d  transcript length against number of exons: 2-D density (as in the pulposeq
##      report, fig-density / fig-density-novel), one row per dataset, one column
##      per category; density scaled to its maximum within each panel
##
## Categories: the reference transcripts kept in the final transcriptome
## (protein-coding, lncRNA; pulposeq's annotated_{protein-coding,lncRNAs} sets)
## and the validated novel transcripts by pulposeq biotype (novel protein-coding,
## novel lncRNA, novel non-coding), in the biotype palette of every figure.
##
## Inputs (read only), for each run folder <project>/{ont_cdna,ont_drs,pb,dong}:
##   transcriptome/curated/metadata/annotated_{protein-coding,lncRNAs}_metadata.csv
##       (transcript_length, num_exons)
##   transcriptome/curated/metadata/annotated_{protein-coding,lncRNAs}_exonlength.csv (width)
##   transcriptome/curated/metadata/novel_transcripts_validated_metadata.csv
##       (len, num_exons, pulposeq_transcript_biotype)
##   transcriptome/curated/annotation/novel_transcripts_validated.gtf (exons)
## Outputs: transcript-structure/results/
##            transcript_structure_per_transcript.csv.gz   (run, category, id, length, exons)
##            transcript_structure_summary.csv             (n, min, quartiles, max per metric)
##            source_data_ed6.rds                          (figure data: a-c summaries, d per transcript)
##          <project>/figures/"Extended data figure 3".{pdf,svg,png}
## Needs: tidyverse, data.table, patchwork, svglite, showtext
## ===========================================================================
suppressPackageStartupMessages({ library(tidyverse); library(data.table); library(patchwork) })
source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
source(file.path(repo, "downstream/longbench/helpers/plot_theme.R"))      # font_family
source(file.path(repo, "downstream/longbench/helpers/figure_helpers.R"))  # save_panel
out_dir <- file.path(repo, "downstream/transcript-structure/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- c(ont_cdna = "LongBench\nONT cDNA", ont_drs = "LongBench\nONT DRS",
          pb = "LongBench\nPacBio", dong = "Dong\nONT cDNA")
cat_levels <- c("protein-coding", "lncRNA", "novel protein-coding", "novel lncRNA", "novel non-coding")
## the biotype palette of every figure (as make_palette in dong/06_isoformic.R)
cat_cols   <- c("protein-coding" = "#A6D854", "lncRNA" = "#E78AC3", "novel protein-coding" = "#238B45",
                "novel lncRNA" = "#AE017E", "novel non-coding" = "#6A51A3")  # biotype palette
novel_lab <- c(novel_protein_coding = "novel protein-coding", novel_lncRNA = "novel lncRNA",
               novel_non_coding = "novel non-coding")

## ------------------------------------------------------------- read -------
tx <- list(); ex <- list()
for (r in names(runs)) {
  md <- file.path(project, r, "transcriptome/curated/metadata")
  an <- file.path(project, r, "transcriptome/curated/annotation")
  for (b in c(`protein-coding` = "protein-coding", lncRNA = "lncRNAs")) {
    cat_b <- if (b == "lncRNAs") "lncRNA" else "protein-coding"
    m <- fread(file.path(md, paste0("annotated_", b, "_metadata.csv")),
               select = c("transcript_id_version", "transcript_length", "num_exons"))
    tx[[length(tx) + 1]] <- data.table(run = r, category = cat_b, id = m$transcript_id_version,
                                       length = m$transcript_length, exons = m$num_exons)
    e <- fread(file.path(md, paste0("annotated_", b, "_exonlength.csv")), select = c("transcript_id_version", "width"))
    ex[[length(ex) + 1]] <- data.table(run = r, category = cat_b, id = e$transcript_id_version, width = e$width)
  }
  nv <- fread(file.path(md, "novel_transcripts_validated_metadata.csv"),
              select = c("qry_id", "len", "num_exons", "pulposeq_transcript_biotype"))
  unk <- setdiff(unique(nv$pulposeq_transcript_biotype), names(novel_lab))
  if (length(unk)) message("  ! ", r, ": novel biotypes not shown: ", paste(unk, collapse = ", "))
  nv <- nv[pulposeq_transcript_biotype %in% names(novel_lab)]
  tx[[length(tx) + 1]] <- data.table(run = r, category = novel_lab[nv$pulposeq_transcript_biotype], id = nv$qry_id,
                                     length = nv$len, exons = nv$num_exons)
  g <- fread(file.path(an, "novel_transcripts_validated.gtf"), header = FALSE, sep = "\t", quote = "",
             col.names = c("chr", "src", "feature", "start", "end", "score", "strand", "frame", "attr"))[feature == "exon"]
  ## the model's own transcript_id: a greedy '.*transcript_id "' would capture the
  ## last attribute ending in transcript_id (e.g. compared_transcript_id)
  g[, id := sub('^(?:.*?; *)?transcript_id "([^"]+)".*$', "\\1", attr, perl = TRUE)]
  if (any(!grepl("^Bambu", g$id))) stop(r, ": novel GTF transcript_id not parsed (", head(g$id[!grepl("^Bambu", g$id)], 1), ")")
  g <- g[id %in% nv$qry_id]
  g[, category := novel_lab[nv$pulposeq_transcript_biotype[match(id, nv$qry_id)]]]
  ex[[length(ex) + 1]] <- data.table(run = r, category = g$category, id = g$id, width = g$end - g$start + 1)
  ## the GTF's exons must agree with the metadata's exon counts
  k <- g[, .N, by = id]; bad <- sum(k$N != nv$num_exons[match(k$id, nv$qry_id)])
  if (bad) message(sprintf("  ! %s: %d novel transcripts with a GTF exon count different from num_exons", r, bad))
  message(sprintf("  %s: %s reference + %s novel transcripts, %s exons", r,
                  format(sum(sapply(tx[length(tx) - 2:1], nrow)), big.mark = ","), format(nrow(nv), big.mark = ","),
                  format(sum(sapply(ex[length(ex) - 2:0], nrow)), big.mark = ",")))
}
tx <- rbindlist(tx); ex <- rbindlist(ex)
lvl <- function(d) d[, `:=`(run = factor(runs[run], levels = runs), category = factor(category, levels = cat_levels))]
lvl(tx); lvl(ex)
## exon length separately for multi- and mono-exonic transcripts, side by side
## within each biotype (a mono-exonic transcript's only exon is the transcript)
ex[, structure := factor(ifelse(.N == 1, "mono-exonic", "multi-exonic"),
                         levels = c("multi-exonic", "mono-exonic")), by = .(run, id)]
tx[, structure := factor(ifelse(exons == 1, "mono-exonic", "multi-exonic"), levels = c("multi-exonic", "mono-exonic"))]
fwrite(tx, file.path(out_dir, "transcript_structure_per_transcript.csv.gz"))

## ---------------------------------------------------------- summaries ------
summ <- function(d, v, metric, by = c("run", "category")) {
  d <- d[, c(mget(by), list(x = as.numeric(get(v))))]
  d[, .(metric = metric, n = .N, min = min(x), q25 = quantile(x, 0.25, names = FALSE),
        median = median(x), q75 = quantile(x, 0.75, names = FALSE), max = max(x)), by = by]
}
sm <- rbind(summ(tx, "length", "transcript length (nt)", by = c("run", "category", "structure")), summ(tx, "exons", "number of exons"),
            summ(ex, "width", "exon length (nt)", by = c("run", "category", "structure")), fill = TRUE)
sm[, run := gsub("\n", " ", run)]
setorder(sm, metric, run, category)
fwrite(sm, file.path(out_dir, "transcript_structure_summary.csv"))
print(sm[, .(metric, run, category, n, median)], nrows = 100)

## ------------------------------------------------------------- plots -------
base <- 8
## boxed facets with grey strips, as in d (theme_bw)
th <- theme_bw(base_size = base, base_family = font_family) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linetype = 2, linewidth = 0.2),
        axis.text = element_text(size = base, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey96"), strip.text = element_text(size = base),
        plot.tag = element_text(size = 12, face = "plain"))
nt_lab <- function(x) ifelse(x >= 1000, paste0(signif(x / 1000, 3), " kb"), paste0(round(x), " nt"))

## two violins, a white notched box, max (diamond) and min (circle)
violins <- function(d, v, ylab, lab_fun) {
  d  <- d[, .(run, category, value = get(v))]
  mm <- rbind(d[, .(stat = "max", value = max(value)), by = .(run, category)],
              d[, .(stat = "min", value = min(value)), by = .(run, category)])
  ggplot(d, aes(category, value, fill = category, colour = category)) +
    geom_violin(width = 1, alpha = 0.5, colour = "white") +
    geom_violin(width = 0.5, alpha = 1, colour = "white") +
    geom_boxplot(width = 0.08, colour = "grey25", fill = "white", linewidth = 0.25,
                 outliers = FALSE, notch = TRUE, notchwidth = 0.2) +
    geom_point(data = mm[stat == "max"], shape = 23, size = 1.6, fill = "white", stroke = 0.6) +
    geom_point(data = mm[stat == "min"], shape = 21, size = 1.2, fill = "white", stroke = 0.6) +
    geom_text(data = mm, aes(label = lab_fun(value)), hjust = -0.35, vjust = 0.5, size = 6 / .pt,
              colour = "black", family = font_family) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_fill_manual(values = cat_cols) + scale_colour_manual(values = cat_cols) +
    facet_wrap(~ run, nrow = 1) +
    labs(x = NULL, y = ylab) + th
}
## a and c: multi-exonic (solid) and mono-exonic (pale) transcripts side by side
## within each biotype, dodged. Groups of < 10 transcripts are drawn as points, not
## violins; a biotype with one structure keeps the same violin width. Extremes are
## unlabelled (see the summary table).
dg <- position_dodge(width = 0.9, preserve = "single")
split_violins <- function(d, ylab, legend = FALSE) {
  d <- copy(d); d[, n := .N, by = .(run, category, structure)]
  mm <- rbind(d[n >= 10, .(stat = "max", value = max(value)), by = .(run, category, structure)],
              d[n >= 10, .(stat = "min", value = min(value)), by = .(run, category, structure)])
  ggplot(d[n >= 10], aes(category, value, fill = category, colour = category, group = interaction(category, structure))) +
    geom_violin(aes(alpha = structure), position = dg, width = 0.85, colour = "white", scale = "width") +
    geom_boxplot(position = dg, width = 0.12, colour = "grey25", fill = "white", linewidth = 0.25,
                 outliers = FALSE, notch = TRUE, notchwidth = 0.2) +
    geom_point(data = mm[stat == "max"], position = dg, shape = 23, size = 1.3, fill = "white", stroke = 0.5) +
    geom_point(data = mm[stat == "min"], position = dg, shape = 21, size = 1.0, fill = "white", stroke = 0.5) +
    geom_point(data = d[n < 10], aes(alpha = structure), position = dg, shape = 16, size = 0.9) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_fill_manual(values = cat_cols, guide = "none") + scale_colour_manual(values = cat_cols, guide = "none") +
    scale_alpha_manual(values = c("multi-exonic" = 1, "mono-exonic" = 0.35), name = NULL,
                       guide = if (legend) guide_legend(override.aes = list(fill = "grey40", colour = NA)) else "none") +
    scale_x_discrete(drop = FALSE) +
    facet_wrap(~ run, nrow = 1) +
    labs(x = NULL, y = ylab) + th +
    theme(legend.position = "top", legend.justification = "left", legend.key.size = unit(3, "mm"),
          legend.margin = margin(0, 0, 0, 0))
}
pa <- split_violins(tx[, .(run, category, structure, value = length)], "Transcript length (nt)", legend = TRUE)
pc <- split_violins(ex[, .(run, category, structure, value = width)], "Exon length (nt)")

## b: median number of exons per biotype, one line per dataset
run_cols <- setNames(c(protocol_cols[["ONT cDNA"]], protocol_cols[["ONT DRS"]], protocol_cols[["PacBio"]], "#27408B"),
                     runs)                                                     # as Extended Data Figs. 1 and 2
med_ex <- tx[, .(median = median(exons)), by = .(run, category)]
pb <- ggplot(med_ex, aes(category, median, colour = run, group = run)) +
  geom_line(linewidth = 0.5, position = position_dodge(width = 0.25)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.25)) +
  scale_colour_manual(values = run_cols, labels = function(x) gsub("\n", " ", x), name = NULL,
                      guide = guide_legend(ncol = 2)) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_pretty(), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Median number of exons") +
  theme_classic(base_size = base, base_family = font_family) +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = 2, linewidth = 0.2),
        axis.text = element_text(size = base, colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(3, "mm"), plot.tag = element_text(size = 12, face = "plain"))
print(dcast(copy(med_ex)[, run := gsub("\n", " ", run)], category ~ run, value.var = "median"))

## d: as the pulposeq report (stat_density2d polygons, log10 axes); density
## scaled to 1 at its maximum in each panel so the colours compare across panels
pd <- ggplot(tx, aes(exons, length)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour_var = "ndensity", bins = 8) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100)) +
  scale_y_log10(breaks = c(300, 1000, 3000, 10000, 30000), labels = scales::label_comma()) +
  scale_fill_viridis_c(name = "Relative\ndensity", limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1)) +
  facet_grid(run ~ category) +
  labs(x = "Number of exons (log10 scale)", y = "Transcript length (nt, log10 scale)") +
  theme_bw(base_size = base, base_family = font_family) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey96"),
        strip.text = element_text(size = base), axis.text = element_text(size = base - 1, colour = "black"),
        legend.title = element_text(size = base), legend.text = element_text(size = base - 1),
        plot.tag = element_text(size = 12, face = "plain"))

## a on top; b (narrow) to the left of c; d below. Blank spacers (untagged)
## separate b from c and the three rows from each other.
ed6 <- (pa / plot_spacer() / ((pb | plot_spacer() | pc) + plot_layout(widths = c(1, 0.12, 3.2))) /
          plot_spacer() / pd) +
  plot_layout(heights = c(1, 0.1, 1.1, 0.1, 2.2)) +
  plot_annotation(tag_levels = "a")
page <- c(260, 257)  # width, height (mm); a / b|c / d plus the two row gaps (d shortened)
save_panel(ed6, "Extended data figure 3", page, file.path(project, "figures"), dpi = 300)

## figure data, kept for sharing: a-c summaries, d per transcript
source(file.path(repo, "downstream/helpers/source_data.R"))
sm_by <- split(sm, sm$metric)
saveRDS(list(ED6a = list(sm_by[["transcript length (nt)"]]),
             ED6b = list(sm_by[["exon length (nt)"]]),
             ED6c = list(as.data.frame(tx[, .(run = gsub("\n", " ", run), category, transcript_id = id,
                                              transcript_length = length, num_exons = exons)]))),
        file.path(out_dir, "source_data_ed6.rds"))
message("Done: ", file.path(project, "figures", "Extended data figure 3"), " and ", out_dir)
