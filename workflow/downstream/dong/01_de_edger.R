## ===========================================================================
## 01_de_edger.R - differential expression, H1975 (TKI_resistant) vs HCC827
##                 (TKI_sensitive), gene and isoform level, edgeR
## Counts: oarfish's, divided by their
## read-to-transcript ambiguity (RTA) overdispersion estimated from the bootstrap
## replicates (downstream/helpers/oarfish_counts.R):
##   transcript level  edgeR::catchOarfish()     -> counts / overdispersion
##   gene level        edgeR::catchOarfishGene() -> counts / gene overdispersion
## Divided counts are not rounded. The expression filter is applied to the
## undivided read counts; edgeR tests the divided counts. Every DE table gains an
## `overdispersion` column. Needs edgeR >= 4.99.12.
## edgeR v4: drop rRNA / processed pseudogenes, keep >= 10 reads in >= 3
## samples, TMM, glmQLFit(robust) -> glmTreat(lfc = 0.5), FDR 0.05; display
## log2FC = predFC(prior.count = 2).
## Outputs used downstream:
##   de_edger/all_results.csv                      volcano + heatmap
##   de_edger/transcript/sample_spearman.csv       sample correlation: Spearman
##                                                 on logCPM (prior.count 2) of the tested isoforms
##   de_edger/transcript/logcpm.rds                heatmap, circos (04), isoformic (06)
##   de_edger/{gene,transcript}/de_TKI_resistant_vs_TKI_sensitive.csv   circos, isoformic
## ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)          # >= 4.0 (normLibSizes, new glmQLFit)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
})

## -------------------------------------------------------------- paths -----
## Inputs : the Dong run folder (pulposeq --outdir, data_dir) + data/samplesheet.csv
## Outputs: downstream/dong/diff-exp/results/ (every output of this script)
source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")   # data_dir, sheet, de_dir
source(file.path(shared_helper_dir, "oarfish_counts.R"))              # read_oarfish()
outroot  <- de_dir

source(file.path(helper_dir, "id_helpers.R"))
source(file.path(helper_dir, "plot_theme.R"))

## ------------------------------------------------------------- settings ---
dataset_lab   <- "H1975 vs HCC827"
levels_to_run <- c("gene", "transcript")
grp_ref       <- "TKI_sensitive"          # denominator
grp_num       <- "TKI_resistant"          # numerator
line_lab      <- c(TKI_sensitive = "HCC827", TKI_resistant = "H1975")

alpha       <- 0.05   # FDR
lfc         <- 0.5    # log2FC threshold tested by glmTreat()
min_count   <- 10
min_samples <- 3      # smallest group has 3 replicates
n_both      <- 2      # Spearman "expressed in both lines": >= min_count in >= n_both samples of each line
filter_mode <- "same"  # "same" = rule above (as DESeq2); "filterByExpr" = edgeR default
fc_prior    <- 2      # prior.count for the display log2FC (predFC) and logCPM
spike_re    <- "^(R[12]_)"
## removed before filtering / normalisation: rRNA is not poly(A) (its level
## tracks library prep), and processed pseudogenes are intronless copies of
## parent mRNAs that take the parent's reads (e.g. RPS28P7 >> RPS28)
drop_biotypes <- c("rRNA", "rRNA_pseudogene", "Mt_rRNA",
                   "processed_pseudogene", "transcribed_processed_pseudogene",
                   "translated_processed_pseudogene")
n_top       <- 10     # heatmap: top n up + top n down
n_panel     <- 5      # top-5 panels
n_dot       <- 15     # dot plot: top n up + top n down

lnc_classes <- c("lncRNA", "novel_lncRNA")
## annotated vs novel lncRNA: isoformic pink / magenta (see cat_cols below)
class_cols  <- c("lncRNA" = "#E78AC3", "novel_lncRNA" = "#AE017E")

## feature categories for the dot plot and the composition plot
cat_levels <- c("Annotated protein-coding", "Annotated lncRNA", "Annotated other",
                "Novel protein-coding", "Novel lncRNA", "Novel non-coding")
## the biotype palette of every figure (as make_palette in 06_isoformic.R);
## "other" in light grey
cat_cols   <- setNames(c("#A6D854", "#E78AC3", "#C9C9C9",
                         "#238B45", "#AE017E", "#6A51A3"),
                       cat_levels)

## DE heatmaps (row z-scores): Moreland "cool-warm" blue / light grey / red,
## interpolated in Lab space so the transition is perceptually even.
## Sample correlation heatmap: YlGnBu (sequential, all values positive), so
## it never looks like a DE heatmap.
hm_cols  <- colorRampPalette(c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2",
                               "#F2CBB7", "#EE8468", "#B40426"), space = "Lab")(100)
cor_cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(50)

## boxed labels with a leader arrow (white box, thin black border)
repel_labels <- function(data, mapping, size = 2.5, ...) {
  geom_label_repel(data = data, mapping = mapping, size = size, colour = "black",
                   fill = "white", label.size = 0.2, label.r = unit(0.05, "lines"),
                   label.padding = unit(0.12, "lines"), box.padding = 0.5,
                   point.padding = 0.2, min.segment.length = 0, max.overlaps = Inf,
                   segment.colour = "black", segment.size = 0.3,
                   arrow = arrow(length = unit(0.012, "npc"), type = "open"),
                   show.legend = FALSE, seed = 1, ...)
}

## group_cols (TKI_sensitive / TKI_resistant) come from helpers/plot_theme.R

coef_name <- paste0("group", grp_num)          # model.matrix(~ group) column
ct_name   <- paste0(grp_num, "_vs_", grp_ref)

## ------------------------------------------------------------ functions ---

## feature ID -> biotype
biotype_map <- function(ref, nov, level) {
  if (level == "transcript") {
    ann <- setNames(as.character(ref$transcript_biotype), ref$transcript_id_version)
    nv  <- setNames(as.character(nov$pulposeq_transcript_biotype), nov$qry_id)
  } else {
    ann <- setNames(as.character(ref$gene_biotype), ref$gene_id_version)
    q   <- as.character(nov$qry_gene_id)
    bt  <- as.character(nov$pulposeq_transcript_biotype)
    bg  <- is_novel(q)
    nv  <- character(0)
    if (any(bg)) {
      pick <- function(x) if ("novel_lncRNA" %in% x) "novel_lncRNA" else
                          if ("novel_protein_coding" %in% x) "novel_protein_coding" else x[1]
      g  <- tapply(bt[bg], q[bg], pick)
      nv <- setNames(as.character(g), names(g))
    }
    extra <- setNames(c(as.character(nov$qry_gene_biotype)[!bg],
                        as.character(nov$compared_gene_biotype)),
                      c(q[!bg], as.character(nov$compared_gene_id)))
    nv <- c(nv, extra[grepl("^ENS", names(extra))])
  }
  m <- c(ann, nv)
  m <- m[!blank(names(m)) & !blank(m)]
  m[!duplicated(names(m))]                 # reference table wins
}

## oarfish counts of the curated set, sequins removed, rows renamed to symbols.
##   counts = divided counts (what edgeR tests), reads = undivided (filter,
##   excluded-biotype totals), od = RTA overdispersion; all in the same row order
read_counts <- function(d) {
  ref <- read.csv(file.path(d, pp$ann_meta))
  nov <- read.csv(file.path(d, pp$novel_meta))
  one_level <- function(lv) {
    x  <- read_oarfish(d, lv, cache_dir = oarfish_cache)
    df <- as.data.frame(x$divided, check.names = FALSE)
    if (lv == "transcript") df$GENEID <- x$gene
    sp <- grepl(spike_re, rownames(df))                       # sequins
    if (lv == "transcript") sp <- sp | grepl(spike_re, df$GENEID)
    df <- df[!sp, , drop = FALSE]
    orig <- rownames(df)
    sw <- swap_ids(df, ref, nov, level = lv)
    sw$GENEID <- NULL
    new <- rownames(sw)
    list(counts  = as.matrix(sw),
         reads   = `rownames<-`(x$counts[orig, , drop = FALSE], new),
         od      = setNames(unname(x$overdispersion[orig]), new),
         id      = setNames(orig, new),
         gid     = if (lv == "transcript") setNames(unname(x$gene[match(orig, rownames(x$counts))]), new),
         biotype = biotype_map(ref, nov, lv))
  }
  list(gene = one_level("gene"), transcript = one_level("transcript"))
}

## sample table from the samplesheet; count columns are "<sample>.bam"
make_coldata <- function(counts, sheet) {
  ss  <- read.csv(sheet)
  smp <- sub("[.]bam$", "", colnames(counts))
  grp <- gsub("-", "_", ss$group[match(smp, ss$sample)])   # TKI-resistant -> TKI_resistant
  if (anyNA(grp) || !all(grp %in% c(grp_ref, grp_num)))
    stop("Count columns not matched to the samplesheet: ",
         paste(smp[is.na(grp) | !grp %in% c(grp_ref, grp_num)], collapse = ", "))
  data.frame(group = factor(grp, levels = c(grp_ref, grp_num)), row.names = smp)
}

## the single contrast: glmQLFit -> glmTreat (H0 |log2FC| <= lfc)
## y: normalised DGEList (no dispersions: glmQLFit v4 estimates its own)
## yd: same with estimateDisp() dispersions, used only by predFC
## ncnt: normalised counts (count scale)
test_contrast <- function(y, yd, design, ncnt, ids, biotype) {
  fit <- glmQLFit(y, design, robust = TRUE)
  tr  <- glmTreat(fit, coef = coef_name, lfc = lfc)
  tt  <- topTags(tr, n = Inf, sort.by = "none")$table

  ## display log2FC: more shrinkage for low counts, used like apeglm
  shr <- predFC(yd, design, prior.count = fc_prior)[, coef_name]

  ## 0 reads in every sample of one group: call is fine, fold change is not
  cnt  <- y$counts
  grp  <- y$samples$group
  zero <- rowSums(cnt[, grp == grp_num, drop = FALSE]) == 0 |
          rowSums(cnt[, grp == grp_ref, drop = FALSE]) == 0

  tab <- tibble(feature     = rownames(tt),
                id          = unname(ids[rownames(tt)]),
                biotype     = unname(biotype[unname(ids[rownames(tt)])]),
                baseMean    = rowMeans(ncnt[rownames(tt), , drop = FALSE]),
                logCPM      = tt$logCPM,
                log2FC      = unname(shr[rownames(tt)]),
                log2FC_test = tt$logFC,
                pvalue      = tt$PValue,
                padj        = tt$FDR,
                zero_in_one_group = unname(zero[rownames(tt)])) %>%
    mutate(biotype = replace_na(biotype, "unknown"),
           novel   = is_novel(id),
           category = feature_category(biotype, novel),
           call = case_when(padj < alpha & log2FC_test > 0 ~ "Up",
                            padj < alpha & log2FC_test < 0 ~ "Down",
                            TRUE                           ~ "NS")) %>%
    arrange(padj)
  list(tab = tab, fit = fit)
}

## BH again within lncRNAs only
lnc_table <- function(tab) {
  tab %>%
    filter(biotype %in% lnc_classes) %>%
    mutate(padj_lnc = p.adjust(pvalue, method = "BH"),
           call_lnc = case_when(padj_lnc < alpha & log2FC_test > 0 ~ "Up",
                                padj_lnc < alpha & log2FC_test < 0 ~ "Down",
                                TRUE                               ~ "NS")) %>%
    arrange(padj_lnc)
}

## annotated x protein-coding / lncRNA / other;
## novel x protein-coding / lncRNA / non-coding (pulposeq classes; every
## novel feature has one, so anything else stops the run)
feature_category <- function(biotype, novel) {
  type <- case_when(biotype %in% c("protein_coding", "novel_protein_coding") ~ "protein-coding",
                    biotype %in% c("lncRNA", "novel_lncRNA")                 ~ "lncRNA",
                    biotype == "novel_non_coding"                            ~ "non-coding",
                    TRUE                                                     ~ "other")
  out <- factor(paste(ifelse(novel, "Novel", "Annotated"), type), levels = cat_levels)
  if (anyNA(out)) stop(sum(is.na(out)), " novel feature(s) without a pulposeq biotype")
  out
}

fmt_p <- function(p) ifelse(p < 1e-300, "< 1e-300", formatC(p, format = "e", digits = 1))

## ---- plots: overall ----

volcano <- function(tab, title) {
  top <- tab %>%
    filter(call != "NS") %>%
    group_by(call) %>%
    slice_min(padj, n = 10, with_ties = FALSE)

  ggplot(tab, aes(log2FC, -log10(pmax(padj, 1e-300)), colour = call)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_hline(yintercept = -log10(alpha), linetype = 2) +
    repel_labels(top, aes(label = feature)) +
    scale_colour_manual(values = c(Up = "darkred", Down = "darkblue", NS = "grey70")) +
    labs(x = "log2 fold change (predFC, display only)", y = "-log10(padj)",
         title = title,
         subtitle = paste0("FDR < ", alpha, ", glmTreat H0: |log2FC| <= ", lfc))
}

top_heatmap <- function(tab, vmat, group, file, title) {
  pick <- function(dirn) tab %>%
    filter(call == dirn) %>%
    arrange(padj, dplyr::desc(abs(log2FC))) %>%
    head(n_top)
  top <- bind_rows(pick("Up"), pick("Down"))

  m <- vmat[top$feature, , drop = FALSE]
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  if (nrow(m) < 2) return(invisible(NULL))

  z   <- t(scale(t(m)))
  lim <- max(abs(z))
  pheatmap(z,
           color  = hm_cols,
           breaks = seq(-lim, lim, length.out = 101),
           cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_col = data.frame(group = group, row.names = colnames(z)),
           annotation_colors = list(group = group_cols),
           border_color = NA, fontfamily = font_family, fontsize = 9, fontsize_row = 8, angle_col = 90,
           main = title, filename = file, width = 6, height = 7)
}

## yd: DGEList with estimateDisp() dispersions; fit: glmQLFit object;
## lcpm: logCPM matrix
qc_plots <- function(yd, fit, lcpm, outdir, title) {
  ## BCV: NB dispersions from estimateDisp (glmQLFit v4 estimates its own)
  y <- yd
  png(file.path(outdir, "bcv.png"), width = 6, height = 5, units = "in", res = 300)
  par(family = font_family)
  plotBCV(yd, main = title)
  dev.off()

  png(file.path(outdir, "ql_dispersion.png"), width = 6, height = 5, units = "in", res = 300)
  par(family = font_family)
  plotQLDisp(fit, main = title)
  dev.off()

  ## MDS (edgeR / limma plotMDS): distance between two samples = root mean
  ## square of the top 500 log2FC between them ("leading log2FC"), picked
  ## separately for each pair; logCPM with prior.count 2
  md <- plotMDS(yd, top = 500, gene.selection = "pairwise", plot = FALSE)
  ve <- if (!is.null(md$var.explained)) round(100 * md$var.explained[1:2], 1) else c(NA, NA)
  pd <- data.frame(dim1 = md$x, dim2 = md$y, name = colnames(yd), group = yd$samples$group)
  p <- ggplot(pd, aes(dim1, dim2, colour = group)) +
    geom_point(size = 3) +
    repel_labels(NULL, aes(label = name), size = 3) +
    scale_colour_manual(values = group_cols, labels = function(x) gsub("_", "-", x),
                        name = NULL) +
    scale_y_continuous(expand = expansion(mult = 0.25)) +
    labs(x = paste0("Leading log2FC dim 1 (", ve[1], "%)"),
         y = paste0("Leading log2FC dim 2 (", ve[2], "%)"),
         title = title, subtitle = "plotMDS, top 500 features per pair (pairwise), logCPM")
  ggsave(file.path(outdir, "mds.png"), p, width = 6.5, height = 5.5, dpi = 300, bg = "white")

  rho <- cor(lcpm, method = "spearman")
  write.csv(rho, file.path(outdir, "sample_spearman.csv"))
  pheatmap(rho,
           color = cor_cols,
           clustering_distance_rows = as.dist(1 - rho),
           clustering_distance_cols = as.dist(1 - rho),
           annotation_col = data.frame(group = y$samples$group, row.names = colnames(rho)),
           annotation_colors = list(group = group_cols),
           border_color = NA, fontfamily = font_family, main = paste("Sample Spearman -", title),
           filename = file.path(outdir, "sample_spearman.png"), width = 6.5, height = 5.5)
}

## ---- plots: lncRNA ----

## y = raw p-value, so lncRNAs and the grey background share one axis
lnc_volcano <- function(full, l, title) {
  ## same layout as the overall volcano, lncRNAs only:
  ## y = -log10(padj_lnc) (BH within lncRNAs), dashed line at alpha;
  ## colour = annotated / novel lncRNA when DE, grey when not
  ## (`full` is kept only so the call in the run block stays unchanged)
  lab_cls <- c(lncRNA = "Annotated lncRNA", novel_lncRNA = "Novel lncRNA")
  l <- l %>%
    mutate(colour = ifelse(call_lnc == "NS", "not DE", lab_cls[biotype]),
           colour = factor(colour, levels = c("not DE", lab_cls))) %>%
    arrange(colour)                                      # DE points on top
  top <- l %>%
    filter(call_lnc != "NS") %>%
    group_by(biotype, call_lnc) %>%
    slice_min(padj_lnc, n = 5, with_ties = FALSE)
  n_de <- table(l$colour)

  ggplot(l, aes(log2FC, -log10(pmax(padj_lnc, 1e-300)), colour = colour)) +
    geom_point(size = 1.1, alpha = 0.7, show.legend = TRUE) +
    geom_hline(yintercept = -log10(alpha), linetype = 2) +
    repel_labels(top, aes(label = feature)) +
    scale_colour_manual(values = c("not DE" = "grey70",
                                   "Annotated lncRNA" = class_cols[["lncRNA"]],
                                   "Novel lncRNA"     = class_cols[["novel_lncRNA"]]),
                        labels = function(x) paste0(x, " (", n_de[x], ")"),
                        name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(x = "log2 fold change (predFC, display only)",
         y = "-log10(padj within lncRNAs)",
         title = title,
         subtitle = paste0("FDR < ", alpha, ", glmTreat H0: |log2FC| <= ", lfc, "; BH within annotated + novel lncRNAs"))
}

lnc_heatmap <- function(l, vmat, group, file, title) {
  top <- l %>%
    filter(call_lnc != "NS", feature %in% rownames(vmat)) %>%
    group_by(call_lnc, biotype) %>%
    slice_min(padj_lnc, n = 5, with_ties = FALSE) %>%   # novel ones always shown
    ungroup()

  m   <- vmat[top$feature, , drop = FALSE]
  ok  <- apply(m, 1, sd) > 0
  m   <- m[ok, , drop = FALSE]
  top <- top[ok, ]
  if (nrow(m) < 2) return(invisible(NULL))

  z   <- t(scale(t(m)))
  lim <- max(abs(z))
  pheatmap(z,
           color  = hm_cols,
           breaks = seq(-lim, lim, length.out = 101),
           cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_row = data.frame(class = top$biotype, row.names = top$feature),
           annotation_col = data.frame(group = group, row.names = colnames(m)),
           annotation_colors = list(class = class_cols, group = group_cols),
           border_color = NA, fontfamily = font_family, fontsize = 9, fontsize_row = 8, angle_col = 90,
           main = title, filename = file, width = 6.5, height = 7)
}

## ---- Spearman: sample pairs ----

## Spearman rho for every pair of samples on log2(normalised count + 1).
## type = "within <line>" (replicate agreement) or "between lines".
## Between-line pairs are oriented x = reference line (HCC827), y = H1975.
sample_pairs <- function(l2, group) {
  s  <- colnames(l2)
  g  <- as.character(group[s])
  pr <- t(combn(seq_along(s), 2))
  i  <- pr[, 1]; j <- pr[, 2]
  flip <- g[i] != g[j] & g[i] == grp_num
  a  <- ifelse(flip, j, i); b <- ifelse(flip, i, j)
  tibble(s1 = s[a], s2 = s[b], g1 = g[a], g2 = g[b]) %>%
    mutate(t1 = ifelse(g1 == g2 & s1 > s2, s2, s1),          # within: sort names
           s2 = ifelse(g1 == g2 & s1 > s2, s1, s2),
           s1 = t1,
           type = ifelse(g1 == g2, paste("within", line_lab[g1]), "between lines"),
           rho  = map2_dbl(s1, s2, ~ cor(l2[, .x], l2[, .y], method = "spearman"))) %>%
    dplyr::select(s1, s2, type, rho)
}

## ---- Spearman scatter: group mean vs group mean ----

## x = mean log2(normalised count + 1) in the reference group, y = same in
## the other group; every tested feature, coloured by DE call; the top n up
## and top n down (by padj) are labelled. The corner label gives the
## reference values needed to read rho:
##   group means, all tested features
##   group means, features expressed in both lines (>= min_count reads in
##     >= n_both samples of EACH line) - removes the on/off features
##   single samples: within-line range (replicate agreement = ceiling) and
##     between-line median (range)
## cnt: raw counts (for the "expressed in both lines" rule)
spearman_scatter <- function(tab, ncnt, cnt, group, file, title, n_lab = 5) {
  group <- setNames(group, colnames(ncnt))
  l2    <- log2(ncnt[tab$feature, , drop = FALSE] + 1)
  rc    <- cnt[tab$feature, , drop = FALSE]
  d  <- tab %>%
    mutate(x = rowMeans(l2[, group == grp_ref, drop = FALSE]),
           y = rowMeans(l2[, group == grp_num, drop = FALSE]),
           both = rowSums(rc[, group == grp_ref, drop = FALSE] >= min_count) >= n_both &
                  rowSums(rc[, group == grp_num, drop = FALSE] >= min_count) >= n_both,
           call = factor(call, levels = c("NS", "Down", "Up"))) %>%
    arrange(call)                                        # DE points on top

  rho_all  <- cor(d$x, d$y, method = "spearman")
  rho_both <- cor(d$x[d$both], d$y[d$both], method = "spearman")
  pr  <- sample_pairs(l2, group)
  win <- range(pr$rho[pr$type != "between lines"])
  btw <- pr$rho[pr$type == "between lines"]

  stats <- tibble(measure = c("group_means_all", "group_means_expressed_in_both",
                              "single_within_min", "single_within_max",
                              "single_between_median", "single_between_min",
                              "single_between_max"),
                  rho = c(rho_all, rho_both, win, median(btw), range(btw)),
                  n_features = c(nrow(d), sum(d$both), rep(nrow(d), 5)))
  write_csv(stats, sub("[.]png$", "_rho.csv", file))

  big <- function(x) formatC(x, format = "d", big.mark = ",")
  lab_txt <- sprintf(paste0("Spearman \u03c1, group means: all tested %.3f (n = %s); ",
                            "expressed in both lines %.3f (n = %s)\n",
                            "Spearman \u03c1, single samples: within line %.3f\u2013%.3f; ",
                            "between lines %.3f (%.3f\u2013%.3f)"),
                     rho_all, big(nrow(d)), rho_both, big(sum(d$both)),
                     win[1], win[2], median(btw), min(btw), max(btw))

  lab <- d %>%
    filter(call != "NS") %>%
    group_by(call) %>%
    slice_min(padj, n = n_lab, with_ties = FALSE) %>%
    ungroup()
  n_de <- table(d$call)

  p <- ggplot(d, aes(x, y, colour = call)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = 2, linewidth = 0.3) +
    geom_point(size = 0.6, alpha = 0.6, stroke = 0) +
    repel_labels(lab, aes(label = feature)) +
    scale_colour_manual(values = c(NS = "grey75", Down = "darkblue", Up = "darkred"),
                        labels = c(NS   = paste0("NS (", n_de[["NS"]], ")"),
                                   Down = paste0("Down (", n_de[["Down"]], ")"),
                                   Up   = paste0("Up (", n_de[["Up"]], ")")),
                        name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    coord_equal() +
    labs(x = paste0("HCC827 (", gsub("_", "-", grp_ref), ")\nmean log2(normalized count + 1)"),
         y = paste0("H1975 (", gsub("_", "-", grp_num), ")\nmean log2(normalized count + 1)"),
         title = title,
         subtitle = paste0(paste0("FDR < ", alpha, ", glmTreat H0: |log2FC| <= ", lfc), "; dashed = identity\n",
                           "expressed in both lines = >= ", min_count, " reads in >= ",
                           n_both, " samples of each line"),
         caption = lab_txt) +
    theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.1),
          plot.caption.position = "plot")
  ggsave(file, p, width = 7.5, height = 7, dpi = 300, bg = "white")
}

## ---- Spearman scatter: replicate vs replicate ----

## one plain scatter per sample pair, log2(normalised count + 1), every
## tested feature; rows: within H1975, within HCC827, then between lines
## (x = HCC827 sample, y = H1975 sample); rho in each panel
spearman_replicates <- function(tab, ncnt, group, file, title) {
  group <- setNames(group, colnames(ncnt))
  l2    <- log2(ncnt[tab$feature, , drop = FALSE] + 1)
  types <- c(paste("within", line_lab[grp_num]), paste("within", line_lab[grp_ref]),
             "between lines")
  pr <- sample_pairs(l2, group) %>%
    mutate(type  = factor(type, levels = types)) %>%
    arrange(type, s2, s1) %>%
    mutate(panel = factor(sprintf("x: %s | y: %s", s1, s2),
                          levels = sprintf("x: %s | y: %s", s1, s2)))
  write_csv(dplyr::select(pr, -panel), sub("[.]png$", ".csv", file))

  pts <- map_dfr(seq_len(nrow(pr)), function(k)
    tibble(panel = pr$panel[k], type = pr$type[k],
           x = l2[, pr$s1[k]], y = l2[, pr$s2[k]]))

  type_cols <- setNames(c(group_cols[[grp_num]], group_cols[[grp_ref]], "grey45"), types)
  n_row <- ceiling(nrow(pr) / 3)

  p <- ggplot(pts, aes(x, y, colour = type)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey30", linetype = 2, linewidth = 0.3) +
    geom_point(size = 0.35, alpha = 0.35, stroke = 0) +
    geom_text(data = pr, aes(x = -Inf, y = Inf, label = sprintf("\u03c1 = %.3f", rho)),
              inherit.aes = FALSE, hjust = -0.1, vjust = 1.4, size = 3) +
    facet_wrap(~ panel, ncol = 3) +
    scale_colour_manual(values = type_cols, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    coord_equal() +
    labs(x = "log2(normalized count + 1), sample on x",
         y = "log2(normalized count + 1), sample on y",
         title = title,
         subtitle = paste0(nrow(l2), " tested features; within line = replicate agreement, ",
                           "between lines = line difference; dashed = identity")) +
    theme(legend.position = "top", strip.text = element_text(size = 7.5))
  ggsave(file, p, width = 8, height = 2.6 * n_row + 1.4, dpi = 300, bg = "white",
         limitsize = FALSE)
}

## ---- top-5 reference + top-5 novel: dumbbell plot ----

## top n per class; q = padj (overall) or padj_lnc (lncRNA)
pick_top <- function(de, n = n_panel) {
  de %>%
    group_by(class) %>%
    arrange(q, dplyr::desc(abs(log2FC)), .by_group = TRUE) %>%
    slice_head(n = n) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

## one row per feature on ONE shared log2 axis: replicates (small points),
## line means (large points) joined by a grey bar. Reference block on top,
## novel block below, each with a small left-aligned heading.
top_panel <- function(top, ncnt, group, row_labs, file, title, subtitle) {
  if (!nrow(top)) return(invisible(NULL))
  group <- setNames(group, colnames(ncnt))
  top <- top %>%
    mutate(class = factor(class, levels = c("Reference", "Novel"))) %>%
    arrange(class, rank)
  f_lev <- rev(top$feature)                                # rank 1 on top
  top <- mutate(top, panel = factor(row_labs[as.integer(class)], levels = row_labs))

  pts <- as_tibble(log2(ncnt[top$feature, , drop = FALSE] + 1), rownames = "feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "value") %>%
    mutate(group = factor(as.character(group[sample]), levels = c(grp_ref, grp_num))) %>%
    left_join(dplyr::select(top, feature, panel), by = "feature") %>%
    mutate(feature = factor(feature, levels = f_lev))
  mns <- pts %>%
    group_by(feature, panel, group) %>%
    summarise(value = mean(value), .groups = "drop")
  seg <- mns %>%
    mutate(group = as.character(group)) %>%
    pivot_wider(names_from = group, values_from = value)

  grp_lab <- function(x) gsub("_", "-", x)
  p <- ggplot() +
    geom_segment(data = seg, aes(x = .data[[grp_ref]], xend = .data[[grp_num]],
                                 y = feature, yend = feature),
                 colour = "grey80", linewidth = 1.6) +
    geom_point(data = pts, aes(value, feature, colour = group), size = 1.1, alpha = 0.8,
               position = position_jitter(width = 0, height = 0.15, seed = 1)) +
    geom_point(data = mns, aes(value, feature, fill = group), shape = 21, size = 3.2,
               colour = "black", stroke = 0.3) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = group_cols, labels = grp_lab, name = NULL) +
    scale_fill_manual(values = group_cols, labels = grp_lab, name = NULL) +
    labs(x = "log2(normalized count + 1)", y = NULL, title = title,
         subtitle = paste0(subtitle, "; small points = replicates, large = line mean")) +
    theme(legend.position = "top",
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold", size = 9),
          panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3))

  ggsave(file, p, width = 7.5, height = 2.2 + 0.3 * nrow(top), dpi = 300, bg = "white")
  write_csv(dplyr::select(top, -panel), sub("[.]png$", ".csv", file))
}

## ---- DE dot plot (after Leung et al. 2024 Nat Commun, Fig. 4a) ----

## Down (left, labels on the left) and Up (right, labels on the right), each
## panel only as wide as its own log2FC range. de: DE rows, q = padj to rank by.
de_dotplot <- function(de, file, title, subtitle, n = n_dot) {
  top <- de %>%
    mutate(direction = ifelse(log2FC_test > 0, "Up", "Down")) %>%
    group_by(direction) %>%
    slice_min(q, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label = paste0(feature, ifelse(zero_in_one_group, " *", "")))
  if (!nrow(top)) return(invisible(NULL))

  ## shared scales so both halves use the same size / colour keys
  size_lim <- range(log10(top$baseMean + 1))
  one_side <- function(dirn) {
    d <- filter(top, direction == dirn) %>%
      mutate(label = fct_reorder(label, abs(log2FC)))     # largest |log2FC| on top
    if (!nrow(d)) return(plot_spacer())
    ggplot(d, aes(log2FC, label)) +
      geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
      geom_segment(aes(x = 0, xend = log2FC, yend = label), colour = "grey80", linewidth = 0.4) +
      geom_point(aes(size = log10(baseMean + 1), fill = category), shape = 21,
                 colour = "black", stroke = 0.3, show.legend = TRUE) +   # keys for absent levels too
      scale_y_discrete(position = if (dirn == "Up") "right" else "left") +
      scale_fill_manual(values = cat_cols, name = NULL, drop = FALSE,
                        limits = intersect(cat_levels, as.character(top$category))) +
      scale_size_continuous(range = c(1.5, 4.5), limits = size_lim,
                            name = "log10(mean\nnorm. count + 1)") +
      guides(fill = guide_legend(override.aes = list(size = 3), order = 1)) +
      labs(x = NULL, y = NULL, subtitle = paste(dirn, "in H1975")) +
      theme(axis.text.y = element_text(size = 7),
            plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 9))
  }

  p <- (one_side("Down") | one_side("Up")) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle,
                    caption = paste0("x = log2 fold change (predFC), H1975 / HCC827",
                                     if (any(top$zero_in_one_group))
                                       "; * 0 reads in every sample of one line" else ""))

  n_row <- max(table(top$direction))
  ggsave(file, p, width = 10, height = 1.9 + 0.2 * n_row, dpi = 300, bg = "white",
         limitsize = FALSE)
  write_csv(top, sub("[.]png$", ".csv", file))
}

## ---- DE composition by feature category ----

## Left: % of tested features called Down (left of 0) and Up (right of 0)
## per category. Right: odds ratio of being DE vs all other tested features
## (two-sided Fisher, 95% CI, log2 axis); filled = BH padj < alpha.
## A category with a zero cell (e.g. every feature DE) has an infinite
## Fisher OR; it is drawn with the Haldane-Anscombe estimate (+0.5 to every
## cell, Woolf CI), marked with a dagger, and the CSV says so (or_corrected = TRUE).
de_composition <- function(tab, file, title) {
  d    <- mutate(tab, de = call != "NS")
  cats <- intersect(cat_levels, unique(as.character(d$category)))

  ## sensitivity checks (CSV only, BH across categories like the main test):
  ##  *_no_onoff: same Fisher test after removing features with 0 reads in
  ##    every sample of one group ("all-zero-in-one-condition")
  ##  *_mh_expr: Mantel-Haenszel OR stratified by expression decile, i.e.
  ##    each category vs other features of similar mean expression
  d$expr_bin <- factor(cut(rank(d$baseMean, ties.method = "first"), 10, labels = FALSE))
  or_test <- function(de, inn, strata = NULL) {
    out <- c(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
    if (!any(inn) || all(inn)) return(out)
    dims <- list(factor(de, c(TRUE, FALSE)), factor(inn, c(TRUE, FALSE)))   # OR = DE odds in / out
    if (!is.null(strata)) dims <- c(dims, list(strata))
    tb <- do.call(table, dims)
    tt <- tryCatch(if (is.null(strata)) fisher.test(tb) else mantelhaen.test(tb),
                   error = function(e) NULL)
    if (is.null(tt)) return(out)
    c(or = unname(tt$estimate), lo = tt$conf.int[1], hi = tt$conf.int[2], p = tt$p.value)
  }
  on <- d$zero_in_one_group %in% TRUE

  enr <- map_dfr(cats, function(cc) {
    inn <- as.character(d$category) == cc
    f0  <- or_test(d$de[!on], inn[!on])
    mh  <- or_test(d$de, inn, d$expr_bin)
    ## 2x2: rows = DE / NS, columns = this category / all others
    k  <- c(de_in  = sum(d$de & inn),  ns_in  = sum(!d$de & inn),
            de_out = sum(d$de & !inn), ns_out = sum(!d$de & !inn))
    ft <- fisher.test(matrix(k[c("de_in", "ns_in", "de_out", "ns_out")], nrow = 2))
    zero <- any(k == 0)
    kh   <- k + 0.5
    or_h <- (kh[["de_in"]] * kh[["ns_out"]]) / (kh[["ns_in"]] * kh[["de_out"]])
    se_h <- sqrt(sum(1 / kh))
    tibble(category = cc, n = sum(inn), n_up = sum(inn & d$call == "Up"),
           n_down = sum(inn & d$call == "Down"), frac_de = k[["de_in"]] / sum(inn),
           odds_ratio = if (zero) or_h else unname(ft$estimate),
           ci_low  = if (zero) exp(log(or_h) - 1.96 * se_h) else ft$conf.int[1],
           ci_high = if (zero) exp(log(or_h) + 1.96 * se_h) else ft$conf.int[2],
           or_corrected = zero, pvalue = ft$p.value,
           median_baseMean = median(d$baseMean[inn]),
           n_onoff = sum(inn & on),
           or_no_onoff = f0[["or"]], ci_low_no_onoff = f0[["lo"]],
           ci_high_no_onoff = f0[["hi"]], pvalue_no_onoff = f0[["p"]],
           or_mh_expr = mh[["or"]], ci_low_mh_expr = mh[["lo"]],
           ci_high_mh_expr = mh[["hi"]], pvalue_mh_expr = mh[["p"]])
  }) %>%
    mutate(padj = p.adjust(pvalue, method = "BH"),
           padj_no_onoff = p.adjust(pvalue_no_onoff, method = "BH"),
           padj_mh_expr  = p.adjust(pvalue_mh_expr, method = "BH"),
           category = factor(category, levels = cat_levels))
  write_csv(enr, sub("[.]png$", ".csv", file))

  big   <- function(x) formatC(x, format = "d", big.mark = ",")
  y_lev <- rev(sprintf("%s (n = %s)", enr$category, big(enr$n)))  # first category on top
  enr   <- mutate(enr, ylab = factor(sprintf("%s (n = %s)", category, big(n)), levels = y_lev))
  y_scale <- scale_y_discrete(limits = y_lev, drop = FALSE,        # same rows in both panels;
                              expand = expansion(add = c(0.6, 1.2)))  # top: room for labels

  bars <- enr %>%
    transmute(ylab, Down = -100 * n_down / n, Up = 100 * n_up / n) %>%
    pivot_longer(c(Down, Up), names_to = "dir", values_to = "x") %>%
    mutate(dir = factor(dir, levels = c("Down", "Up")))

  p1 <- ggplot(bars, aes(x, ylab, fill = dir)) +
    geom_col(width = 0.65) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_text(data = filter(bars, abs(x) >= 0.5),
              aes(label = sprintf("%.0f%%", abs(x)), hjust = ifelse(x < 0, 1.2, -0.2)),
              size = 2.6) +
    y_scale +
    scale_fill_manual(values = c(Down = "darkblue", Up = "darkred"), name = NULL) +
    scale_x_continuous(labels = function(v) paste0(abs(v), "%"),
                       expand = expansion(mult = 0.18)) +
    labs(x = "% of tested features called DE", y = NULL) +
    theme(legend.position = "top")

  ## right panel: one point per category. OR = (odds that a feature of this
  ## category is DE) / (same odds for all other tested features).
  ## OR > 1 (right of the dashed line) = category enriched for DE; OR < 1 =
  ## depleted. Bar = 95% CI; filled = BH padj < alpha; a dagger after padj
  ## marks a Haldane-Anscombe OR (zero cell in the 2x2 table).
  enr <- mutate(enr,
                sig_lab = factor(ifelse(padj < alpha, "sig", "ns"), levels = c("sig", "ns")),
                p_lab   = paste0("padj ", fmt_p(padj), ifelse(or_corrected, " \u2020", "")))
  ok    <- is.finite(enr$ci_low) & is.finite(enr$ci_high) & enr$ci_low > 0
  top_y <- length(y_lev) + 0.8                                     # direction labels row
  p2 <- ggplot(enr, aes(odds_ratio, ylab)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.3) +
    geom_linerange(data = enr[ok, ], aes(xmin = ci_low, xmax = ci_high),
                   colour = "grey30", linewidth = 0.5) +
    geom_point(aes(fill = sig_lab), shape = 21, size = 3, colour = "black", stroke = 0.5,
               show.legend = TRUE) +                        # keep the open key if all are sig
    geom_text(aes(x = Inf, label = p_lab), hjust = -0.1, size = 2.5) +
    annotate("text", x = 1 / 1.15, y = top_y, label = "\u2190 depleted", hjust = 1,
             size = 2.6, colour = "grey30") +
    annotate("text", x = 1.15, y = top_y, label = "enriched \u2192", hjust = 0,
             size = 2.6, colour = "grey30") +
    y_scale +
    scale_x_continuous(trans = "log2") +
    expand_limits(x = 1) +
    scale_fill_manual(values = c(sig = "black", ns = "white"),
                      labels = c(sig = paste0("< ", alpha), ns = paste0("\u2265 ", alpha)),
                      name = "BH padj", drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = "odds ratio of being DE vs all other features (log2 axis, 95% CI)", y = NULL) +
    theme(legend.position = "top",
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 60, 5.5, 5.5))

  cap <- if (any(enr$or_corrected))
    "\u2020 zero count in the 2x2 table: Haldane-Anscombe OR (+0.5 per cell, Woolf 95% CI)" else NULL
  p <- p1 + p2 + plot_layout(widths = c(1.3, 1)) +
    plot_annotation(title = title, caption = cap,
                    subtitle = "OR vs all other tested features (Fisher, two-sided), BH across categories")
  ggsave(file, p, width = 10.5, height = 2.5 + 0.45 * length(cats), dpi = 300, bg = "white")
}

## ----------------------------------------------------------------- run ----
dat <- read_counts(data_dir)

## the counts edgeR is given (divided, sequins removed, original IDs), so
## 06_isoformic.R reads exactly this matrix
for (lv in levels_to_run) {
  dir.create(file.path(outroot, lv), recursive = TRUE, showWarnings = FALSE)
  m  <- dat[[lv]]$counts
  id <- unname(dat[[lv]]$id[rownames(m)])
  ec <- if (lv == "transcript") {
    data.frame(TXNAME = id, GENEID = unname(dat[[lv]]$gid[rownames(m)]), m, check.names = FALSE)
  } else {
    data.frame(GENEID = id, m, check.names = FALSE)
  }
  write_tsv(ec, file.path(outroot, lv, "edger_input_counts.tsv"))
}
all_res <- list(); all_lnc <- list()

for (lv in levels_to_run) {
  cts  <- dat[[lv]]$counts                      # divided counts (edgeR input)
  rds  <- dat[[lv]]$reads                       # undivided reads (filter)
  bt   <- unname(dat[[lv]]$biotype[dat[[lv]]$id[rownames(cts)]])
  drop <- !is.na(bt) & bt %in% drop_biotypes
  excluded <- tibble(feature = rownames(cts)[drop], biotype = bt[drop],
                     total_count = rowSums(rds)[drop])
  cts  <- cts[!drop, , drop = FALSE]
  rds  <- rds[!drop, , drop = FALSE]
  message(lv, ": removed ", sum(drop), " rRNA / processed-pseudogene features")
  cd  <- make_coldata(cts, sheet)
  colnames(cts) <- colnames(rds) <- rownames(cd)

  y      <- DGEList(counts = cts, group = cd$group)
  design <- model.matrix(~ group, data = cd)
  stopifnot(coef_name %in% colnames(design))
  ## filter on the read counts (undivided), then test the divided counts
  keep <- if (filter_mode == "filterByExpr") {
    filterByExpr(DGEList(counts = rds, group = cd$group), design = design)
  } else {
    rowSums(rds >= min_count) >= min_samples
  }
  y  <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])   # TMM
  yd <- estimateDisp(y, design, robust = TRUE)           # BCV plot + predFC only

  outdir <- file.path(outroot, lv)
  lncdir <- file.path(outdir, "lncrna")
  dir.create(lncdir, recursive = TRUE, showWarnings = FALSE)
  write_csv(excluded, file.path(outdir, "excluded_biotypes.csv"))
  ttl <- paste(dataset_lab, lv)

  ## normalised counts on the count scale + logCPM for heatmaps / PCA
  eff  <- y$samples$lib.size * y$samples$norm.factors
  ncnt <- cpm(y) * mean(eff) / 1e6
  vmat <- cpm(y, log = TRUE, prior.count = fc_prior)
  grp  <- setNames(y$samples$group, colnames(y))
  saveRDS(list(mat = vmat, norm = ncnt,
               id = unname(dat[[lv]]$id[rownames(y)]), group = y$samples$group),
          file.path(outdir, "logcpm.rds"))

  ## ---- overall DE ----
  res <- test_contrast(y, yd, design, ncnt, dat[[lv]]$id, dat[[lv]]$biotype)
  tab <- res$tab %>% mutate(overdispersion = unname(dat[[lv]]$od[feature]),
                            .after = baseMean)
  qc_plots(yd, res$fit, vmat, outdir, ttl)
  write_csv(tab, file.path(outdir, paste0("de_", ct_name, ".csv")))
  ggsave(file.path(outdir, paste0("volcano_", ct_name, ".png")),
         volcano(tab, paste(ttl, ct_name)),
         width = 7, height = 5.5, dpi = 300, bg = "white")
  top_heatmap(tab, vmat, y$samples$group,
              file.path(outdir, paste0("heatmap_top_", ct_name, ".png")),
              paste0(ttl, " - ", ct_name, "\ntop ", n_top, " up + top ", n_top,
                     " down, logCPM row z-score"))

  ## ---- lncRNA DE ----
  lnc <- lnc_table(tab)
  write_csv(lnc, file.path(lncdir, paste0("lncrna_", ct_name, ".csv")))
  ggsave(file.path(lncdir, paste0("volcano_lncrna_", ct_name, ".png")),
         lnc_volcano(tab, lnc, paste("lncRNAs -", ttl, ct_name)),
         width = 7, height = 5.5, dpi = 300, bg = "white")
  lnc_heatmap(lnc, vmat, y$samples$group,
              file.path(lncdir, paste0("heatmap_lncrna_", ct_name, ".png")),
              paste("Top DE lncRNAs -", ttl, "\nlogCPM, row z-score"))

  ## ---- Spearman scatter ----
  spearman_scatter(tab, ncnt, y$counts, y$samples$group,
                   file.path(outdir, paste0("spearman_scatter_", ct_name, ".png")),
                   paste(ttl, "- group means"))
  spearman_replicates(tab, ncnt, y$samples$group,
                      file.path(outdir, paste0("spearman_replicates_", ct_name, ".png")),
                      paste(ttl, "- replicate vs replicate"))

  ## ---- top-5 panels ----
  de_all <- tab %>%
    filter(call != "NS") %>%
    mutate(q = padj, class = ifelse(novel, "Novel", "Reference"))
  top_panel(pick_top(de_all), ncnt, grp,
            c(paste("Top", n_panel, "reference DE", lv, "(annotated)"),
              paste("Top", n_panel, "novel DE", lv)),
            file.path(outdir, paste0("top", n_panel, "_ref_novel_", ct_name, ".png")),
            paste(ttl, "-", ct_name, "- top DE, reference vs novel"),
            paste0("ranked by padj (H0 |log2FC| <= ", lfc, "); label log2FC = predFC"))

  de_lnc <- lnc %>%
    filter(call_lnc != "NS") %>%
    mutate(q = padj_lnc, class = ifelse(biotype == "novel_lncRNA", "Novel", "Reference"))
  top_panel(pick_top(de_lnc), ncnt, grp,
            c(paste("Top", n_panel, "annotated lncRNA DE", lv),
              paste("Top", n_panel, "novel lncRNA DE", lv)),
            file.path(lncdir, paste0("top", n_panel, "_lncrna_ref_novel_", ct_name, ".png")),
            paste(ttl, "-", ct_name, "- top DE lncRNAs, annotated vs novel"),
            "ranked by padj within lncRNAs (padj_lnc); label log2FC = predFC")

  ## ---- DE dot plots ----
  de_dotplot(de_all,
             file.path(outdir, paste0("dotplot_top_", ct_name, ".png")),
             paste(ttl, "-", ct_name, "- top DE"),
             paste0("top ", n_dot, " Up + top ", n_dot, " Down by padj (H0 |log2FC| <= ", lfc, ")"))
  de_dotplot(de_lnc,
             file.path(lncdir, paste0("dotplot_top_lncrna_", ct_name, ".png")),
             paste(ttl, "-", ct_name, "- top DE lncRNAs"),
             paste0("top ", n_dot, " Up + top ", n_dot, " Down by padj within lncRNAs"))

  ## ---- DE composition by category ----
  de_composition(tab, file.path(outdir, paste0("de_composition_", ct_name, ".png")),
                 paste(ttl, "-", ct_name, "- DE by feature category"))

  ## ---- log ----
  message(sprintf("%-10s tested %6d  up %5d  down %5d | lncRNA tested %5d  up %4d  down %4d",
                  lv, nrow(tab), sum(tab$call == "Up"), sum(tab$call == "Down"),
                  nrow(lnc), sum(lnc$call_lnc == "Up"), sum(lnc$call_lnc == "Down")))
  all_res[[lv]] <- mutate(tab, level = lv, contrast = ct_name)
  all_lnc[[lv]] <- mutate(lnc, level = lv, contrast = ct_name)
}

write_csv(bind_rows(all_res), file.path(outroot, "all_results.csv"))
lnc_summary <- bind_rows(all_lnc) %>%
  dplyr::count(level, contrast, biotype, call_lnc) %>%
  pivot_wider(names_from = call_lnc, values_from = n, values_fill = 0)
write_csv(lnc_summary, file.path(outroot, "lncrna_summary.csv"))
print(lnc_summary, n = Inf)

message("Done: ", outroot)
