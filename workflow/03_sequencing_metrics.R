## ===========================================================================
## 03_sequencing_metrics.R - sequencing, orientation and alignment metrics per
## library (LongBench ONT cDNA, ONT DRS, PacBio; Dong ONT cDNA), from the
## NanoComp and Restrander outputs pulposeq writes. Extended Data figure
## (sequencing metrics) + the tables behind Supplementary Table 1.
##
##   a  reads per library at each stage the pipeline reports: raw reads,
##      reads oriented by Restrander (ONT cDNA only), alignment records (BAM)
##   b  Restrander orientation per library (ONT cDNA): +, - and undetermined
##   c  read length of the raw reads: median and N50 per library
##   d  median read quality (Phred) of the raw reads per library
##   e  alignment: median percent identity and share of bases aligned
##      (total bases aligned / total bases of the aligned reads) per library
##
## Inputs (read only): <project>/<run>/nanocomp/{raw_reads,rest_reads,mapping}/AllNanoStats.txt
##                     <project>/<run>/restrander/*.restrander_stats.json
## Outputs: downstream/sequencing-metrics/results/
##            sequencing_metrics_per_library.csv  (NanoStats general summary, one row per library and stage)
##            restrander_orientation.csv          (one row per library)
##          <project>/figures/"Extended data figure 2".{pdf,svg,png}
## ===========================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
source(file.path(repo, "downstream/longbench/helpers/plot_theme.R"))      # font_family, protocol_cols
source(file.path(repo, "downstream/longbench/helpers/figure_helpers.R"))  # save_panel
out_dir <- file.path(repo, "downstream/sequencing-metrics/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- c(ont_cdna = "LongBench ONT cDNA", ont_drs = "LongBench ONT DRS",
          pb = "LongBench PacBio", dong = "Dong ONT cDNA")
run_cols <- c("LongBench ONT cDNA" = protocol_cols[["ONT cDNA"]], "LongBench ONT DRS" = protocol_cols[["ONT DRS"]],
              "LongBench PacBio" = protocol_cols[["PacBio"]], "Dong ONT cDNA" = "#27408B")   # as Extended Data Fig. 1
stages <- c(raw_reads = "Raw reads", rest_reads = "Restranded reads", mapping = "Alignments (BAM)")

## ---- NanoStats "General summary" block: one row per library ----------------
nanostats <- function(f) {
  l <- readLines(f, warn = FALSE)
  i <- which(grepl("^General summary:", l)); j <- which(!grepl(":", l))[1] - 1
  hdr <- strsplit(trimws(sub("^General summary:", "", l[i])), "\\s+")[[1]]
  rows <- l[(i + 1):(if (is.na(j)) length(l) else j)]
  rows <- rows[grepl("^[A-Za-z].*:\\s", rows) & !grepl("^Number, percentage|^Top 5", rows)]
  m <- do.call(rbind, lapply(rows, function(r) {
    k <- sub(":.*$", "", r); v <- strsplit(trimws(sub("^[^:]*:", "", r)), "\\s+")[[1]]
    c(k, v[seq_along(hdr)])
  }))
  as_tibble(m, .name_repair = "minimal") %>% setNames(c("metric", hdr)) %>%
    pivot_longer(-metric, names_to = "library", values_to = "value") %>%
    mutate(value = as.numeric(gsub(",", "", value)),
           library = sub("\\.restranded$", "", sub("\\.(fastq|fq)(\\.gz)?$|\\.bam$|_sorted$", "", library))) %>%
    pivot_wider(names_from = metric, values_from = value)
}
ns <- map_dfr(names(runs), function(r) map_dfr(names(stages), function(st) {
  f <- file.path(project, r, "nanocomp", st, "AllNanoStats.txt")
  if (!file.exists(f)) return(NULL)
  nanostats(f) %>% mutate(run = runs[[r]], stage = stages[[st]], .before = 1)
}))
## the same library is named differently by stage (raw H146_ONT, aligned H146_ONTc):
## `sample` drops the protocol suffix so a library can be followed across stages
sample_of <- function(x) sub("_(ONT_DRS|ONTc|ONTd|ONT|PB)$", "", x)
ns <- ns %>% mutate(dataset = ifelse(grepl("^Dong", run), "Dong", "LongBench"), .before = 1) %>%
  mutate(sample = sample_of(library), .after = library)
## NanoComp prints "Fraction of bases aligned" to one decimal (0.9 for every ONT
## cDNA library); recompute it from the two totals NanoComp reports for the
## alignments, as panel e does, so Supplementary Table 1a carries the real value
if (all(c("Total bases aligned", "Total bases") %in% names(ns)))
  ns <- ns %>% mutate(`Fraction of bases aligned` = round(`Total bases aligned` / `Total bases`, 4))
write_csv(ns, file.path(out_dir, "sequencing_metrics_per_library.csv"))

rs0 <- map_dfr(c("ont_cdna", "dong"), function(r) {
  fs <- list.files(file.path(project, r, "restrander"), "restrander_stats\\.json$", full.names = TRUE)
  map_dfr(fs, function(f) {
    j <- jsonlite::fromJSON(f)$stats
    tibble(run = runs[[r]], library = sub("\\.restrander_stats\\.json$", "", basename(f)),
           total_reads = j$totalReads, plus = j$strandStats[["+"]], minus = j$strandStats[["-"]],
           undetermined = j$strandStats[["?"]],
           `TSO-TSO` = j$artefactStats[["TSO-TSO"]], `RTP-RTP` = j$artefactStats[["RTP-RTP"]])
  })
}) %>% mutate(oriented_fraction = (plus + minus) / total_reads)
rs <- rs0 %>% mutate(sample = sample_of(library), .after = library)
write_csv(rs, file.path(out_dir, "restrander_orientation.csv"))

## libraries missing at a stage (e.g. a restranded file NanoComp did not get)
miss <- ns %>% distinct(run, stage, sample) %>% count(run, sample) %>%
  left_join(ns %>% distinct(run, stage) %>% count(run, name = "n_stages"), by = "run") %>% filter(n < n_stages)
if (nrow(miss)) message("  ! libraries absent from a NanoComp stage: ",
                        paste(sprintf("%s %s", miss$run, miss$sample), collapse = "; "))

## ---------------------------------------------------------------- figure ---
theme_set(theme_bw(base_family = font_family, base_size = 8) +
  theme(text = element_text(face = "plain"), plot.title = element_text(size = 8.5),
        plot.title.position = "plot", strip.background = element_rect(fill = "grey94"),
        strip.text = element_text(size = 7.5), panel.grid.minor = element_blank(),
        legend.key.size = unit(3, "mm"), plot.margin = margin(1, 2, 1, 1, "mm")))
run_lv <- unname(runs)
## run labels under c-e: one line each, set at 45 degrees as in a
short  <- c("LongBench ONT cDNA" = "ONT cDNA", "LongBench ONT DRS" = "ONT DRS",
            "LongBench PacBio" = "PacBio", "Dong ONT cDNA" = "ONT cDNA (Dong)")
x45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
two_line <- c("LongBench ONT cDNA" = "LongBench\nONT cDNA", "Dong ONT cDNA" = "Dong\nONT cDNA")
nsf <- ns %>% mutate(run = factor(run, levels = run_lv), stage = factor(stage, levels = stages))
jit <- position_jitter(width = 0.15, height = 0, seed = 1)
dots <- function(d, y, ylab, title) ggplot(d, aes(run, {{ y }}, colour = run)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey30") +
  geom_point(position = jit, size = 1.1, alpha = 0.85) +
  scale_colour_manual(values = run_cols, guide = "none") +
  scale_x_discrete(labels = short) + labs(x = NULL, y = ylab, title = title) +
  theme(panel.grid.major.x = element_blank()) + x45

## a: reads per library and stage (one line per library)
pA <- ggplot(nsf, aes(stage, `Number of reads` / 1e6, group = sample, colour = run)) +
  geom_line(linewidth = 0.3, alpha = 0.6) + geom_point(size = 1, alpha = 0.85) +
  facet_wrap(~ run, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = run_cols, guide = "none") +
  scale_x_discrete(labels = c("Raw reads" = "Raw", "Restranded reads" = "Oriented",
                              "Alignments (BAM)" = "Aligned")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Reads per library (millions)",
       title = "Reads per library: raw, oriented by Restrander (ONT cDNA) and aligned (BAM records); one line per library")
raw <- nsf %>% filter(stage == "Raw reads")
## c: read length (median and N50, raw reads)
pB <- raw %>% select(run, library, Median = `Median read length`, N50 = `Read length N50`) %>%
  pivot_longer(c(Median, N50), names_to = "metric") %>%
  ggplot(aes(run, value / 1000, colour = run, shape = metric)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1), size = 1.1) +
  scale_colour_manual(values = run_cols, guide = "none") +
  scale_shape_manual(values = c(Median = 16, N50 = 2), name = NULL) +
  scale_x_discrete(labels = short) + scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Read length (kb)", title = "Read length, raw reads") +
  theme(panel.grid.major.x = element_blank(), legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_blank()) + x45
## d: read quality (raw)
pC <- dots(raw, `Median read quality`, "Median read quality (Phred)", "Read quality, raw reads")
## e: alignment identity and share of bases aligned
aln <- nsf %>% filter(stage == "Alignments (BAM)") %>%
  transmute(run, library, `Identity, median (%)` = `Median percent identity`,
            `Bases aligned (%)` = 100 * `Total bases aligned` / `Total bases`) %>%
  pivot_longer(-c(run, library), names_to = "metric") %>%
  mutate(metric = factor(metric, levels = c("Identity, median (%)", "Bases aligned (%)")))
pD <- dots(aln, value, NULL, "Alignment (minimap2)") + facet_wrap(~ metric, nrow = 1, scales = "free_y")
## b: Restrander orientation (drawn beside a)
pE <- rs %>% mutate(run = factor(run, levels = run_lv)) %>%
  pivot_longer(c(plus, minus, undetermined), names_to = "strand") %>%
  mutate(strand = factor(strand, levels = c("plus", "minus", "undetermined"),
                         labels = c("+ kept", "\u2212 reversed", "undetermined, removed")),
         share = value / total_reads) %>%
  ## vertical bars: one column per library, datasets side by side
  ggplot(aes(sample, share * 100, fill = strand)) +
  geom_col(width = 0.75) +
  facet_grid(. ~ run, scales = "free_x", space = "free_x",
             labeller = as_labeller(c("LongBench ONT cDNA" = "LongBench ONT cDNA", "Dong ONT cDNA" = "Dong ONT cDNA"))) +
  scale_fill_manual(values = c("grey35", "grey65", "#E69F00"), name = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Reads (%)", title = "Restrander orientation (ONT cDNA)") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        panel.grid.major.x = element_blank(), legend.margin = margin(0, 0, 0, 0))

## top row: a | b (orientation), same height and width; below: length,
## quality, alignment. Letters follow reading order.
pA <- pA + labs(title = "Reads per library (aligned = BAM records)") +
  facet_wrap(~ run, nrow = 1, scales = "free_x", labeller = as_labeller(c("LongBench ONT cDNA" = "ONT cDNA",
             "LongBench ONT DRS" = "ONT DRS", "LongBench PacBio" = "PacBio", "Dong ONT cDNA" = "ONT cDNA\n(Dong)"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
h_top <- 70; h_m <- 65   # bottom row taller: room for the 45-degree run labels
fig_w <- 180   # mm
p <- ((pA | pE) + plot_layout(widths = c(1.35, 1))) / ((pB | pC | pD) + plot_layout(widths = c(1, 1, 1.7))) +
  plot_layout(heights = c(h_top, h_m)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "plain", family = font_family))
save_panel(p, "Extended data figure 2", c(fig_w, h_top + h_m), file.path(project, "figures"))
message("Done: ", out_dir)
