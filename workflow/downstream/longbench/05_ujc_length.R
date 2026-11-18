## ===========================================================================
## 05_ujc_length.R - cross-protocol novel transcripts, stratified by length  (Fig. 1j)
## The spike-ins showed a length window per protocol:
##   PacBio   under-captures < 1 kb (Kinnex size selection)
##   ONT cDNA under-captures > 2 kb
##   ONT DRS  flattest, drops only > ~5 kb
## This script asks whether protocol-specific novel transcripts follow the
## same windows.
##
## Inputs: unique_intron_chains/results/ujc_table.tsv (04_ujc_report.R),
##         <run>/transcriptome/curated/metadata/novel_transcripts_validated_metadata.csv
## Multi-exon UJCs only (mono-exon models are never matched across protocols).
## Length = spliced transcript length (`len` / `spliced_length_*`).
##
## Outputs: unique_intron_chains/results/length/
##   1. novel_length_by_protocol.csv   length distribution per protocol (all / novel lncRNA)
##   2. ujc_category_by_length.csv     sharing category composition per length bin
##   3. corroboration_by_length.csv    % of protocol X's novel transcripts whose UJC
##                                     is also reported by Y, by X's own length (Fig. 1j)
## NOTE: a shared UJC is the same intron chain, not the same transcript model;
##       ends and lengths can differ between protocols. No consolidated IDs.
## ===========================================================================

suppressPackageStartupMessages(library(tidyverse))

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
indir  <- ujc_dir                       # ujc_table.tsv from 04_ujc_report.R
outdir <- file.path(indir, "length")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prot  <- c(ONTd = "ONT DRS", ONTc = "ONT cDNA", PB = "PacBio")
files <- c(ONTd = "ont_drs", ONTc = "ont_cdna", PB = "pb")     # run folders
bins  <- c(0, 500, 1000, 2000, 5000, Inf)
labs  <- c("<0.5 kb", "0.5-1 kb", "1-2 kb", "2-5 kb", ">5 kb")
len_bin <- function(x) cut(x, bins, labels = labs, right = FALSE)  # [0.5, 1) kb: a 500-nt model falls in 0.5-1 kb

## ------------------------------------ 1. length of novel transcripts ------
nov <- imap_dfr(files, function(f, k)
  read_csv(file.path(run_dir[[f]], pp$novel_meta),
           show_col_types = FALSE) %>%
    transmute(protocol = prot[[k]], id = qry_id, len,
              biotype = pulposeq_transcript_biotype))

nov2 <- bind_rows(mutate(nov, set = "all novel transcripts"),
                  nov %>% filter(biotype == "novel_lncRNA") %>% mutate(set = "novel lncRNAs")) %>%
  mutate(protocol = factor(protocol, levels = prot))

len_tab <- nov2 %>%
  mutate(bin = len_bin(len)) %>%
  group_by(set, protocol) %>%
  mutate(n_total = n(), median_len = median(len)) %>%
  group_by(set, protocol, n_total, median_len, bin) %>%
  summarise(pct = round(100 * n() / dplyr::first(n_total), 1), .groups = "drop") %>%
  pivot_wider(names_from = bin, values_from = pct, values_fill = 0)
write_csv(len_tab, file.path(outdir, "novel_length_by_protocol.csv"))
print(len_tab, width = Inf)

## ------------------------------------------ 2. UJC sharing by length ------
u <- read_tsv(file.path(indir, "ujc_table.tsv"), show_col_types = FALSE) %>%
  filter(structure == "multi-exon") %>%
  mutate(across(starts_with("present_"), ~ as.logical(as.character(.x))))
sl    <- as.matrix(u[, paste0("spliced_length_", names(prot))])
u$len <- apply(sl, 1, median, na.rm = TRUE)            # representative length of the UJC
u$bin <- len_bin(u$len)

cat_tab <- u %>% dplyr::count(bin, category) %>%
  group_by(bin) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
write_csv(cat_tab, file.path(outdir, "ujc_category_by_length.csv"))

## ---------------------------------- 3. corroboration by own length --------
## of protocol X's novel transcripts, the % whose UJC is also reported by Y
## (by X's spliced length); "not reported" != "not expressed" in Y
corr <- map_dfr(names(prot), function(x) map_dfr(setdiff(names(prot), x), function(y) {
  ux <- u[u[[paste0("present_", x)]] %in% TRUE, ]
  tibble(from = prot[[x]], found_in = prot[[y]],
         len = ux[[paste0("spliced_length_", x)]],
         found = ux[[paste0("present_", y)]] %in% TRUE)
})) %>%
  mutate(bin = len_bin(len))

corr_tab <- corr %>%
  group_by(from, found_in, bin) %>%
  summarise(n = n(), pct_found = round(100 * mean(found), 1), .groups = "drop")
write_csv(corr_tab, file.path(outdir, "corroboration_by_length.csv"))
print(corr_tab, n = Inf)

message("Done: ", outdir)
