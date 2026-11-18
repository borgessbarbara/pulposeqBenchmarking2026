## ---------------------------------------------------------------------------
## 05_cis_lncRNA_top10.R - alternative cis-pair circos + heatmap (Dong):
## the top 10 novel lncRNAs in potential cis-acting pairs with REFERENCE genes
##
## Runs 04_cis_circos.R with a variant (same data, statistics, drawing and print
## sizes; nothing in 04's own results changes):
##   pairs   novel transcript = "Novel lncRNA"; partner = an annotated (reference)
##           transcript of any biotype of ANOTHER gene (no same-gene pairs); both DE (04's de_filter "both"); within
##           1 Mb, |rho| > 0.5 (04's settings); one partner per novel lncRNA (the
##           strongest), ranked by |rho| (2 decimals), then by the lncRNA's DE FDR; top 10
##   heatmap every selected pair (any reference biotype)
## Outputs: correlation-cis-1Mb/results_lncRNA_top10/
##   paper_cis_circos.{pdf,svg,png}, paper_cis_heatmap.{pdf,svg,png},
##   circos_cis_de_pairs.csv (the 10 pairs), and 04's other files for this selection
## Run: Rscript 05_cis_lncRNA_top10.R  (from this folder)
## ---------------------------------------------------------------------------

source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/dong/00_paths.R")   # cis_dir
cis_variant <- list(
  outdir    = file.path(dirname(cis_dir), "results_lncRNA_top10"),
  max_pairs = 10,
  ## with 3 + 3 samples |rho| saturates (most pairs separating the two lines
  ## have |rho| = 1), so ties are broken by the novel lncRNA's DE evidence
  ## (smallest FDR first); one partner per lncRNA, then the top 10
  filter    = function(sel) {
    ## no same-gene pairs: a novel isoform with a reference isoform of its own
    ## gene is not a cis pair between loci
    sel <- sel[as.character(sel$novel_category) == "Novel lncRNA" & !sel$partner_novel & !sel$same_gene, ]
    sel[order(-round(abs(sel$rho), 2), sel$novel_padj), ]
  },
  hm_keep   = c("Annotated protein-coding", "Annotated lncRNA", "Annotated other"),
  title     = "Novel lncRNAs in cis-acting pairs"
)
source(file.path(dg_dir, "04_cis_circos.R"))
