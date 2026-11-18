# Unique Intron Chains - Aggregating novel transcripts across protocols

Novel transcripts of the three LongBench protocols grouped by unique intron chain (UIC), to see
which novel structures several protocols agree on. No counts are modelled, so these scripts are
the pre-integration ones with the new input paths
(`transcriptome/curated/annotation/novel_transcripts_validated.gtf`,
`transcriptome/curated/metadata/novel_transcripts_validated_metadata.csv`).

| Script | What it does | Outputs (`results/`) |
|---|---|---|
| `04_ujc_report.R` | novel transcripts grouped by intron chain; sharing across protocols, class-code concordance, transcript-end differences | `ujc_table.tsv`, `ujc_pair_deltas.tsv`, `summary_*.tsv` |
| `05_ujc_length.R` | sharing by transcript length (needs 04) | `length/` |
| `07_supp_tables_ujc.R` | supplementary tables ST2-ST7 (needs 04) | `supp_tables/` |

The novel set itself comes from the oarfish runs: the expression filter used oarfish counts
(>= 5 summed across samples), so it can differ from the pre-integration set.
