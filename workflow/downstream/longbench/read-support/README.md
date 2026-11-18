# Read-level support per transcript biotype (Fig. 1e)

`06_read_support.R`: for every curated transcript, whether Bambu has at least one full-length
read and at least one uniquely assigned read, per biotype and protocol.

- Full-length and unique support are Bambu's (`BambuOutput_{fullLength,unique}Counts_transcript_validated.txt`).
- The curated set is chosen on oarfish counts in these runs, so transcripts kept on oarfish
  counts alone have no Bambu reads and fall in "neither".
- Outputs (`results/`): `read_support_<run>.csv`, `read_support_all_runs_{long,wide}.csv`,
  `read_support_tables.xlsx`.
