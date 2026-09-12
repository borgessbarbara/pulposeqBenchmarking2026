# pulposeqBenchmarking2026

Collection of the [integrativebioinformatics/pulposeq](https://github.com/integrativebioinformatics/pulposeq.git) workflow integration and benchmarking for isoform-level analysis of bulk long-read RNA-seq public datasets.

## Datasets under benchmarking

### ENCODE: Rush Alzheimer’s Disease Study (PacBio Sequel II)
Data chosen corresponds to tissue samples from the dorsolateral prefrontal cortex of female donors.

Access the RUSH AD datasets at the [ENCODE website](https://www.encodeproject.org/brain-matrix/?type=Experiment&status=released&internal_tags=RushAD).

**Alzheimer’s Cognitive Impairment** (90+ years)

ENCFF708BOPO, ENCFF785KVJ, ENCFF446EFU and ENCFF156TTD

**No cognitive impairment**

- 79 years: ENCFF838DFB
- 88 years: ENCFF206TQZ
- 90+ years: ENCFF827DUW and ENCFF260AWP

The ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57–74 (2012). https://doi.org/10.1038/nature11247

### Lung Adenocarcinoma cell lines H1975 and HCC827 (PromethION ONT cDNA)

BioProject:  PRNA723287

- H1975 replicates: SRR14286054, SRR14286055 and SRR14286056.
- HCC827 replicates: SRR14286063, SRR14286064 and SRR14286065.

Dong, X., Du, M.R.M., Gouil, Q. et al. Benchmarking long-read RNA-sequencing analysis tools using in silico mixtures. Nat Methods 20, 1810–1821 (2023). https://doi.org/10.1038/s41592-023-02026-3

### LongBench
From the [LongBench](https://github.com/mritchielab/LongBench.io) sequencing of LUAD, SCLC-A and SCLC-P cell lines, we evaluated pipeline runs across the bulk PacBio, ONT cDNA, ONT direct RNA-seq datasets.

Sequencing data are available from Gene Expression Omnibus under accession number GSE303762 and [AWS S3](https://registry.opendata.aws/longbench/).

You can directly browse and download individual files using the [LongBench data browser](https://mritchielab.github.io/LongBench.io/browse-data/).

You, Y. et al. Benchmarking long-read RNA-sequencing technologies with LongBench: a cross-platform reference dataset profiling cancer cell lines with bulk and single-cell approaches. bioRxiv; 2025.  [https://doi.org/10.1101/2025.09.11.675724]().

