#!/usr/bin/bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz

gunzip gencode.v50.primary_assembly.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

cat gencode.v50.primary_assembly.annotation.gtf sequins_mix_A_B.gtf > gencode.v50.primary_assembly.annotation.dong.sequins.gtf
cat GRCh38.primary_assembly.genome.fa sequins_mix_A_B.fasta > GRCh38.primary_assembly.genome.dong.sequins.fa

cat gencode.v50.primary_assembly.annotation.gtf sequins_sirv_ercc_longbench.gtf > gencode.v50.primary_assembly.annotation.longbench.spikeins.gtf
cat GRCh38.primary_assembly.genome.fa sequins_sirv_ercc_longbench.fasta > GRCh38.primary_assembly.genome.longbench.spikeins.fa