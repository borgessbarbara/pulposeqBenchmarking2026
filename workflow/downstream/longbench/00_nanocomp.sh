#!/bin/bash
#---------------Script SBATCH - NPAD ----------------
#SBATCH -J nanocomp_pb_ont_drs
#SBATCH -p amd-512
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=mbarborgess@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 8:00:00
#SBATCH -o nanocomp_pb_ont_%j.err.out
#SBATCH -e nanocomp_pb_ont_%j.err.out

NanoComp --fastq H1975_ONT.fastq.gz H1975_ONT_DRS.fastq.gz H1975_PB.fastq.gz \
    HCC827_ONT.fastq.gz HCC827_ONT_DRS.fastq.gz HCC827_PB.fastq.gz \
    H2228_ONT.fastq.gz H2228_ONT_DRS.fastq.gz H2228_PB.fastq.gz \
    H526_ONT.fastq.gz  H526_ONT_DRS.fastq.gz H526_PB.fastq.gz \
    H211_ONT.fastq.gz H211_ONT_DRS.fastq.gz H211_PB.fastq.gz \
    SHP77_ONT.fastq.gz SHP77_ONT_DRS.fastq.gz SHP77_PB.fastq.gz \
    H146_ONT.fastq.gz H146_ONT_DRS.fastq.gz H146_PB.fastq.gz \
    H69_ONT.fastq.gz H69_ONT_DRS.fastq.gz H69_PB.fastq.gz \
    -t 24 -o nanocomp/ --plot 'violin' --format 'svg' \
    --names LUAD_ONTc_H1975 LUAD_ONTd_H1975 LUAD_PB_H1975 \
    LUAD_ONTc_HCC827 LUAD_ONTd_HCC827 LUAD_PB_HCC827 \
    LUAD_ONTc_H2228 LUAD_ONTd_H2228 LUAD_PB_H2228 \
    SCLC-P_ONTc_H526 SCLC-P_ONTd_H526 SCLC-P_PB_H526 \
    SCLC-P_ONTc_H211 SCLC-P_ONTd_H211 SCLC-P_PB_H211 \
    SCLC-A_ONTc_SHP77 SCLC-A_ONTd_SHP77 SCLC-A_PB_SHP77 \
    SCLC-A_ONTc_H146 SCLC-A_ONTd_H146 SCLC-A_PB_H146 \
    SCLC-A_ONTc_H69 SCLC-A_ONTd_H69 SCLC-A_PB_H69
