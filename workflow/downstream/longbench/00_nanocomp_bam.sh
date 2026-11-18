#!/bin/bash
#SBATCH -J nanocomp_bam_pb_ont_drs
#SBATCH -p amd-512
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-user=mbarborgess@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 8:00:00
#SBATCH -o nanocomp_pb_ont_%j.err.out
#SBATCH -e nanocomp_pb_ont_%j.err.out

BASE=/home/mbbdsantana/scratch/paper
OUTBASE=$BASE/mapping_stats
WORKDIR=$OUTBASE/bam_input

# Stage BAM/BAI symlinks in a single directory so NanoComp receives bare filenames
mkdir -p "$WORKDIR"
ln -sf "$BASE"/ont/results/minimap2/*.bam "$BASE"/ont/results/minimap2/*.bam.bai "$WORKDIR"/
ln -sf "$BASE"/ont_drs/results/minimap2/*.bam "$BASE"/ont_drs/results/minimap2/*.bam.bai "$WORKDIR"/
ln -sf "$BASE"/pb/results/minimap2/*.bam "$BASE"/pb/results/minimap2/*.bam.bai "$WORKDIR"/
cd "$WORKDIR" || exit 1

NanoComp --bam H1975_ONTc.bam H1975_ONTd.bam H1975_PB.bam \
    HCC827_ONTc.bam HCC827_ONTd.bam HCC827_PB.bam \
    H2228_ONTc.bam H2228_ONTd.bam H2228_PB.bam \
    H526_ONTc.bam H526_ONTd.bam H526_PB.bam \
    H211_ONTc.bam H211_ONTd.bam H211_PB.bam \
    SHP77_ONTc.bam SHP77_ONTd.bam SHP77_PB.bam \
    H146_ONTc.bam H146_ONTd.bam H146_PB.bam \
    H69_ONTc.bam H69_ONTd.bam H69_PB.bam \
    -t 24 -o "$OUTBASE"/nanocomp_violin/ --plot 'violin' --format 'svg' \
    --names LUAD_ONTc_H1975 LUAD_ONTd_H1975 LUAD_PB_H1975 \
    LUAD_ONTc_HCC827 LUAD_ONTd_HCC827 LUAD_PB_HCC827 \
    LUAD_ONTc_H2228 LUAD_ONTd_H2228 LUAD_PB_H2228 \
    SCLC-P_ONTc_H526 SCLC-P_ONTd_H526 SCLC-P_PB_H526 \
    SCLC-P_ONTc_H211 SCLC-P_ONTd_H211 SCLC-P_PB_H211 \
    SCLC-A_ONTc_SHP77 SCLC-A_ONTd_SHP77 SCLC-A_PB_SHP77 \
    SCLC-A_ONTc_H146 SCLC-A_ONTd_H146 SCLC-A_PB_H146 \
    SCLC-A_ONTc_H69 SCLC-A_ONTd_H69 SCLC-A_PB_H69

