#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=20G

# load CellRanger for generate ref genome
module load cellranger/7.1.0

cd /projects/e30836/project/group7/ref_genome/hg38

cellranger mkref \
    --genome=cellranger_GRCh38 \
    --fasta=GRCh38.p14.genome.fa \
    --genes=gencode.v47.annotation.filtered.gtf \
    --ref-version=v47