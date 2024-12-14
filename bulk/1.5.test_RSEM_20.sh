#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=30G

# load RSEM
module load rsem/1.3.0 

mkdir -p /projects/e30836/project/group7/bulk/expression/test

# Loop through SRA IDs
rsem-calculate-expression -p 10 \
    --bam /projects/e30836/project/group7/bulk/aligned/SRR13118820_Aligned.toTranscriptome.out.bam \
    /projects/e30836/hw1/hg38.index/STAR/gencode.v28 \
    /projects/e30836/project/group7/bulk/expression/test/SRR13118820
