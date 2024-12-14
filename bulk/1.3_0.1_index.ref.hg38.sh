#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=15G

# load STAR v2.7
module load STAR/2.7.5

# change directory
cd /projects/e30836/project/group7/ref_genome/hg38

STAR --runThreadN 10 \
     --runMode genomeGenerate \
     --genomeDir star_index/ \
     --genomeFastaFiles GRCh38.p14.genome.fa \
     --sjdbGTFfile gencode.v47.annotation.gtf \
     --sjdbOverhang 99