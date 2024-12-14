#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=20G

# load STAR v2.6
module load STAR/2.6.0
 
# make your working directory
mkdir -p /projects/e30836/project/group7/bulk/aligned

# A command you actually want to execute:
STAR --runThreadN 10\
    --quantMode TranscriptomeSAM \
    --genomeDir /projects/e30836/hw1/hg38.index/STAR \
    --readFilesIn "/projects/e30836/project/group7/bulk/trimmed/SRR13118820_1_trimmed.fq" \
    --outFileNamePrefix "/projects/e30836/project/group7/bulk/aligned/SRR13118820_"
done