#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4096

# load CellRanger for generate ref genome
module load cellranger/7.1.0

# change working directory
cd /projects/e30836/project/group7/ref_genome/hg38/

# make gtf of hg38 for cellranger
cellranger mkgtf gencode.v47.annotation.gtf  gencode.v47.annotation.filtered.gtf --attribute=gene_biotype:protein_coding
