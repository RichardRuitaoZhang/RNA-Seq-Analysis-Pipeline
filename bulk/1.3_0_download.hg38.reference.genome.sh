#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4096

# make directory
mkdir -p /projects/e30836/project/group7/ref_genome/hg38
# change directory
cd /projects/e30836/project/group7/ref_genome/hg38

# download genome fasta file and unzip
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.p14.genome.fa.gz
gunzip GRCh38.p14.genome.fa.gz

# download gtf annotation file and unzip
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
gunzip gencode.v47.annotation.gtf.gz
