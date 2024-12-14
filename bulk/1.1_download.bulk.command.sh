#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4096

# Load sra toolkit modules to download datasets
module load sratoolkit/3.0.0

# Download bulk RNA-Seq datasets and convert into fastq files
mkdir -p /projects/e30836/project/group7/bulk/data

for sra_id in $(cat sra_ids_bulk.txt); do


    fastq-dump -I --split-files $sra_id -O /projects/e30836/project/group7/bulk/data
done