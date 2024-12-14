#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4096

# load trim-galore
module load TrimGalore/0.6.10 

# run fastQC
mkdir -p /projects/e30836/project/group7/bulk/trimmed

for sra_id in $(cat sra_ids_bulk.txt); do

    trim_galore --quality 20 --length 20 --fastqc "/projects/e30836/project/group7/bulk/data/${sra_id}_1.fastq" -o /projects/e30836/project/group7/bulk/trimmed
done