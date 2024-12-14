#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4096

# load fastQC
module load fastqc/0.12.0

# run fastQC
mkdir -p /projects/e30836/project/group7/bulk/qc_reports

for sra_id in $(cat sra_ids_bulk.txt); do


     fastqc "/projects/e30836/project/group7/bulk/data/${sra_id}_1.fastq" -o /projects/e30836/project/group7/bulk/qc_reports
done