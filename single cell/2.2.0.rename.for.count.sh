#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=4096

# Download single-cell RNA-Seq datasets and convert into fastq files
mkdir -p /projects/e30836/project/group7/sc/rename

# rename the fastq files to meet cellranger readin
for sra_id in $(cat sra_ids_sc.txt); do

    cp /projects/b1080/zrt/bme311/project/sc/data/data/${sra_id}_1.fastq /projects/b1080/zrt/bme311/project/sc/data/${sra_id}_S1_L001_I1_001.fastq.gz
    cp /projects/b1080/zrt/bme311/project/sc/data/data/${sra_id}_2.fastq /projects/b1080/zrt/bme311/project/sc/data/${sra_id}_S1_L001_R1_001.fastq.gz
    cp /projects/b1080/zrt/bme311/project/sc/data/data/${sra_id}_3.fastq /projects/b1080/zrt/bme311/project/sc/data/${sra_id}_S1_L001_R2_001.fastq.gz
done