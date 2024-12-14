#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4096

# load modules you need to use
module load samtools/1.6

# make directory
mkdir -p /projects/e30836/project/group7/bulk/bam
mkdir -p /projects/e30836/project/group7/bulk/bam/sorted

# convert sam to bam
for sra_id in $(cat sra_ids_bulk.txt); do
    samtools view -b /projects/e30836/project/group7/bulk/aligned/${sra_id}_Aligned.out.sam > /projects/e30836/project/group7/bulk/bam/${sra_id}_Aligned.out.bam
done

# sort
for sra_id in $(cat sra_ids_bulk.txt); do
    samtools sort /projects/e30836/project/group7/bulk/bam/${sra_id}_Aligned.out.bam > /projects/e30836/project/group7/bulk/bam/sorted/${sra_id}_Aligned.out.sort.bam
done

# index
for sra_id in $(cat sra_ids_bulk.txt); do
    samtools index /projects/e30836/project/group7/bulk/bam/sorted/${sra_id}_Aligned.out.sort.bam
done