#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=20G

# load RSEM
module load rsem/1.3.0 

# change working dir
cd /projects/e30836/project/group7/bulk/expression

# construct expression matrix
# it can be done by gene output
head -n 9 sra_ids_bulk.txt | xargs -I {} echo {}.genes.results | xargs rsem-generate-data-matrix > pc9_expression_matrix.txt  # generate matrix for PC9 cells
