#!/bin/bash
#SBATCH -A e30836
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=20G

# load CellRanger for generate ref genome
module load cellranger/7.1.0

# Set paths
REFERENCE="/projects/e30836/project/group7/ref_genome/hg38/cellranger_GRCh38/"
DATA_DIR="/projects/b1080/zrt/bme311/project/sc/data/"
OUTPUT_DIR="/projects/e30836/project/group7/sc/cellranger"  # Define output directory
SAMPLE_LIST="/projects/e30836/project/group7/scripts/sra_ids_sc.txt"  

# make output dir
mkdir -p $OUTPUT_DIR

# Run cellranger count
for sra_id in $(cat $SAMPLE_LIST); do
  # Create individual output directory
  mkdir -p "$OUTPUT_DIR/$sra_id"

  # Change to the output directory for the sample
  cd "$OUTPUT_DIR/$sra_id" || exit 1


  cellranger count \
    --id="$sra_id" \
    --transcriptome="$REFERENCE" \
    --fastqs="$DATA_DIR" \
    --sample="$sra_id" \
    --localcores=8 \

  # Return to the script's original directory
  cd - || exit 1
done

