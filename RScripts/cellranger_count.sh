#!/bin/bash
# Bash script for running cellranger count

cellranger count \
  --id=group7_sc_output \
  --transcriptome=/projects/e30836/project/group7/ref_genome/hg38 \
  --fastqs=/projects/e30836/project/group7/sc/fastqs \
  --sample=sample_name \
  --localcores=8 \
  --localmem=64

