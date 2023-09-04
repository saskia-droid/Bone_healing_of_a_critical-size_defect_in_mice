#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --output=/home/sperretgentil/Bone/p03_Franziska_TruSeq/outputs_and_errors/%j.out
#SBATCH --error=/home/sperretgentil/Bone/p03_Franziska_TruSeq/outputs_and_errors/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=30:00:00
#SBATCH --mail-user=saskia.perret-gentil@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

# add software to environment
module add UHTS/Quality_control/fastqc/0.11.9

READ_DIR=/data/projects/p639_rnaseq_hofstetter/03_Franziska_TruSeq/raw
OUT_DIR=/home/sperretgentil/Bone/p03_Franziska_TruSeq/results/reads_quality

fastqc --outdir $OUT_DIR $READ_DIR/*.fastq.gz --threads 1
