#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=35:00:00
#SBATCH --job-name=gene_counts_STAR
#SBATCH --mail-user=saskia.perret-gentil@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/sperretgentil/Bone/p03_Franziska_TruSeq/outputs_and_errors/%j.out
#SBATCH --error=/home/sperretgentil/Bone/p03_Franziska_TruSeq/outputs_and_errors/%j.err
#SBATCH --array=1-96

module load UHTS/Aligner/STAR/2.7.10a_alpha_220818

read_1=/data/projects/p639_rnaseq_hofstetter/03_Franziska_TruSeq/raw/mRNA_${SLURM_ARRAY_TASK_ID}_R1.fastq.gz
read_2=/data/projects/p639_rnaseq_hofstetter/03_Franziska_TruSeq/raw/mRNA_${SLURM_ARRAY_TASK_ID}_R2.fastq.gz

INDEX_DIR=/data/users/sperretgentil/Franziska_data/Mus_musculus_GRCm39_index_STAR

OutFilePathPrefix=/data/users/sperretgentil/Franziska_data/TruSeq_mapping_counts/S${SLURM_ARRAY_TASK_ID}_

STAR --quantMode GeneCounts --genomeDir ${INDEX_DIR} --runThreadN 8 --readFilesIn ${read_1} ${read_2} --readFilesCommand zcat --outFileNamePrefix ${OutFilePathPrefix} --outSAMtype BAM SortedByCoordinate
