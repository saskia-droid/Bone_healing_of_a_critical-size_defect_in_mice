#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=03:00:00
#SBATCH --job-name=index_STAR
#SBATCH --mail-user=saskia.perret-gentil@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/sperretgentil/Bone/p03_Franziska_TruSeq/outputs_and_errors/%j.out
#SBATCH --error=/home/sperretgentil/Bone/p03_Franziska_TruSeq/outputs_and_errors/%j.err

module load UHTS/Aligner/STAR/2.7.10a_alpha_220818

FASTA_FILE=/data/references/Mus_musculus/Ensembl/GRCm39/Sequence/Mus_musculus.GRCm39.dna.primary_assembly.fa
GTF_FILE=/data/references/Mus_musculus/Ensembl/GRCm39/Annotation/Genes/build107/Mus_musculus.GRCm39.107.gtf
OUTPUT_DIR=/data/users/sperretgentil/Franziska_data/Mus_musculus_GRCm39_index_STAR

STAR --runMode genomeGenerate --genomeDir ${OUTPUT_DIR} --genomeFastaFiles ${FASTA_FILE} --sjdbGTFfile ${GTF_FILE} --runThreadN 8

