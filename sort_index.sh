#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=16g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

sample=$(echo $1)

module load samtools
samtools sort -o ${sample}_sorted.bam ${sample}Aligned.out.bam

samtools index ${sample}_sorted.bam

rm ${sample}Aligned.out.bam
