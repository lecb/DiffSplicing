#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=40:00:00
#SBATCH --mem=8g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

#sbatch --array=1-17 regtools.sh

module load biobuilds

module load conda
source activate /home/egs1g09/.conda/envs/LECB #installation of regtools

while read -r line; do
    sample=$(echo "$line"| cut -f1)
    bam=$(echo "$line" | cut -f2)
    cd $sample
    echo Converting $sample to $sample.junc
    samtools index $bam
    regtools junctions extract -s 2 -a 8 -m 50 -M 500000 $bam -o ${sample}.junc
    echo ${sample}.junc >>../DRAGON_junctiles_${SLURM_ARRAY_TASK_ID}.txt
    cd ../
done < chunks${SLURM_ARRAY_TASK_ID}

conda deactivate
