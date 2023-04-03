#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mem=50g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

sample=$(echo $1)

module load STAR

STAR --genomeDir /mainfs/hgig/private/RNAseq/STAR_index\
 --readFilesIn ${sample}_1.fq.gz,${sample}_2.fq.gz\
  --readFilesCommand zcat \
  --runThreadN 8\
   --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${sample}
