#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=16g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

module load biobuilds

module load conda
source activate /home/egs1g09/.conda/envs/LC

Rscript /home/egs1g09/leafcutter/scripts/leafcutter_ds.R -e /mainfs/hgig/private/RNAseq/ACCORD/annotation_codes/gencode_hg38/gencode_hg38_all_exons.txt.gz --num_threads 4 DRAGON_LC_perind_numers.counts.gz  samples_with_o2_data_passed2
