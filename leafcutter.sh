#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=16g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

module load biobuilds

module load conda
source activate /home/egs1g09/.conda/envs/LC

#intron clustering
python /home/egs1g09/leafcutter/clustering/leafcutter_cluster_regtools.py -j DRAGON_junctiles_paths.txt -m 50 -o DRAGON_LC -l 500000
