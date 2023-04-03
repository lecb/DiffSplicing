#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mem=32g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

ls -1 *R1*.gz | xargs gunzip -c | gzip > ${1}_1.fq.gz
ls -1 *R2*.gz | xargs gunzip -c | gzip > ${1}_2.fq.gz
