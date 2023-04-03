#################################################################################
#### Differential splicing workflow for RNA-seq DRAGON data; Ellie Seaby '23 ####
#################################################################################

#Move to DRAGON working directory
cd /mainfs/hgig/private/RNAseq/DRAGON

#NB: There is a back up of all DRAGON fastqs in /research/HGIGData/RNASEQ/DRAGON/
#Paths to all fastqs are here: /mainfs/hgig/private/RNAseq/DRAGON/fastqs_DRAGON

cp /mainfs/hgig/private/RNAseq/DRAGON/fastqs_DRAGON sample.map

#copy all fastqs into new directories
while read -r line; do
    sample=$(echo "$line"| cut -f1)
    file=$(echo "$line" | cut -f2)
    mkdir $sample
    cd $sample
    \cp $file .
    cd ../
done < sample.map

#Make unique sample list
cut -f1 sample.map | sort | uniq > uniq_samples

#################################################################################

while read -r line; do  \cp -f concat_fastqs.sh ./$line; cd $line; sbatch concat_fastqs.sh $line; cd /mainfs/hgig/private/RNAseq/DRAGON; done < uniq_samples

#################################################################################

#concat_fastqs.sh

#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mem=32g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

ls -1 *R1*.gz | xargs gunzip -c | gzip > ${1}_1.fq.gz
ls -1 *R2*.gz | xargs gunzip -c | gzip > ${1}_2.fq.gz

#################################################################################

#run_STAR_align.sh
while read -r line; do \cp -f STAR_align.sh ./$line; cd $line; sbatch STAR_align.sh $line; cd /mainfs/hgig/private/RNAseq/DRAGON; done < uniq_samples

#################################################################################

#star_align.sh

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

#################################################################################

#run_sort_index.sh
while read -r line; do \cp sort_index.sh ./$line; cd $line; sbatch sort_index.sh $line; cd /mainfs/hgig/private/RNAseq/DRAGON; done < uniq_samples

#################################################################################

#sort_index.sh

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

#################################################################################

cd /mainfs/hgig/private/RNAseq/DRAGON

find . -name \*sorted.bam > sorted_bams
cut -f2 -d '/' sorted_bams > samples
cut -f3 -d'/' sorted_bams > tmp
#cut -f2,3 -d'_' tmp > samples
paste samples tmp > samples_sorted_bams

#extract only those with O2 data
grep -wf samples_o2_passed samples_sorted_bams >samples_sorted_bams_o2

#split into chunks
split -l 7 --numeric-suffixes samples_sorted_bams_o2 chunks

#rename chunks0 as last chunk

sbatch --array=1-17 regtools.sh

#################################################################################

#regtools.sh
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

###############################################################################

#concat junctile paths, add paths and save as DRAGON_junctiles_paths.txt
cat DRAGON_junctiles*txt > all_junctiles
cut -d'.' -f1 all_junctiles > a
sed 's:$:/:' a >a2
paste a2 all_junctiles > tmp
perl -i -pe 's/[ \t]+//' tmp
sed 's:^:/mainfs/hgig/private/RNAseq/DRAGON/:' tmp > DRAGON_junctiles_paths.txt
rm tmp a a2

###############################################################################

#Installling leafcutter is a HUGE pain  - please use following instructions for conda install

conda env create -n LC
conda activate LC
conda install -y r-base=3.6
conda install -c r r-rstan
conda install -c bio
conda samtools
conda install -c bio
conda regtools
conda install -c conda-forge r-devtools
conda install -c bioconda bioconductor-biobase
conda install -c bioconda bioconductor-dirichletmultinomial
conda install r-recommended r-irkernel
conda install -c conda-forge r-r.utils
conda install -c conda-forge r-hmisc
conda install -c conda-forge r-domc
conda install -c conda-forge r-optparse
conda install -c conda-forge r-intervals
conda install -c conda-forge r-shinycssloaders
conda install -c r r-bh
conda install -c bio
conda bioconductor-biobaseR -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter")'

sbatch leafcutter.sh

###############################################################################
#leafcutter.sh
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

#Differential intron excision analysis # need to add in groups file
# . e.g.
# Sample_001-125027_D1	SoC
#Sample_002-125027_D5	SoC
#Sample_003-125028_D1	SoC
#Sample_004-125028_D5	SoC
#Sample_006-125029_D1	S003
#Sample_008-125025_D1	S003
#Sample_009-125025_D5	S003
#Sample_010-125026_D1	S003

###############################################################################
#leafcutterR.sh

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

###############################################################################
#LC_cluster.sh

#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=16g
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

module load biobuilds

module load conda
source activate /home/egs1g09/.conda/envs/LC

Rscript /home/egs1g09/leafcutter/scripts/ds_plots.R -e /mainfs/hgig/private/RNAseq/ACCORD/annotation_codes/gencode_hg38/gencode_hg38_all_exons.txt.gz DRAGON_LC_perind_numers.counts.gz samples_with_o2_data_passed2  significant_results.txt -f 0.05
