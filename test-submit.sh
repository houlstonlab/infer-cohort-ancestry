#!/bin/bash

#SBATCH -o test/test.out
#SBATCH -e test/test.err
#SBATCH -J test
#SBATCH -p master-worker
#SBATCH -t 120:00:00

# Setup test directory
mkdir -p test/ test/input
cd test/

# # Download test data
# URL="https://figshare.com/ndownloader/files"

# wget -c $URL/52216769 -O input/vcf-ancestry.tar.gz

# # Unzip the files
# tar -xzvf input/vcf-ancestry.tar.gz -C input/
# cp /data/reference-data/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta input/assembly38.fasta
# cp /data/reference-data/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai input/assembly38.fasta.fai 

# Run nextflow
module load Nextflow

# nextflow run houlstonlab/infer-cohort-ancestry -r main \
nextflow run ../main.nf \
    --output_dir ./results/ \
    -profile local,test \
    -resume

# usage: nextflow run [ local_dir/main.nf | git_url ]  
# These are the required arguments:
#     -r            {main,dev} to run specific branch
#     -profile      {local,cluster} to run using differens resources
#     -params-file  params.json to pass parameters to the pipeline
#     -resume       To resume the pipeline from the last checkpoint

mv .nextflow.log nextflow.log
