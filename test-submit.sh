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

# wget -c $URL/50487621 -O input/pheno.variants.vcf.gz
# wget -c $URL/50487624 -O input/pheno.variants.vcf.gz.tbi
# wget -c $URL/50385591 -O input/pheno.cases.txt
# wget -c $URL/50810508 -O input/cohorts_info.csv

# wget -c $URL/50791842 -O input/dbsnp.146.vcf.gz
# wget -c $URL/50791845 -O input/dbsnp.146.vcf.gz.tbi

# wget -c $URL/50793063 -O input/assembly38.fasta
# wget -c $URL/50793066 -O input/assembly38.fasta.fai
# wget -c $URL/50812212 -O input/ld_regions.txt

# wget -c $URL/50813238 -O input/populations.txt
# wget -c $URL/50813415 -O input/populations_id.txt
# wget -c $URL/50813337 -O input/populations_info.txt
# wget -c $URL/50813235 -O input/clusters.txt

# # Copy as reference files
# cp input/pheno.variants.vcf.gz input/ref.variants.vcf.gz
# cp input/pheno.variants.vcf.gz.tbi input/ref.variants.vcf.gz.tbi
# cp input/pheno.cases.txt input/ref.cases.txt

# Run nextflow
module load Nextflow

# nextflow run houlstonlab/infer-cohort-ancestry -r main \
nextflow run ../main.nf \
    --output_dir ./results/ \
    -profile local,gha \
    -resume

# usage: nextflow run [ local_dir/main.nf | git_url ]  
# These are the required arguments:
#     -r            {main,dev} to run specific branch
#     -profile      {local,cluster} to run using differens resources
#     -params-file  params.json to pass parameters to the pipeline
#     -resume       To resume the pipeline from the last checkpoint

mv .nextflow.log nextflow.log
