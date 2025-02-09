### Introduction

This workflow infers the genetic ancestry of a cohort of individuals. 
The workflow selects informative variants from the cohort and a reference 
dataset, and calculates the prinicipal components representing the variance in
genetic backgrounds. The workflow is designed and tested on a cohort generated 
from the 1KG project.

### Usage

Different versions of the workflow can be called using `-r` and output directed
to `--output_dir`. The typical command looks like the following. 

```bash
nextflow run houlstonlab/infer-cohort-ancestry \
    -r main \
    --output_dir results/ \
    --cohorts input/cohort.info.csv \
    --fasta input/Homo_sapiens_assembly38.fasta \
    --dbsnp input/dbsnp_146.hg38.vcf.gz{,.tbi} \
    --ld_regions input/high_ld_regions.txt
```

### Inputs & Parameters

- `cohorts`: a csv file with 6 columns: 
    - `cohort`: cohort name
    - `type`: one of `'cases'` or `'reference'`
    - `size`: the number of individuals in the cohort
    - `vars_file`: a vcf file
    - `vars_index`: an index file
    - `population`: a tab delimiter file. with 4 columns: family id, sample id,
    population, and super population assignments. For samples whoes ancestry is
    to be inferred, the two last columns should be assinged as `NA`.

- References:
    - `fasta`: assembly fasta file
    - `dbsnp`: dbsnp VCF file and index
    - `ld_regions`: bed file of genomic regions to exclude

- Parameters
    - `AF`: allele frequency. Default is 0
    - `HWE`: HWE test p-value threshold. Default is 0
    - `F_MISSING`: the fraction of missing. Default is 1
    - `window`: the size of the window in bp to consider for pruning SNPs in LD.
    Default is 5. 
    - `step`: the size of the step in bp. Default is 5 
    - `rsquared`: the R^2 value threshold. Default is 0.00001
    - `relatedness`: a cutoff to remove related individuals. Default is 0.000025
    - `N_VARS`: the number of SNPs to use. Default is 5000.
    - `N_DIMS`: the number of dimenstion to use. Default is 2
    - `modes`: one or more of `'clusters,noclusters,mds'`. The first two are the
    PCA analysis within or without clusters/populations.

### Output

- `cohorts/` : subsets of the VCF files
- `updated/`  : updated VCF files with variant IDs
- `filtered/` : filtered VCF files
- `plinked/`  : the selected variants in plink format
- `removed/`  : ambigious calls removed 
- `updated/`  : updated alleles as in dbsnp
- `fixed/`    : fixed reference allele
- `combined/` : combined plink files
- `merged/`   : merged cases and references in one plink file
- `pruned/`   : pruned variants
- `scaled/`   : principal components or MDS
- `assigned/` : assinged ancestry
- `plots/`    : dimension reduction plots
