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
    --vcf "input/*.variants.vcf.gz{,.tbi}"
    --cases "input/*.cases.txt"
    --cohorts_info "input/cohorts_info.csv"
    --fasta "input/assembly38.fasta{,.fai}"
    --dbsnp "input/dbsnp.146.vcf.gz{,.tbi}"
    --ld_regions "input/ld_regions.txt"
    --populations "input/populations.txt"
    --populations_id "input/populations_id.txt"
    --populations_info "input/populations_info.txt"
    --clusters "input/clusters.txt"
```

### Inputs & Parameters

- Variants:
    - `vcf`: VCF files containing the variants of the cases and references
    - `cases`: a text file with cases per line to include from the VCFs
    - `cohorts_info`: a csv with cohort information. Expected to have three 
    columns `cohort`, `type` and `size`

- References:
    - `fasta`: assembly fasta file and index
    - `dbsnp`: dbsnp VCF file
    - `ld_regions`: bed file of genomic regions to exclude

- Population info
    - `populations`: a list of the different populations, one per line
    - `populations_id`: a tsv file with three columns: `FID`, `IID` and `Pop`
    - `populations_info`: a tsv file with further description of the populations
    - `clusters`: a tsv file with three columns: cluster, sample name, population

- Parameters
    - `AF`: allele frequency. Default is 0
    - `HWE`: HWE test p-value threshold. Default is 0
    - `F_MISSING`: the fraction of missing. Default is 1
    - `window`: the size of the window in bp to consider for pruning SNPs in LD 
    - `step`: the size of the step in bp 
    - `rsquared`: the R^2 value threshold
    
### Output

- `cohorts/`: subsets of the VCF files
- `updated/`: updated VCF files with variant IDs
- `filtered/`: filtered VCF files
- `plinked/`: the selected variants in plink format
- `merged/`: merged cases and references in one plink file
- `pruned/`: pruned variants
- `pca/`: principal components
- `plots/`: principal components plots
