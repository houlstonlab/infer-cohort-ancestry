process SUBSET {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/cohorts", mode: 'copy')

    input:
    tuple val(cohort), path(vcf_in), path(index_in),
          val(type), val(size), path(cases), 
          val(chrom), path(coordinates)

    output:
    tuple val(cohort), val(type), val(chrom),
          path("${cohort}.${type}.${chrom}.snps.vcf.gz"),
          path("${cohort}.${type}.${chrom}.snps.vcf.gz.tbi")
 
    script:
    """
    #!/bin/bash
    # Extract cases
    cat ${cases} | awk '{ print \$2 }' | uniq > cases.txt

    # Extract coordinates
    cat ${coordinates} | sort | uniq -c | awk '{ if(\$1 > 2) print \$2"\t"\$3}' > regions.txt
    
    # Extract SNPs
    bcftools view -S cases.txt -R regions.txt ${vcf_in} | \
    bcftools norm -d none | \
    bcftools view -v snps -m2 -M2 | \
    bcftools +setGT -- -t . -n 0 | \
    bcftools view --threads ${task.cpu} -Oz -o ${cohort}.${type}.${chrom}.snps.vcf.gz
    tabix ${cohort}.${type}.${chrom}.snps.vcf.gz
    """
}
