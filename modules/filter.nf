process FILTER {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/filtered", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom),
          path(vcf_in), path(index_in) 

    output:
    tuple val(cohort), val(type), val(chrom),
          path("${cohort}.${type}.${chrom}.filtered.vcf.gz"),
          path("${cohort}.${type}.${chrom}.filtered.vcf.gz.tbi")

    script:
    """
    #!/bin/bash
    bcftools view ${vcf_in} | \
    bcftools filter -e 'ID="."' | \
    bcftools norm -d none | \
    bcftools view -m2 -M2 | \
    bcftools filter -i 'AF > ${params.AF} && HWE > ${params.HWE} && F_MISSING < ${params.F_MISSING}' | \
    bcftools view --threads ${task.cpu} -Oz -o ${cohort}.${type}.${chrom}.filtered.vcf.gz

    tabix ${cohort}.${type}.${chrom}.filtered.vcf.gz
    """
}
