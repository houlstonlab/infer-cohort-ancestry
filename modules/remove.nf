process REMOVE {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/removed", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom),
          path(vcf_in), path(index_in)

    output:
    tuple val(cohort), val(type), val(chrom),
          path("${cohort}.${type}.${chrom}.removed.vcf.gz"),
          path("${cohort}.${type}.${chrom}.removed.vcf.gz.tbi")
     
    script:
    """
    #!/bin/bash
    # Remove A/T and G/C SNPs
    bcftools view \
        -e '(REF="A" & ALT="T") || (REF="G" & ALT="C")' \
        ${vcf_in} \
        --threads ${task.cpu} \
        -Oz -o ${cohort}.${type}.${chrom}.removed.vcf.gz
    
    tabix ${cohort}.${type}.${chrom}.removed.vcf.gz
    """
}
