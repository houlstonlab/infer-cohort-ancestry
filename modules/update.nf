process UPDATE {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/updated", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path(vcf_in), path(index_in), 
          val(dbsnp), path(dbsnp_in), path(dbsnp_index)

    output:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path("${cohort}.${type}.${chrom}.updated.vcf.gz"),
          path("${cohort}.${type}.${chrom}.updated.vcf.gz.tbi")
     
    script:
    """
    #!/bin/bash
    bcftools annotate \
        -a ${dbsnp_in} \
        -c ID \
        ${vcf_in} \
        --threads ${task.cpu} \
        -Oz -o ${cohort}.${type}.${chrom}.updated.vcf.gz
    
    tabix ${cohort}.${type}.${chrom}.updated.vcf.gz
    n_vars=\$(bcftools index -n ${cohort}.${type}.${chrom}.updated.vcf.gz)
    """
}
