process UPDATE {
    tag "${cohort}:${type}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/updated", mode: 'copy')

    input:
    tuple val(cohort), val(type),
          path(vcf_in), path(index_in), 
          val(dbsnp), path(dbsnp_in), path(dbsnp_index)

    output:
    tuple val(cohort), val(type),
          path("${cohort}.${type}.updated.vcf.gz"),
          path("${cohort}.${type}.updated.vcf.gz.tbi")
     
    script:
    """
    #!/bin/bash
    bcftools annotate \
        -a ${dbsnp_in} \
        -c ID \
        ${vcf_in} \
        --threads ${task.cpu} \
        -Oz -o ${cohort}.${type}.updated.vcf.gz
    
    tabix ${cohort}.${type}.updated.vcf.gz
    """
}
