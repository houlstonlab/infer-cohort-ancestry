process CONVERT {
    tag "${cohort}:${type}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/plinked", mode: 'copy')

    input:
    tuple val(cohort), val(type),
          path(vcf_in), path(index_in)

    output:
    tuple val(cohort), val(type),
          path("${cohort}.${type}.bim"),
          path("${cohort}.${type}.bed"),
          path("${cohort}.${type}.fam"),
          path("${cohort}.${type}.nosex"),
          path("${cohort}.${type}.log")

    script:
    """
    #!/bin/bash
    plink \
        --vcf ${vcf_in} \
        --make-bed \
        --const-fid 0 \
        --out ${cohort}.${type}
    """
}