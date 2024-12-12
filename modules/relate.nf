process RELATE {
    tag "${cohort}:${type}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/relatedness", mode: 'copy')

    input:
    tuple val(cohort), val(type),
          path(bim), path(bed), path(fam), path(nosex), path(log)

    output:
    tuple val(cohort), val(type),
          path("${cohort}.${type}.genome")

    script:
    """
    #!/bin/bash
    plink \
        --bfile ${bim.baseName} \
        --genome \
        --out ${cohort}.${type}
    """
}
