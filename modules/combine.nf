process COMBINE {
    tag "${cohort}:${type}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/combined", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom),
          path(bim), path(bed), path(fam), path(nosex), path(in), path(out), path(log)

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
    # Create a list of all files
    echo "${bed.join('\n')}" > bed.txt
    echo "${bim.join('\n')}" > bim.txt
    echo "${fam.join('\n')}" > fam.txt
    paste -d ' ' bed.txt bim.txt fam.txt | sort -V > allfiles.txt

    plink \
        --make-bed \
        --nonfounders \
        --merge-list allfiles.txt \
        --out ${cohort}.${type}
    """
}
