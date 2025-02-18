process CONVERT {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/plinked", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path(vcf_in), path(index_in)

    output:
    tuple val(cohort), val(type), val(chrom),
          path("${cohort}.${type}.${chrom}.bim"),
          path("${cohort}.${type}.${chrom}.bed"),
          path("${cohort}.${type}.${chrom}.fam"),
          path("${cohort}.${type}.${chrom}.nosex"),
          path("${cohort}.${type}.${chrom}.log")

    script:
    """
    #!/bin/bash
    echo '.' > tmp.exclude
    plink \
        --vcf ${vcf_in} \
        --make-bed \
        --const-fid 0 \
        --exclude tmp.exclude \
        --fill-missing-a2 \
        --out ${cohort}.${type}.${chrom}
    """
}