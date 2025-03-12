process FILTER {
    tag "${ref}:${cohort}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/filtered", mode: 'copy')

    input:
    tuple val(ref), val(cohort),
          path(bim), path(bed), path(fam), path(nosex), path(log),
          path(pop),
          path(snplist)

    output:
    tuple val(ref), val(cohort),
          path("${ref}.${cohort}.filtered.bim"),
          path("${ref}.${cohort}.filtered.bed"),
          path("${ref}.${cohort}.filtered.fam"),
          path("${ref}.${cohort}.filtered.nosex"),
          path("${ref}.${cohort}.filtered.log"),
          path(pop)

    script:
    """
    #!/bin/bash
    # Filter variants
    plink --bfile ${bim.baseName} \
        --maf ${params.MAF} \
        --hwe ${params.HWE} \
        --geno ${params.F_MISSING} \
        --extract ${snplist} \
        --make-bed \
        --out ${ref}.${cohort}.filtered
    """
}
