process FILTER {
    tag "${ref}:${cohort}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/filtered", mode: 'copy')

    input:
    tuple val(ref), val(cohort),
          path(bim), path(bed), path(fam), path(nosex), path(log),
          path(pop)

    output:
    tuple val(ref), val(cohort),
          path("${ref}.${cohort}.bim"),
          path("${ref}.${cohort}.bed"),
          path("${ref}.${cohort}.fam"),
          path("${ref}.${cohort}.nosex"),
          path("${ref}.${cohort}.log"),
          path(pop)

    script:
    """
    #!/bin/bash
    # Filter variants
    plink --bfile ${bim.baseName} \
        --maf ${params.MAF} \
        --hwe ${params.HWE} \
        --geno ${params.F_MISSING} \
        --make-bed \
        --out filtered
    
    # Select N_VARS random variants
    RANDOM=42; shuf -n ${params.N_VARS} filtered.bim | \
    cut -f 2 | \
    uniq \
    > ${ref}.${cohort}.variants.txt

    plink --bfile filtered \
        --extract ${ref}.${cohort}.variants.txt \
        --make-bed \
        --out ${ref}.${cohort}
    """
}
