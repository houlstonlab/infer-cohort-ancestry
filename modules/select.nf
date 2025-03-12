process SELECT {
    tag "${ref}:${cohort}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/selected", mode: 'copy')

    input:
    tuple val(ref), val(cohort),
          path(bim), path(bed), path(fam), path(nosex), path(log),
          path(pop)

    output:
    tuple val(ref), val(cohort),
          path("${ref}.${cohort}.selected.bim"),
          path("${ref}.${cohort}.selected.bed"),
          path("${ref}.${cohort}.selected.fam"),
          path("${ref}.${cohort}.selected.nosex"),
          path("${ref}.${cohort}.selected.log"),
          path(pop)

    script:
    """
    #!/bin/bash
    # Filter variants
    plink --bfile ${bim.baseName} \
        --write-snplist \
        --out selected
    
    # Select N_VARS random variants
    RANDOM=42; shuf -n ${params.N_VARS} selected.snplist > ${ref}.${cohort}.variants.txt

    plink --bfile ${bim.baseName} \
        --extract ${ref}.${cohort}.variants.txt \
        --make-bed \
        --out ${ref}.${cohort}.selected
    """
}
