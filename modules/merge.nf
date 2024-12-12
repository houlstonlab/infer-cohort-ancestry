process MERGE {
    tag "${ref}:${cohort}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/merged", mode: 'copy')

    input:
    tuple val(cohort), val(cohort_type),
          path(cohort_bim), path(cohort_bed), path(cohort_fam), path(cohort_nosex), path(cohort_log),
    	  val(ref), val(ref_type),
          path(ref_bim), path(ref_bed), path(ref_fam), path(ref_nosex), path(ref_log)

    output:
    tuple val(ref), val(cohort),
          path("${ref}.${cohort}.bim"),
          path("${ref}.${cohort}.bed"),
          path("${ref}.${cohort}.fam"),
          path("${ref}.${cohort}.nosex"),
          path("${ref}.${cohort}.log")

    script:
    """
    #!/bin/bash
    # Merge the reference and cohort files, keeping only the intersection variants
    plink --bfile ${ref_bim.baseName} \
        --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
        --make-bed \
        --out ${ref}.${cohort}
    """
}
