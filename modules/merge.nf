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
    # Attempt to merge and identify problematic SNPs
    plink --bfile ${ref_bim.baseName} \
        --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
        --merge-mode 6 \
        --out merge_attempt || true
    
    # Remove problematic SNPs from both datasets
    plink --bfile ${ref_bim.baseName} --exclude merge_attempt.missnp --make-bed --out ref_cleaned
    plink --bfile ${cohort_bim.baseName} --exclude merge_attempt.missnp --make-bed --out cases_cleaned

    # Attempt to merge with removing problematic SNPs
    plink --bfile ref_cleaned \
        --bmerge cases_cleaned.bed cases_cleaned.bim cases_cleaned.fam \
        --make-bed \
        --maf ${params.AF} \
        --hwe ${params.HWE} \
        --mind ${params.F_MISSING} \
        --out ${ref}.${cohort}
    """
}
