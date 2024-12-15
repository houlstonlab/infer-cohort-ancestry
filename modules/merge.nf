process MERGE {
    tag "${ref}:${cohort}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/merged", mode: 'copy')

    input:
    tuple val(cohort), val(cohort_type),
          path(cohort_bim), path(cohort_bed), path(cohort_fam), path(cohort_nosex), path(cohort_log), path(cohort_pop),
    	  val(ref), val(ref_type),
          path(ref_bim), path(ref_bed), path(ref_fam), path(ref_nosex), path(ref_log), path(ref_pop)

    output:
    tuple val(ref), val(cohort),
          path("${ref}.${cohort}.bim"),
          path("${ref}.${cohort}.bed"),
          path("${ref}.${cohort}.fam"),
          path("${ref}.${cohort}.nosex"),
          path("${ref}.${cohort}.log"),
          path("${ref}.${cohort}.pop")

    script:
    """
    #!/bin/bash
    # Return population file
    cat ${ref_pop} ${cohort_pop} > ${ref}.${cohort}.pop
    cat ${ref}.${cohort}.pop | awk '{ print "0", \$2, \$1, \$2}' > famids.txt
    
    # Attempt to merge and identify problematic SNPs
    plink --bfile ${ref_bim.baseName} \
        --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
        --merge-mode 6 \
        --out merge_attempt || true
    
    if [ -f merge_attempt.missnp ]; then
        # Remove problematic SNPs from both datasets
        plink --bfile ${ref_bim.baseName} --exclude merge_attempt.missnp --make-bed --out ref_cleaned
        plink --bfile ${cohort_bim.baseName} --exclude merge_attempt.missnp --make-bed --out cases_cleaned

        # Attempt to merge with removing problematic SNPs
        plink --bfile ref_cleaned \
            --bmerge cases_cleaned.bed cases_cleaned.bim cases_cleaned.fam \
            --make-bed \
            --update-ids famids.txt \
            --maf ${params.AF} \
            --hwe ${params.HWE} \
            --mind ${params.F_MISSING} \
            --out ${ref}.${cohort}
    else
        plink --bfile ${ref_bim.baseName} \
            --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
            --make-bed \
            --update-ids famids.txt \
            --maf ${params.AF} \
            --hwe ${params.HWE} \
            --mind ${params.F_MISSING} \
            --out ${ref}.${cohort}
    fi
    """
}
