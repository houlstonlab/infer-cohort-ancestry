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
    cat \
        <(cat ${ref_pop}    | awk '{print \$0, "${ref_type}"}') \
        <(cat ${cohort_pop} | awk '{print \$0, "${cohort_type}"}') \
        > ${ref}.${cohort}.pop
    
    # Attempt to merge and identify problematic SNPs
    plink --bfile ${ref_bim.baseName} \
        --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
        --out merge_attempt || true
    
    if [ -f merge_attempt.missnp ]; then
        # Remove problematic SNPs from both datasets
        plink --bfile ${ref_bim.baseName} --exclude merge_attempt.missnp --write-snplist --out ref_cleaned
        plink --bfile ${cohort_bim.baseName} --exclude merge_attempt.missnp --write-snplist --out cases_cleaned
        
        # Common SNPs from both datasets
        comm -12 ref_cleaned.snplist cases_cleaned.snplist > common_snps.txt

        # Attempt to merge with removing problematic SNPs
        plink --bfile ${ref_bim.baseName} \
            --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
            --extract common_snps.txt \
            --make-bed \
            --out ${ref}.${cohort}
    else
        # Common SNPs from both datasets
        plink --bfile ${ref_bim.baseName} --write-snplist --out ref_cleaned
        plink --bfile ${cohort_bim.baseName} --write-snplist --out cases_cleaned

        comm -12 ref_cleaned.snplist cases_cleaned.snplist > common_snps.txt

        plink --bfile ${ref_bim.baseName} \
            --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
            --make-bed \
            --out ${ref}.${cohort}
    fi
    """
}
