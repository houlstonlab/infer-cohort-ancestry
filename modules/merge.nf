process MERGE {
    tag "${ref}:${cohort}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/merged", mode: 'copy')

    input:
    tuple val(cohort), val(cohort_type),
          path(cohort_bim), path(cohort_bed), path(cohort_fam), path(cohort_nosex),
          path(cohort_in), path(cohort_out), path(cohort_log)
    tuple val(ref), val(ref_type),
          path(ref_bim), path(ref_bed), path(ref_fam), path(ref_nosex),
          path(ref_in), path(ref_out), path(ref_log)

    output:
    tuple val(ref), val(cohort),
          path("${ref}.snplist"),
          path("${cohort}.snplist"),
          path("${ref}.${cohort}.snplist"),
          path("${ref}.${cohort}.bim"),
          path("${ref}.${cohort}.bed"),
          path("${ref}.${cohort}.fam"),
          path("${ref}.${cohort}.nosex"),
          path("${ref}.${cohort}.log")

    script:
    """
    #!/bin/bash
    # Write the list of SNPs in the reference and cohort files to text files
    plink --bfile ${ref_bim.baseName} --write-snplist | sort -V > ${ref}.snplist
    plink --bfile ${cohort_bim.baseName} --write-snplist | sort -V > ${cohort}.snplist
    
    # Find the intersection of the lists of SNPs
    grep -wFf ${ref}.snplist ${cohort}.snplist > ${ref}.${cohort}.snplist

    # Merge the reference and cohort files, keeping only the intersection variants
    plink --bfile ${ref_bim.baseName} \
        --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} \
        --make-bed \
        --out ${ref}.${cohort}
    """
}
