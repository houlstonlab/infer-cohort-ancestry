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
    
    # Extract snp lists, missnps and flipped snps
    plink --bfile ${ref_bim.baseName} --write-snplist --out refsnps || true
    plink --bfile ${cohort_bim.baseName} --write-snplist --out refsnps || true
    plink --bfile ${ref_bim.baseName} --bmerge ${cohort_bed} ${cohort_bim} ${cohort_fam} --write-snplist --out merged || true
    
    # Get common snps
    comm -12 <(sort refsnps.snplist) <(sort refsnps.snplist) | \
    comm -3 - <(sort merged.missnp) \
    > common_snps.txt

    # Merge
    plink --bfile ${ref_bim.baseName} \
        --extract common_snps.txt \
        --make-bed \
        --out ref

    plink --bfile ${cohort_bim.baseName} \
        --extract common_snps.txt \
        --make-bed \
        --out cohort
    
    plink --bfile ref \
        --bmerge cohort.bed cohort.bim cohort.fam \
        --extract common_snps.txt \
        --make-bed \
        --out ${ref}.${cohort}
    """
}
