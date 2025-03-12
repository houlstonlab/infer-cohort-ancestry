process PRUNE {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'
    label 'plink'

    publishDir("${params.output_dir}/pruned", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom),
          path(bim), path(bed), path(fam), path(nosex), path(log),
          path(ld_regions)
          
    output:
    tuple val(cohort), val(type), val(chrom),
          path("${cohort}.${type}.${chrom}.prune.bim"),
          path("${cohort}.${type}.${chrom}.prune.bed"),
          path("${cohort}.${type}.${chrom}.prune.fam"),
          path("${cohort}.${type}.${chrom}.prune.nosex"),
        //   path("${cohort}.${type}.${chrom}.prune.in"),
        //   path("${cohort}.${type}.${chrom}.prune.out"),
          path("${cohort}.${type}.${chrom}.prune.log")

    script:
    """
    #!/bin/bash        
    plink --bfile ${bim.baseName} \
        --exclude range ${ld_regions} \
        --indep-pairwise ${params.window} ${params.step} ${params.rsquared} \
        --const-fid 0 \
        --out plink_tmp
    
    plink --bfile ${bim.baseName} \
        --extract plink_tmp.prune.in \
        --make-bed \
        --const-fid 0 \
        --out ${cohort}.${type}.${chrom}.prune
    
    # Explicitly rename output files
    mv plink_tmp.prune.in ${cohort}.${type}.${chrom}.prune.in
    mv plink_tmp.prune.out ${cohort}.${type}.${chrom}.prune.out
    """
}
