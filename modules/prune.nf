process PRUNE {
    tag "${cohort}:${type}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/pruned", mode: 'copy')

    input:
    tuple val(cohort), val(type),
          path(bim), path(bed), path(fam), path(nosex), path(log),
          path(ld_regions)
          
    output:
    tuple val(cohort), val(type),
          path("${cohort}.${type}.prune.bim"),
          path("${cohort}.${type}.prune.bed"),
          path("${cohort}.${type}.prune.fam"),
          path("${cohort}.${type}.prune.nosex"),
          path("${cohort}.${type}.prune.in"),
          path("${cohort}.${type}.prune.out"),
          path("${cohort}.${type}.prune.log")

    script:
    """
    #!/bin/bash        
    plink --bfile ${bim.baseName} \
        --exclude range ${ld_regions} \
        --indep-pairwise ${params.window} ${params.step} ${params.rsquared} \
        --rel-cutoff ${params.relatedness} \
        --const-fid 0 \
        --out plink_tmp
    
    plink --bfile ${bim.baseName} \
        --extract plink_tmp.prune.in \
        --make-bed \
        --const-fid 0 \
        --out ${cohort}.${type}.prune
    
    # Explicitly rename output files
    mv plink_tmp.prune.in ${cohort}.${type}.prune.in
    mv plink_tmp.prune.out ${cohort}.${type}.prune.out
    """
}
