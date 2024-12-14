process PCA {
    tag "${ref}:${cohort}:${mode}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/pca", mode: 'copy')

    input:
    tuple val(ref), val(cohort),
          path(bim), path(bed), path(fam), path(nosex), path(log), 
          val(mode), path(populations), path(clusters)

    output:
    tuple val(ref), val(cohort), val(mode),
          path("${ref}.${cohort}.${mode}.variants.txt"),
          path("${ref}.${cohort}.${mode}.eigenvec"),
          path("${ref}.${cohort}.${mode}.eigenval"),
          path("${ref}.${cohort}.${mode}.clst"),
          path("${ref}.${cohort}.${mode}.log")

    script:
    if (mode == 'clusters') {
        """
        #!/bin/bash
        # Select random variants
        RANDOM=42; shuf -n ${params.N_VARS} ${bim} | \
        cut -f 2 \
        > ${ref}.${cohort}.${mode}.variants.txt
        
        # Perform PCA
        plink --bfile ${bim.baseName} \
            --extract ${ref}.${cohort}.${mode}.variants.txt \
            --pca \
            --pca-clusters ${clusters} \
            --within ${populations} \
            --write-cluster \
            --out ${ref}.${cohort}.${mode}
        """
    } else if (mode == 'noclusters') {
        """
        #!/bin/bash
        # Select random variants
        RANDOM=42; shuf -n ${params.N_VARS} ${bim} | \
        cut -f 2 \
        > ${ref}.${cohort}.${mode}.variants.txt
        
        # Perform PCA
        plink --bfile ${bim.baseName} \
            --extract ${ref}.${cohort}.${mode}.variants.txt \
            --pca \
            --out ${ref}.${cohort}.${mode}
        touch ${ref}.${cohort}.${mode}.clst
        """
    } else {
        throw new Exception("Unknown mode: ${mode}")
    }
}
