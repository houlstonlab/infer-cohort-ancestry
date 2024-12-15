process PCA {
    tag "${ref}:${cohort}:${mode}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/pca", mode: 'copy')

    input:
    tuple val(ref), val(cohort),
          path(bim), path(bed), path(fam), path(nosex), path(log), path(pop),
          val(mode)

    output:
    tuple val(ref), val(cohort), val(mode),
          path("${ref}.${cohort}.${mode}.variants.txt"),
          path("${ref}.${cohort}.${mode}.eigenvec"),
          path("${ref}.${cohort}.${mode}.eigenval"),
          path("${ref}.${cohort}.${mode}.clst"),
          path("${ref}.${cohort}.${mode}.log"),
          path("${ref}.${cohort}.${mode}.pop")

    script:
    if (mode == 'clusters') {
        """
        #!/bin/bash
        cat ${pop} > ${ref}.${cohort}.${mode}.pop
        cat ${pop} | awk '{ print \$1, \$2, \$3}' > populations.txt
        cat ${pop} | awk '{if (\$3 != 'NA'); print \$3}' | sort -u > clusters.txt
           
        # Select random variants
        RANDOM=42; shuf -n ${params.N_VARS} ${bim} | \
        cut -f 2 \
        > ${ref}.${cohort}.${mode}.variants.txt
        
        # Perform PCA
        plink --bfile ${bim.baseName} \
            --extract ${ref}.${cohort}.${mode}.variants.txt \
            --pca \
            --pca-clusters clusters.txt \
            --within populations.txt \
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
                
        cat ${pop} > ${ref}.${cohort}.${mode}.pop
        touch ${ref}.${cohort}.${mode}.clst
        """
    } else {
        throw new Exception("Unknown mode: ${mode}")
    }
}
