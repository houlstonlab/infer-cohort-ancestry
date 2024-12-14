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
          path("${ref}.${cohort}.${mode}.eigenvec"),
          path("${ref}.${cohort}.${mode}.eigenval")

    script:
    if (mode == 'clusters') {
        """
        #!/bin/bash
        plink --bfile ${bim.baseName} \
            --maf ${params.AF} \
            --hwe ${params.HWE} \
            --mind ${params.F_MISSING} \
            --pca \
            --within ${clusters} \
            --pca-clusters ${populations} \
            --write-cluster \
            --out ${ref}.${cohort}.${mode}
        """
    } else if (mode == 'noclusters') {
        """
        #!/bin/bash
        plink --bfile ${bim.baseName} \
            --maf ${params.AF} \
            --hwe ${params.HWE} \
            --mind ${params.F_MISSING} \
            --pca \
            --out ${ref}.${cohort}.${mode}
        """
    } else {
        throw new Exception("Unknown mode: ${mode}")
    }
}
