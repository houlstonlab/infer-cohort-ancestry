process PLOT {
    tag "${ref}:${cohort}:${mode}"

    label 'simple'

    container = params.rocker

    publishDir("${params.output_dir}/plots", mode: 'copy')

    input:
    tuple val(ref), val(cohort), val(mode),
          path(file), path(log), path(pop)
    
    output:
    tuple val(ref), val(cohort), val(mode),
          path("${ref}.${cohort}.${mode}.png")
    
    script:
    """
    #!/bin/bash
    plot_pca.R ${ref} ${cohort} ${mode} ${file} ${pop}
    """
}