process FIX {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'
    label 'bcftools'

    publishDir("${params.output_dir}/fixed", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path(vcf_in), path(index_in), 
          val(fasta), path(fasta_in), path(fasta_index)

    output:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path("${cohort}.${type}.${chrom}.fixed.vcf.gz"),
          path("${cohort}.${type}.${chrom}.fixed.vcf.gz.tbi")
     
    script:
    """
    #!/bin/bash
    # Fix VCF file
    bcftools +fixref \
        ${vcf_in} \
        -Oz -o ${cohort}.${type}.${chrom}.fixed.vcf.gz \
        -- -d \
        -f ${fasta_in} \
        -m flip

    tabix ${cohort}.${type}.${chrom}.fixed.vcf.gz
    n_vars=\$(bcftools index -n ${cohort}.${type}.${chrom}.fixed.vcf.gz)
    """
}
