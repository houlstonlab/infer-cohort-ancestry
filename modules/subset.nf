process SUBSET {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/cohorts", mode: 'copy')

    input:
    // tuple val(chrom), path(coordinates),
    tuple val(cohort), val(type), val(size), 
		  path(vcf_in), path(index_in),
          path(population),
          val(chrom)

    output:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path("${cohort}.${type}.${chrom}.snps.vcf.gz"),
          path("${cohort}.${type}.${chrom}.snps.vcf.gz.tbi")
          
 
    script:
    """
    #!/bin/bash
	# Extract SNPs
    bcftools view \
        ${vcf_in} \
        -S <(awk '{ print \$2 }' ${population}) \
        -r ${chrom} \
        -g het \
        -v snps -m2 -M2 \
        --threads ${task.cpu} \
        -Oz -o ${cohort}.${type}.${chrom}.snps.vcf.gz
    tabix ${cohort}.${type}.${chrom}.snps.vcf.gz
    n_vars=\$(bcftools index -n ${cohort}.${type}.${chrom}.snps.vcf.gz)
    """
}

    // bcftools norm -d none | \
    // bcftools view -v snps -m2 -M2 | \
