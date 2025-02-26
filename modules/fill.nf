process FILL {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/filled", mode: 'copy')

    input:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path(vcf_in), path(index_in)

    output:
    tuple val(cohort), val(type), val(chrom), env(n_vars),
          path("${cohort}.${type}.${chrom}.filled.vcf.gz"),
          path("${cohort}.${type}.${chrom}.filled.vcf.gz.tbi")
     
    script:
    """
    #!/bin/bash
    # Fill VCF file
	bcftools view  -i 'FILTER="PASS"' ${vcf_in} | \
    bcftools norm -m -any | \
    bcftools +setGT -- -t . -n 0 | \
    bcftools +fill-tags -- -t all | \
    bcftools view -e 'HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
    bcftools +setGT -- -t q -n 0 -i 'FMT/GQ < ${params.GQ} | FMT/DP < ${params.DP} | VAF < ${params.VAF}' | \
    bcftools +fill-tags -- -t all | \
    bcftools view -e 'MAF < ${params.MAF}' | \
 	bcftools view -g het --threads ${task.cpu} -Oz -o ${cohort}.${type}.${chrom}.filled.vcf.gz
    tabix ${cohort}.${type}.${chrom}.filled.vcf.gz
    n_vars=\$(bcftools index -n ${cohort}.${type}.${chrom}.filled.vcf.gz)
	"""
}
