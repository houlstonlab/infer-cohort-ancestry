process COORDINATES {
    tag "${key}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/coordinates", mode: 'copy')

    input:
    tuple val(key), path(vcf), path(index), val(chrom)

    output:
    tuple val(chrom), path("${key}.${chrom}.snps.txt")
 
    script:
    """
    #!/bin/bash
    bcftools view -v snps -m2 -M2 -r ${chrom} ${params.common ? "-i 'COMMON=1'" : ""} ${vcf} | \
    bcftools query -f '%CHROM\t%POS\n' | \
    uniq > ${key}.${chrom}.snps.txt
    """
}
