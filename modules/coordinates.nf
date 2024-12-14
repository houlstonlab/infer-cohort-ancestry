process COORDINATES {
    tag "${key}.${chrom}"

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
    bcftools query -r chr${chrom} -f '%CHROM\t%POS\n' ${vcf} | uniq > ${key}.${chrom}.snps.txt
    """
}
