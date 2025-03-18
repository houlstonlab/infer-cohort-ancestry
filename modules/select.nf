process SELECT {
    tag "${key}:${chrom}"

    label 'simple'
    label 'bcftools'

    publishDir("${params.output_dir}/selected", mode: 'copy')

    input:
    tuple val(key), path(file), path(index), 
          val(chrom)

    output:
    tuple val(key), val(chrom),
          path("${key}.${chrom}.vcf.gz"),
          path("${key}.${chrom}.vcf.gz.tbi"),
          path("${key}.${chrom}.txt")

    script:
    """
    #!/bin/bash
    # Select variants
    bcftools view \
    -r ${chrom} \
    -i 'VC=="SNV" && SAO!=2 && COMMON==1 && G5==1 && KGPhase3==1' \
    ${file} \
    --threads ${task.cpus} \
    -Oz -o ${key}.${chrom}.vcf.gz

    # Index
    tabix ${key}.${chrom}.vcf.gz

    # Create BED file
    bcftools query -f '%CHROM\t%POS\t%ID\n' ${key}.${chrom}.vcf.gz > ${key}.${chrom}.txt
    """
}
