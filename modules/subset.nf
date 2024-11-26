process SUBSET {
    tag "${cohort}:${type}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/cohorts", mode: 'copy')

    input:
    tuple val(cohort), path(vcf_in), path(index_in),
          val(type), val(size), path(cases), 
          val(fasta), path(fasta_in), path(fasta_index)

    output:
    tuple val(cohort), val(type),
          path("${cohort}.${type}.snps.vcf.gz"),
          path("${cohort}.${type}.snps.vcf.gz.tbi")
 
    script:
    if (type == 'references') {
        """
        #!/bin/bash
        bcftools view -v snps -S ${cases} ${vcf_in} | \
        bcftools +fill-tags -- -t all | \
        bcftools view --threads ${task.cpu} -Oz -o ${cohort}.${type}.snps.vcf.gz
        
        tabix ${cohort}.${type}.snps.vcf.gz
        """
    } else if (type == 'cases') {
        """
        #!/bin/bash
        bcftools view -v snps -S ${cases} ${vcf_in} | \
        bcftools filter -i 'FILTER="PASS"' | \
        bcftools annotate -x INFO/CSQ | \
        bcftools +setGT -- -t '.' -n '0/0' | \
        bcftools +fixploidy -- | \
        bcftools +fill-tags -- -t all | \
        bcftools view --threads ${task.cpu} -Oz -o ${cohort}.${type}.snps.vcf.gz
        
        tabix ${cohort}.${type}.snps.vcf.gz
        """
    } else {
        printl "Unknown cohort type: ${type}"
    }
}
