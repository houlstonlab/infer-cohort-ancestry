process SUBSET {
    tag "${cohort}:${type}:${chrom}"

    label 'simple'

    container = params.bcftools

    publishDir("${params.output_dir}/cohorts", mode: 'copy')

    input:
    tuple val(cohort), path(vcf_in), path(index_in),
          val(type), val(size), path(cases), 
          val(fasta), path(fasta_in), path(fasta_index),
          val(chrom), path(coordinates)

    output:
    tuple val(cohort), val(type), val(chrom),
          path("${cohort}.${type}.${chrom}.snps.vcf.gz"),
          path("${cohort}.${type}.${chrom}.snps.vcf.gz.tbi")
 
    script:
    if (type == 'references') {
        """
        #!/bin/bash
        echo -e "${chrom}\tchr${chrom}" > rename_chrs.txt
        cat ${coordinates} | sort | uniq -d | head -10000 > regions.txt

        bcftools annotate -R regions.txt -S ${cases} --rename-chrs rename_chrs.txt ${vcf_in} | \
        bcftools norm -d none | \
        bcftools view -v snps -m2 -M2 | \
        bcftools view --threads ${task.cpu} -Oz -o ${cohort}.${type}.${chrom}.snps.vcf.gz
        
        tabix ${cohort}.${type}.${chrom}.snps.vcf.gz
        """
    } else if (type == 'cases') {
        """
        #!/bin/bash
        echo -e "${chrom}\tchr${chrom}" > rename_chrs.txt
        cat ${coordinates} | sort | uniq -d | head -10000 > regions.txt
        
        bcftools annotate -R regions.txt -S ${cases} --rename-chrs rename_chrs.txt ${vcf_in} | \
        bcftools view -v snps | \
        bcftools filter -i 'FILTER="PASS"' | \
        bcftools annotate -x INFO/CSQ | \
        bcftools +setGT -- -t '.' -n '0/0' | \
        bcftools +fixploidy -- | \
        bcftools +fill-tags -- -t all | \
        bcftools view --threads ${task.cpu} -Oz -o ${cohort}.${type}.${chrom}.snps.vcf.gz
        
        tabix ${cohort}.${type}.${chrom}.snps.vcf.gz
        """
    } else {
        printl "Unknown cohort type: ${type}"
    }
}
