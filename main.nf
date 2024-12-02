#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { SUBSET }  from './modules/subset.nf'
include { UPDATE }  from './modules/update.nf'
include { FILTER }  from './modules/filter.nf'
include { CONVERT } from './modules/convert.nf'
include { RELATE } from './modules/relate.nf'
include { PRUNE }   from './modules/prune.nf'
include { MERGE }   from './modules/merge.nf'
include { PCA }     from './modules/pca.nf'
include { PLOT }    from './modules/plot.nf'

// Define input channels
variants_ch = Channel.fromFilePairs(params.vcf, flat: true)
cases_ch = Channel.fromPath(params.cases)
    | map { [it.simpleName, it] }
cohorts_info_ch = Channel.fromPath(params.cohorts_info)
        | splitCsv(header: true, sep: ',')
        | map { row -> [ row.cohort, row.type, row.size ] }

// Reference files
dbsnp       = Channel.fromFilePairs(params.dbsnp, flat: true)
fasta       = Channel.fromFilePairs(params.fasta, flat: true)
ld_regions  = Channel.fromPath(params.ld_regions)

populations         = Channel.fromPath(params.populations)
populations_id      = Channel.fromPath(params.populations_id)
populations_info    = Channel.fromPath(params.populations_info)

clusters    = Channel.fromPath(params.clusters)
modes       = Channel.of('clusters', 'noclusters')

// worflow
workflow {
    // Prepare SNPs
    // Subset, update, filter, convert, and prune
    variants_ch 
        | combine(cohorts_info_ch, by: 0) 
        | combine(cases_ch, by: 0) 
        | combine(fasta)
        | SUBSET
        | combine(dbsnp)
        | UPDATE
        | FILTER
        | CONVERT
        | combine(ld_regions)
        | PRUNE
        | branch {
            ref     : it[1] == 'references'
            cohort  : it[1] != 'references'
        }
        | set { snps }

    // Merge cohorts, and run PCA
    MERGE( snps.cohort, snps.ref )
        | combine(modes)
        | combine(populations)
        | combine(clusters)
        | PCA
        | combine(populations_id)
        | combine(populations_info)
        | PLOT

    // Calculate relationship matrix
    CONVERT.out | RELATE
}
