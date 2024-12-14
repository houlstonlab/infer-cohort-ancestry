#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { UPDATE }      from './modules/update.nf'
include { FIX }         from './modules/fix.nf'
include { CONVERT }     from './modules/convert.nf'
include { PRUNE }       from './modules/prune.nf'
include { COMBINE }     from './modules/combine.nf'
include { MERGE }       from './modules/merge.nf'
include { PCA }         from './modules/pca.nf'
include { PLOT }        from './modules/plot.nf'

// Define input channels
chroms_ch =  Channel.of (1..22)

variants_ch = Channel.fromFilePairs(params.vcf, flat: true)
cases_ch = Channel.fromPath(params.cases)
    | map { [it.simpleName, it] }

cohorts_info_ch = Channel.fromPath(params.cohorts_info)
    | splitCsv(header: true, sep: ',')
    | map { row -> [ row.cohort, row.type, row.size ] }

// Reference files
dbsnp               = Channel.fromFilePairs(params.dbsnp, flat: true)
fasta               = Channel.fromFilePairs(params.fasta, flat: true)
ld_regions          = Channel.fromPath(params.ld_regions)
populations         = Channel.fromPath(params.populations)
populations_id      = Channel.fromPath(params.populations_id)
populations_info    = Channel.fromPath(params.populations_info)
clusters            = Channel.fromPath(params.clusters)
modes_ch            = Channel.of(params.modes.split(','))

// worflow
workflow {
    // Prepare SNPs
    variants_ch
        | concat(dbsnp)
        | combine(chroms_ch)
        | COORDINATES
        | groupTuple(by: 0)
        | set { coordinates }

    // Subset, update, filter, convert, and prune
    variants_ch 
        | combine(cohorts_info_ch, by: 0) 
        | combine(cases_ch, by: 0) 
        | combine(coordinates)
        | SUBSET
        | combine(dbsnp)
        | UPDATE
        | combine(fasta)
        | FIX
        | CONVERT
        | combine(ld_regions)
        | PRUNE
        | groupTuple(by: [0,1])
        | COMBINE
        | branch {
            ref     : it[1] == 'references'
            cohort  : it[1] == 'cases'
        }
        | set { snps }

    // Merge cohorts, and run PCA
    snps.cohort
        | combine(snps.ref)
        | MERGE
        | combine(modes_ch)
        | combine(populations)
        | combine(clusters)
        | PCA
        | combine(populations_id)
        | combine(populations_info)
        | PLOT
}
