#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { REMOVE }      from './modules/remove.nf'
include { UPDATE }      from './modules/update.nf'
include { FIX }         from './modules/fix.nf'
include { CONVERT }     from './modules/convert.nf'
include { PRUNE }       from './modules/prune.nf'
include { COMBINE }     from './modules/combine.nf'
include { MERGE }       from './modules/merge.nf'
include { FILTER }      from './modules/filter.nf'
include { SCALE }       from './modules/scale.nf'
include { ASSIGN }      from './modules/assign.nf'
include { PLOT }        from './modules/plot.nf'

// Define input channels
// cohorts_info_ch = Channel.fromPath(params.cohorts)
//     | splitCsv(header: true, sep: ',')
//     | map { row -> [ row.cohort, row.type, row.size ] }

variants_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    // | map { row -> [ row.cohort, file(row.vars_file), file(row.vars_index) ] }
    | map { row -> [ 
        row.cohort, row.type, row.size,
        file(row.vars_file), file(row.vars_index),
        file(row.population)
     ] }

population_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    | map { row -> [ row.cohort, file(row.population) ] }

dbsnp       = Channel.fromFilePairs(params.dbsnp, flat: true)
fasta       = Channel.fromFilePairs(params.fasta, flat: true)
ld_regions  = Channel.fromPath(params.ld_regions)
chroms_ch   = Channel.of (1..22) | map { "chr$it" }
modes_ch    = Channel.of(params.modes.split(','))

// worflow
workflow {
    // Prepare SNPs
    // variants_ch
    //     | concat(dbsnp)
    dbsnp
        | combine(chroms_ch)
        | COORDINATES
        | set { coordinates }

    // // Subset, update, filter, convert, and prune
    variants_ch 
    //     | combine(cohorts_info_ch, by: 0) 
    //     | combine(population_ch, by: 0) 
        | combine(coordinates)
        | SUBSET
        | filter { it[3].toInteger() > 0 }
        | REMOVE
        | filter { it[3].toInteger() > 0 }
        | combine(dbsnp)
        | UPDATE
        | filter { it[3].toInteger() > 0 }
        | combine(fasta)
        | FIX
        | filter { it[3].toInteger() > 0 }
        | CONVERT
        | combine(ld_regions)
        | PRUNE
        | groupTuple(by: [0,1])
        | combine(population_ch, by: 0) 
        | COMBINE
        | branch {
            ref     : it[1] == 'references'
            cohort  : it[1] == 'cases'
        }
        | set { snps }

    // Merge cohorts, filter, scale and plot
    snps.cohort
        | combine(snps.ref)
        | MERGE
        | FILTER
        | combine(modes_ch)
        | SCALE
        | ASSIGN
        | PLOT
}
