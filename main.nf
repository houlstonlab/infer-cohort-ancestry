#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { REMOVE }      from './modules/remove.nf'
include { UPDATE }      from './modules/update.nf'
include { FIX }         from './modules/fix.nf'
include { FILL }        from './modules/fill.nf'
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
chroms_ch   = Channel.of (1..2) | map { "chr$it" }
modes_ch    = Channel.of(params.modes.split(','))

// worflow
workflow {
    // Get coordinates
    dbsnp
        | combine(chroms_ch)
        | COORDINATES
        | combine(variants_ch)
        | SUBSET
        | filter { it[3].toInteger() > 0 }
        | ( params.remove ? REMOVE : map {it} )
        | filter { it[3].toInteger() > 0 }
        | ( params.update ? combine(dbsnp) : map {it} )
        | ( params.update ? UPDATE : map {it} )
        | filter { it[3].toInteger() > 0 }
        | ( params.fix    ? combine(fasta) : map {it} )
        | ( params.fix    ? FIX    : map {it} )
        | filter { it[3].toInteger() > 0 }
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    // Fill cases
    snps.cases
        | ( params.fill ? FILL : map {it} )
        | set { cases }

    // Convert    
    snps.references
        | concat(cases)
        | filter { it[3].toInteger() > 0 }
        | CONVERT
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    // Prune
    snps.cases
        | ( params.prune ? combine(ld_regions) : map {it} )
        | ( params.prune ? PRUNE : map {it} )
        | set { cases }
    
    // Combine
    snps.references
        | concat(cases)
        | groupTuple(by: [0,1])
        | combine(population_ch, by: 0) 
        | COMBINE
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    // Merge, filter, scale, assign and plot
    snps.cases
        | combine(snps.references)
        | MERGE
        | FILTER
        | combine(modes_ch)
        | SCALE
        | ASSIGN
        | PLOT
}
