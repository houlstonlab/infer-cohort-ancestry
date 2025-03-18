#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { SUBSET }      from './modules/subset.nf'
include { REMOVE }      from './modules/remove.nf'
include { FIX }         from './modules/fix.nf'
include { FILL }        from './modules/fill.nf'
include { CONVERT }     from './modules/convert.nf'
include { PRUNE }       from './modules/prune.nf'
include { COMBINE }     from './modules/combine.nf'
include { MERGE }       from './modules/merge.nf'
include { FILTER }      from './modules/filter.nf'
include { SELECT }      from './modules/select.nf'
include { SCALE }       from './modules/scale.nf'
include { ASSIGN }      from './modules/assign.nf'
include { PLOT }        from './modules/plot.nf'

// Define input channels
variants_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    | map { row -> [ 
        row.cohort, row.type, row.size,
        file(row.vars_file), file(row.vars_index),
        file(row.population)
     ] }

population_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    | map { row -> [ row.cohort, file(row.population) ] }

dbsnp       = Channel.fromFilePairs(params.dbsnp, flat: true) | map { ['dbsnp', it[1], it[2]] }
fasta       = Channel.fromFilePairs(params.fasta, flat: true)
ld_regions  = Channel.fromPath(params.ld_regions)
chroms_ch   = Channel.of (1..22) | map { "chr$it" }
modes_ch    = Channel.of(params.modes.split(','))

// worflow
workflow {
    // Select dbsnp variants and subset Cohorts
    // Optional: remove, and fix
    dbsnp
        | combine(chroms_ch)
        | SELECT
        | map { it.last() }
        | splitText(by: params.chunk, file: true) \
        | map { it -> 
            // Split the file name and extract the first par
            def chrom = it.baseName.split("\\.")[1]
            def chunk = it.baseName.split("\\.")[2]
            return [chrom, chunk, it]
        }
        | combine(variants_ch)
        | SUBSET
        | filter { it.last().toInteger() > 0 }
        | ( params.remove ? REMOVE : map {it} )
        | filter { it.last().toInteger() > 0 }
        | ( params.fix    ? combine(fasta) : map {it} )
        | ( params.fix    ? FIX    : map {it} )
        | filter { it.last().toInteger() > 0 }
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    // Fill the study SNPs
    // Along with the reference, convert and combine
    snps.cases
        | ( params.fill && map{ it[1] } == 'cases' ? FILL : map {it} )
        | filter { it.last().toInteger() > 0 }
        | concat(snps.references)
        | CONVERT
        | groupTuple(by: [0,1])
        | combine(population_ch, by: 0) 
        | COMBINE
        | branch {
            references : it[1] == 'references'
            cases      : it[1] == 'cases'
        }
        | set { snps }

    // Prune the study SNPs
    // Merge with the rerences, filter, scale, assign and plot
    snps.cases
        | ( params.prune ? combine(ld_regions) : map {it} )
        | ( params.prune ? PRUNE : map {it} )
        | combine(snps.references)
        | MERGE
        | FILTER
        | combine(modes_ch)
        | SCALE
        | ASSIGN
        | PLOT
}
