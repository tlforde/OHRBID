#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import subworkflows
include {bactocap} from './workflows/bactocap.nf'

// main workflow
workflow {

    // add a trailing slash if it was not originally provided to --inputDir
    inputDir_amended = "${params.inputDir}".replaceFirst(/$/, "/") 

    inDir = inputDir_amended

    pattern = params.pattern
    reads = inDir + "**/" + pattern

    Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
	   .ifEmpty { error "cannot find any reads matching ${pattern} in ${inDir}" }
	   .set{ input_files }

    // create channel for the bowtie index and the reference
    bowtieDir = Channel.fromPath( "${params.refDir}/*.bt2" )

    refFastaDir = Channel.fromPath( "${params.refDir}/${params.refGenomeName}" )

    genBankDir = Channel.fromPath( "${params.refDir}/${params.genBank}" )

    discRegCustDir = Channel.fromPath( "${params.refDir}/${params.discRegFileCustom}" )

    // main workflow
    main:

      // BACTOCAP SUB-WORKFLOW
      bactocap(input_files, bowtieDir, refFastaDir, genBankDir, discRegCustDir)
}
