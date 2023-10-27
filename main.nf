#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import subworkflows
include {bactocap} from './workflows/bactocap.nf'
include {bactocap1} from './workflows/bactocapSaveSpace.nf'
include {bactocap2} from './workflows/bactocapSaveSpace.nf'
include {bactocap3} from './workflows/bactocapSaveSpace.nf'
include {bactocap1Repeat} from './workflows/bactocapSaveSpace.nf'

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

      if (params.saveSpace) {

        bactocap1(input_files, bowtieDir)

        picardBam = bactocap1.out.picardBam
        picardMetrics =  bactocap1.out.picardMetrics

        bactocap2(refFastaDir, genBankDir, discRegCustDir, picardBam, picardMetrics)

        trigger = bactocap2.out.bactocap1Repeat_trigger

        bactocap1Repeat(input_files, bowtieDir, trigger)

        picardBam2 = bactocap1Repeat.out.picardBam
        panSNPs = bactocap2.out.panSNPs
        panSNPsIntervals = bactocap2.out.panSNPsIntervals
        
        bactocap3(refFastaDir, picardBam2, panSNPs, panSNPsIntervals)

      }

      else {
      
        bactocap(input_files, bowtieDir, refFastaDir, genBankDir, discRegCustDir)

      }
}
