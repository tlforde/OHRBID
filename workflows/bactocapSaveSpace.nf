// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {createStatsFile} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {trimmomatic} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {bowtie2} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {sam2bam} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {indexbam} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {picard} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {bedtools} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {computeStatistics} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {mpileupVarscan} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {discardRegions} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {extractPanSNPs} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {bamreadcount} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {makeVCtable} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {makeVCtableNoMonomorphic} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {makeVCflagTable} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {removeSelectedRegions} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {generateAlignmentTable} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {filterPanSNPs} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {filterPanSNPsNoCustom} from '../modules/bactocapSaveSpaceModules.nf' params(params)
include {extractMutatedGenes} from '../modules/bactocapSaveSpaceModules.nf' params(params)

// define workflow component
workflow bactocap1 {

    take:

      input_files
      bowtieDir

    main:

      trimmomatic(input_files)

      bowtie2(trimmomatic.out.trim_paired, bowtieDir.toList())

      sam2bam(bowtie2.out.bowtie2_out)

      indexbam(sam2bam.out.sam2bam_out)

      picard(sam2bam.out.sam2bam_out.join(indexbam.out.indexbam_out, by: 0))

    emit:

      picardBam = picard.out.picard_out
      picardMetrics = picard.out.picard_metrics

}

workflow bactocap2 {

    take:

      refFasta
      genBank
      discRegFileCustom
      picardBam
      picardMetrics

    main:

      discardRegions(genBank.toList())

      bedtools(picardBam)

      computeStatistics(bedtools.out.bedtools_out.join(picardMetrics, by: 0))

      createStatsFile(computeStatistics.out.summary_stats.collect())

      mpileupVarscan(picardBam, refFasta.toList())

      extractPanSNPs(mpileupVarscan.out.varscan_out.collect())

      extractMutatedGenes(genBank.toList(), extractPanSNPs.out.extpansnps_out)

      if ( params.discRegFileCustom != null ) {

          filterPanSNPs(extractPanSNPs.out.extpansnps_out, discardRegions.out.discardregions_out, discRegFileCustom.toList())

      }

      else {

          filterPanSNPsNoCustom(extractPanSNPs.out.extpansnps_out, discardRegions.out.discardregions_out)

      }

    emit:

      panSNPs = extractPanSNPs.out.extpansnps_out
      panSNPsIntervals = extractPanSNPs.out.extpansnps_intervals

}

workflow bactocap3 {

    take:

      refFasta
      picardBam
      panSNPs
      panSNPsIntervals

    main:

      bamreadcount(picardBam, panSNPsIntervals, refFasta.toList())

      makeVCtable(panSNPs, bamreadcount.out.bamreadcount_out.collect())

      if ( params.noMonomorphic == "yes" ) {

          makeVCtableNoMonomorphic(makeVCtable.out.makeVCtable_out)

          makeVCflagTable(makeVCtableNoMonomorphic.out.makeVCtableNoMono_out)

          removeSelectedRegions(makeVCtableNoMonomorphic.out.makeVCtableNoMono_out, makeVCflagTable.out.makeVCflagTable_probinter)

          generateAlignmentTable(removeSelectedRegions.out.rmSelectRegions_out)

      }

      else {

          makeVCflagTable(makeVCtable.out.makeVCtable_out)

          removeSelectedRegions(makeVCtable.out.makeVCtable_out, makeVCflagTable.out.makeVCflagTable_probinter)

          generateAlignmentTable(removeSelectedRegions.out.rmSelectRegions_out)

      }
}
