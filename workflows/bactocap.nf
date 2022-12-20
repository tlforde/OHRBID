// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {createStatsFile} from '../modules/bactocapModules.nf' params(params)
include {trimmomatic} from '../modules/bactocapModules.nf' params(params)
include {bowtie2} from '../modules/bactocapModules.nf' params(params)
include {sam2bam} from '../modules/bactocapModules.nf' params(params)
include {indexbam} from '../modules/bactocapModules.nf' params(params)
include {picard} from '../modules/bactocapModules.nf' params(params)
include {bedtools} from '../modules/bactocapModules.nf' params(params)
include {computeStatistics} from '../modules/bactocapModules.nf' params(params)
include {mpileupVarscan} from '../modules/bactocapModules.nf' params(params)
include {discardRegions} from '../modules/bactocapModules.nf' params(params)
include {extractPanSNPs} from '../modules/bactocapModules.nf' params(params)
include {bamreadcount} from '../modules/bactocapModules.nf' params(params)
include {makeVCtable} from '../modules/bactocapModules.nf' params(params)
include {makeVCtableNoMonomorphic} from '../modules/bactocapModules.nf' params(params)
include {makeVCflagTable} from '../modules/bactocapModules.nf' params(params)
include {removeSelectedRegions} from '../modules/bactocapModules.nf' params(params)
include {generateAlignmentTable} from '../modules/bactocapModules.nf' params(params)
include {filterPanSNPs} from '../modules/bactocapModules.nf' params(params)
include {filterPanSNPsNoCustom} from '../modules/bactocapModules.nf' params(params)
include {extractMutatedGenes} from '../modules/bactocapModules.nf' params(params)

// define workflow component
workflow bactocap {

    take:

      input_files
      bowtieDir
      refFasta
      genBank
      discRegFileCustom

    main:

      trimmomatic(input_files)

      bowtie2(trimmomatic.out.trim_paired, bowtieDir.toList())

      sam2bam(bowtie2.out.bowtie2_out)

      indexbam(sam2bam.out.sam2bam_out)

      picard(sam2bam.out.sam2bam_out.join(indexbam.out.indexbam_out, by: 0))

      bedtools(picard.out.picard_out)

      computeStatistics(bedtools.out.bedtools_out.join(picard.out.picard_metrics, by: 0))

      createStatsFile(computeStatistics.out.summary_stats.collect())

      mpileupVarscan(picard.out.picard_out, refFasta.toList())

      discardRegions(genBank.toList())

      extractPanSNPs(mpileupVarscan.out.varscan_out.collect())

      bamreadcount(picard.out.picard_out, extractPanSNPs.out.extpansnps_intervals, refFasta.toList())

      if ( params.discRegFileCustom != null ) {

          filterPanSNPs(extractPanSNPs.out.extpansnps_out, discardRegions.out.discardregions_out, discRegFileCustom.toList())

      }

      else {

          filterPanSNPsNoCustom(extractPanSNPs.out.extpansnps_out, discardRegions.out.discardregions_out)

      }

      extractMutatedGenes(genBank.toList(), extractPanSNPs.out.extpansnps_out)

      makeVCtable(extractPanSNPs.out.extpansnps_out, bamreadcount.out.bamreadcount_out.collect())

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
