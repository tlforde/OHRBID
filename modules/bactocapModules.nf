
process trimmomatic {
    /**
    * Trim adapters from reads (https://github.com/usadellab/Trimmomatic)
    * @input tuple sample_id, path(forward), path(reverse)
    * @output trim_paired tuple sample_id, path("*_fp.fq"), path("*_rp.fq")
    * @output trim_unpaired tuple sample_id, path("*_fu.fq.gz"), path("*_ru.fq.gz")
    */

    tag { sample_id }
    label 'low_cpu'
    label 'low_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*_{fp,rp}.fq", mode: 'copy'

    input:
    tuple val(sample_id), path(forward), path(reverse)

    output:
    tuple val(sample_id), path("*_fp.fq"), path("*_rp.fq"), emit: trim_paired optional true
    tuple val(sample_id), path("*_fu.fq.gz"), path("*_ru.fq.gz"), emit: trim_unpaired optional true
    path "${sample_id}.err", emit: trim_error optional true

    script:
    error_log = "${sample_id}.err"

    """
    if [[ \$(zcat ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      echo "error: Sample ${sample_id} is empty" > ${error_log}
      exit 0
    else
      trimmomatic PE ${forward} ${reverse} -threads ${task.cpus} -phred33 ${sample_id}_fp.fq ${sample_id}_fu.fq.gz ${sample_id}_rp.fq ${sample_id}_ru.fq.gz SLIDINGWINDOW:4:${params.qThreshold} MINLEN:${params.minLen}
    fi
    """
}

process bowtie2 {
    /**
    * Align sequence reads with bowtie2 (https://github.com/BenLangmead/bowtie2)
    * @input tuple sample_id, path(forward), path(reverse)
    * @input index
    * @output bowtie2_out tuple sample_id, path("*.sam")
    */

    tag { sample_id }
    label 'medium_cpu'
    label 'medium_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*.sam", mode: 'copy'

    input:
    tuple val(sample_id), path(forward), path(reverse)
    path(index)

    output:
    tuple val(sample_id), path("*.sam"), emit: bowtie2_out

    script:
    sam = "${sample_id}.sam"

    """
    bowtie2 -p ${task.cpus} -x ./${params.bowtieIndexName} -1 $forward -2 $reverse --maxins ${params.maxins} > ${sam}
    """
}

process sam2bam {
    /**
    * Convert sam to bam with samtools (https://github.com/samtools/samtools)
    * @input tuple sample_id, path(sam)
    * @output sam2bam_out tuple sample_id, path("*.bam")
    */

    tag { sample_id }

    label 'medium_cpu'
    label 'low_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("*.bam"), emit: sam2bam_out

    script:
    bam = "${sample_id}.bam"

    """
    samtools sort -@ ${task.cpus} -o ${bam} ${sam}
    """
}

process indexbam {
    /**
    * Index bam with samtools (https://github.com/samtools/samtools)
    * @input tuple sample_id, path(bam)
    * @output sam2bam_out tuple sample_id, path("*.bam")
    */

    tag { sample_id }

    label 'medium_cpu'
    label 'low_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*bai", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("*.bai"), emit: indexbam_out

    script:
    """
    samtools index ${bam} -@ ${task.cpus}
    """
}


process picard {
    /**
    * Identify duplicate reads with picard (https://github.com/broadinstitute/picard)
    * @input tuple sample_id, path(bam)
    * @output picard_out tuple sample_id, path("*.bam")
    */

    tag { sample_id }

    label 'medium_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*{.bam,.bai,_metrics}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*_noDuplicateReads.bam"), path("*_noDuplicateReads.bai"), emit: picard_out
    tuple val(sample_id), path("*_metrics"), emit: picard_metrics

    script:
    bam_nodup = "${sample_id}_noDuplicateReads.bam"
    metrics = "${sample_id}_noDuplicateReads.dup_metrics"

    """
    picard MarkDuplicates INPUT=${bam} OUTPUT=${bam_nodup} METRICS_FILE=${metrics} CREATE_INDEX=true REMOVE_DUPLICATES=true
    """
}

process bedtools {
    /**
    * Compute histogram of feature coverage with bedtools (https://github.com/arq5x/bedtools2)
    * @input tuple sample_id, path(bam)
    * @output bedtools_out tuple sample_id, path("*.genomecov.hist")
    */

    tag { sample_id }

    label 'low_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*.genomecov.hist", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.genomecov.hist"), emit: bedtools_out

    script:
    hist = "${sample_id}.genomecov.hist"

    """
    bedtools genomecov -ibam ${bam} > ${hist}
    """
}

process computeStatistics {
    /**
    * Run python script to compute statistics and draw histogram
    * @input tuple sample_id, path(hist)
    * @output stats_out tuple sample_id, path("*.genomestats.txt")
    * @output hist_out tuple sample_id, path("*_coverage_histogram.png")
    */

    tag { sample_id }

    label 'low_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*{.txt,.png}", mode: 'copy'

    input:
    tuple val(sample_id), path(histogram), path(picardMetrics)

    output:
    tuple val(sample_id), path("*.genomestats.txt"), emit: stats_out
    tuple val(sample_id), path("*_coverage_histogram.png"), emit: stats_png
    path("*.summarystats.txt"), emit: summary_stats 

    script:
    genomeStats = "${sample_id}.genomestats.txt"
    hist = "${sample_id}_coverage_histogram.png"

    """
    python3 ${baseDir}/bin/compute_statistics.py -n ${sample_id} -c ${params.covThreshold} -b ${histogram} -p ${picardMetrics}
    """
}

process createStatsFile {
    /**
    * Run python script to create summary statistics file
    * @input none
    * @output summary_out path("summary_statistics.txt")
    */

    label 'low_memory'

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*_statistics.txt", mode: 'copy'

    input:
    path(statsFiles)

    output:
    path("summary_statistics.txt", emit: summary_out)

    script:

    """
    python3 ${baseDir}/bin/create_summary_stats.py -c ${params.covThreshold}


    cat *.summarystats.txt >> summary_statistics.txt
    """
}

process mpileupVarscan {
    /**
    * Call vars with mpileup and pass to varscan (https://github.com/dkoboldt/varscan)
    * @input tuple sample_id, path(bam), path(refSeq)
    * @output varscan_out tuple sample_id, path("*.varscan")
    */

    tag { sample_id }

    label 'medium_memory'

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*.varscan", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(refFasta)

    output:
    path("*.varscan", emit: varscan_out)

    script:
    varscan = "${sample_id}.varscan"

    """
    samtools mpileup -f ${params.refGenomeName} ${bam} | varscan pileup2snp --min-coverage ${params.covThresVarscan} --min-reads2 ${params.minReads} --min-var-freq ${params.minVarFreq} --min-avg-qual ${params.minAvgQual} > ${varscan}
    """
}

process discardRegions {
    /**
    * Run python script to get intervals to be discarded
    * @input path(genbank)
    * @output discardregions_out tuple path("discarded_regions.txt")
    */

    label 'low_memory'

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.txt", mode: 'copy'

    input:
    path(genbank)

    output:
    path("discarded_regions.txt", emit: discardregions_out)

    script:

    """
    python3 ${baseDir}/bin/discard_regions.py -g ${genbank}
    """
}

process extractPanSNPs {
    /**
    * Run python script to extract pan snps
    * @input path(varscan)
    * @output extractpansnps_out tuple path("panSNPs.txt"), path("panSNPs_intervals.txt")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.txt", mode: 'copy'

    input:
    path(varscan)

    output:
    path("panSNPs.txt", emit: extpansnps_out)
    path("panSNPs_intervals.txt", emit: extpansnps_intervals)

    script:

    """
    python3 ${baseDir}/bin/extract_panSNPs.py
    awk '{printf "%s\\t%s\\t%s\\n", \$1, \$2, \$2}' panSNPs.txt > panSNPs_intervals.txt
    """
}

process bamreadcount {
    /**
    * Run bam-readcount
    * @input tuple sample_id, path(bam), path(bai)
    * @input path(snp)
    * @input path(refseq)
    * @output bamreadcount_out tuple sample_id, path("*.bamreadcount")
    */

    tag { sample_id }

    publishDir "${params.outputDir}/${sample_id}/${task.process.replaceAll(":", "_")}", pattern: "*.bamreadcount", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(snpintervals)
    path(refseq)

    output:
    path("*.bamreadcount", emit: bamreadcount_out)

    script:
    readcount = "${sample_id}.bamreadcount"

    """
    bam-readcount -w 1 -b ${params.minBaseQual} -l ${snpintervals} -f ${params.refGenomeName} ${bam} > ${readcount}
    """
}

process makeVCtable {
    /**
    * Run python script to get generate VCtable
    * @input path(snps)
    * @output makeVCtable_out path("VC_table.dat")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.dat", mode: 'copy'

    input:
    path(snps)
    path(bamreadcount)

    output:
    path("VC_table.dat", emit: makeVCtable_out)

    script:

    """
    python3 ${baseDir}/bin/makeVCtable.py ${snps} .
    """
}

process makeVCtableNoMonomorphic {
    /**
    * Run python script to generate variant call table (tab-delimited) with monomorphic sites excluded
    * @input path(vctable)
    * @output makeVCtable_out path("VC_table.dat")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*{.dat,.txt}", mode: 'copy'

    input:
    path(vctable)

    output:
    path("monomorphic_sites.txt", emit: makeVCtableNoMono_sites) 
    path("VC_table_noMonomorphic.dat", emit: makeVCtableNoMono_out)

    script:

    """
    python3 ${baseDir}/bin/makeVCtableNoMonomorphic.py ${vctable}
    """
}

process makeVCflagTable {
    /**
    * Run python script to generate flag table (tab-delimited)
    * @input path(vctable)
    * @output makeVCflagTable_out path("VC_flag_table.dat")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.dat", mode: 'copy'

    input:
    path(vctable)

    output:
    path("VC_flag_table.dat", emit: makeVCflagTable_out)
    path("VC_flag_table_problematic_sites_intervals.txt", emit: makeVCflagTable_probinter)
    path("VC_flag_table_refined.dat", emit: makeVCflagTable_refine)

    script:

    """
    python3 ${baseDir}/bin/makeVCflagTable.py ${vctable} ${params.covThresVC} ${params.baseFreqThres}

    awk -F 'below' '{if (NF-1 < 2) {print \$0}}' VC_flag_table.dat > VC_flag_table_refined.dat
    awk -F 'below' '{if (NF-1 >= 2) {print \$0}}' VC_flag_table.dat > VC_flag_table_problematic_sites.dat

    awk '{printf "%s\\t%s\\n",\$2,\$2}' VC_flag_table_problematic_sites.dat > VC_flag_table_problematic_sites_intervals.txt
    """
}

process removeSelectedRegions {
    /**
    * Run python script to remove selected regions
    * @input path(vctable)
    * @input path(probintertxt)
    * @output makeVCtable_out path("VC_table_filtered.dat")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.dat", mode: 'copy'

    input:
    path(vctable)
    path(probintertxt)

    output:
    path("VC_table_filtered.dat", emit: rmSelectRegions_out)

    script:

    """
    python3 ${baseDir}/bin/removeSelectedRegions.py ${vctable} ${probintertxt}
    """
}


process generateAlignmentTable {
    /**
    * Run python script to generate alignment table
    * @input path(vctablefilter)
    * @output makeVCtable_out path("alignment_table.dat")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.dat", mode: 'copy'

    input:
    path(vctablefilter)

    output:
    path("alignment_table.dat", emit: genaligntable_out)

    script:

    """
    python3 ${baseDir}/bin/generateAlignmentTable.py ${vctablefilter} ${params.covThresVC} ${params.baseFreqThres}
    """
}

process filterPanSNPs {
    /**
    * Run python script to filter pan snps
    * @input path(snpintervals)
    * @input path(discardedregions)
    * @input path(probregions)
    * @output filterpansnps_out path("panSNPs_intervals_refined.txt")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.txt", mode: 'copy'

    input:
    path(snpintervals)
    path(discardedregions)
    path(probregions)

    output:
    path("panSNPs_intervals_refined.txt", emit: filterpansnps_intervals)

    script:

    """
    python3 ${baseDir}/bin/filter_panSNPs.py -i ${snpintervals} -d ${discardedregions} -c ${params.discRegFileCustom}
    """
}

process filterPanSNPsNoCustom {
    /**
    * Run python script to filter pan snps
    * @input path(snpintervals)
    * @input path(discardedregions)
    * @input path(probregions)
    * @output filterpansnps_out path("panSNPs_intervals_refined.txt")
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.txt", mode: 'copy'

    input:
    path(snpintervals)
    path(discardedregions)

    output:
    path("panSNPs_intervals_refined.txt", emit: filterpansnps_intervals)

    script:

    """
    python3 ${baseDir}/bin/filter_panSNPs.py -i ${snpintervals} -d ${discardedregions}
    """
}

process extractMutatedGenes {
    /**
    * Run python script to extract mutated genes
    * @input path(snpfile)
    * @input path(genbank)
    * @output NEEDS TO BE DEFINED
    */

    publishDir "${params.outputDir}/${task.process.replaceAll(":", "_")}", pattern: "*.txt", mode: 'copy'

    input:
    path(snpfile)
    path(genbank)

    script:

    """
    python3 ${baseDir}/bin/extract_mutated_genes.py -s ${snpfile} -g ${genbank}
    """
}
