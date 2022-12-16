#!/usr/bin/env python3

import argparse
import re
import traceback
import matplotlib
matplotlib.use('agg')# to avoid using problematic tkinter
import matplotlib.pyplot as plt
import numpy

def compute_statistics(sampleName, coverageThreshold, coverageHistogram, picardMetrics):
    try:
        binPattern = re.compile(r'genome\t(\d+)\t(\d+)\t.+')
        coverageValues = []
        counts = []
        coveragesCumulative = []
        countsOverCovThreshold = 0
        totalCoverage = 0
        inHistFile = open(coverageHistogram, 'r')
        for line in inHistFile:
            data = binPattern.findall(line)
            if data:
                coverageValues.append(int(data[0][0]))
                counts.append(int(data[0][1]))
                if int(data[0][0]) >= coverageThreshold:
                    countsOverCovThreshold += int(data[0][1])
                if int(data[0][0]) > 0:
                    coveragesCumulative.extend([int(data[0][0]) for i in range(int(data[0][1]))])
        metricsPattern = re.compile(r'.+?\t(\d+)\t(\d+)\t\d+\t(\d+)\t\d+\t\d+\t\d+\t(.+?)\t\d+\n')
        unpairedReadCount = 0
        pairedReadCount = 0
        unmappedReadCount = 0
        duplicatesFraction = 0.0
        inMetricsFile = open(picardMetrics, 'r')
        for line in inMetricsFile:
            data = metricsPattern.findall(line)
            if data:
                unpairedReadCount = int(data[0][0])
                pairedReadCount = int(data[0][1])
                unmappedReadCount = int(data[0][2])
                duplicatesFraction = float(data[0][3])
                break
        draw_histogram(sampleName, coverageValues, counts)
        outFile = open('%s.genomestats.txt' % sampleName, 'w')
        outFile.write('genome fraction with coverage >=%sX: %s%%\t(%s out of %s bases)\ngenome fraction with coverage >0X: ' \
            '%s%%\t(%s out of %s bases)\navg coverage of covered regions (with coverage >0X): %.1fX\n' \
            'median coverage of covered regions: %sX\nmax coverage: %sX\nstandard deviation: %.1fX\n\n' \
            'number of unpaired reads: %s\nnumber of read pairs: %s\nnumber of unmapped reads: %s\n' \
            'total number of reads: %s\nproportion of mapped reads: %.1f%%\nfraction of read duplications: %.2f%%\n' % (
                coverageThreshold,
                round((countsOverCovThreshold/sum(counts))*100,1),
                countsOverCovThreshold,
                sum(counts),
                round((len(coveragesCumulative)/sum(counts))*100,1),
                len(coveragesCumulative),
                sum(counts),
                numpy.mean(coveragesCumulative),
                numpy.median(coveragesCumulative),
                str(coveragesCumulative[-1]),
                numpy.std(coveragesCumulative),
                unpairedReadCount,
                pairedReadCount,
                unmappedReadCount,
                (unpairedReadCount+(pairedReadCount*2)+unmappedReadCount),
                round(((unpairedReadCount+(pairedReadCount*2))/(unpairedReadCount+(pairedReadCount*2)+unmappedReadCount))*100,1),
                round(duplicatesFraction,2)))
        summaryFile = open('%s.summarystats.txt' % sampleName, 'w')
        summaryFile.write('%s\t%s\t%s\t%.1f\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%.1f\t%.2f\n' % (
            sampleName,
            round((countsOverCovThreshold/sum(counts))*100,1),
            round((len(coveragesCumulative)/sum(counts))*100,1),
            numpy.mean(coveragesCumulative),
            numpy.median(coveragesCumulative),
            str(coveragesCumulative[-1]),
            numpy.std(coveragesCumulative),
            unpairedReadCount,
            pairedReadCount,
            unmappedReadCount,
            (unpairedReadCount+(pairedReadCount*2)+unmappedReadCount),
            round(((unpairedReadCount+(pairedReadCount*2))/(unpairedReadCount+(pairedReadCount*2)+unmappedReadCount))*100,1),
            round(duplicatesFraction,2)))
        return True
    except:
        traceback.print_exc()
        return False
    finally:
        inHistFile.close()
        inMetricsFile.close()
        outFile.close()

def draw_histogram(sampleName, xAxis, yAxis):
    try:
        plt.xlabel('X coverage')
        plt.ylabel('n bases')
        plt.hist(xAxis, xAxis, weights=yAxis)
        plt.savefig('%s_coverage_histogram.png' % sampleName, dpi=400)
        plt.close()
        return True
    except:
        traceback.print_exc()
        return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--sample-name", dest="sampleName", required=True,
                        help="Sample Name")
    parser.add_argument("-c", "--cov-threshold", dest="covThreshold", required=True,
                        help="Coverage Threshold")
    parser.add_argument("-b", "--bedtools-histogram", dest="covHistogram", required=True,
                        help="Genome coverage histogram file from bedtools")
    parser.add_argument("-p", "--picard-metrics", dest="picardMetrics", required=True,
                        help="Picard metrics file")

    args = parser.parse_args()
    sampleName = args.sampleName
    covThreshold = int(args.covThreshold)
    covHistogram = args.covHistogram
    picardMetrics = args.picardMetrics

    compute_statistics(sampleName, covThreshold, covHistogram, picardMetrics)
