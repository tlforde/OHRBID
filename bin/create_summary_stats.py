#!/usr/bin/env python3

import argparse

def create_file(coverageThreshold):

    with open('summary_statistics.txt', 'w') as summaryStatsFile:
        summaryStatsFile.write('#columns represent following attributes:\n#sample name; genome fraction with coverage >=%sX; ' \
                               'genome fraction with coverage >0X; average coverage of covered regions (with coverage >0X); ' \
                               'median coverage of covered regions; maximum coverage; coverage standard deviation; number of unpaired reads; ' \
                               'number of read pairs; number of unmapped reads; total number of reads; proportion of mapped reads; ' \
                               'fraction of read duplications\n' % coverageThreshold)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cov-threshold", dest="covThreshold", required=True,
                        help="Coverage Threshold")

    args = parser.parse_args()
    covThreshold = args.covThreshold

    create_file(covThreshold)
