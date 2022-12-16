#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def get_intervals_to_be_discarded(gbFile, omitMolsList):
    try:
        with open('discarded_regions.txt', 'w') as outFile:
            for record in SeqIO.parse(open(gbFile, 'r'), 'genbank'):
                for feature in record.features:
                    if feature.type in omitMolsList:
                        outFile.write('%s\t%s\t%s\n' % (record.id, feature.location.start, feature.location.end))
    except:
        traceback.print_exc()
        return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gen-bank", dest="genBank", required=True,
                        help="GenBank File")

    args = parser.parse_args()
    genBank = args.genBank

    get_intervals_to_be_discarded(genBank, ['rRNA', 'tRNA'])
