#!/usr/bin/env python3

import re
import sys
import argparse
from Bio import SeqIO

def extract_mutated_genes(snpInFile, gbFile):
    try:
        with open(snpInFile, 'r') as inSnpFile:
            snpPattern = re.compile(r'(.+?)\t(\d+)\t.+')
            snps = []
            for line in inSnpFile:
                data = snpPattern.findall(line)
                if not data:
                    sys.stderr.write('ERROR!')
                    return False
                snps.append([str(data[0][0]), int(data[0][1])])

            for record in SeqIO.parse(open(gbFile, 'r'), 'genbank'):
                for feature in record.features:
                    for snp in snps:
                        if snp[0] == record.id and snp[1] >= feature.location.start+1 and snp[1] <= feature.location.end:
                            if feature.type == 'CDS':
                                if 'EC_number' in feature.qualifiers:
                                    mutatedGenes.append({'id': record.id, 'locus_tag': feature.qualifiers['locus_tag'][0], 'start': feature.location.start+1, 'end': feature.location.end, 'strand': feature.location.strand, 'seq': feature.location.extract(record).seq, 'ec': feature.qualifiers['EC_number'], 'product': feature.qualifiers['product']})
                                else:
                                    mutatedGenes.append({'id': record.id, 'locus_tag': feature.qualifiers['locus_tag'][0], 'start': feature.location.start+1, 'end': feature.location.end, 'strand': feature.location.strand, 'seq': feature.location.extract(record).seq, 'ec': '-', 'product': feature.qualifiers['product']})
                            elif feature.type == 'misc_binding':
                                mutatedGenes.append({'id': record.id, 'locus_tag': 'misc_binding: bound_moiety=%s' % feature.qualifiers['bound_moiety'], 'start': feature.location.start+1, 'end': feature.location.end, 'strand': feature.location.strand, 'seq': feature.location.extract(record).seq, 'ec': '-', 'product': 'NA'})
                            elif feature.type == 'regulatory':
                                mutatedGenes.append({'id': record.id, 'locus_tag': 'regulatory: bound_moiety=%s' % feature.qualifiers['bound_moiety'], 'start': feature.location.start+1, 'end': feature.location.end, 'strand': feature.location.strand, 'seq': feature.location.extract(record).seq, 'ec': '-', 'product': 'NA'})
                            elif feature.type == 'misc_feature':
                                mutatedGenes.append({'id': record.id, 'locus_tag': feature.qualifiers['note'], 'start': feature.location.start+1, 'end': feature.location.end, 'strand': feature.location.strand, 'seq': feature.location.extract(record).seq, 'ec': '-', 'product': 'NA'})
    except:
        traceback.print_exc()
        return False


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gen-bank", dest="genBank", required=True,
                        help="GenBank File")
    parser.add_argument("-s", "--snp-file", dest="snpFile", required=True,
                        help="SNP File")

    args = parser.parse_args()
    genBank = args.genBank
    snpFile = args.snpFile

    extract_mutated_genes(snpFile, genBank)
