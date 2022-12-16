#!/usr/bin/env python3

import re
import glob

def extract_panSNPs():
    positionPattern = re.compile(r'(.+?)\t(\d+?)\t(.)\t.+')# reference contig, position, reference base
    SNPs = []
    for filename in glob.iglob('*varscan'):
        inFile = open(filename, 'r')
        for line in inFile:
            snpData = positionPattern.findall(line)
            isInSNPs = False
            if snpData:
                for snp in SNPs:
                    if int(snpData[0][1]) == snp[1]:
                        isInSNPs = True
                        break
                if not isInSNPs:
                    SNPs.append([snpData[0][0], int(snpData[0][1]), snpData[0][2]])
        inFile.close()
    SNPs.sort(key=lambda x: x[1], reverse=False)

    outFile = open('panSNPs.txt', 'w')
    for snp in SNPs:
        outFile.write("%s\t%s\t%s\n" % (str(snp[0]), str(snp[1]), str(snp[1])))
    outFile.close()


if __name__ == "__main__":

    extract_panSNPs()
