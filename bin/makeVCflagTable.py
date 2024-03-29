#!/usr/bin/env python

#created by Matej Medvecky

#VC_table.dat must be tab-limited table generated by either makeVCtable.py or makeVCtableNoMonomorphic.py

import sys
import re

def make_flag_table():
    if len(sys.argv) != 4:
        print("USAGE: %s <VC_table.dat> <coverage_threshold> <frequency threshold>\nCAUTION: Coverage threshold must be an integer value.\nCAUTION: Frequency threshold must be flat lower that 1.0. (e.g. frequency of 20%% must be specified as 0.2)." % sys.argv[0])
        return False

    with open(sys.argv[1], 'r') as snpFile, open('VC_flag_table.dat', 'w') as outFile:
        if float(sys.argv[3]) > 1.0:
            print("ERROR: Frequency threshold is higher than 1.0 (i.e. 100%).")
            return False
        else:
            freqThreshold = float(sys.argv[3])
        headerPattern = re.compile(r'(^chromosome\tposition\tref allele\t.+)')
        headerData = headerPattern.findall(snpFile.readline())
        if not headerData:
            print("ERROR: Correct header is missing! Are you sure you input table generated by makeVCtable.py script?")
            return False
        else:
            outFile.write("%s\n" % headerData[0])
        snps = []
        for line in snpFile:
            snps.append(list(line.rstrip().split("\t")))
        if not snps:
            print("ERROR: Something went wrong, list of snps is empty! Are you sure you input table generated by makeVCtable.py script?")
            return False

        snpPattern = re.compile(r'^(.+?):(\d+)\(.+?:(\d+)\(.+?\/(\d+)\(.+?\/(\d+)\(.+?\/(\d+)\(.+')
        for snp in snps:
            outFile.write("%s\t%d\t%s" % (str(snp[0]), int(snp[1]), str(snp[2])))
            for currSample in snp[3:]:
                snpData = snpPattern.findall(currSample)
                if snpData:
                    outFile.write("\t%s:" % str(snpData[0][0]))
                    isFirstCaution = True
                    if str(snpData[0][0]) == 'A':
                        if int(snpData[0][1]) > 0 and (float(snpData[0][2])/int(snpData[0][1])) < freqThreshold:
                            outFile.write(" below FREQUENCY THRESHOLD (%.2f)" % (float(snpData[0][2])/int(snpData[0][1])))
                            isFirstCaution = False
                    elif str(snpData[0][0]) == 'C':
                        if int(snpData[0][1]) > 0 and (float(snpData[0][3])/int(snpData[0][1])) < freqThreshold:
                            outFile.write(" below FREQUENCY THRESHOLD (%.2f)" % (float(snpData[0][3])/int(snpData[0][1])))
                            isFirstCaution = False
                    elif str(snpData[0][0])== 'G':
                        if int(snpData[0][1]) > 0 and (float(snpData[0][4])/int(snpData[0][1])) < freqThreshold:
                            outFile.write(" below FREQUENCY THRESHOLD (%.2f)" % (float(snpData[0][4])/int(snpData[0][1])))
                            isFirstCaution = False
                    elif str(snpData[0][0]) == 'T':
                        if int(snpData[0][1]) > 0 and (float(snpData[0][5])/int(snpData[0][1])) < freqThreshold:
                            outFile.write(" below FREQUENCY THRESHOLD (%.2f)" % (float(snpData[0][5])/int(snpData[0][1])))
                            isFirstCaution = False
                    if int(snpData[0][1]) < int(sys.argv[2]):
                        if isFirstCaution:
                            outFile.write(" below COVERAGE THRESHOLD (%d)" % int(snpData[0][1]))
                            isFirstCaution = False
                        else:
                            outFile.write("; below COVERAGE THRESHOLD (%d)" % int(snpData[0][1]))
                    if str(snpData[0][0]) == "AMBIGUOUS":
                        if isFirstCaution:
                            outFile.write(" AMBIGUOUS")
                            isFirstCaution = False
                        else:
                            outFile.write("; AMBIGUOUS")
                    if str(snpData[0][0]) == "INDEL":
                        if isFirstCaution:
                            outFile.write(" INDEL")
                            isFirstCaution = False
                        else:
                            outFile.write("; INDEL")
                    if isFirstCaution:
                        outFile.write(" OK")
                else:
                    print("ERROR: Wrong snp format! Are you sure you input table generated by makeVCtable.py script?")
                    return False
            outFile.write("\n")
make_flag_table()
