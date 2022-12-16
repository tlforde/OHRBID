#!/usr/bin/env python

#created by Matej Medvecky

#each line in snp_file.txt must consist of following three tab-delimited info: chromosome position referece_allele
#bam-readcount files must have 'bam-readcount' suffix and must be derived from the same reference genome as snp_file!

import sys
import re
import glob

def make_variant_calling_table():
    if len(sys.argv) != 3:
        print("USAGE: %s <snp_file.txt> </path/to/dir/with/*bam-readcount/files>" % sys.argv[0])
        return False

    with open(sys.argv[1], 'r') as snpFile, open('VC_table.dat', 'w') as outFile:
        snpPattern = re.compile(r'(.+?)\t(\d+?)\t(.)\n')
        snps = []
        for line in snpFile:
            snpData = snpPattern.findall(line)
            if snpData:
                snps.append([snpData[0][1], snpData[0][0], snpData[0][2], []])
            else:
                print("ERROR: SNP entry has wrong format!\nLine: %s\nExiting.." % line)
                return False

        outFile.write("chromosome\tposition\tref allele")
        namePattern = re.compile(r'([^\/]+)\.bam-readcount$')
        brPattern = re.compile(r'(.+?)\t(\d+)\t(.)\t(\d+)\t.+?\tA:(\d+):.+?:.+?:.+?:(\d+):(\d+).+?\tC:(\d+):.+?:.+?:.+?:(\d+):(\d+).+?\tG:(\d+):.+?:.+?:.+?:(\d+):(\d+).+?\tT:(\d+):.+?:.+?:.+?:(\d+):(\d+).+?\tN:(\d+):.+?:.+?:.+?:(\d+):(\d+).+')#pattern for matching read count matrics for particular bases in bam-readcount output
        for filename in glob.glob(sys.argv[2] + '/*bam-readcount'):
            name = namePattern.findall(filename)
            if name:
                outFile.write("\t%s" % name[0])
            else:
                print("ERROR: Could not extract name of the file %s. Exiting.." % filename)
                return False
            try:
                f = open(filename, 'r')
                for snp in snps:
                    f.seek(0)
                    for line in f:
                        brData = brPattern.findall(line)
                        if brData:
                            isInBrData = False
                            doesMatchPosition = False
                            if snp[0] == brData[0][1]:
                                doesMatchPosition = True
                                if snp[1] == brData[0][0]:
                                    if snp[2] == brData[0][2]:
                                        snp[3].append([brData[0][3], brData[0][4], brData[0][5], brData[0][6], brData[0][7], brData[0][8], brData[0][9], brData[0][10], brData[0][11], brData[0][12], brData[0][13], brData[0][14], brData[0][15], brData[0][16], brData[0][17], brData[0][18]])
                                        isInBrData = True
                                        break
                                    else:
                                        print("ERROR: Reference allele in SNP_file is differing from reference allele in bam-readcount file %s! Exiting.." % filename)
                                        return False
                                else:
                                    continue
                        else:
                            print("ERROR: bam-readcount entry has wrong format!\nLine: %s\nExiting.." % line)
                            return False
                    if doesMatchPosition and not isInBrData:
                        print("CAUTION: SNP with position %d was found in %s file only on different chromosome and therefore SNP was not called. Possible reference naming problem?" % (snp[0], filename))
                    if not isInBrData:
                        snp[3].append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            except IOError:
                print("ERROR: Could not open %s file. Exiting...\n" % f)
                return False
            finally:
                f.close()
        outFile.write("\n")

        for snp in snps:
            outFile.write("%s\t%d\t%s" % (snp[1], int(snp[0]), snp[2]))
            for sampleSnp in snp[3]:
                if int(sampleSnp[0]) > 0:
                    if ((int(sampleSnp[1])+int(sampleSnp[4])+int(sampleSnp[7])+int(sampleSnp[10])+int(sampleSnp[13]))/float(sampleSnp[0])) < 0.5:
                        outFile.write("\tINDEL:")#write INDEL if more than 50% of mapped reads have indel at given site
                    elif int(sampleSnp[1]) > int(sampleSnp[4]) and int(sampleSnp[1]) > int(sampleSnp[7]) and int(sampleSnp[1]) > int(sampleSnp[10]) and int(sampleSnp[1]) > int(sampleSnp[13]):
                        outFile.write("\tA:")
                    elif int(sampleSnp[4]) > int(sampleSnp[1]) and int(sampleSnp[4]) > int(sampleSnp[7]) and int(sampleSnp[4]) > int(sampleSnp[10]) and int(sampleSnp[4]) > int(sampleSnp[13]):
                        outFile.write("\tC:")
                    elif int(sampleSnp[7]) > int(sampleSnp[1]) and int(sampleSnp[7]) > int(sampleSnp[4]) and int(sampleSnp[7]) > int(sampleSnp[10]) and int(sampleSnp[7]) > int(sampleSnp[13]):
                        outFile.write("\tG:")
                    elif int(sampleSnp[10]) > int(sampleSnp[1]) and int(sampleSnp[10]) > int(sampleSnp[4]) and int(sampleSnp[10]) > int(sampleSnp[7]) and int(sampleSnp[10]) > int(sampleSnp[13]):
                        outFile.write("\tT:")
                    elif int(sampleSnp[13]) > int(sampleSnp[1]) and int(sampleSnp[13]) > int(sampleSnp[4]) and int(sampleSnp[13]) > int(sampleSnp[7]) and int(sampleSnp[13]) > int(sampleSnp[10]):
                        outFile.write("\tN:")
                    else:
                        outFile.write("\tAMBIGUOUS:")
                elif int(sampleSnp[0]) == 0:
                    outFile.write("\t-:")
                outFile.write("%d(%d,%d):%d(%d,%d)/%d(%d,%d)/%d(%d,%d)/%d(%d,%d)/%d(%d,%d)/%d" % (int(sampleSnp[0]), (int(sampleSnp[2])+int(sampleSnp[5])+int(sampleSnp[8])+int(sampleSnp[11])+int(sampleSnp[14])), (int(sampleSnp[3])+int(sampleSnp[6])+int(sampleSnp[9])+int(sampleSnp[12])+int(sampleSnp[15])), int(sampleSnp[1]), int(sampleSnp[2]), int(sampleSnp[3]), int(sampleSnp[4]), int(sampleSnp[5]), int(sampleSnp[6]), int(sampleSnp[7]), int(sampleSnp[8]), int(sampleSnp[9]), int(sampleSnp[10]), int(sampleSnp[11]), int(sampleSnp[12]), int(sampleSnp[13]), int(sampleSnp[14]), int(sampleSnp[15]), (int(sampleSnp[0])-(int(sampleSnp[1])+int(sampleSnp[4])+int(sampleSnp[7])+int(sampleSnp[10])+int(sampleSnp[13])))))
            outFile.write("\n")

make_variant_calling_table()
