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

    outFile = open('panSNPs_intervals.txt', 'w')
    for snp in SNPs:
        outFile.write("%s\t%s\t%s\n" % (str(snp[0]), str(snp[1]), str(snp[1])))
    outFile.close()


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


def filter_panSNPs(snpInFile, discRegFileGenerated, discRegFileCustom):
    try:
        with open(snpInFile, 'r') as inSnpFile, open('panSNPs_intervals_refined.txt', 'w') as outSnpFile:
            discRegions = []
            discRegsPattern = re.compile(r'(.+?)\t(\d+)\t(\d+)\n')
            if discRegFileGenerated:
                inDiscRegFileGenerated = open(discRegFileGenerated, 'r')
                for line in inDiscRegFileGenerated:
                    data = discRegsPattern.findall(line)
                    if not data:
                        sys.stderr.write('ERROR1!')
                        inDiscRegFileGenerated.close()
                        return False
                    discRegions.append([data[0][0], data[0][1], data[0][2]])
                inDiscRegFileGenerated.close()
            if discRegFileCustom:
                inDiscRegFileCustom = open(discRegFileCustom, 'r')
                for line in inDiscRegFileCustom:
                    data = discRegsPattern.findall(line)
                    if not data:
                        sys.stderr.write('ERROR2!')
                        inDiscRegFileCustom.close()
                        return False
                    discRegions.append([data[0][0], data[0][1], data[0][2]])
                inDiscRegFileCustom.close()
            discRegions.sort(key=lambda x:int(x[1]), reverse=False)

            snpPattern = re.compile(r'(.+?)\t(\d+).+')
            for line in inSnpFile:
                isDiscarded = False
                data = snpPattern.findall(line)
                if not data:
                    sys.stderr.write('ERROR3!')
                    return False
                for region in discRegions:
                    if str(region[0]) == str(data[0][0]) and int(region[1]) <= int(data[0][1]) <= int(region[2]):
                        isDiscarded = True
                        break
                    elif str(region[0]) == str(data[0][0]) and int(region[1]) > int(data[0][1]):
                            break
                if not isDiscarded:
                    outSnpFile.write(line)
    except:
        traceback.print_exc()
        return False


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gen-bank", dest="genBank", required=True,
                        help="GenBank File")

    args = parser.parse_args()
    genBank = args.genBank

    extract_panSNPs()
    get_intervals_to_be_discarded(genBank, ['rRNA', 'tRNA'])
    filter_panSNPs('panSNPs_intervals.txt', 'discarded_regions.txt', null)
