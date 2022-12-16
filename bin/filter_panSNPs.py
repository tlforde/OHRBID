#!/usr/bin/env python3

import re
import sys
import argparse

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
    parser.add_argument("-i", "--pansnps-intervals", dest="panSNPsIntervals", required=True,
                        help="pan SNPs intervals txt file")
    parser.add_argument("-d", "--discard-regions", dest="discardRegions", required=True,
                        help="txt file from discardRegions")
    parser.add_argument("-c", "--custom-regions", dest="customRegions", default=None, type=str,
                        help="Custom discard regions file")

    args = parser.parse_args()
    panSNPsIntervals = args.panSNPsIntervals
    discardRegions = args.discardRegions
    customRegions = args.customRegions

    filter_panSNPs(panSNPsIntervals, discardRegions, customRegions)
