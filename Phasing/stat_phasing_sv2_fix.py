

import os
import sys
import re
import csv
from sys import argv

input1 = open(argv[1], 'r') #nanomonsv.supporting_read_phased_summary.txt
out1 = open(argv[2], 'w') #nanomonsv.supporting_read_phased_summary2.txt
phaseblock = open(argv[3], 'r') #germlinesnvs.t_readphased.vcf.blocklist

phaseblock_dict = {}

### make phaseblock dict
header1 = next(phaseblock)

for line1 in csv.reader(phaseblock,delimiter='\t'):
    ps_name1 = 'PS_' + line1[2]
    chr1 = str(line1[1])
    start1 = str(line1[3])
    end1 = str(line1[4])
    phaseblock_dict[ps_name1] = chr1 + '&' + start1 + '&' + end1

### add SV info into dicts
header1 = next(csv.reader(input1, delimiter = '\t'))
out1.write('\t'.join(header1[:-1]) + '\t' + 'Phaseblock_breakpoint1' + '\t' + 'Haplotype_breakpoint1'+ '\t' + 'Phaseblock_breakpoint2' + '\t' + 'Haplotype_breakpoint2' + '\n')

for line1 in csv.reader(input1,delimiter='\t'):
    chr1_sv = line1[header1.index('SV_Key')].split(',')[0]
    pos1_sv = int(line1[header1.index('SV_Key')].split(',')[1])
    chr2_sv = line1[header1.index('SV_Key')].split(',')[3]
    pos2_sv = int(line1[header1.index('SV_Key')].split(',')[4])
    stat1 = line1[header1.index('Phasing_summary(PS;read_ratio(1:2)|... | Unphased_read)')].split('|')
    n_ps1 = 0
    ps_dict = {}
    flag1 = ''
    ps_list = []
    breakpoint1 = 'Unphased'
    breakpoint2 = 'Unphased'
    haplotype1 = 'NA'
    haplotype2 = 'NA'

    for ps1 in stat1:
        if ps1.find('PS_') != -1:
            ps_name1 = ps1.split(';')[0]
            n_ps1 += 1
            ratio1 = ps1.split(';')[1].split(':')
            chr1 = phaseblock_dict[ps_name1].split('&')[0]
            start1 = str(phaseblock_dict[ps_name1].split('&')[1])
            end1 = str(phaseblock_dict[ps_name1].split('&')[2])

            if max(int(ratio1[0]),int(ratio1[1])) >=3 and max(int(ratio1[0]),int(ratio1[1]))/(int(ratio1[0])+int(ratio1[1])) >= 0.7:
                if chr1_sv == chr1 and int(start1) <= pos1_sv <= int(end1):
                    breakpoint1 = ps_name1
                    if int(ratio1[0]) > int(ratio1[1]):
                        haplotype1 = '1'
                    elif int(ratio1[1]) > int(ratio1[0]):
                        haplotype1 = '2'
                if chr2_sv == chr1 and int(start1) <= pos2_sv <= int(end1):
                    breakpoint2 = ps_name1
                    if int(ratio1[0]) > int(ratio1[1]):
                        haplotype2 = '1'
                    elif int(ratio1[1]) > int(ratio1[0]):
                        haplotype2 = '2'

    out1.write('\t'.join(line1[:-1]) + '\t' + str(breakpoint1) + '\t' + str(haplotype1) + '\t' + str(breakpoint2) + '\t' + str(haplotype2) +'\n')


input1.close()
out1.close()
