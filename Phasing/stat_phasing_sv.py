

import os
import sys
import re
import csv
from sys import argv

input1 = open(argv[1], 'r') ##nanomonsv.supporting_read_phased.txt
input2 = open(argv[2], 'r') ##nanomonsv.supporting_read_supplementary_phased.txt
out1 = open(argv[3], 'w') ##nanomonsv.phasing_summary.txt

chr_list = ['chr' + str(i) for i in range(1,23)]

svtype_dict = {}
unphased_dict = {}
ps_dict = {}
hp_dict = {}

### add SV info into dicts
header1 = next(csv.reader(input1, delimiter = '\t'))
for line1 in csv.reader(input1,delimiter='\t'):
    key1 = line1[header1.index('SV_Key')]
    chr1 = key1.split(',')[0]
    chr2 = key1.split(',')[3]
    svtype1 = line1[header1.index('SV_Type')]
    ps1 = str(line1[header1.index('PS_FLAG')])
    hp1 = str(line1[header1.index('HP_FLAG')])

    if chr1 in chr_list and chr2 in chr_list:
        if key1 not in svtype_dict:
            svtype_dict[key1] = svtype1
            ps_dict[key1] = set()
            hp_dict[key1] = {}
            unphased_dict[key1] = 0

        if ps1 != 'NA' and hp1 != 'NA':
            ps_dict[key1].add(ps1)
            if ps1 not in hp_dict[key1]:
                hp_dict[key1][ps1] = []
            hp_dict[key1][ps1].append(hp1)

        else:
            unphased_dict[key1] += 1


header1 = next(csv.reader(input2, delimiter = '\t'))
for line1 in csv.reader(input2,delimiter='\t'):
    key1 = line1[header1.index('SV_Key')]
    chr1 = key1.split(',')[0]
    chr2 = key1.split(',')[3]
    svtype1 = line1[header1.index('SV_Type')]
    ps1 = str(line1[header1.index('PS_FLAG')])
    hp1 = str(line1[header1.index('HP_FLAG')])

    if chr1 in chr_list and chr2 in chr_list:
        if key1 not in svtype_dict:
            svtype_dict[key1] = svtype1
            ps_dict[key1] = set()
            hp_dict[key1] = {}
            unphased_dict[key1] = 0

        if ps1 != 'NA' and hp1 != 'NA':
            ps_dict[key1].add(ps1)
            if ps1 not in hp_dict[key1]:
                hp_dict[key1][ps1] = []
            hp_dict[key1][ps1].append(hp1)

        else:
            unphased_dict[key1] += 1 


### make summary for each SV
out1.write('SV_Key' + '\t' + 'SV_Type' + '\t' + 'Phasing_summary(PS;read_ratio(1:2)|... | Unphased_read)' + '\n')

for key1 in svtype_dict:
  # readname_list = readname_dict[key1]
    svtype = svtype_dict[key1]
    ps_set = ps_dict[key1]
    unphased = str(unphased_dict[key1])
    ps_ratio = ''

    if len(ps_set) == 0:
        ps_ratio = str('Unphased;' + unphased)
    else:
        for ps1 in ps_set:
            ps_ratio += str('PS_' + ps1 + ';'+ str(hp_dict[key1][ps1].count('1')) + ':' + str(hp_dict[key1][ps1].count('2')) + '|')
        ps_ratio += str('Unphased;' + unphased)
    
    out1.write(key1 + '\t' + svtype + '\t' + ps_ratio + '\n')


input1.close()
input2.close()
out1.close()
