

import os
import sys
import re
import csv
from sys import argv

svtype = open(argv[1], 'r') #nanomonsv.result.filt.svtype.txt
support_read = open(argv[2], 'r') #nanomonsv.supporting_read.txt
phasing = open(argv[3], 'r') #phased_read.tsv
phase_block = open(argv[4], 'r') #germlinesnvs.t_readphased.vcf.blocklist
out1 = open(argv[5], 'w') #nanomonsv.supporting_read_phased.txt


### make svtype_dict
svtype_dict = {} #{'chr,XXXXX,+,chr,YYYYY':'Deletion'}
header1 = next(csv.reader(svtype, delimiter = '\t'))
for line1 in csv.reader(svtype,delimiter='\t'):
    key1 = ','.join(line1[0:7])
    if line1[header1.index('Is_filter')] == 'PASS':
        svtype_dict[key1] = line1[header1.index('SV_Type')]

### make phase_dict
phase_dict = {} #{readname:PS_tag&HP_tag}
header1 = next(csv.reader(phasing, delimiter = '\t'))
for line1 in csv.reader(phasing,delimiter='\t'):
    readname1 = line1[header1.index('Read_name')]
    ps_tag1 = line1[header1.index('PS_tag')]
    hp_tag1 = line1[header1.index('HP_tag')]
    
    ps_tag1 = str(line1[6])
    hp_tag1 = str(line1[7])
    phase_dict[readname1] = ps_tag1 + '&' + hp_tag1

### make block_dict
block_dict = {} #{PS_tag:chr&start&end}
header1 = next(phase_block)
for line1 in csv.reader(phase_block,delimiter='\t'):
    chr1 = line1[1]
    ps_tag1 = str(line1[2])
    start1 = str(line1[3])
    end1 = str(line1[4])
    block_dict[ps_tag1] = chr1 + '&' + start1 + '&' + end1

### add svtype and phasing information for SV supporting reads
out1.write('SV_Key' + '\t' + 'Readname' + '\t' + 'Info_1' + '\t' + 'Info_2' + '\t' + 'SV_Type' + '\t' + 'PS_FLAG' + '\t' + 'HP_FLAG' + '\n')
for line1 in csv.reader(support_read,delimiter = '\t'):
    svtype_flag = 'NA'
    ps_flag = 'NA'
    hp_flag = 'NA'

    key1 = line1[0]
    chr1 = key1.split(',')[0]
    pos1 = int(key1.split(',')[1])
    chr2 = key1.split(',')[3]
    pos2 = int(key1.split(',')[4])

    readname1 = line1[1]
    if key1 in svtype_dict:
        svtype_flag = svtype_dict[key1]

        if readname1 in phase_dict:
            ps_flag = phase_dict[readname1].split('&')[0]
            hp_flag = phase_dict[readname1].split('&')[1]

        if ps_flag != 'NA':
            ps_chr1 = block_dict[ps_flag].split('&')[0]
            ps_start1 = int(block_dict[ps_flag].split('&')[1])
            ps_end1 = int(block_dict[ps_flag].split('&')[2])

            if (chr1 == ps_chr1 and ps_start1 <= pos1 <= ps_end1) or (chr2 == ps_chr1 and ps_start1 <= pos2 <= ps_end1):
                out1.write('\t'.join(line1) + '\t' + svtype_flag + '\t' +  ps_flag + '\t' + hp_flag + '\n')
            else:
                out1.write('\t'.join(line1) + '\t' + svtype_flag + '\t' + 'NA' + '\t' + 'NA' + '\n')
        else:
                out1.write('\t'.join(line1) + '\t' + svtype_flag + '\t' + '\t' + 'NA' + '\t' + 'NA' + '\n')


svtype.close()
support_read.close()
phasing.close()
out1.close()
