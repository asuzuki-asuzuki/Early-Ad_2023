

import os
import sys
import re
import csv
from sys import argv

input1 = open(argv[1], 'r') #nanomonsv.supporting_read_filt.txt
input2 = open(argv[2], 'r') #germlinesnps_pos_list.txt
input3 = open(argv[3], 'r') #samtools_mpileup_list_supplementary.txt
phase_block = open(argv[4], 'r') #blocklist
#svtype = open(argv[5], 'r') #nanomonsv.result.filt.svtype.txt
out1 = open(argv[6], 'w') #nanomonsv.supporting_read_supplementary.txt

readname_dict = {}
read_hp_dict = {}
snp_dict = {}

### make svtype_dict
svtype_dict = {} #{'chr,XXXXXX,+,chr,YYYYYY':'Deletion'}
with open(argv[5], 'r')as svtype:
    for num,line1 in enumerate(svtype):
        if num != 0:
            line1_list=line1[0:len(line1)-1].split("\t")
            key1 = (line1_list[0]+","+line1_list[1]+","+line1_list[2]+","+line1_list[3]+","+line1_list[4]+","+line1_list[5]+","+line1_list[6])
            if line1_list[13] == 'PASS':
                svtype_dict[key1] = line1_list[14]


### add SV info into dict
for line1 in csv.reader(input1,delimiter='\t'):
    key1 = line1[0]
    readname1 = line1[1]

    if readname1 not in readname_dict:
        readname_dict[readname1] = []
        readname_dict[readname1].append(key1)
    if readname1 not in read_hp_dict:
        read_hp_dict[readname1] = {}

### add germline SNPs and haplotype into dict
for line1 in csv.reader(input2,delimiter='\t'):
    chr1 = line1[0]
    pos1 = str(line1[1])
    ref1 = str(line1[2]).upper()
    alt1 = str(line1[3]).upper()
    info_list = line1[4].split(':')
    hp1_tmp = line1[5].split(':')[info_list.index('GT')]
    ps1 = line1[5].split(':')[info_list.index('PS')]

    if hp1_tmp == '1|0':
       hp_alt = '1'
       hp_ref = '2'
    elif hp1_tmp == '0|1':
        hp_ref = '1'
        hp_alt = '2'
    else:
        hp1 = ''
        hp2 = ''
    snp_dict[chr1 + '_' + pos1 + '_' + alt1] = chr1 + '_' + ps1 + '_' + hp_alt
    snp_dict[chr1 + '_' + pos1 + '_' + ref1] = chr1 + '_' + ps1 + '_' + hp_ref


### make block_dict
block_dict = {} #{PS_tag:chr&start&end}
header1 = next(phase_block)
for line1 in csv.reader(phase_block,delimiter='\t'):
    chr1 = line1[1]
    ps_tag1 = str(line1[2])
    start1 = str(line1[3])
    end1 = str(line1[4])
    block_dict[ps_tag1] = chr1 + '&' + start1 + '&' + end1

### add samtools_mpileup_info
for line1 in csv.reader(input3, delimiter='\t'):
    chr1 = line1[0]
    pos1 = str(line1[1])

    base_raw = str.upper(line1[4])
    read_list = line1[6].split(',')

    symbol_list = ['^']
    for symbol1 in symbol_list:
        while base_raw.find(symbol1) != -1:
            base_tmp_list = list(base_raw)
            base_tmp_list[base_tmp_list.index(symbol1):base_tmp_list.index(symbol1)+2] = ''
            base_raw = ''.join(base_tmp_list)

    base_raw = ''.join(list(base_raw.replace('$', '')))

    base_tmp = ''
    m = re.finditer(r'\d{1,10}', base_raw)
    n = 0
    for i in m:
        base_tmp += base_raw[n:int(i.start())-1]
        n = int(i.start()) + int(i.group()) + len(str(i.group()))
    base_tmp += base_raw[n:]

    base_list = list(base_tmp)

    if len(base_list) != len(read_list):
        print(line1)
        print(''.join(base_list) + '___' + str(len(read_list)))

    for i in range(0, len(read_list)):
        read1 = read_list[i]
        base1 = base_list[i].upper()
    
        if chr1 + '_' + pos1 + '_' + base1 in snp_dict:
            info1 = snp_dict[chr1 + '_' + pos1 + '_' + base1]
            if info1 not in read_hp_dict[read1]:
                read_hp_dict[read1][info1] = 0
            read_hp_dict[read1][info1] += 1


### make summary for each SV
out1.write('SV_Key' + '\t' + 'Readname' + '\t' + 'Info_1' + '\t' + 'SV_Type' + '\t' + 'PS_FLAG' + '\t' 'HP_FLAG' + '\t' + 'PS_Start' + '\t' + 'PS_End' + '\n')

for read1 in read_hp_dict:
    d = read_hp_dict[read1]
    max_kv_list = [kv for kv in d.items() if kv[1] == max(d.values())]

    if len(max_kv_list) == 1:
        i = max_kv_list[0]
        chr_tag1 = i[0].split('_')[0]
        ps_tag1 = str(i[0].split('_')[1])
        hp_tag1 = str(i[0].split('_')[2])
        n_snp1 = int(d[i[0]])
        if hp_tag1 == '1':
            i_tmp = chr_tag1 + '_' + ps_tag1 + '_' + '2'
        else:
            i_tmp = chr_tag1 + '_' + ps_tag1 + '_' + '1'

        if i_tmp in d:
            n_snp2 = int(d[i_tmp])
        else:
            n_snp2 = 0

        if n_snp1/(n_snp1 + n_snp2) >= 0.7 and n_snp1 >= 2:
            ps_chr1 = block_dict[ps_tag1].split('&')[0]
            ps_start1 = int(block_dict[ps_tag1].split('&')[1])
            ps_end1 = int(block_dict[ps_tag1].split('&')[2])

            for key1 in readname_dict[read1]:
                try:
                    sv_type1 = svtype_dict[key1]
                    chr1 = key1.split(',')[0]
                    pos1 = int(key1.split(',')[1])
                    chr2 = key1.split(',')[3]
                    pos2 = int(key1.split(',')[4])

                    if (chr1 == ps_chr1 and ps_start1 <= pos1 <= ps_end1) or (chr2 == ps_chr1 and ps_start1 <= pos2 <= ps_end1):
                        result_list = [key1, read1, 'Supplementary', sv_type1, ps_tag1, hp_tag1, str(ps_start1), str(ps_end1)]
                        out1.write('\t'.join(result_list) + '\n')
                except:
                    print(svtype_dict)
                    print(key1)
                    sv_type1 = svtype_dict[key1]


input1.close()
input2.close()
input3.close()
phase_block.close()
out1.close()
