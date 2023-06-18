
import os
import sys
import re
import pysam
import datetime
from sys import argv

print('extract_phased_read was started_' + str(datetime.datetime.now()))

### when using bam files, change 'r' to 'rb'!
input1 = pysam.AlignmentFile(argv[1], 'rb')
### when using bam files, change 'w' to 'wb'!
out1 = open(argv[2], 'w')
out1.write('Read_name' + '\t' + 'FLAG' + '\t' + 'Chromosome' + '\t' + 'Position' + '\t' + 'MAPQ' + '\t' + 'Read_length' + '\t' + 'PS_tag' + '\t' + 'HP_tag' + '\n')

for read1 in input1:
    tag_list = read1.tags
    ps_tag1 = ''
    hp_tag1 = ''
    read_name1 = ''
    for tag1 in tag_list:
        if tag1[0] == 'PS':
            ps_tag1 = str(tag1[1])
        elif tag1[0] == 'HP':
            hp_tag1 = str(tag1[1])

    if ps_tag1 != '' and hp_tag1 != '':
        read_name1 = read1.query_name
        flag1 = str(read1.flag)
        chr1 = str(read1.reference_name)
        mapq1 = str(read1.mapq)
        pos1 = str(read1.reference_start + 1)
        length1 = str(read1.infer_read_length())
        out1.write(read_name1 + '\t' + flag1 + '\t' + chr1 + '\t' + pos1 + '\t' + mapq1 + '\t' + length1 + '\t' + ps_tag1 + '\t' + hp_tag1 + '\n')

input1.close()
out1.close()
print('extract_phased_read was successfully finished!' + str(datetime.datetime.now()))
