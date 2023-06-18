
import os
import sys
import re
import csv
import datetime
from sys import argv
import pysam

print('extracting SV reads for nanomonsv was started_' + str(datetime.datetime.now()))

# input file (bam)
input_bam = pysam.AlignmentFile(argv[1], 'rb')
# nanomonsv.supporting_read.txt
input_nanomonsv = open(argv[2],'r')
# output file (bam)
out_primary_bam = pysam.AlignmentFile(argv[3], "wb", template = input_bam)
out_supplementary_bam = pysam.AlignmentFile(argv[4], "wb", template = input_bam)


csv.field_size_limit(1000000000)
readname_dict = {}

# extract readname from input_readname
for line1 in csv.reader(input_nanomonsv, delimiter = '\t'):
    key1 = line1[0]
    readname1 = line1[1]

    if readname1 not in readname_dict:
        readname_dict[readname1] = []
    readname_dict[readname1].append(key1)

# extract SV supporting reads from bam
for read1 in input_bam:
    readname1 = str(read1.query_name)
    
    if readname1 in readname_dict:
        sv_list = readname_dict[readname1]
        read_chr1 = read1.reference_name
        start1 = int(read1.reference_start)
        end1 = int(read1.reference_end)

        for sv1 in sv_list:
            chr1 = sv1.split(',')[0]
            pos1 = int(sv1.split(',')[1])
            chr2 = sv1.split(',')[3]
            pos2 = int(sv1.split(',')[4])

            if (chr1 == read_chr1 and start1 <= pos1 <= end1) or (chr2 == read_chr1 and start1 <= pos2 <= end1):
                if read1.flag == 0 or read1.flag == 16:
                    out_primary_bam.write(read1)
                elif str(read1.is_supplementary) == 'True':
                    out_supplementary_bam.write(read1)


input_bam.close()
input_nanomonsv.close()
out_primary_bam.close()
out_supplementary_bam.close()

print('extracting SV reads for nanomonsv was successfully finished_' + str(datetime.datetime.now()))


