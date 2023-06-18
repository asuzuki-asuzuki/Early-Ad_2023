
import os
import sys
import re
import csv
import datetime
from sys import argv
import pysam

print('extracting uniq reads was started_' + str(datetime.datetime.now()))

### when using bam files, change 'r' to 'rb'!
input1 = pysam.AlignmentFile(argv[1], 'rb')
### when using bam files, change 'w' to 'wb'!
out1 = pysam.AlignmentFile(argv[2], 'wb', template = input1)

readname_dict = {}

for read1 in input1:
    name1 = read1.query_name
    if name1 not in readname_dict:
        readname_dict[name1] = 0
    readname_dict[name1] += 1


input1.close()

input1 = pysam.AlignmentFile(argv[1], 'rb')


### extract uniq reads from supplementary alignment
for read1 in input1:
    name1 = read1.query_name

    if readname_dict[name1] == 1:
        out1.write(read1)

input1.close()
out1.close()

print('extracting uniq reads was successfully finished_' + str(datetime.datetime.now()))

