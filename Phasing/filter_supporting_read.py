

import os
import sys
import re
import csv
from sys import argv

input1 = open(argv[1], 'r') #nanomonsv.result.filt.svtype.txt
input2 = open(argv[2], 'r') #nanomonsv.supporting_read.txt
out1 = open(argv[3], 'w') #nanomonsv.supporting_read_filt.txt

sv_set = set()

### make sv_set for passed SVs
header1 = next(csv.reader(input1, delimiter = '\t'))
for line1 in csv.reader(input1,delimiter='\t'):
    if line1[header1.index('Is_filter')] == 'PASS':
        sv_id1 = ','.join(line1[0:7])
        sv_set.add(sv_id1)
### filter supporting read results
with open(argv[3], 'w') as out1:
    for line2 in csv.reader(input2,delimiter='\t'):
        line2_str=(line2[0]+","+line2[1]+","+line2[2]+","+line2[3]+","+line2[4]+","+line2[5]+","+line2[6])
        if line2_str in sv_set:
            out1.write(line2_str+"\t"+line2[8]+"\t"+line2[9]+"\t"+line2[10]+"\n")

input1.close()
input2.close()
out1.close()

