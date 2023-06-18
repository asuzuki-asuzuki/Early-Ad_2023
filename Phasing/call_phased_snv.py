'''
Call and count phased SNVs using by pileup file
'''

import argparse
import sys

class pileup_info:
    def __init__(self):
        self.prev_chrom = ''
        self.prev_pos = 0
        self.prev_block_num = 0
        self.prev_haplo_info = []
        self.chrom = ''
        self.prev_pos = 0
        self.block_num = 0
        self.haplo_info = []

    def update_prev_info(self):
        self.prev_chrom = self.chrom
        self.prev_pos = self.pos
        self.prev_block_num = self.block_num
        self.prev_haplo_info = self.haplo_info

    def update_current_info(self, chrom, pos, block_num, haplo_info):
        self.chrom = chrom
        self.pos = pos
        self.block_num = block_num
        self.haplo_info = haplo_info

class snv_info:
    def __init__(self):
        self.ref = ''
        self.alt = ''
        self.vaf = 0.0

    def insert(self, ref, alt, vaf):
        self.ref = ref
        self.alt = alt
        self.vaf = vaf

def load_snv(ifile:str) -> dict:

    f = open(ifile, 'r', encoding="utf-8")

    out = dict()

    for line in f:
        if line[0] == '#': continue

        items = line.rstrip('\n').split('\t')
        if len(items[3]) != 1: continue
        if len(items[4]) != 1: continue

        info = snv_info()
        r, v = items[-1].split(':')[1].split(',')
        info.insert(items[3], items[4], int(v) / (int(r) + int(v))) 

        key = items[0] + ':' + items[1]

        out[key] = info

    f.close()

    count_SNV = dict()

    for key in out:
        chrom = key.split(':')[0]
        if chrom in ['chrEBV', 'chrM']: continue
        if chrom in count_SNV:
            count_SNV[chrom] += 1
        else:
            count_SNV[chrom] = 1

    total = 0

    for num in count_SNV.values():
        total += num - 1
    print("Total SNV pairs:", total, sep='\t', file=sys.stderr) 

    return out

def count_phased_snv(ifile:str, snv_dict:dict) -> None:

    f = open(ifile, 'r')
    info = pileup_info()
    print('Prev_info', 'curr_info', '1_0|0', '1_0|1', '1_1|0', '1_1|1', '2_0|0', '2_0|1', '2_1|0', '2_1|1', sep='\t')
    for i, line in enumerate(f):
        items = line.rstrip('\n').split('\t')

        chrom = items[0]
        pos = items[1]
        num = int(items[3])
        p_base = items[4]
        rnames = items[6].split(',')
        haplo_num = items[7].split(',')
        block_num = items[8].split(',')

        if num == 0: continue

        # Parse p_base
        j = 0
        base = list()
        while j < len(p_base):
            if p_base[j] in ['*', '#']:
                base.append(p_base[j])
                j += 1
            elif p_base[j] == '^':
                j += 2
            elif p_base[j] in ['+', '-']:
                l = j + 1
                n = ''
                while p_base[l] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    n += p_base[l]
                    l += 1
                j = l + int(n)
            elif p_base[j] in ['>', '<', '$']:
                j += 1
            else:
                base.append(p_base[j])
                j += 1

        if len(base) != num:
            print('Some errors occur in parsing p_base', chrom, pos, base, num, sep='\t', file=sys.stderr)
            sys.exit(1)

        h_info = list()
        b_num = 0
        for k in range(0, num):
            if haplo_num[k] == '*': continue
            if b_num == 0:
                b_num = block_num[k]
            h_info.append((rnames[k], base[k], haplo_num[k]))

        info.update_current_info(chrom, pos, b_num, h_info)

        if i == 0 or info.prev_block_num != info.block_num or len(h_info) == 0:
            info.update_prev_info()
            continue

        tmp_dict = {'1':dict(), '2':dict()}

       
        prev_snv = snv_dict[info.prev_chrom + ':' + str(info.prev_pos)]
        print(info.prev_chrom, info.prev_pos, prev_snv.ref, prev_snv.alt, prev_snv.vaf, sep=',', end='\t')

        curr_snv = snv_dict[info.chrom + ':' + str(info.pos)]
        print(info.chrom, info.pos, curr_snv.ref, curr_snv.alt, curr_snv.vaf, sep=',', end='\t')

        for item in info.prev_haplo_info:
            if item[1].upper() == prev_snv.ref.upper():
                tmp_dict[item[2]][item[0]] = [0]
            elif item[1].upper() == prev_snv.alt.upper():
                tmp_dict[item[2]][item[0]] = [1]

        for item in info.haplo_info:
            if item[0] not in tmp_dict[item[2]]: continue
            if item[1].upper() == curr_snv.ref.upper():
                tmp_dict[item[2]][item[0]].append(0)
            elif item[1].upper() == curr_snv.alt.upper():
                tmp_dict[item[2]][item[0]].append(1)

        for hp in tmp_dict:
            out = [0, 0, 0, 0]
            '''
            Index    Mean
              0       0|0
              1       0|1
              2       1|0
              3       1|1
            '''
            for rname in tmp_dict[hp]:
                if tmp_dict[hp][rname] == [0, 0]: out[0] += 1
                elif tmp_dict[hp][rname] == [0, 1]: out[1] += 1
                elif tmp_dict[hp][rname] == [1, 0]: out[2] += 1
                elif tmp_dict[hp][rname] == [1, 1]: out[3] += 1
            
            if hp == '1':
                print(out[0], out[1], out[2], out[3], sep='\t', end='\t')
            elif hp == '2':
                print(out[0], out[1], out[2], out[3], sep='\t')


        info.update_prev_info()

    f.close()


def main():
    snv_file = sys.argv[1]
    pile_file = sys.argv[2]

    snv_dict = load_snv(snv_file)
    count_phased_snv(pile_file, snv_dict)

if __name__ == '__main__':
    main()
