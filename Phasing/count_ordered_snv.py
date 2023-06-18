'''
Count ordered SNV pairs
'''

import sys

def count(ifile:str, threshold:int) -> None:
    f = open(ifile, 'r')

    cnt = 0
    for i, line in enumerate(f):
        if i == 0: continue

        items = line.rstrip('\n').split('\t')


        if int(items[5]) > threshold and int(items[4]) > threshold and int(items[3]) == 0:
            cnt += 1
            print(items[0].split(',')[0], items[0].split(',')[1], items[1].split(',')[1], sep='\t')
            continue

        if int(items[5]) > threshold and int(items[3]) > threshold and int(items[4]) == 0:
            cnt += 1
            print(items[0].split(',')[0], items[0].split(',')[1], items[1].split(',')[1], sep='\t')
            continue

        if int(items[9]) > threshold and int(items[7]) > threshold and int(items[8]) == 0:
            cnt += 1
            print(items[0].split(',')[0], items[0].split(',')[1], items[1].split(',')[1], sep='\t')
            continue

        if int(items[9]) > threshold and int(items[8]) > threshold and int(items[7]) == 0:
            cnt += 1
            print(items[0].split(',')[0], items[0].split(',')[1], items[1].split(',')[1], sep='\t')
            continue

    print("# of ordered SNV pairs", cnt, sep='\t', file=sys.stderr)

    f.close()

def main():
    count(sys.argv[1], int(sys.argv[2]))

if __name__ == '__main__':
    main()
