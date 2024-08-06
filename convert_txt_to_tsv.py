import os
import sys

fname = f'{os.getcwd()}/{sys.argv[1]}'

'''
this script is essentially just a hacky tool to convert my copy-pasted
lcms results from the txt file into something that can be easily copy-pasted
into google sheets
'''

results = []
with open(fname) as f:
    split_lines = [l.strip().split() for l in f.readlines() if l.strip()]
    for idx,l in enumerate(split_lines):
        if len(l) > 2:
            if l[1] == 'Name':
                tb = float(split_lines[idx+1][8])
                caff = float(split_lines[idx+2][8])

                results.append((tb, caff))
                
                # cga = float(split_lines[idx+3][8])
                # results.append((tb, caff, cga))

for tb, caff in results:
    print(f'{tb}\t{caff}')

# for tb, caff, cga in results:
#     print(f'{tb}\t{caff}\t{cga}')
