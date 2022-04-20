#!/usr/bin/env python3

### Input will be piped in from zcat (as if from stdin), and the output file
### should be supplied as an argument, like such:
### zcat .../SRRXXXXXXX.fastq.gz | ./count_barcodes.py .../SRRXXXXXXX.pickle

import sys
import fileinput
from collections import Counter
#import matplotlib.pyplot as plt
import pickle

if len(sys.argv) != 2:
    raise ValueError('Single argument should be output pickle filepath')

bc_counts = Counter()
i = 0
for line in fileinput.input('-'):
    i += 1
    # try:
    #     line = input()
    # except EOFError:
    #     break # no more input, on to the next step

    # only the second of every 4 lines has the actual sequence
    if i % 4 != 2:
        continue # only the second of every 4 lines has the actual sequence
    barcode = line[:19]
    bc_counts[barcode] += 1

    # if i >= 4000:
    #     break

with open(sys.argv[1], 'wb') as f:
    pickle.dump(bc_counts, f)

#plt.hist(bc_counts.values())
#plt.show()
