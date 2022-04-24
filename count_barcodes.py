#!/usr/bin/env python3

### Input will be piped in from zcat (as if from stdin), and the output file
### should be supplied as an argument, like such:
### zcat .../SRRXXXXXXX.fastq.gz | ./count_barcodes.py .../SRRXXXXXXX.pickle

import sys
import fileinput
from collections import Counter
import pickle

if len(sys.argv) != 2:
    raise ValueError('Single argument should be output pickle filepath')

bc_counts = Counter() # dict with def. value 0
i = 0
for line in fileinput.input('-'):
    # only the second of every 4 lines has the actual sequence
    i += 1
    if i % 4 != 2:
        continue

    barcode = line[:19] # barcodes are the first 19 bp of the seq
    bc_counts[barcode] += 1

with open(sys.argv[1], 'wb') as f:
    pickle.dump(bc_counts, f)
