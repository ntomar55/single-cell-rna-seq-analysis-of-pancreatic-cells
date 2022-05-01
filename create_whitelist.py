#!/usr/bin/env python3

import pickle
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

pickle_dir="/projectnb2/bf528/users/saxophone/data_p4/fastq_bc/"

#barcode_sets = []
bc_counts = []
whitelist = set()

samples = ['SRR3879604', 'SRR3879605', 'SRR3879606']
for i in range(3):
    with open(pickle_dir + samples[i] + '.pickle', 'rb') as p:
        bc_counts = pickle.load(p)
    # Filter barcodes with less than 10^4.75 reads
    cutoff = 10 ** 4.75
    ind_whitelist = set()
    for bc in bc_counts:
        if bc_counts[bc] > cutoff:
            whitelist.add(bc)
            ind_whitelist.add(bc)
    with open("/projectnb2/bf528/users/saxophone/data_p4/" + samples[i]
              + "_whitelist.txt", 'w') as f:
        for bc in ind_whitelist:
            f.write(bc + '\n')

with open("/projectnb2/bf528/users/saxophone/data_p4/whitelist.txt", 'w') as f:
    for bc in whitelist:
        f.write(bc + '\n')
