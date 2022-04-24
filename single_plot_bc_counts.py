#!/usr/bin/env python3

import pickle
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

pickle_dir="/projectnb2/bf528/users/saxophone/data_p4/fastq_bc/"

#barcode_sets = []
bc_counts = []
lumped = Counter()

samples = ['SRR3879604', 'SRR3879605', 'SRR3879606']
for i in range(3):
    sample = samples[i]
    with open(pickle_dir + sample + '.pickle', 'rb') as p:
        bc_counts.append(pickle.load(p))
    lumped += bc_counts[i]
    (n, bins, patches) = plt.hist(np.log10(list(bc_counts.values())),
                                  cumulative=-1, density=False,
                                  histtype='step', bins='fd', label=sample)
    #print(f'{n=}, {bins=}, {patches=}')
    #plt.plot(bins[:-1], n, label='Test redraw')
    #break

plt.ylim(0,2000)
plt.xlim(left=4)
plt.xlabel('Filter cutoff')
plt.ylabel('# Barcodes included')
plt.legend(loc='best')
plt.show()

# count_freq = Counter()
# for count in lumped.values():
#     count_freq[count] += 1

# plt.loglog(count_freq.keys(), count_freq.values(), '.')
# plt.show()
