#!/usr/bin/env python3

import pickle
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

fig, axs = plt.subplots(4, 1, sharex=True)
pickle_dir="/projectnb2/bf528/users/saxophone/data_p4/fastq_bc/"

#barcode_sets = []
lumped = Counter()

samples = ['SRR3879604', 'SRR3879605', 'SRR3879606']
for i in range(3):
    sample = samples[i]
    ax = axs[i]
    with open(pickle_dir + sample + '.pickle', 'rb') as p:
        bc_counts = pickle.load(p)
    lumped += bc_counts
    ax.hist(np.log10(list(bc_counts.values())), density=False, cumulative=True,
            histtype='step', bins='fd')
    #ax.hist(np.log10(list(bc_counts.values())), density=True, bins=np.log10([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 20.5, 30.5, 40.5, 50.5, 100.5, 200.5, 500.5, 1000.5, 2000.5, 5000.5, 10_000.5, 100_000.5, 1_000_000.5]))


axs[3].hist(np.log10(list(lumped.values())), density=True, cumulative=True,
            histtype='step', bins='fd')
plt.show()

# count_freq = Counter()
# for count in lumped.values():
#     count_freq[count] += 1

# plt.loglog(count_freq.keys(), count_freq.values(), '.')
# plt.show()
