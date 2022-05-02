#!/usr/bin/env python3

import pickle
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

pickle_dir="/projectnb2/bf528/users/saxophone/data_p4/fastq_bc/"
figure_dir="/projectnb2/bf528/users/saxophone/data_p4/figures/"

bc_counts = []

samples = ['SRR3879604', 'SRR3879605', 'SRR3879606']
for i in range(3):
    sample = samples[i]
    with open(pickle_dir + sample + '.pickle', 'rb') as p:
        bc_counts.append(pickle.load(p))
    (n, bins, patches) = plt.hist(np.log10(list(bc_counts[i].values())),
                                  cumulative=True, density=True,
                                  histtype='step', bins='fd', label=sample)
    patches[0].set_xy(patches[0].get_xy()[:-1])

plt.xlabel('log 10 barcode count')
plt.ylabel('Cumulative distribution')
plt.legend(loc='lower right')
plt.savefig(figure_dir + 'cum_dist.png', dpi=300)
plt.show()

for i in range(3):
    (n, bins, patches) = plt.hist(np.log10(list(bc_counts[i].values())),
                                  cumulative=-1, density=False,
                                  histtype='step', bins='fd', label=sample)

plt.ylim(0,2000)
plt.xlim(left=4)
plt.xlabel('Filter cutoff (log 10 count)')
plt.ylabel('# Barcodes included')
plt.legend(loc='upper right')
plt.savefig(figure_dir + 'filter_sel.png', dpi=300)
plt.show()
