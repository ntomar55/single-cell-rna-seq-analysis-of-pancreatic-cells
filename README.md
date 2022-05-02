# Project Description

We attempt to replicate the results of Baron et al., processing single-cell
RNA-seq data to identify commmon and uncommon cell types in the human
pancreas.

# Contributors

Daniel Gealow, Nikita Tomar, and David Lenci

# Repository Contents

`Project 4 Report.pdf`: Our final report.

`count_barcodes.qsub`: A shell script (to be submitted to the qsub queue)
that runs three instances of zcat in parallel, each gradually piping the
contents of one of the sample barcode fastq files to an instance of 
`count_barcodes.py`.

`count_barcodes.py`: Recieves a fastq file as piped input and counts the
number of occure barcodes (first 19 bp in the sequence) into a `Counter()`
dictionary, which is then saved to a pickle file.

`plot_bc_counts.py`: Plots the distribution of barcode counts in each
pickle file in two figures to help determine an appropriate filter
cutoff.

`create_whitelist.py`: Determines the set of barcodes that appear
at least 10^4.5 times in any of the pickled counters, and writes them
to whitelist files. (The combined `whitelist.txt` file is the one that we
actually use in our further analysis).

`generate_index.qsub`: Runs `salmon index` to generate an index from
the gencode v40 human reference transcriptome.

`run_alevin.qsub`: Runs `salmon alevin` to generate the UMI count matrix.
Requires the barcode and read 2 files for each of the three samples,
the `whitelist.txt` file, a transcript-to-gene mapping file (`t2g_map.tsv`),
and the index created by `generate_index.qsub`.

`Programmer.R`: Processes the avelin data, filters out low quality genes,
reduces dimension and performs clustering on them.

`analyst_main.R`: Contains the code for identifying potential cell markers,
labeling clusters, and then generating the clustered heatmap.
