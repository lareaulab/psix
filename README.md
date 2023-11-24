# Psix 0.11.0

Psix is a computational tool for identifying cell-state associated alternative splicing events in single cell RNA-seq (scRNA-seq) data, described in 
[Buen Abad Najar et al, 2022](https://genome.cshlp.org/content/32/7/1385).

This branch contains a new version pf Psix currently under development. By compiling functions with [numba](https://numba.pydata.org/) and building lookup tables, Psix 0.11.0 is much faster and efficient than the original version of Psix. Other features include built-in functions to transform splice junction read counts into mRNA count estimates. In the near future we will implement built-in dimentionality reduction and visualization functions. The goal is to reduce the required four inputs of Psix to a minimum of one input (splice junction read counts). 

## Installation

Psix 0.11.0 is a Python module and it requires Python version 3.7 or higher. Psix 0.11.0 is installed directly from github using the following command:

```
pip install git+https://github.com/lareaulab/psix.git@psix_numba
```

Missing package dependencies will be automatically installed. Psix 0.11.0 is in working condition, although it is still under development and testing. Please make sure to raise an [issue](https://github.com/lareaulab/psix/issues) if you run into any problems.

## Building lookup tables

Most of the time consuming steps in Psix consist of thousands of numerical probability calculations that sometimes require hundreds of thousands of steps. In practice, the majority of these calculations are repeated most of the time. By pre-computing these probabilities and storing them in lookup tables, we can skip the vast majority of the time consuming calculations and speed Psix up significantly.

We can create a directory with lookup tables using the following command:

```
psix.psix.build_lookup(out_dir = 'lookup/') # This function will create a directory named "lookup"
```

This step usually takes just a couple of minutes, and lookup tables typically take 1.5Mb of storage.

## Creating a Psix object

In Psix 0.11.0 you can create a Psix object using a matrix of PSI values for each exon in each cell, and a corresponding matrix of splice junction read counts:

```
psix_object = psix.Psix(psi_table = 'path/to/psi.tab.gz',
                        counts_per_event = 'path/to/read_counts.tab.gz',
                        counts_type = 'reads', # Make sure to specify that these are read counts
                       )
```

The counts_per_reads table aggregates all the informative read counts for a given cassette exon. For example, if `exon_j` in `cell_i` has `5` reads covering each of the two splice junctions consistent with exon inclusion (`10` reads total), and `5` reads covering the splice junction consistent with exon exclusion, the PSI will be `(5+5)/(5*2 + 5 + 5) = 0.5` and the counts_per_event will be `5*3 = 15`.

You can also provide mRNA counts instead of read counts to the `counts_per_event` table. Just make sure to pass the argument `counts_type = 'mrna'`

You can create a Psix object from splice junction tables and a TPM table as described [here](https://github.com/lareaulab/psix#creating-a-psix-object-with-smart-seq2-data).

## Running Psix using lookup tables
You can run Psix using the following command:
```
psix_object.run_psix(latent='data/pc2_rd.renamed.tab.gz',
                    lookup = 'lookup/' # specify a pre-built lookup tables directory
                    )
```
You can specify the location of the lookup tables directory using the argument `lookup`. This way you only need to build the lookup tables once. Like in the original version of Psix, you can further increase speed by parallelizing the processes using the argument `n_jobs=5` (or any number of jobs). 

You can run Psix 0.11.0 without lookup tables by using the argument `no_lookup=True`. This will run each probability calculation anew for every exon in every cell, like in the earlier version of Psix. These functions have been improved using numba, which means that there will still be a significant improvement on speed over the original version of Psix. Using lookup tables, however, is the fastest approach. 
