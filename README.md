# Psix

Psix is a computational tool for identifying cassette exons with informative biological variation in a single-cell dataset.

Inspired in auto-correlation approaches, Psix will tell you if an exon's splicing is significantly associated with a cell-metric that shows the relationships between single cells. In practice, this map could be a low-dimensional representation of the gene expression of a single cell population. 

[Coverage dependent biases](https://elifesciences.org/articles/54603) add unwanted technical variation to splicing observations in single cells. Psix uses a probabilistic approach to determined if the observed variance of an exon accross a phenptypic landscape is biological, or if it's the result of these biases.

ADD LINKS TO EXAMPLES, INCLUDING FOR THE PAPER, HERE

## Installation

Psix is a Python module and it requires Python version 3.6 or higher. Psix is installed directly from github using the following command:

```
pip install git+https://github.com/lareaulab/psix.git
```

Missing package dependencies will be automatically installed.

## Getting started

This section was written with smart-seq2 data in mind, but running Psix on UMI data is not much different. For the specifics on running Psix on UMI-based scRNA-seq data, go to **Running Psix in UMI data**.

Psix requires three inputs from the user:
* A directory containing SJ.out.tab files from STAR. This is used to calculate the exon's observed $\hat{\Psi}$.
* A matrix of gene expression in transcripts per million (TPM; for smart-seq2 data only). This is used to estimate the captured mRNA molecules per observation.
* A low-dimensional cell space; e.g., a PCA projection of normalized gene expression. This is used to define cell neighborhoods.

### Mapping and preprocessing scRNA-seq data

#### SJ.out.tab files

We recommend mapping raw scRNA-seq reads using STAR version $\geq 2.5.3a$. Psix uses the ```SJ.out.tab``` files from the STAR aligner. Individual files from each single cell should be stored in the same directory with the following naming format: ```cellID.SJ.out.tab```. The files can be gzipped or uncompressed. If you are using STARsolo, go to **Running Psix with STARsolo**.

#### TPM matrix for smart-seq2 only

The TPM matrix of **gene** expression can be obtained running different methods. We use RSEM version $\geq 1.2.31$ because it can be run using STAR as the aligner. Other methods such as Kallisto can also be used to generate this matrix. This is only required for smart-seq2 data. Go to **Running Psix in UMI data** to see .

#### Low-dimensional cell space

A low-dimensional cell space is provided by the user, since Psix does not perform dimensionality reduction. In principle Psix can run with any metric. However, we recommend using interpretable dimensionality reduction methods such as PCA over the normalized gene expression, while avoiding non-interpretable methods such as tSNE except for visualization. 

For small smart-seq2 datasets (fewer than 5000 cells), we recommend using SCONE to select the best normalization method before applying a linear dimensionality reduction such as PCA. Alternatively, other methods such as ZINB-Wave can be used on this data. For larger datasets, we recommend using the latent space of scVI as the low-dimensional cell space.

### Running Psix

You can import psix to python simply by running:

```import psix```
