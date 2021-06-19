# Psix

Psix is a computational tool for identifying cassette exons with informative biological variation in a single-cell dataset.

Inspired in auto-correlation approaches, Psix will tell you if an exon's splicing is significantly associated with a cell-metric that shows the relationships between single cells. In practice, this map could be a low-dimensional representation of the gene expression of a single cell population. 

[Coverage dependent biases](https://elifesciences.org/articles/54603) add unwanted technical variation to splicing observations in single cells. Psix uses a probabilistic approach to determined if the observed variance of an exon accross a phenptypic landscape is biological, or if it's the result of these biases.

ADD LINKS TO EXAMPLES HERE

## Installation

Psix requires Python version 3.6 or higher. Psix is installed directly from github using the following command:

```
pip install git+https://github.com/lareaulab/psix.git
```

Missing package dependencies will be automatically installed.

## Getting started

### Mapping and preprocessing scRNA-seq data

We recommend mapping raw scRNA-seq reads using STAR version $\geq 2.5.2$. Psix uses the ```SJ.out.tab``` files from the STAR aligner. Individual files from each single cell should be stored in the same directory with the following naming format: ```cellID.SJ.out.tab```. The files can be gzipped or uncompressed. If you are using STARsolo, go to the **Running Psix with STARsolo section**.

