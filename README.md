# Psix

Psix is a computational tool for identifying cassette exons with informative biological variation in a single-cell dataset.

Inspired in auto-correlation approaches, Psix will tell you if an exon's splicing is significantly associated with a cell-metric that shows the relationships between single cells. In practice, this map could be a low-dimensional representation of the gene expression of a single cell population.

A characteristic feature of Psix is that it is aware of [coverage dependent biases](https://elifesciences.org/articles/54603) that add unwanted technical variation to splicing observations in single cells. As a result, it's precision is superior to other approaches that might be mislead by such effects.
