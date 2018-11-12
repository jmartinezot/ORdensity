# findDE

In this repository is located the R package that implements the statistical method presented in the paper [Identification of differentially expressed genes by means of outlier detection](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2318-8).

## Introduction

An important issue in microarray data is to select, from thousands of genes, a small number of informative differentially expressed (DE) genes which may be key elements for a disease. If each gene is analyzed individually, there is a big number of hypotheses to test and a multiple comparison correction method must be used. Consequently, the resulting cut-off value may be too small. Moreover, an important issue is the selection’s replicability of the DE genes. We present a new method, called ORdensity, to obtain a reproducible selection of DE genes. It takes into account the relation between all genes and it is not a gene-by-gene approach, unlike the usually applied techniques to DE gene selection.

The proposed method returns three measures, related to the concepts of outlier and density of false positives in a neighbourhood, which allow us to identify the DE genes with high classification accuracy. To assess the performance of ORdensity, we used simulated microarray data and four real microarray cancer data sets. The results indicated that the method correctly detects the DE genes; it is competitive with other well accepted methods; the list of DE genes that it obtains is useful for the correct classification or diagnosis of new future samples and, in general, it is more stable than other procedures.

## Installation

To install the package from this repository, just run the following code

```
library('devtools')
install_github('jmartinezot/findDE')
```

This package requires the ```cluster``` library to be installed; otherwise it will automatically install and load it.

To start working with it, just load it in the R enviroment with the following command

```
library('findDE')
```

## Example

There is a example dataframe called ```example``` shipped with the package.
```
x <- example[, 3:32]
y <- example[, 33:62]
positive <- as.matrix(x)
negative <- as.matrix(y)
myfindDE <- new("findDE", positive = positive, negative = negative)
```
```
summaries(myfindDE)
plotFPvsOR(myfindDE)
plotclusters(myfindDE)
clusplotk(myfindDE, 2)
clusplotk(myfindDE, 4)
```
