# findDE

In this repository is located the R package that implements the statistical method presented in the paper [Identification of differentially expressed genes by means of outlier detection](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2318-8).

## Introduction

An important issue in microarray data is to select, from thousands of genes, a small number of informative differentially expressed (DE) genes which may be key elements for a disease. If each gene is analyzed individually, there is a big number of hypotheses to test and a multiple comparison correction method must be used. Consequently, the resulting cut-off value may be too small. Moreover, an important issue is the selection’s replicability of the DE genes. We present a new method, called ORdensity, to obtain a reproducible selection of DE genes. It takes into account the relation between all genes and it is not a gene-by-gene approach, unlike the usually applied techniques to DE gene selection.

The proposed method returns three measures, related to the concepts of outlier and density of false positives in a neighbourhood, which allow us to identify the DE genes with high classification accuracy. The first measure is a statistic called OR, introduced in \[1\]\[2\]

## References

\[1\] [Arenas C, Toma C, Cormand B, Irigoien I. Identifying extreme observations, outliers and noise in clinical and genetic data. Current Bioinformatics 2017;12(2):101–17.](http://www.eurekaselect.com/142998/article)

\[2\] [Arenas C, Irigoien I, Mestres F, Toma C, Cormand B. Extreme observations in biomedical data. In: Ainsbury EA, Calle ML, Cardis E, et al., editors. Extended Abstracts Fall 2015. Trends in Mathematics vol 7. Birkhäuser, Cham: Springer; 2017. p. 3–8.] (https://link.springer.com/chapter/10.1007/978-3-319-55639-0_1)

## Installation

To install the package from this repository, just run the following code

```
library('devtools')
install_github('jmartinezot/findDE')
```

This package requires the ```cluster``` library to be installed; otherwise it will automatically install and load it. Likewise, the ```foreach``` library is used for parallelization.

To start working with the package, just load it in the R enviroment with the following command

```
library('findDE')
```

## Example

There is a example dataframe called ```example``` shipped with the package. It contains 1000 observations of 62 variables. Each row correspond to a gene and contains 62 values: DEgen, gap and the values for the gene expression in 30 positive cases and in 30 negative cases.

First, let us extract the positive and negative cases from the ```example``` database.

```
x <- example[, 3:32]
y <- example[, 33:62]
positive <- as.matrix(x)
negative <- as.matrix(y)
```
To create an S4 object to perform the analysis, follow this syntaxis 

```
myfindDE <- new("findDE", positive = positive, negative = negative)
```
By default, no parallelizing is enabled. To enable it, just run instead

```
myfindDE <- new("findDE", positive = positive, negative = negative, parallel = TRUE)
```
There is also a ```verbose``` option, but currently is only intended for information to developers

```
myfindDE <- new("findDE", positive = positive, negative = negative, parallel = TRUE, verbose = TRUE)
```
When the object is created its statistics are computed, so afterwards it is possible to extract the information and generate several plots with the following instructions

```
summaries(myfindDE)
```
```
$summaryOR
  [1] 175.043223 172.287790 155.536259 152.547051 149.653354 148.872937 144.807897
  [8] 143.422005 134.877088 130.618427 120.829166 104.994030 103.662961  98.361779
 [15]  93.452659  93.163103  92.594912  92.201463  91.677326  88.445307  88.439498
 [22]  84.472748  84.100832  77.322660  76.090339  74.376304  73.184401  69.740607
 [29]  68.716381  65.358008  64.385061  64.092504  62.752272  59.933443  58.385493
 [36]  54.302410  52.915735  51.702722  50.880954  48.721447  48.442793  46.001939
 [43]  45.903855  45.506047  44.777559  44.274336  44.218182  43.849973  42.720397
 [50]  42.017388  41.820112  40.453281  40.134270  39.639383  38.792334  38.549919
 [57]  38.549461  35.147004  35.045910  34.617445  34.412237  54.068005  34.258620
 [64]  31.738086  31.395715  30.785716  35.751730  33.959616  38.606681  29.471714
 [71]  19.429654  24.244624  23.047096  18.742760  20.625201  24.793410  23.117921
 [78]  26.192692  22.766881  18.374832  18.328136  18.911661  16.859134  14.632071
 [85]  13.488933  14.734285  17.892081  17.597946  11.898146  10.339745  14.767943
 [92]  11.400038   9.955954  14.577235  11.009624   8.452682   9.085678   7.992158
 [99]   9.543673   8.201157  14.155428   7.385797   8.108665   8.133859   8.115440
[106]   7.659571   7.613712

$summarymeanFP
  [1] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
 [17] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
 [33] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
 [49] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.01 0.01 0.01
 [65] 0.02 0.02 0.03 0.03 0.06 0.08 0.15 0.16 0.19 0.23 0.25 0.28 0.28 0.29 0.29 0.30
 [81] 0.30 0.34 0.36 0.38 0.42 0.43 0.46 0.48 0.61 0.61 0.62 0.62 0.63 0.64 0.74 0.79
 [97] 0.82 0.83 0.85 0.88 0.91 0.98 0.99 1.00 1.00 1.00 1.00

$summarydFP
  [1] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
  [8] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
 [15] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
 [22] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
 [29] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
 [36] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
 [43] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
 [50] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
 [57] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.02250759 0.03491742
 [64] 0.03257241 0.05862612 0.07784708 0.09040119 0.09969090 0.16041671 0.23734397
 [71] 0.52936672 0.46940196 0.64966166 0.88534314 0.80602538 0.74330105 0.84999062
 [78] 0.93827338 0.89804348 0.67547883 1.17118130 1.22784138 1.43636830 1.36470675
 [85] 1.65248914 1.67501787 1.36187497 1.79681251 2.38536830 4.00322958 2.69880220
 [92] 2.73011412 3.55633124 2.86300369 4.40165453 5.81985305 5.54480306 6.84330731
 [99] 5.76309784 7.52559437 3.21192579 4.29292572 5.70884616 3.84099601 9.13942664
[106] 5.33755737 8.53078041
```

```
plotFPvsOR(myfindDE)
```
```
plotclusters(myfindDE)
```
```
clusplotk(myfindDE, 2)
```
```
clusplotk(myfindDE, 4)
```
