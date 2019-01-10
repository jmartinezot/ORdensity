# findDE

In this repository is located the R package that implements the statistical method presented in the paper [Identification of differentially expressed genes by means of outlier detection](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2318-8).

## Introduction

An important issue in microarray data is to select, from thousands of genes, a small number of informative differentially expressed (DE) genes which may be key elements for a disease. If each gene is analyzed individually, there is a big number of hypotheses to test and a multiple comparison correction method must be used. Consequently, the resulting cut-off value may be too small. Moreover, an important issue is the selection’s replicability of the DE genes. We present a new method, called ORdensity, to obtain a reproducible selection of DE genes. It takes into account the relation between all genes and it is not a gene-by-gene approach, unlike the usually applied techniques to DE gene selection.

The proposed method returns three measures, related to the concepts of outlier and density of false positives in a neighbourhood, which allow us to identify the DE genes with high classification accuracy. The first measure is a statistic called OR, introduced in \[1\]\[2\]

## References

\[1\] Arenas C, Toma C, Cormand B, Irigoien I. [Identifying extreme observations, outliers and noise in clinical and genetic data.](http://www.eurekaselect.com/142998/article) Current Bioinformatics 2017;12(2):101–17.

\[2\] Arenas C, Irigoien I, Mestres F, Toma C, Cormand B. [Extreme observations in biomedical data.](https://link.springer.com/chapter/10.1007/978-3-319-55639-0_1) In: Ainsbury EA, Calle ML, Cardis E, et al., editors. Extended Abstracts Fall 2015. Trends in Mathematics vol 7. Birkhäuser, Cham: Springer; 2017. p. 3–8.

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

There is a example dataframe called ```example``` shipped with the package. This data is the result of a simulation of 100 differentially expressed genes in a pool of 1000 genes. It contains 1000 observations of 62 variables. Each row correspond to a gene and contains 62 values: DEgen, gap and the values for the gene expression in 30 positive cases and in 30 negative cases. The DEgen field value is 1 for differentially expressed genes and 0 for those which are not.

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

If the researcher just wants to extract the differentially expressed genes detected by the ORdensity method, a call to findDEgenes will return a list with the clusters found, along with their mean value of the OR statistic. Higher OR values mean higher probability of true differentially expressed.

For example, after running this code

```
result <- findDEgenes(myfindDE)
```

and being told that the optimal clustering consists of just two clusters, 

```
The ORdensity method has found that the optimal clustering of the data consists of 2 clusters
```

we could then take a look to the genes corresponding to the generated clusters

```
> result$DFgenes
[[1]]
[[1]]$cluster_number
[1] 1

[[1]]$genes
 [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  20  21  22  23  24  25  26  27  28  29  30  32  34  35  36  37  39  40 42  46  47  48  49  50  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  76  79  80  81  82  83 84  86  87  88  89  90  91  92  93  95  96  97  98  99 100

[[1]]$meanOR
[1] 62.46812


[[2]]
[[2]]$cluster_number
[1] 2

[[2]]$genes
 [1]  18  19  31  33  38  41  43  44  45  51  75  78  85  94 104 277 399 598 618 651 670 946

[[2]]$meanOR
[1] 10.64626
```
The clusters are ordered in decreasing order according to the value of the mean of the OR statistic. We see that the mean is higher in the first cluster (62.46812) than in the second one (10.64626), which means that the first cluster is more likely composed of true differentially expressed genes, and the second one to be composed of false positives. With more clusters, the last ones are likely false negatives.

We could also check the values associated to the individual genes in the clusters

```
result$clusters
[[1]]
       id        OR     DifExp minFP meanFP maxFP    density    radius
 [1,]  62 175.04322 -0.8237232   0.0   0.00   0.0 0.00000000 1.1259372
 [2,]  50 172.28779 -0.8237232   0.0   0.00   0.0 0.00000000 1.0924395
 [3,]  61 155.53626 -0.8237232   0.0   0.00   0.0 0.00000000 0.8990361
 [4,]  70 152.54705 -0.8237232   0.0   0.00   0.0 0.00000000 0.8642152
 [5,]   2 149.65335 -0.8237232   0.0   0.00   0.0 0.00000000 0.8167136
 [6,]  68 148.87294 -0.8237232   0.0   0.00   0.0 0.00000000 0.7960463
 [7,]  52 134.87709 -0.8237232   0.0   0.00   0.0 0.00000000 0.6182741
 [8,]  34  77.32266 -0.8237232   0.0   0.00   0.0 0.00000000 0.3456795
 [9,]  28  62.75227 -0.8237232   0.0   0.00   0.0 0.00000000 0.4647984
[10,]  40 104.99403 -0.8137232   0.0   0.01   0.1 0.01991604 0.5371505
[11,]  82  93.45266 -0.8137232   0.0   0.01   0.1 0.02892853 0.4128200
[12,]  15  93.16310 -0.8137232   0.0   0.01   0.1 0.02956274 0.4129561
[13,]  64  92.59491 -0.8137232   0.0   0.01   0.1 0.02606554 0.4483877
[14,]  73  91.67733 -0.8137232   0.0   0.01   0.1 0.02812384 0.4228071
[15,]  67  88.44531 -0.8137232   0.0   0.01   0.1 0.03110092 0.3420220
[16,]  29  88.43950 -0.8137232   0.0   0.01   0.1 0.02927928 0.3963377
[17,]  81  54.30241 -0.8137232   0.0   0.01   0.1 0.02837142 0.3602596
[18,]  25  51.70272 -0.8137232   0.0   0.01   0.1 0.03502176 0.2983137
[19,]  90  42.72040 -0.8137232   0.0   0.01   0.1 0.04340459 0.2356958
[20,]   9  42.01739 -0.8137232   0.0   0.01   0.1 0.04779317 0.2128020
[21,]  83  39.63938 -0.8137232   0.0   0.01   0.1 0.04505324 0.2469755
[22,]  76  38.54992 -0.8137232   0.0   0.01   0.1 0.04492580 0.2384937
[23,]  14  38.54946 -0.8137232   0.0   0.01   0.1 0.04340459 0.2312692
[24,]  53  69.74061 -0.8037232   0.0   0.02   0.1 0.04959866 0.4111845
[25,]  92  45.90385 -0.8037232   0.0   0.02   0.1 0.04838004 0.4275156
[26,]  42  44.77756 -0.8037232   0.0   0.02   0.1 0.09793269 0.2112018
[27,]  39  44.21818 -0.8037232   0.0   0.02   0.1 0.06279485 0.3293868
[28,]  16  43.84997 -0.8037232   0.0   0.02   0.1 0.09271376 0.2449293
[29,]   6  84.47275 -0.7937232   0.0   0.03   0.1 0.08493890 0.3846254
[30,]  23  76.09034 -0.7937232   0.0   0.03   0.1 0.06951063 0.4469373
[31,]  71  74.37630 -0.7937232   0.0   0.03   0.1 0.09291337 0.3240597
[32,]   3  46.00194 -0.7937232   0.0   0.03   0.1 0.12042185 0.2728672
[33,]   1  40.13427 -0.7937232   0.0   0.03   0.1 0.11677728 0.2616333
[34,]   4  38.79233 -0.7937232   0.0   0.03   0.2 0.11428770 0.2821381
[35,]  26  58.38549 -0.7837232   0.0   0.04   0.1 0.12398654 0.3442278
[36,]  21  50.88095 -0.7837232   0.0   0.04   0.2 0.12891887 0.3228514
[37,]  99  48.72145 -0.7837232   0.0   0.04   0.2 0.13486838 0.3370886
[38,]  97  48.44279 -0.7837232   0.0   0.04   0.2 0.11703689 0.3451836
[39,]  58  45.50605 -0.7837232   0.0   0.04   0.1 0.14878155 0.2813442
[40,]  87  44.27434 -0.7837232   0.0   0.04   0.1 0.16164926 0.2795572
[41,]  48  41.82011 -0.7837232   0.0   0.04   0.1 0.16978212 0.2589286
[42,]  22  40.45328 -0.7837232   0.0   0.04   0.3 0.16928717 0.2458142
[43,]  36 120.82917 -0.7737232   0.0   0.05   0.2 0.06519188 0.8265846
[44,]  93 103.66296 -0.7737232   0.0   0.05   0.2 0.08567963 0.5959258
[45,]  11  68.71638 -0.7737232   0.0   0.05   0.1 0.12380840 0.4120751
[46,]  89  64.38506 -0.7737232   0.0   0.05   0.1 0.11819244 0.4255925
[47,]  91  52.91573 -0.7737232   0.0   0.05   0.2 0.16534165 0.3091681
[48,]   7 143.42201 -0.7637232   0.0   0.06   0.2 0.05724804 1.0933673
[49,]  10 130.61843 -0.7637232   0.0   0.06   0.2 0.06702330 0.9488431
[50,]  17  84.10083 -0.7637232   0.0   0.06   0.2 0.11512295 0.5528341
[51,]   5  73.18440 -0.7637232   0.0   0.06   0.2 0.12805340 0.4756794
[52,]  37  64.09250 -0.7637232   0.0   0.06   0.2 0.14079606 0.4318049
[53,]  88  35.14700 -0.7637232   0.0   0.06   0.1 0.26485308 0.2356897
[54,]  80  34.41224 -0.7637232   0.0   0.06   0.2 0.20050219 0.3058787
[55,]  32 144.80790 -0.7537232   0.0   0.07   0.2 0.06450527 1.1129392
[56,]  65  98.36178 -0.7537232   0.0   0.07   0.2 0.10961285 0.6438543
[57,]  63  65.35801 -0.7537232   0.0   0.07   0.2 0.16773789 0.4298988
[58,]   8  35.75173 -0.7537232   0.0   0.07   0.2 0.20765235 0.3446413
[59,]  79  35.04591 -0.7537232   0.0   0.07   0.2 0.24553401 0.2939147
[60,]  24  92.20146 -0.7437232   0.0   0.08   0.2 0.12637911 0.6530435
[61,]  54  59.93344 -0.7437232   0.0   0.08   0.2 0.23161371 0.3659407
[62,]  56  34.25862 -0.7437232   0.0   0.08   0.2 0.27959278 0.2866394
[63,] 100  30.78572 -0.7437232   0.0   0.08   0.2 0.32405298 0.2544298
[64,]  47  38.60668 -0.7337232   0.0   0.09   0.2 0.24010010 0.3758630
[65,]  35  34.61744 -0.7337232   0.0   0.09   0.2 0.30338254 0.2994968
[66,]  84  31.73809 -0.7337232   0.0   0.09   0.2 0.29753366 0.3054070
[67,]  86  31.39571 -0.7237232   0.0   0.10   0.3 0.30289156 0.3380744
[68,]  72  54.06801 -0.7137232   0.0   0.11   0.3 0.25201319 0.4447310
[69,]  13  33.95962 -0.7037232   0.0   0.12   0.2 0.39522229 0.3078496
[70,]  20  29.47171 -0.7037232   0.0   0.12   0.2 0.36064999 0.3376266
[71,]  95  19.42965 -0.6637232   0.1   0.16   0.3 0.56731613 0.2828725
[72,]  66  24.24462 -0.6337232   0.0   0.19   0.4 0.66112753 0.3234741
[73,]  59  18.74276 -0.5737232   0.2   0.25   0.4 0.99452059 0.2602576
[74,]  98  18.32814 -0.5737232   0.1   0.25   0.4 1.00183786 0.2584957
[75,]  60  23.04710 -0.5637232   0.1   0.26   0.4 0.98442393 0.2805144
[76,]  55  26.19269 -0.5237232   0.1   0.30   0.6 0.97959618 0.3133187
[77,]  27  24.79341 -0.5137232   0.2   0.31   0.4 0.86305155 0.3698262
[78,]  12  20.62520 -0.5137232   0.2   0.31   0.5 1.02266734 0.3066674
[79,]  69  22.76688 -0.5037232   0.2   0.32   0.4 1.01699885 0.3229340
[80,]  49  23.11792 -0.4837232   0.2   0.34   0.5 1.04442334 0.3276343
[81,]  96  16.85913 -0.4537232   0.3   0.37   0.4 1.45362576 0.2590543
[82,]  57  18.91166 -0.4437232   0.2   0.38   0.5 1.43391223 0.2730387
[83,]  30  18.37483 -0.4337232   0.3   0.39   0.5 0.91488614 0.4300349
[84,]  46  14.63207 -0.3937232   0.3   0.43   0.6 1.59270035 0.2732312
[85,]  74  17.89208 -0.3537232   0.4   0.47   0.6 1.42662586 0.3349624

[[2]]
       id        OR       DifExp minFP meanFP maxFP  density    radius
 [1,]  31 17.597946 -0.303723229   0.4   0.52   0.6 2.148361 0.2501579
 [2,]  85 13.488933 -0.293723229   0.4   0.53   0.7 2.220875 0.2417213
 [3,]  45 14.734285 -0.283723229   0.4   0.54   0.7 2.356679 0.2337810
 [4,]  33 11.898146 -0.193723229   0.3   0.63   0.7 2.493585 0.2560273
 [5,]  78 10.339745 -0.193723229   0.5   0.63   0.7 4.131262 0.1540423
 [6,]  43 14.577235 -0.173723229   0.4   0.65   0.8 3.002335 0.2194208
 [7,]  75 14.767943 -0.153723229   0.6   0.67   0.8 3.309754 0.2070907
 [8,]  18  9.955954 -0.153723229   0.6   0.67   0.7 4.037818 0.1670731
 [9,]  19 11.400038 -0.123723229   0.7   0.70   0.7 3.408621 0.2076696
[10,]  51 11.009624 -0.093723229   0.6   0.73   0.8 4.575104 0.1629221
[11,]  94  8.452682 -0.023723229   0.8   0.80   0.8 6.077158 0.1320430
[12,] 399  7.992158 -0.003723229   0.8   0.82   0.9 6.843839 0.1203501
[13,] 104  9.543673  0.016276771   0.7   0.84   0.9 6.096918 0.1430609
[14,]  44  9.085678  0.036276771   0.8   0.86   0.9 6.395257 0.1364570
[15,]  38  8.201157  0.046276771   0.8   0.87   0.9 7.140990 0.1237924
[16,]  41 14.155428  0.116276771   0.9   0.94   1.0 3.502890 0.2737847
[17,] 598  7.385797  0.166276771   0.9   0.99   1.0 4.342387 0.2294194
[18,] 277  8.133859  0.176276771   1.0   1.00   1.0 3.576610 0.2821911
[19,] 618  8.115440  0.176276771   1.0   1.00   1.0 9.172236 0.1107120
[20,] 946  8.108665  0.176276771   1.0   1.00   1.0 5.777176 0.1742042
[21,] 651  7.659571  0.176276771   1.0   1.00   1.0 5.163137 0.1952410
[22,] 670  7.613712  0.176276771   1.0   1.00   1.0 8.051950 0.1256902
```

As a rule of thumb, differentially expressed genes are expected to present high values of OR and low values of meanFP and density. We could also analyze each gene individually inside each cluster. The motivation of the clustering is to distinguish those false positives that score high in OR and low in meanFP and density, but are similar to other known false positives obtained by boostrapping. The procedure is detailed in the paper referenced above. 

If the researcher is interested in a more thorough analysis, other functions are at their service.

A summary of OR, meanFP and dFP (density) can be obtained.

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
A plot with a representation of the potential genes based on OR (vertical axis), FP (horizontal axis) and dFP (size of the circle is inversely proportional to its value) can also be obtained. The plot is similar to Fig.3b in the paper.

```
plotFPvsOR(myfindDE)
```
The plot of k values against the silhouette measure is also provided.

```
plotclusters(myfindDE)
```
It is also possible to see a graphic representation of the clustering projected onto the first two principal components

```
clusplotk(myfindDE, k = 2)
```
```
clusplotk(myfindDE, k = 4)
```
