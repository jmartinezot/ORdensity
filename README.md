# ORdensity

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
install_github('jmartinezot/ORdensity')
```

This package requires the ```cluster``` library to be installed; otherwise it will automatically install and load it. Likewise, the ```foreach``` library is used for parallelization.

To start working with the package, just load it in the R enviroment with the following command

```
library('ORdensity')
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
myORdensity <- new("ORdensity", positive = positive, negative = negative)
```
By default, no parallelizing is enabled. To enable it, just run instead

```
myORdensity <- new("ORdensity", positive = positive, negative = negative, parallel = TRUE)
```
It is also possible to enable or disable replicability, and to pass the seed to the pseudorandom number generator. The default values are 

```
myORdensity <- new("ORdensity", positive = positive, negative = negative, replicable = TRUE, seed = 0)
```
with the function using the given seed to set the random generator. If replicable = FALSE, no seed is used.

If the researcher just wants to extract the differentially expressed genes detected by the ORdensity method, a call to findDEgenes will return a list with the clusters found, along with their mean value of the OR statistic. Higher OR values mean higher probability of true differentially expressed.

For example, after running this code

```
result <- findDEgenes(myORdensity)
```

and being told that the optimal clustering consists of just two clusters, 

```
The ORdensity method has found that the optimal clustering of the data consists of 2 clusters
```

we could then take a look to the genes corresponding to the generated clusters

```
> result
[[1]]
[[1]]$cluster_number
[1] 1

[[1]]$genes
 [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  20  21  22  23  24  25  26  27  28  29  30  32  34  35  36  37  39  40  42  46  47  48  49  50  52  53  54  55  56
[47]  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  76  79  80  81  82  83  84  86  87  88  89  90  91  92  93  95  96  97  98  99 100

[[1]]$meanOR
[1] 62.99879


[[2]]
[[2]]$cluster_number
[1] 2

[[2]]$genes
 [1]  18  19  31  33  38  41  43  44  45  51  74  75  78  85  94 104 277 399 598 618 651 670 946

[[2]]$meanOR
[1] 10.96129

```
The clusters are ordered in decreasing order according to the value of the mean of the OR statistic. We see that the mean is higher in the first cluster (62.99879) than in the second one (10.96129), which means that the first cluster is more likely composed of true differentially expressed genes, and the second one to be composed of false positives. With more clusters, the last ones are likely false negatives.

We could also check a more detailed summary of the object.

```
summary(myORdensity)
```
The output would be the following

```
This is the proposed clustering made by the ORdensity method
For the computation of FP and dFP a total of 10 neighbours have been taken into account
The expected number of false positives neighbours is 8.237232 
The ORdensity method has found that the optimal clustering of the data consists of 2 clusters

$neighbours
[1] 10

$expectedFalsePositiveNeighbours
[1] 8.237232

$clusters
$clusters[[1]]
       id        OR  FP        dFP
 [1,]  62 175.04322 0.0  0.0000000
 [2,]  50 172.28779 0.0  0.0000000
 [3,]  61 155.53626 0.0  0.0000000
 [4,]  70 152.54705 0.0  0.0000000
 [5,]   2 149.65335 0.0  0.0000000
 [6,]  68 148.87294 0.0  0.0000000
 [7,]  32 144.80790 0.0  0.0000000
 [8,]   7 143.42201 0.0  0.0000000
 [9,]  52 134.87709 0.0  0.0000000
[10,]  10 130.61843 0.0  0.0000000
[11,]  36 120.82917 0.0  0.0000000
[12,]  40 104.99403 0.0  0.0000000
[13,]  93 103.66296 0.0  0.0000000
[14,]  65  98.36178 0.0  0.0000000
[15,]  82  93.45266 0.0  0.0000000
[16,]  15  93.16310 0.0  0.0000000
[17,]  64  92.59491 0.0  0.0000000
[18,]  24  92.20146 0.0  0.0000000
[19,]  73  91.67733 0.0  0.0000000
[20,]  67  88.44531 0.0  0.0000000
[21,]  29  88.43950 0.0  0.0000000
[22,]   6  84.47275 0.0  0.0000000
[23,]  17  84.10083 0.0  0.0000000
[24,]  34  77.32266 0.0  0.0000000
[25,]  23  76.09034 0.0  0.0000000
[26,]  71  74.37630 0.0  0.0000000
[27,]   5  73.18440 0.0  0.0000000
[28,]  53  69.74061 0.0  0.0000000
[29,]  11  68.71638 0.0  0.0000000
[30,]  89  64.38506 0.0  0.0000000
[31,]  37  64.09250 0.0  0.0000000
[32,]  28  62.75227 0.0  0.0000000
[33,]  26  58.38549 0.0  0.0000000
[34,]  81  54.30241 0.0  0.0000000
[35,]  91  52.91573 0.0  0.0000000
[36,]  25  51.70272 0.0  0.0000000
[37,]  21  50.88095 0.0  0.0000000
[38,]  99  48.72145 0.0  0.0000000
[39,]  97  48.44279 0.0  0.0000000
[40,]   3  46.00194 0.0  0.0000000
[41,]  92  45.90385 0.0  0.0000000
[42,]  58  45.50605 0.0  0.0000000
[43,]  42  44.77756 0.0  0.0000000
[44,]  87  44.27434 0.0  0.0000000
[45,]  39  44.21818 0.0  0.0000000
[46,]  16  43.84997 0.0  0.0000000
[47,]  90  42.72040 0.0  0.0000000
[48,]   9  42.01739 0.0  0.0000000
[49,]  48  41.82011 0.0  0.0000000
[50,]  22  40.45328 0.0  0.0000000
[51,]   1  40.13427 0.0  0.0000000
[52,]  83  39.63938 0.0  0.0000000
[53,]  76  38.54992 0.0  0.0000000
[54,]  14  38.54946 0.0  0.0000000
[55,]   8  35.75173 0.0  0.0000000
[56,]  80  34.41224 0.0  0.0000000
[57,]  13  33.95962 0.0  0.0000000
[58,]  63  65.35801 0.1  0.2393753
[59,]  54  59.93344 0.1  0.2820190
[60,]  47  38.60668 0.1  0.2667662
[61,]  79  35.04591 0.1  0.3444442
[62,]  86  31.39571 0.1  0.2931306
[63,]  72  54.06801 0.2  0.4522637
[64,]   4  38.79233 0.2  0.7073464
[65,]  35  34.61744 0.2  0.6733581
[66,]  88  35.14700 0.3  1.4112193
[67,] 100  30.78572 0.3  1.1781804
[68,]  56  34.25862 0.4  1.3966969
[69,]  84  31.73809 0.6  1.9543445
[70,]  20  29.47171 1.1  3.2552493
[71,]  66  24.24462 1.4  4.0342057
[72,]  95  19.42965 1.6  5.6854818
[73,]  60  23.04710 1.8  6.0885440
[74,]  12  20.62520 2.2  6.9215094
[75,]  59  18.74276 2.5 10.1692216
[76,]  27  24.79341 2.6  6.9045174
[77,]  69  22.76688 2.7  8.3615886
[78,]  98  18.32814 2.8 11.2156424
[79,]  49  23.11792 3.0  9.1467849
[80,]  55  26.19269 3.3 10.7893489
[81,]  30  18.37483 3.5  7.9978672
[82,]  57  18.91166 3.6 12.7433027
[83,]  96  16.85913 3.8 14.9634715
[84,]  46  14.63207 3.8 13.7055582

$clusters[[2]]
       id        OR   FP      dFP
 [1,]  74 17.892081  4.9 14.84705
 [2,]  31 17.597946  5.1 20.46774
 [3,]  45 14.734285  5.2 22.29119
 [4,]  85 13.488933  5.2 21.39663
 [5,]  33 11.898146  5.9 22.35787
 [6,]  78 10.339745  6.0 38.41902
 [7,]  43 14.577235  6.3 28.06259
 [8,]  19 11.400038  6.6 31.42799
 [9,]  18  9.955954  6.6 38.49118
[10,]  75 14.767943  6.9 33.21473
[11,]  51 11.009624  7.4 45.27537
[12,]  94  8.452682  8.0 60.76025
[13,]  44  9.085678  8.3 57.26337
[14,] 399  7.992158  8.3 68.18916
[15,] 104  9.543673  8.5 58.97102
[16,]  38  8.201157  8.9 76.21379
[17,]  41 14.155428  9.5 33.61409
[18,] 598  7.385797  9.8 45.59191
[19,] 277  8.133859 10.0 37.81764
[20,] 618  8.115440 10.0 85.55168
[21,] 946  8.108665 10.0 58.13935
[22,] 651  7.659571 10.0 56.17629
[23,] 670  7.613712 10.0 84.82799

```


As a rule of thumb, differentially expressed genes are expected to present high values of OR and low values of meanFP and density. We could also analyze each gene individually inside each cluster. The motivation of the clustering is to distinguish those false positives that score high in OR and low in meanFP and density, but are similar to other known false positives obtained by boostrapping. The procedure is detailed in the paper referenced above. 

If the researcher is interested in a more thorough analysis, other functions are at their service.

The data before being clustered can be obtained with the following function

```
preclusteredData(myORdensity)
Column S denotes the cases when FP=0
Column F denotes the cases when FP < expectedFalsePositives
     id         OR   FP        dFP S F
1    62 175.043223  0.0  0.0000000 S F
2    50 172.287790  0.0  0.0000000 S F
3    61 155.536259  0.0  0.0000000 S F
4    70 152.547051  0.0  0.0000000 S F
5     2 149.653354  0.0  0.0000000 S F
6    68 148.872937  0.0  0.0000000 S F
7    32 144.807897  0.0  0.0000000 S F
8     7 143.422005  0.0  0.0000000 S F
9    52 134.877088  0.0  0.0000000 S F
10   10 130.618427  0.0  0.0000000 S F
11   36 120.829166  0.0  0.0000000 S F
12   40 104.994030  0.0  0.0000000 S F
13   93 103.662961  0.0  0.0000000 S F
14   65  98.361779  0.0  0.0000000 S F
15   82  93.452659  0.0  0.0000000 S F
16   15  93.163103  0.0  0.0000000 S F
17   64  92.594912  0.0  0.0000000 S F
18   24  92.201463  0.0  0.0000000 S F
19   73  91.677326  0.0  0.0000000 S F
20   67  88.445307  0.0  0.0000000 S F
21   29  88.439498  0.0  0.0000000 S F
22    6  84.472748  0.0  0.0000000 S F
23   17  84.100832  0.0  0.0000000 S F
24   34  77.322660  0.0  0.0000000 S F
25   23  76.090339  0.0  0.0000000 S F
26   71  74.376304  0.0  0.0000000 S F
27    5  73.184401  0.0  0.0000000 S F
28   53  69.740607  0.0  0.0000000 S F
29   11  68.716381  0.0  0.0000000 S F
30   89  64.385061  0.0  0.0000000 S F
31   37  64.092504  0.0  0.0000000 S F
32   28  62.752272  0.0  0.0000000 S F
33   26  58.385493  0.0  0.0000000 S F
34   81  54.302410  0.0  0.0000000 S F
35   91  52.915735  0.0  0.0000000 S F
36   25  51.702722  0.0  0.0000000 S F
37   21  50.880954  0.0  0.0000000 S F
38   99  48.721447  0.0  0.0000000 S F
39   97  48.442793  0.0  0.0000000 S F
40    3  46.001939  0.0  0.0000000 S F
41   92  45.903855  0.0  0.0000000 S F
42   58  45.506047  0.0  0.0000000 S F
43   42  44.777559  0.0  0.0000000 S F
44   87  44.274336  0.0  0.0000000 S F
45   39  44.218182  0.0  0.0000000 S F
46   16  43.849973  0.0  0.0000000 S F
47   90  42.720397  0.0  0.0000000 S F
48    9  42.017388  0.0  0.0000000 S F
49   48  41.820112  0.0  0.0000000 S F
50   22  40.453281  0.0  0.0000000 S F
51    1  40.134270  0.0  0.0000000 S F
52   83  39.639383  0.0  0.0000000 S F
53   76  38.549919  0.0  0.0000000 S F
54   14  38.549461  0.0  0.0000000 S F
55    8  35.751730  0.0  0.0000000 S F
56   80  34.412237  0.0  0.0000000 S F
57   13  33.959616  0.0  0.0000000 S F
58   63  65.358008  0.1  0.2393753   F
59   54  59.933443  0.1  0.2820190   F
60   47  38.606681  0.1  0.2667662   F
61   79  35.045910  0.1  0.3444442   F
62   86  31.395715  0.1  0.2931306   F
63   72  54.068005  0.2  0.4522637   F
64    4  38.792334  0.2  0.7073464   F
65   35  34.617445  0.2  0.6733581   F
66   88  35.147004  0.3  1.4112193   F
67  100  30.785716  0.3  1.1781804   F
68   56  34.258620  0.4  1.3966969   F
69   84  31.738086  0.6  1.9543445   F
70   20  29.471714  1.1  3.2552493   F
71   66  24.244624  1.4  4.0342057   F
72   95  19.429654  1.6  5.6854818   F
73   60  23.047096  1.8  6.0885440   F
74   12  20.625201  2.2  6.9215094   F
75   59  18.742760  2.5 10.1692216   F
76   27  24.793410  2.6  6.9045174   F
77   69  22.766881  2.7  8.3615886   F
78   98  18.328136  2.8 11.2156424   F
79   49  23.117921  3.0  9.1467849   F
80   55  26.192692  3.3 10.7893489   F
81   30  18.374832  3.5  7.9978672   F
82   57  18.911661  3.6 12.7433027   F
83   96  16.859134  3.8 14.9634715   F
84   46  14.632071  3.8 13.7055582   F
85   74  17.892081  4.9 14.8470517   F
86   31  17.597946  5.1 20.4677432   F
87   45  14.734285  5.2 22.2911883   F
88   85  13.488933  5.2 21.3966265   F
89   33  11.898146  5.9 22.3578662   F
90   78  10.339745  6.0 38.4190204   F
91   43  14.577235  6.3 28.0625912   F
92   19  11.400038  6.6 31.4279856   F
93   18   9.955954  6.6 38.4911843   F
94   75  14.767943  6.9 33.2147288   F
95   51  11.009624  7.4 45.2753660   F
96   94   8.452682  8.0 60.7602533   F
97   44   9.085678  8.3 57.2633727    
98  399   7.992158  8.3 68.1891610    
99  104   9.543673  8.5 58.9710191    
100  38   8.201157  8.9 76.2137941    
101  41  14.155428  9.5 33.6140872    
102 598   7.385797  9.8 45.5919137    
103 277   8.133859 10.0 37.8176444    
104 618   8.115440 10.0 85.5516822    
105 946   8.108665 10.0 58.1393519    
106 651   7.659571 10.0 56.1762919    
107 670   7.613712 10.0 84.8279861 
```

A plot with a representation of the potential genes based on OR (vertical axis), FP (horizontal axis) and dFP (size of the circle is inversely proportional to its value) can also be obtained. The plot is similar to Fig.3b in the paper.

```
plotFPvsOR(myORdensity)
```

![plot1](/images/plotFPvsOR.png)

By default, the number of clusters computed by the ORdensity method is used. Other values for the number of clusters can be specified.

```
plotFPvsOR(myORdensity, k = 5)
```

![plot1](/images/plotFPvsOR5.png)

The plot of k values against the silhouette measure is also provided.

```
silhouetteAnalysis(myORdensity)
```

![plot2](/images/silhouetteAnalysis.png)

It is also possible to see a graphic representation of the clustering projected onto the first two principal components

```
clusplotk(myORdensity)
```
![plot3](/images/clusplotk.png)

Other number of clusters can also be checked

```
clusplotk(myORdensity, k = 4)
```
![plot4](/images/clusplotk4.png)

