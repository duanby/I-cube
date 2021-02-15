# Interactive identification of individual treatment effect
This code is an accompaniment to the paper titled "Interactive identification of individual treatment effect with FDR control", and produce all plots in the paper.

## Overview
We propose methods to identify subjects with positive effects, with a finite-sample error control on the proportion of false identification, in a randomized experiment.

To reproduce the experiments:
```
$Rscript code/reproduce.R
```
To generate the figures:
```
$Rscript code/plots.R
```


## Dependencies
The code was tested using R (version 3.6.0) and the following packages:
magrittr, splines, ggplot2, randomForest, caTools, dplyr, tidyr, doParallel, foreach
