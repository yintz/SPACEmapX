<img src="https://www.nds.ox.ac.uk/images/logos/secondary-logo" height="75" /> <img src="https://www.nds.ox.ac.uk/images/logos/primary-logo" height="75"/> 

<a href="https://twitter.com/intent/follow?screen_name=Wien_Yin">
<img src="https://img.shields.io/twitter/follow/Wien_Yin?style=social&logo=X",alt="follow on Twitter"></a>
<a href="https://twitter.com/intent/follow?screen_name=OxPCaBiol">
<img src="https://img.shields.io/twitter/follow/OxPCaBiol?style=social&logo=X",alt="follow on Twitter"></a>

[![](https://img.shields.io/badge/SPACEmapX-version0.99-blue.svg)](https://github.com/yintz/SPACEmapX/releases)[![](https://img.shields.io/github/last-commit/yintz/SPACEmapX.svg)](https://github.com/yintz/SPACEmapX/commits/main)


# SPACEMapX
This package provides guidance for analysis of spatially resolved tissue (Visium ST), focusing on clonal dynamics derived from infer copy number status.

Our particular goal has been identification of the "lethal clone", defined as that which metasises from primary cancer tissue to lymph nodes. however, this package can be used for any analysis of clonal evolution and/or somatic mosaicism in heterogeneous tissue.  

If there is anything you need to assist with running this package, please feel free to contact us via GitHub "issue" messaging

## Functions

This packages contains 4 major functions.
1 Infer Copy Number Variation based on transcriptomics data.
2 Spatial data wrapper
3 Heatmap Plots
4 Dendrogram selector


## Pre-requirement
This package is based on InferCNV. SPACEmapX is designed as a helper package for it.


#### Software Requirements
JAGS
R libraries: 
graphics, grDevices, RColorBrewer, gplots, futile.logger, stats, utils, methods, ape, Matrix, fastcluster, dplyr, HiddenMarkov, ggplot2, edgeR, coin, caTools, digest, reshape, rjags, fitdistrplus, future, foreach, doParallel, BiocGenerics, SummarizedExperiment, SingleCellExperiment, tidyr, parallel, coda, gridExtra, argparse


# tutorial 
https://github.com/yintz/SPACEmapX/wiki/SPACEmapX-tutorial


## Installation
package can be installed through GitHub using;
``` r
install.packages("remotes")
remotes::install_github("yintz/SPACEmapX")
```

## Functions
``` r
ShowTwig(TrigID)
ShowTwigSectionName(TrigID)
```



`SpaceMapX` is an R package that collects useful tools for Spatially Resolved Transcriptomics data analysis and visualization.


# 10x ST Data

This package is designed to work with Visium Fresh Forzen, FFPE V1,V2 data only.

If you need the ST test data, this is the data we used. 
[MendeleyRepository](https://data.mendeley.com/v1/datasets/svw96g68dv/draft?a=3f263217-2bd3-4a3c-8125-8c517c3a9e29).
for more ST data, you can also go to 10X website to download.


# Funding 
