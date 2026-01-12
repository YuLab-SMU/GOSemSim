# GOSemSim: GO semantic similarity measurement

[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#since)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![platform](http://www.bioconductor.org/shields/availability/devel/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#archives)
[![codecov](https://codecov.io/gh/YuLab-SMU/GOSemSim/branch/master/graph/badge.svg)](https://codecov.io/gh/YuLab-SMU/GOSemSim/)

[![](https://img.shields.io/badge/release%20version-2.34.0-green.svg)](https://www.bioconductor.org/packages/GOSemSim)
[![](https://img.shields.io/badge/devel%20version-2.35.0-green.svg)](https://github.com/YuLab-SMU/GOSemSim)
[![](https://img.shields.io/badge/download-1061088/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/GOSemSim)
[![](https://img.shields.io/badge/download-21399/month-blue.svg)](https://bioconductor.org/packages/stats/bioc/GOSemSim)

The semantic comparisons of Gene Ontology (GO) annotations provide
quantitative ways to compute similarities between genes and gene groups,
and have became important basis for many bioinformatics analysis
approaches. GOSemSim is an R package for semantic similarity computation
among GO terms, sets of GO terms, gene products and gene clusters.
GOSemSim implemented five methods proposed by Resnik, Schlicker, Jiang,
Lin and Wang respectively.

## :writing_hand: Authors

Guangchuang YU <https://yulab-smu.top>

School of Basic Medical Sciences, Southern Medical University

Learn more at <https://yulab-smu.top/contribution-knowledge-mining/>.

If you use [GOSemSim](http://bioconductor.org/packages/GOSemSim) in
published research, please cite:

- **Yu G**. [Gene Ontology Semantic Similarity Analysis Using
  GOSemSim](http://dx.doi.org/10.1007/978-1-0716-0301-7_11). In:
  Kidder B. (eds) Stem Cell Transcriptional Networks. ***Methods in
  Molecular Biology***, 2020, 2117:207-215. Humana, New York, NY.
- **Yu G**<sup>\#</sup>, Li F<sup>\#</sup>, Qin Y, Bo X<sup>\*</sup>, Wu
  Y and Wang S<sup>\*</sup>. [GOSemSim: an R package for measuring
  semantic similarity among GO terms and gene
  products](http://dx.doi.org/10.1093/bioinformatics/btq064).
  ***Bioinformatics***. 2010, 26(7):976-978.

## :arrow_double_down: Installation

Get the released version from Bioconductor:

``` r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("GOSemSim")
```

Or the development version from github:

``` r
## install.packages("remotes")
remotes::install_github("YuLab-SMU/GOSemSim")
```

## :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
