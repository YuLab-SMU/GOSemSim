# GOSemSim: GO semantic similarity measurement

[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#since)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![platform](http://www.bioconductor.org/shields/availability/devel/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#archives)
[![codecov](https://codecov.io/gh/GuangchuangYu/GOSemSim/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/GOSemSim/)

[![](https://img.shields.io/badge/release%20version-2.12.0-green.svg)](https://www.bioconductor.org/packages/GOSemSim)
[![](https://img.shields.io/badge/devel%20version-2.13.1-green.svg)](https://github.com/guangchuangyu/GOSemSim)
[![Linux/Mac Travis Build
Status](https://img.shields.io/travis/GuangchuangYu/GOSemSim/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/GOSemSim)
[![AppVeyor Build
Status](https://img.shields.io/appveyor/ci/Guangchuangyu/GOSemSim/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/GOSemSim)

[![](https://img.shields.io/badge/download-176874/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/GOSemSim)
[![](https://img.shields.io/badge/download-5604/month-blue.svg)](https://bioconductor.org/packages/stats/bioc/GOSemSim)
[![download](http://www.bioconductor.org/shields/downloads/release/GOSemSim.svg)](https://bioconductor.org/packages/stats/bioc/GOSemSim)

The semantic comparisons of Gene Ontology (GO) annotations provide
quantitative ways to compute similarities between genes and gene groups,
and have became important basis for many bioinformatics analysis
approaches. GOSemSim is an R package for semantic similarity computation
among GO terms, sets of GO terms, gene products and gene clusters.
GOSemSim implemented five methods proposed by Resnik, Schlicker, Jiang,
Lin and Wang respectively.

[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=GOSemSim&url=http://bioinformatics.oxfordjournals.org/content/26/7/976&screen_name=guangchuangyu)
[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/GuangchuangYu)
[![](https://img.shields.io/badge/follow%20me%20on-WeChat-green.svg)](https://guangchuangyu.github.io/blog_images/biobabble.jpg)

## :writing\_hand: Authors

Guangchuang YU <https://guangchuangyu.github.io>

School of Basic Medical Sciences, Southern Medical University

If you use [GOSemSim](http://bioconductor.org/packages/GOSemSim) in
published research, please cite:

  - **Yu G**. [Gene Ontology Semantic Similarity Analysis Using
    GOSemSim](http://dx.doi.org/10.1007/978-1-0716-0301-7_11). In:
    Kidder B. (eds) Stem Cell Transcriptional Networks. ***Methods in
    Molecular Biology***, 2020, 2117:207-215. Humana, New York, NY.
  - **Yu G**<sup>\#</sup>, Li F<sup>\#</sup>, Qin Y, Bo X<sup>\*</sup>,
    Wu Y and Wang S<sup>\*</sup>. [GOSemSim: an R package for measuring
    semantic similarity among GO terms and gene
    products](http://dx.doi.org/10.1093/bioinformatics/btq064).
    ***Bioinformatics***. 2010, 26(7):976-978.

## :arrow\_double\_down: Installation

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
remotes::install_github("GuangchuangYu/GOSemSim")
```

## :sparkling\_heart: Contributing

We welcome any contributions\! By participating in this project you
agree to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
