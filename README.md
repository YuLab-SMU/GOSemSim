GOSemSim: GO semantic similarity measurement
============================================

[![releaseVersion](https://img.shields.io/badge/release%20version-1.30.3-green.svg?style=flat)](https://bioconductor.org/packages/GOSemSim) [![develVersion](https://img.shields.io/badge/devel%20version-1.99.4-green.svg?style=flat)](https://github.com/GuangchuangYu/GOSemSim) [![Bioc](http://www.bioconductor.org/shields/years-in-bioc/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#since) [![total](https://img.shields.io/badge/downloads-56966/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim) [![month](https://img.shields.io/badge/downloads-2125/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![codecov](https://codecov.io/gh/GuangchuangYu/GOSemSim/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/GOSemSim/) [![Last-changedate](https://img.shields.io/badge/last%20change-2016--09--26-green.svg)](https://github.com/GuangchuangYu/GOSemSim/commits/master) [![commit](http://www.bioconductor.org/shields/commits/bioc/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#svn_source) [![GitHub forks](https://img.shields.io/github/forks/GuangchuangYu/GOSemSim.svg)](https://github.com/GuangchuangYu/GOSemSim/network) [![GitHub stars](https://img.shields.io/github/stars/GuangchuangYu/GOSemSim.svg)](https://github.com/GuangchuangYu/GOSemSim/stargazers)

[![platform](http://www.bioconductor.org/shields/availability/devel/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#archives) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/GOSemSim.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GOSemSim/) [![Linux/Mac Travis Build Status](https://img.shields.io/travis/GuangchuangYu/GOSemSim/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/GOSemSim) [![AppVeyor Build Status](https://img.shields.io/appveyor/ci/Guangchuangyu/GOSemSim/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/GOSemSim) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-green.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-gosemsim/README.html)

The semantic comparisons of Gene Ontology (GO) annotations provide quantitative ways to compute similarities between genes and gene groups, and have became important basis for many bioinformatics analysis approaches. `GOSemSim` is an R package for semantic similarity computation among GO terms, sets of GO terms, gene products and gene clusters. `GOSemSim` implemented five methods proposed by *Resnik*, *Schlicker*, *Jiang*, *Lin* and *Wang* respectively.

[![Twitter](https://img.shields.io/twitter/url/https/github.com/GuangchuangYu/GOSemSim.svg?style=social)](https://twitter.com/intent/tweet?hashtags=GOSemSim&url=http://bioinformatics.oxfordjournals.org/content/26/7/976&screen_name=guangchuangyu)

------------------------------------------------------------------------

Please cite the following article when using `GOSemSim`:

**Yu G**<sup>†</sup>, Li F<sup>†</sup>, Qin Y, Bo X<sup>\*</sup>, Wu Y and Wang S<sup>\*</sup>. GOSemSim: an R package for measuring semantic similarity among GO terms and gene products. ***Bioinformatics***. 2010, 26(7):976-978.

[![doi](https://img.shields.io/badge/doi-10.1093/bioinformatics/btq064-green.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btq064) [![citation](https://img.shields.io/badge/cited%20by-217-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=9484177541993722322) [![Altmetric](https://img.shields.io/badge/Altmetric-6-green.svg?style=flat)](https://www.altmetric.com/details/100979)

------------------------------------------------------------------------

For details, please visit our project website, <https://guangchuangyu.github.io/GOSemSim>.

-   [Documentation](https://guangchuangyu.github.io/GOSemSim/documentation/)
-   [Featured Articles](https://guangchuangyu.github.io/GOSemSim/featuredArticles/)
-   [Feedback](https://guangchuangyu.github.io/GOSemSim/#feedback)

### Citation

[![citation](https://img.shields.io/badge/cited%20by-217-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=9484177541993722322)

       +-+--------+-------+-------+--------+-------+-------+---+
       |                                           *           |
    40 +                          *        *               *   +
       |                                                       |
       |                  *                                    |
    30 +                                                       +
       |                                                       |
       |                                                       |
    20 +                                                       +
       |                                                       |
       |                                                       |
    10 +                                                       +
       | *        *                                            |
       +-+--------+-------+-------+--------+-------+-------+---+
       2010     2011    2012    2013     2014    2015    2016   

### Download stats

[![download](http://www.bioconductor.org/shields/downloads/GOSemSim.svg)](https://bioconductor.org/packages/stats/bioc/GOSemSim/) [![total](https://img.shields.io/badge/downloads-56966/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim) [![month](https://img.shields.io/badge/downloads-2125/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim)

         +-------------+----------------------+---------------------+---------------------+------------+
         |                                                                                       *     |
         |                                                                                        *    |
         |                                                                                     *       |
    2000 +                                                                                             +
         |                                                                                             |
         |                                                                                   *         |
         |                                                                                      *      |
    1500 +                                                                                    *        +
         |                                                                                 *           |
         |                                                                              ***            |
         |                                                                        **        *          |
         |                                                                                             |
    1000 +                                                                   *  **  ****               +
         |                                                              *   * *                        |
         |                                                            ** * *   *                       |
         |                                  *             ** **   ****    *                            |
     500 +                                                  *                                          +
         |                  *                    *     **   *  ***                                     |
         |               *   *                  * **     *                                             |
         |       *** **** **  * ****  *** ** * *    ***                                                |
         |     **   *          *    **   *    *                                                        |
       0 +   **                                                                                        +
         +-------------+----------------------+---------------------+---------------------+------------+
                     2010                   2012                  2014                  2016
