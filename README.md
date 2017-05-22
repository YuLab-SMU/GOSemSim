GOSemSim: GO semantic similarity measurement
============================================

[![](https://img.shields.io/badge/release%20version-2.2.0-green.svg?style=flat)](https://bioconductor.org/packages/GOSemSim) [![](https://img.shields.io/badge/devel%20version-2.3.0-green.svg?style=flat)](https://github.com/guangchuangyu/GOSemSim) [![Bioc](http://www.bioconductor.org/shields/years-in-bioc/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#since) [![](https://img.shields.io/badge/download-48972/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim) [![](https://img.shields.io/badge/download-1985/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![codecov](https://codecov.io/gh/GuangchuangYu/GOSemSim/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/GOSemSim/) [![Last-changedate](https://img.shields.io/badge/last%20change-2017--05--22-green.svg)](https://github.com/GuangchuangYu/GOSemSim/commits/master) [![commit](http://www.bioconductor.org/shields/commits/bioc/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#svn_source) [![GitHub forks](https://img.shields.io/github/forks/GuangchuangYu/GOSemSim.svg)](https://github.com/GuangchuangYu/GOSemSim/network) [![GitHub stars](https://img.shields.io/github/stars/GuangchuangYu/GOSemSim.svg)](https://github.com/GuangchuangYu/GOSemSim/stargazers)

[![platform](http://www.bioconductor.org/shields/availability/devel/GOSemSim.svg)](https://www.bioconductor.org/packages/devel/bioc/html/GOSemSim.html#archives) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/GOSemSim.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GOSemSim/) [![Linux/Mac Travis Build Status](https://img.shields.io/travis/GuangchuangYu/GOSemSim/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/GOSemSim) [![AppVeyor Build Status](https://img.shields.io/appveyor/ci/Guangchuangyu/GOSemSim/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/GOSemSim)

The semantic comparisons of Gene Ontology (GO) annotations provide quantitative ways to compute similarities between genes and gene groups, and have became important basis for many bioinformatics analysis approaches. `GOSemSim` is an R package for semantic similarity computation among GO terms, sets of GO terms, gene products and gene clusters. `GOSemSim` implemented five methods proposed by *Resnik*, *Schlicker*, *Jiang*, *Lin* and *Wang* respectively.

For details, please visit our project website, <https://guangchuangyu.github.io/GOSemSim>.

-   [Documentation](https://guangchuangyu.github.io/GOSemSim/documentation/)
-   [Featured Articles](https://guangchuangyu.github.io/GOSemSim/featuredArticles/)
-   [Feedback](https://guangchuangyu.github.io/GOSemSim/#feedback)

[![Twitter](https://img.shields.io/twitter/url/https/github.com/GuangchuangYu/GOSemSim.svg?style=social)](https://twitter.com/intent/tweet?hashtags=GOSemSim&url=http://bioinformatics.oxfordjournals.org/content/26/7/976&screen_name=guangchuangyu)

------------------------------------------------------------------------

Please cite the following article when using `GOSemSim`:

**Yu G**<sup>†</sup>, Li F<sup>†</sup>, Qin Y, Bo X<sup>\*</sup>, Wu Y and Wang S<sup>\*</sup>. GOSemSim: an R package for measuring semantic similarity among GO terms and gene products. ***Bioinformatics***. 2010, 26(7):976-978.

[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btq064-green.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btq064) [![](https://img.shields.io/badge/Altmetric-18-green.svg?style=flat)](https://www.altmetric.com/details/100979)

------------------------------------------------------------------------

### Citation

[![citation](https://img.shields.io/badge/cited%20by-275-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=9484177541993722322) [![](https://img.shields.io/badge/cited%20in%20Web%20of%20Science%20Core%20Collection--green.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000276045800023) [![](https://img.shields.io/badge/ESI-Highly%20Cited%20Paper-green.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000276045800023)

    70 +-+--------------+-------------+-------------+----------+
       |                                            *          |
    60 +                                                       +
       |                                                       |
    50 +                                     *                 +
       |                              *                        |
    40 +                       *                               +
       |                *                                      |
    30 +                                                       +
       |                                                   *   |
    20 +                                                       +
       |                                                       |
    10 + *       *                                             +
       +-+--------------+-------------+-------------+----------+
       2010           2012          2014          2016          

### Download stats

[![download](http://www.bioconductor.org/shields/downloads/GOSemSim.svg)](https://bioconductor.org/packages/stats/bioc/GOSemSim/) [![](https://img.shields.io/badge/download-48972/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim) [![](https://img.shields.io/badge/download-1985/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim)

         +-------------+-------------------+--------------------+-------------------+------------------+
         |                                                                                       *     |
         |                                                                                        *    |
    3000 +                                                                                             +
         |                                                                                             |
         |                                                                                   * *       |
         |                                                                                      *      |
         |                                                                                   **        |
         |                                                                                **           |
    2000 +                                                                              *              +
         |                                                                                  *          |
         |                                                                             *               |
         |                                                                             * *             |
         |                                                                         * *                 |
         |                                                                  **    * * *                |
         |                                                              *        *                     |
    1000 +                                                         *     * *  ****                     +
         |                                                       ** ****  *                            |
         |                               *            **  *   ***    *                                 |
         |                                    *     *   ** ***                                         |
         |          * *******         *      **** *  *                                                 |
         |      **** * *     ********* ** ***    * *                                                   |
       0 +   ***                                                                                       +
         +-------------+-------------------+--------------------+-------------------+------------------+
                     2010                2012                 2014                2016
