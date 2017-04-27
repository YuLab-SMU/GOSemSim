---
html_preview: False
output:
  md_document:
    variant: markdown
---

GOSemSim: GO semantic similarity measurement
============================================

<!-- AddToAny BEGIN -->
<div class="a2a_kit a2a_kit_size_32 a2a_default_style">

<a class="a2a_dd" href="//www.addtoany.com/share"></a>
<a class="a2a_button_facebook"></a> <a class="a2a_button_twitter"></a>
<a class="a2a_button_google_plus"></a>
<a class="a2a_button_pinterest"></a> <a class="a2a_button_reddit"></a>
<a class="a2a_button_sina_weibo"></a> <a class="a2a_button_wechat"></a>
<a class="a2a_button_douban"></a>

</div>

<script async src="//static.addtoany.com/menu/page.js"></script>
<!-- AddToAny END -->
<link rel="stylesheet" href="https://guangchuangyu.github.io/css/font-awesome.min.css">
<link rel="stylesheet" href="https://guangchuangyu.github.io/css/academicons.min.css">

[![](https://img.shields.io/badge/release%20version-2.2.0-blue.svg?style=flat)](https://bioconductor.org/packages/GOSemSim)
[![](https://img.shields.io/badge/devel%20version-2.3.0-blue.svg?style=flat)](https://github.com/guangchuangyu/GOSemSim)
[![](https://img.shields.io/badge/download-47079/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim)
[![](https://img.shields.io/badge/download-2208/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/GOSemSim)

The semantic comparisons of Gene Ontology (GO) annotations provide
quantitative ways to compute similarities between genes and gene groups,
and have became important basis for many bioinformatics analysis
approaches. `GOSemSim` is an R package for semantic similarity
computation among GO terms, sets of GO terms, gene products and gene
clusters. `GOSemSim` implemented five methods proposed by *Resnik*,
*Schlicker*, *Jiang*, *Lin* and *Wang* respectively.

`GOSemSim` is released within the
[Bioconductor](https://bioconductor.org/packages/GOSemSim) project and
the source code is hosted on
<a href="https://github.com/GuangchuangYu/GOSemSim"><i class="fa fa-github fa-lg"></i>
GitHub</a>.

<i class="fa fa-user"></i> Author
---------------------------------

Guangchuang Yu, School of Public Health, The University of Hong Kong.

<a href="https://twitter.com/guangchuangyu"><i class="fa fa-twitter fa-3x"></i></a>
<a href="https://guangchuangyu.github.io/blog_images/biobabble.jpg"><i class="fa fa-wechat fa-3x"></i></a>
<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=Guangchuang+Yu[Author+-+Full]"><i class="ai ai-pubmed ai-3x"></i></a>
<a href="https://scholar.google.com.hk/citations?user=DO5oG40AAAAJ&hl=en"><i class="ai ai-google-scholar ai-3x"></i></a>
<a href="https://orcid.org/0000-0002-6485-8781"><i class="ai ai-orcid ai-3x"></i></a>
<a href="https://impactstory.org/u/0000-0002-6485-8781"><i class="ai ai-impactstory ai-3x"></i></a>

<i class="fa fa-book"></i> Citation
-----------------------------------

Please cite the following article when using `GOSemSim`:

[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btq064-blue.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btq064)
[![](https://img.shields.io/badge/Altmetric-18-blue.svg?style=flat)](https://www.altmetric.com/details/100979)
[![citation](https://img.shields.io/badge/cited%20by-273-blue.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=9484177541993722322)
[![](https://img.shields.io/badge/cited%20in%20Web%20of%20Science%20Core%20Collection--blue.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000276045800023)
[![](https://img.shields.io/badge/ESI-Highly%20Cited%20Paper-blue.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000276045800023)

**Yu G**<sup>†</sup>, Li F<sup>†</sup>, Qin Y, Bo X<sup>\*</sup>, Wu Y
and Wang S<sup>\*</sup>. GOSemSim: an R package for measuring semantic
similarity among GO terms and gene products. ***Bioinformatics***. 2010,
26(7):976-978.

<i class="fa fa-pencil"></i> Featured Articles
----------------------------------------------

<img src="https://guangchuangyu.github.io/featured_img/GOSemSim/2014PNAS.png" width="650">

<i class="fa fa-hand-o-right"></i> Find out more on
<i class="fa fa-pencil"></i> [Featured
Articles](https://guangchuangyu.github.io/GOSemSim/featuredArticles/).

<i class="fa fa-download"></i> Installation
-------------------------------------------

Install `GOSemSim` is easy, follow the guide in the [Bioconductor
page](https://bioconductor.org/packages/GOSemSim/):

``` {.r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("GOSemSim")
```

<i class="fa fa-cogs"></i> Overview
-----------------------------------

#### <i class="fa fa-angle-double-right"></i> Methods

-   Information content based methods proposed by *Resnik*, *Schlicker*,
    *Jiang* and *Lin*
-   Graph structure based method proposed by *Wang*

#### <i class="fa fa-angle-double-right"></i> Combine methods for aggregating multiple GO terms

-   max
-   avg
-   rcmax
-   BMA

#### <i class="fa fa-angle-double-right"></i> Functions

-   goSim and mgoSim for measureing semantic similarity among GO terms
-   geneSim and mgeneSim for measureing semantic similarity among genes
-   clusterSim and mclusterSim for measureing semantic similarity among
    gene clusters

#### <i class="fa fa-angle-double-right"></i> Supported organisms

-   20 species that have OrgDb available in Bioconductor
-   Many other species with e GO annotation query online via
    [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/))

<i class="fa fa-hand-o-right"></i> Find out details and examples on
<i class="fa fa-book"></i>
[Documentation](https://guangchuangyu.github.io/GOSemSim/documentation/).

<i class="fa fa-wrench"></i> Related Tools
------------------------------------------

<ul class="fa-ul">
    <li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/clusterProfiler">clusterProfiler</a> for Ontologies/pathways enrichment analyses</li>
    <li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/DOSE">DOSE</a> for Disease Ontology Semantic and Enrichment analyses</li>
    <li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/meshes">meshes</a> for MeSH Enrichment and Semantic analysis</li>

</ul>
<i class="fa fa-hand-o-right"></i> Find out more on
[projects](https://guangchuangyu.github.io/#projects).

<i class="fa fa-code-fork"></i> Projects that depend on *GOSemSim*
------------------------------------------------------------------

#### <i class="fa fa-angle-double-right"></i> CRAN packages

-   [BiSEp](https://cran.r-project.org/package=BiSEp): Toolkit to
    Identify Candidate Synthetic Lethality
-   [LANDD](https://cran.r-project.org/package=LANDD): Liquid
    Association for Network Dynamics Detection
-   [ppiPre](https://cran.r-project.org/package=ppiPre): Predict
    Protein-Protein Interactions Based on Functional and Topological
    Similarities

#### <i class="fa fa-angle-double-right"></i> Bioconductor packages

-   [BioMedR](https://www.bioconductor.org/packages/BioMedR): Generating
    Various Molecular Representations for Chemicals, Proteins, DNAs/RNAs
    and Their Interactions
-   [clusterProfiler](https://www.bioconductor.org/packages/clusterProfiler):
    statistical analysis and visualization of functional profiles for
    genes and gene clusters
-   [DOSE](https://www.bioconductor.org/packages/DOSE): Disease Ontology
    Semantic and Enrichment analysis
-   [meshes](https://www.bioconductor.org/packages/meshes): MeSH
    Enrichment and Semantic analyses
-   [Rcpi](https://www.bioconductor.org/packages/Rcpi): Molecular
    Informatics Toolkit for Compound-Protein Interaction in Drug
    Discovery
-   [tRanslatome](https://www.bioconductor.org/packages/tRanslatome):
    Comparison between multiple levels of gene expression

<i class="fa fa-comment"></i> Feedback
--------------------------------------

<ul class="fa-ul">
    <li><i class="fa-li fa fa-hand-o-right"></i> Please make sure you have followed <a href="https://guangchuangyu.github.io/2016/07/how-to-bug-author/"><strong>the important guide</strong></a> before posting any issue/question</li>
    <li><i class="fa-li fa fa-bug"></i> For bugs or feature requests, please post to <i class="fa fa-github-alt"></i> [github issue](https://github.com/GuangchuangYu/GOSemSim/issues)</li>
    <li><i class="fa-li fa fa-question"></i>  For user questions, please post to [Bioconductor support site](https://support.bioconductor.org/) and [Biostars](https://www.biostars.org/). We are following every post tagged with **GOSemSim**</li>
    <li><i class="fa-li fa fa-commenting"></i> Join the group chat on <a href="https://twitter.com/hashtag/GOSemSim"><i class="fa fa-twitter fa-lg"></i></a> and <a href="http://huati.weibo.com/k/GOSemSim"><i class="fa fa-weibo fa-lg"></i></a></li>

</ul>
