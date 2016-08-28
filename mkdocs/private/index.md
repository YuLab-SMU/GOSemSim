<!-- addtoany:= -->

<!-- release:=GOSemSim -->
<!-- devel:=GOSemSim -->
<!-- download:=GOSemSim:=total -->
<!-- download:=GOSemSim:=month -->

The semantic comparisons of Gene Ontology (GO) annotations provide quantitative ways to compute similarities between genes and gene groups, and have became important basis for many bioinformatics analysis approaches. `GOSemSim` is an R package for semantic similarity computation among GO terms, sets of GO terms, gene products and gene clusters. `GOSemSim` implemented five methods proposed by _Resnik_, _Schlicker_, _Jiang_, _Lin_ and _Wang_ respectively.


`GOSemSim` is released within the [Bioconductor](https://bioconductor.org/packages/GOSemSim) project and the source code is hosted on <a href="https://github.com/GuangchuangYu/GOSemSim"><i class="fa fa-github fa-lg"></i> GitHub</a>.

## <i class="fa fa-user"></i> Author

Guangchuang Yu, School of Public Health, The University of Hong Kong.

## <i class="fa fa-book"></i> Citation

<!-- doi:=10.1093/bioinformatics/btq064 -->
<!-- citation:=tuHXwOkdijsC:=9484177541993722322 -->
<!-- altmetric:=100979 -->

Please cite the following article when using `GOSemSim`:

__Yu G__<sup>†</sup>, Li F<sup>†</sup>, Qin Y, Bo X<sup>\*</sup>, Wu Y and Wang S<sup>\*</sup>. 
GOSemSim: an R package for measuring semantic similarity among GO terms and gene products.
__*Bioinformatics*__. 2010, 26(7):976-978. 


## <i class="fa fa-pencil"></i> Featured Articles

<img src="featured_img/2014PNAS.png" width="650">

<i class="fa fa-hand-o-right"></i> Find out more on <i class="fa fa-pencil"></i> [Featured Articles](https://guangchuangyu.github.io/GOSemSim/featuredArticles/).

## <i class="fa fa-download"></i> Installation

Install `GOSemSim` is easy, follow the guide in the [Bioconductor page](https://bioconductor.org/packages/GOSemSim/):

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("GOSemSim")
```

## <i class="fa fa-cogs"></i> Overview

#### <i class="fa fa-angle-double-right"></i> Methods

+ Information content based methods proposed by _Resnik_, _Schlicker_, _Jiang_ and _Lin_
+ Graph structure based method proposed by _Wang_

#### <i class="fa fa-angle-double-right"></i> Combine methods for aggregating multiple GO terms

+ max
+ avg
+ rcmax
+ BMA

#### <i class="fa fa-angle-double-right"></i> Functions

+ goSim and mgoSim for measureing semantic similarity among GO terms
+ geneSim and mgeneSim for measureing semantic similarity among genes
+ clusterSim and mclusterSim for measureing semantic similarity among gene clusters

#### <i class="fa fa-angle-double-right"></i> Supported organisms

+ 20 species that have OrgDb available in Bioconductor
+ Many other species with e GO annotation query online via [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/))

<i class="fa fa-hand-o-right"></i> Find out details and examples on <i class="fa fa-book"></i> [Documentation](https://guangchuangyu.github.io/GOSemSim/documentation/).

## <i class="fa fa-wrench"></i> Related Tools

<ul class="fa-ul">
	<li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/clusterProfiler">clusterProfiler</a> for Ontologies/pathways enrichment analyses</li>
	<li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/DOSE">DOSE</a> for Disease Ontology Semantic and Enrichment analyses</li>
</ul>

<i class="fa fa-hand-o-right"></i> Find out more on [projects](https://guangchuangyu.github.io/#projects).

## <i class="fa fa-code-fork"></i> Projects that depend on GOSemSim

#### <i class="fa fa-angle-double-right"></i> CRAN packages

+ [BiSEp](https://cran.r-project.org/package=BiSEp): Toolkit to Identify Candidate Synthetic Lethality
+ [LANDD](https://cran.r-project.org/package=LANDD): Liquid Association for Network Dynamics Detection 
+ [ppiPre](https://cran.r-project.org/package=ppiPre): Predict Protein-Protein Interactions Based on Functional and Topological Similarities
+ [protr](https://cran.r-project.org/package=protr): Generating Various Numerical Representation Schemes of Protein Sequence

#### <i class="fa fa-angle-double-right"></i> Bioconductor packages

+ [DOSE](https://www.bioconductor.org/packages/DOSE/): Disease Ontology Semantic and Enrichment analysis
+ [Rcpi](https://www.bioconductor.org/packages/Rcpi/): Toolkit for Compound-Protein Interaction in Drug Discovery
+ [tRanslatome](https://www.bioconductor.org/packages/tRanslatome/): Comparison between multiple levels of gene expression


<!--
<i class="fa fa-hand-o-right"></i> Find out more on <i class="fa fa-github-alt"></i> [github](http://scisoft-net-map.isri.cmu.edu/application/GOSemSim/gitprojects).
-->


## <i class="fa fa-comment"></i> Feedback
<ul class="fa-ul">
	<li><i class="fa-li fa fa-hand-o-right"></i> Please make sure you [follow the guide](https://guangchuangyu.github.io/2016/07/how-to-bug-author/) before posting any issue/question</li>
	<li><i class="fa-li fa fa-bug"></i> For bugs or feature requests, please post to <i class="fa fa-github-alt"></i> [github issue](https://github.com/GuangchuangYu/GOSemSim/issues)</li>
	<li><i class="fa-li fa fa-question"></i>  For user questions, please post to [Bioconductor support site](https://support.bioconductor.org/) and [Biostars](https://www.biostars.org/). We are following every post tagged with **GOSemSim**</li>
	<li><i class="fa-li fa fa-commenting"></i> Join the group chat on <a href="https://twitter.com/hashtag/GOSemSim"><i class="fa fa-twitter fa-lg"></i></a> and <a href="http://huati.weibo.com/k/GOSemSim"><i class="fa fa-weibo fa-lg"></i></a></li>
</ul>
