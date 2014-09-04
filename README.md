#  GOSemSim

Implemented five methods proposed by Resnik, Schlicker, Jiang, Lin and Wang respectively for estimating GO semantic similarities. Support many species, including Anopheles, Arabidopsis, Bovine, Canine, Chicken, Chimp, Coelicolor, E coli strain K12 and Sakai, Fly, Human, Malaria, Mouse, Pig, Rhesus, Rat, Worm, Xenopus, Yeast, and Zebrafish.

## Authors ##

Guangchuang YU, School of Public Health, The University of Hong Kong [http://ygc.name](http://ygc.name)

## Citation ##

Please cite the following article when using `GOSemSim`:

```
Yu G, Li F, Qin Y, Bo X, Wu Y and Wang. S (2010). 
GOSemSim: an R package for measuring semantic similarity among GO terms and gene products.
Bioinformatics, 26(7), pp. 976-978. 
```

URL: [http://bioinformatics.oxfordjournals.org/content/26/7/976.full](http://bioinformatics.oxfordjournals.org/content/26/7/976.full)

## License ##

All source code is copyright, under the Artistic-2.0 License.
For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Installation ##

To install:
 * the latest released version:
   `biocLite("GOSemSim")`
 * the latest development version:
   `install_github("GuangchuangYu/GOSemSim")`

Find out more at [http://www.bioconductor.org/packages/release/bioc/html/GOSemSim.html](http://www.bioconductor.org/packages/release/bioc/html/GOSemSim.html) and check out the [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/GOSemSim.pdf).

## Proper use of GOSemSim ##

I am very glad that many people find GOSemSim useful and [GOSemSim](http://bioinformatics.oxfordjournals.org/content/26/7/976.full) has been cited by 114 (by google scholar, Aug, 2014). 

There are two R packages [BiSEp](http://cran.r-project.org/web/packages/BiSEp/index.html) and [tRanslatome](http://www.bioconductor.org/packages/release/bioc/html/tRanslatome.html) depend on `GOSemSim` and three R packages [clusterProfiler](http://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [DOSE](http://www.bioconductor.org/packages/release/bioc/html/DOSE.html) and [Rcpi](http://www.bioconductor.org/packages/release/bioc/html/Rcpi.html) import `GOSemSim`.

[SemDist](http://www.bioconductor.org/packages/devel/bioc/html/SemDist.html) package copy some of the source code from `GOSemSim` with acknowledging within source code and document.

[ppiPre](http://cran.r-project.org/web/packages/ppiPre/index.html) package copy many source code from `GOSemSim` without any acknowledgement in souce code or document and did not cited `GOSemSim` in their [publication](http://www.biomedcentral.com/1752-0509/7/S2/S8). This violates the restriction of open source license.

For R developers, if you found functions provided in `GOSemSim` useful, please depends or imports `GOSemSim`.
If you would like to copy and paste source code, you should acknowledge the source code was copied/derived from `GOSemSim` authored by Guangchuang Yu <guangchuangyu@gmail.com> within source code, add `GOSemSim` in Suggests field and also includes the following reference in the man files for functions that copied/derived from `GOSemSim` and cited in vignettes.

```
\references{
  Yu et al. (2010) GOSemSim: an R package for measuring
  semantic similarity among GO terms and gene products
  \emph{Bioinformatics} (Oxford, England), 26:7 976--978,
  April 2010. ISSN 1367-4803
  \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
  PMID: 20179076
}
```

You are welcome to use `GOSemSim` in the way you like, but please cite it and give it the proper credit. I hope you can understand.


## Bugs/Feature requests ##

 - If you have any, [let me know](https://github.com/GuangchuangYu/GOSemSim/issues). Thx!


