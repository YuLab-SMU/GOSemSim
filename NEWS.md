# GOSemSim 2.29.1.001

+ update `buildGOmap()` parameter to consistent with `enricher()` and `GSEA()` (2024-02-06, Tue, #47)

# GOSemSim 2.29.1

+ extend `godata()` to support passing a data.frame (can be output of `read.gaf()` or `read.blast2go()`) to 'annoDb' (2023-01-16, Tue)
+ deprecate 'OrgDb' and introduce new parameter 'annoDb' in `godata()` 
+ standardize the output of `read.gaf()` and `read.blast2go()`
+ optimize `buildGOmap()`

# GOSemSim 2.28.0

+ Bioconductor RELEASE_3_18 (2023-10-25, Wed)

# GOSemSim 2.27.3

+ use `check_installed()` to check package dependency (2023-09-12, Tue, #43)

# GOSemSim 2.27.2

+ `read.blast2go()` to parse 'blast2go' result (2023-07-10, Mon, #41, #42)
+ move `buildGOmap()` and `read.gaf()` from 'clusterProfiler' (2023-07-10, Mon)

# GOSemSim 2.27.1

+ semantic similarity measurement support for MPO (2023-04-06, Thu)
+ TCSS semantic similarity measurement support for DO and MPO (2023-04-06, Thu)

# GOSemSim 2.26.0

+ Bioconductor RELEASE_3_17 (2023-05-03, Wed)

# GOSemSim 2.24.0

+ Bioconductor RELEASE_3_16 (2022-11-02, Wed)

# GOSemSim 2.23.1

+ Replacing DO.db with HDO.db (2022-07-29, Mon)

# GOSemSim 2.22.0

+ Bioconductor 3.15 release

# GOSemSim 2.21.1

+ Avoid eval-parse in `load_OrgDb()` (2022-01-10, Mon)

# GOSemSim 2.20

+ Bioconductor 3.14 release

# GOSemSim 2.19.1

+ TCSS method (@qibaiqi, #35; 2021-08-02, Mon)

# GOSemSim 2.18.0

+ Bioconductor 3.13 release

# GOSemSim 2.17.1

+ bug fixed according to the update of GO.db (2020-10-29, Thu)
  - <https://github.com/YuLab-SMU/GOSemSim/issues/32>
  
# GOSemSim 2.16.0

+ Bioconductor 3.12 release (2020-10-28, Wed)

# GOSemSim 2.15.2

+ new site, <https://yulab-smu.top/biomedical-knowledge-mining-book/> for documentation (2020-09-04, Fri)
+ update vignette
+ update `data/gotbl`

# GOSemSim 2.15.1

+ bug fixed of IC method when input IDs contain invalid terms. (2020-07-25, Sat)

# GOSemSim 2.14.0

+ Bioconductor 3.11 release


# GOSemSim 2.13.1

+ add new citation (2020-03-19, Thu)
+ fixed compiling error due to the change of Rcpp 
  - <https://github.com/YuLab-SMU/GOSemSim/issues/27>

# GOSemSim 2.7.1

+ `mgeneSim` and `mclusterSim` now always return matrix (2018-08-08, Wed)
    - <https://www.biostars.org/p/330642/#331598>

# GOSemSim 2.5.1

+ return NA for deprecated IDs (2018-01-09, Fri)
    - <https://support.bioconductor.org/p/105822/#105840>
