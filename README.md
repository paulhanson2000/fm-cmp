# Goals:
+ Reproduce fine mapping results from [DIAMANTE paper](https://doi.org/10.1038/s41588-022-01058-3) ([results to reproduce](https://kp4cd.org/index.php/node/869))
+ Try various fine-mapping methods, see how they perform on this multi-ancestry data
  - [CAVIAR](https://github.com/fhormoz/caviar) ([paper1](https://doi.org/10.1534/genetics.114.167908), [paper2](https://doi.org/10.1016/j.ajhg.2016.10.003))
  - [CAVIARBF](https://bitbucket.org/Wenan/caviarbf/src/master/) ([paper](https://doi.org/10.1534/genetics.116.188953))
  - [FINEMAP](http://www.christianbenner.com/) ([paper1](https://doi.org/10.1093/bioinformatics/btw018), [paper2](https://doi.org/10.1016/j.ajhg.2017.08.012), [paper3](https://doi.org/10.1101/318618))
  - [GUESSFM](https://github.com/chr1swallace/GUESSFM) ([paper](https://doi.org/10.1371/journal.pgen.1005272))
  - [MA-FOCUS](https://www.mancusolab.com/ma-focus) ([paper](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00306-8))
  - ~~MANTRA ([paper](https://doi.org/10.1002/gepi.20630))~~
  - [mJAM](https://github.com/USCbiostats/hJAM) ([paper](https://doi.org/10.1101/2022.12.22.521659))
  - [MsCAVIAR](https://github.com/nlapier2/MsCAVIAR) ([paper](https://doi.org/10.1371/journal.pgen.1009733))
  - [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0) ([paper1](https://doi.org/10.1371/journal.pgen.1004722), [paper2](https://www.cell.com/ajhg/fulltext/S0002-9297(15)00243-8), [paper3](https://doi.org/10.1093%2Fbioinformatics%2Fbtw615))
  - ~~[RIVERA](https://github.com/yueli-compbio/RiVIERA) ([paper](https://doi.org/10.1093/nar/gkw627))~~
  - [SuSiEx](https://github.com/getian107/SuSiEx) ([SuSiEx paper](https://doi.org/10.1101/2023.01.07.23284293), [original SuSie paper](https://doi.org/10.1111/rssb.12388) and [susieR package](https://github.com/stephenslab/susieR))
    - [SuSiE-inf and FINEMAP-inf](https://github.com/FinucaneLab/fine-mapping-inf) ([paper](https://doi.org/10.1101/2022.10.21.513123))
+ Try incorporating functional annotations if possible, if not already built-in to the method (e.g. [fGWAS](https://github.com/joepickrell/fgwas), [PolyFun](https://github.com/omerwe/polyfun))
+ Try using the most suitable multi-ancestry methods on traits other than T2D (e.g. HbA1C, fasting insulin, MAGIC fasting glucose...)

# How to run
This pipeline only works on Linux (for now?).\
If you are on a compute cluster using LMod, load the required modules:
```{bash}
source get_modules_FOR_COMPUTE_CLUSTER_PEOPLE.sh
```
If you are not on a cluster, you will need to get the following dependencies yourself:
+ [R](https://cran.r-project.org/mirrors.html)
+ [GSL](https://www.gnu.org/software/gsl/) (for MsCAVIAR)
+ [LAPACK](https://www.netlib.org/lapack/#_software) (for MsCAVIAR)

You will also need to install [Quarto](https://quarto.org/docs/get-started/) to render certain plots.

To download data and the required R packages, run the following. Only needs to be done once.
```{bash}
./get_data
Rscript get_r_pkgs.R
```
<WIP>
