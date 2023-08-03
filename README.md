# How to run TL;DR
Linux only.
```{bash}
# Download this repository, and go inside.
git clone git@github.com:paulhanson2000/fm-cmp.git
cd fm-cmp

# If you're on a compute cluster using LMod, run this every time. Tested on Compute Canada & BU SCC.
# Otherwise, on your personal computer make sure you have the dependencies.
source init/get_lmod_modules.sh

# Need only run once
init/get_data.sh
init/get_tools.sh
Rscript init/get_r_pkgs.R

# Run pipeline (need Quarto, get it from quarto.org)
quarto render

# TODO: Redbean stuff to look at results
```
To use this pipeline on your own data, simply edit the `config/` files!

# Goals:
+ Reproduce fine mapping results from [2022 DIAMANTE paper](https://doi.org/10.1038/s41588-022-01058-3) ([results to reproduce](https://kp4cd.org/index.php/node/869))
+ Try various fine-mapping methods, see how they perform on this multi-ancestry data
  - [CAVIAR](https://github.com/fhormoz/caviar) ([paper1](https://doi.org/10.1534/genetics.114.167908), [paper2](https://doi.org/10.1016/j.ajhg.2016.10.003))
  - [CAVIARBF](https://bitbucket.org/Wenan/caviarbf/src/master/) ([paper](https://doi.org/10.1534/genetics.116.188953))
  - [FINEMAP](http://www.christianbenner.com/) ([paper1](https://doi.org/10.1093/bioinformatics/btw018), [paper2](https://doi.org/10.1016/j.ajhg.2017.08.012), [paper3](https://doi.org/10.1101/318618))
  - [GUESSFM](https://github.com/chr1swallace/GUESSFM) ([paper](https://doi.org/10.1371/journal.pgen.1005272))
  - [MA-FOCUS](https://www.mancusolab.com/ma-focus) ([paper](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00306-8))
  - ~~MANTRA ([paper](https://doi.org/10.1002/gepi.20630))~~
  - [MESuSiE](https://github.com/borangao/meSuSie)
  - [mJAM](https://github.com/USCbiostats/hJAM) ([paper](https://doi.org/10.1101/2022.12.22.521659))
  - [MsCAVIAR](https://github.com/nlapier2/MsCAVIAR) ([paper](https://doi.org/10.1371/journal.pgen.1009733))
  - [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0) ([paper1](https://doi.org/10.1371/journal.pgen.1004722), [paper2](https://www.cell.com/ajhg/fulltext/S0002-9297(15)00243-8), [paper3](https://doi.org/10.1093%2Fbioinformatics%2Fbtw615))
  - ~~[RIVERA](https://github.com/yueli-compbio/RiVIERA) ([paper](https://doi.org/10.1093/nar/gkw627))~~
  - [SparsePro](https://github.com/zhwm/SparsePro) ([paper](https://doi.org/10.1101/2021.10.04.463133))
  - [SuSiEx](https://github.com/getian107/SuSiEx) ([SuSiEx paper](https://doi.org/10.1101/2023.01.07.23284293), [original SuSie paper](https://doi.org/10.1111/rssb.12388) and [susieR package](https://github.com/stephenslab/susieR))
    - [SuSiE-inf and FINEMAP-inf](https://github.com/FinucaneLab/fine-mapping-inf) ([paper](https://doi.org/10.1101/2022.10.21.513123))
+ Incorporating functional annotations if possible, if not already built-in to the method (e.g. [fgwas](https://github.com/joepickrell/fgwas), [PolyFun](https://github.com/omerwe/polyfun))
+ Try multiple (multi-ancestry) reference panels such as 1000 Genomes, TOPMed.
+ Try other traits than T2D (e.g. HbA1C, fasting insulin, MAGIC fasting glucose...)

# Dependencies
If you are not on a cluster, you will need to get the following dependencies yourself:
+ [R](https://cran.r-project.org/mirrors.html) (tested on: version 4.2.2 (works), 4.3.1 (fails))
+ [Python](https://www.python.org/) (tested on: version 3.11.2 (works), 3.8.10 (fails))
  - [virtualenv](https://virtualenv.pypa.io/en/latest/installation.html)
+ [bcftools](http://www.htslib.org/download/)
+ [Boost](https://www.boost.org/) (for fgwas)
+ [GSL](https://www.gnu.org/software/gsl/) (for MsCAVIAR)
+ [LAPACK](https://www.netlib.org/lapack/#_software) (for MsCAVIAR)

# TODO
Document how to get a GitHub token to be allowed to download more stuff, may be required for the R packages hosted on GitHub
Document in a miscellaneous section or s/t: Tip: LD data is big. If you are on a Linux compute cluster, you recommended to make the in/ld directory a symbolic link to your cluster's scratch space.
Add more guidance about how to use your own data, more than "just edit the configs". Mention all files needed and have example code blocks of what each file should look like.
