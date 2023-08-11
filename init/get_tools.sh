#!/bin/bash

git submodule update --init --remote # For tools hosted on GitHub, download their source

# Detect if on Compute Canada. If so, will use pre-installed Python things 
PIPFLAGS=""
uname -n | grep calculquebec
if [ $? -eq 0 ]; then PIPFLAGS="--no-index"; fi

# Exit if anything fails
set -e 

# fgwas
cd third_party/fgwas/
if ! [ -f src/fgwas ]; then
  autoreconf -f -i # Because https://stackoverflow.com/a/33286344
  ./configure
  sed -i -e 's/iter <5000/iter<1000/' -e 's/iter > 4999/iter>999/' src/SNPs.cpp # NOTE: cursed editing of fgwas's source to increase speed by reducing the number of iterations before it concludes a model doesn't converge.
  make
fi
cd -

# FINEMAP
cd third_party/
if ! [ -f finemap/finemap_v1.4.2_x86_64 ]; then
  curl -O http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
  tar -zxf finemap_v1.4.2_x86_64.tgz
  rm       finemap_v1.4.2_x86_64.tgz
  mv       finemap_v1.4.2_x86_64 finemap
fi
cd -

# liftOver
mkdir -p third_party/liftover/
cd       third_party/liftover/
if ! [ -f liftOver ]; then
  curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
  chmod +x liftOver
  curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz # TODO: instead of hardcoded downloads for this specific project, detect which chain files are necessary based on config files, check they're valid and download them automatically. Could be done later, in an R script, if there is internet.
  curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
  curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi
cd -

# MR-MEGA
mkdir -p third_party/mr_mega/
cd       third_party/mr_mega/
if ! [ -f MR-MEGA ]; then
  curl -O https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip
  unzip MR-MEGA_v0.2.zip
  rm    MR-MEGA_v0.2.zip
  make
fi
cd -

# MsCAVIAR
# If you get errors while compiling, you might be missing these dependencies:
  # FlexiBLAS: https://gitlab.mpi-magdeburg.mpg.de/software/flexiblas-release
  # GSL: https://www.gnu.org/software/gsl/
cd third_party/MsCAVIAR/
if ! [ -f MsCAVIAR ]; then
  set +e
  make
  if [ $? -ne 0 ]; then # Didn't work, try Flexiblas.
    sed -i 's/-llapack -lblas/-lflexiblas/' Makefile # Edit makefile to use flexiblas (lol, cursed)
    make
    sed -i 's/-lflexiblas/-llapack -lblas/' Makefile # Undo, in case this script needs to be run again.
    set -e
  fi
fi
cd -

# PAINTOR
cd third_party/PAINTOR_V3.0/
if ! [ -f PAINTOR ]; then
  2to3 -wn PAINTOR_Utilities/AnnotateLocus.py
  chmod +x install.sh
  ./install.sh
fi
cd -

# TODO: major improvements to make for this install:
  # Use --no-index (in PIPFLAGS variable)
  # Maybe create new pyvenvs per-job, as apparently this can be faster? (See CCDB Python docs)
# PolyFun
cd third_party/polyfun/
if ! [ -d polyfun_py_env ]; then virtualenv polyfun_py_env; fi
source polyfun_py_env/bin/activate
set +e
python -c "import rpy2" # TODO: not a foolproof check but w/e
already_installed=$?
set -e
if [  $already_installed -ne 0 ]; then
  python -m pip install --upgrade pip
  python -m pip install numpy==1.25.1
  python -m pip install scipy==1.11.1
  python -m pip install scikit-learn==1.3.0
  python -m pip install pandas==2.0.3
  python -m pip install tqdm==4.65.0
  python -m pip install pyarrow==11.0.0
  python -m pip install bitarray==2.7.6
  python -m pip install networkx==3.1
  python -m pip install pandas-plink==2.2.9
  python -m pip install rpy2==3.5.7

  # TODO: still not working on BU SCC due to rpy2 issues
  set +e
  fix_rpy2.sh # Fix rpy2 bug on BU SCC :V
  set -e
fi
deactivate
cd -

# SuSiEx
cd third_party/SuSiEx/
if ! [ -f bin/SuSiEx ]; then
  cd src
  make all
fi
cd -

# CAVIARBF -- is an R package
# mJAM -- is an R package
