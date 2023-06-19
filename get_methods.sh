#!/bin/bash
set -e # exit if anything fails

git submodule update --remote

# fGWAS
cd ./third_party/fgwas
if ! [ -f src/fgwas ]; then
  ./configure
  make
fi
cd -

# FINEMAP
cd ./third_party/
if ! [ -f finemap/finemap_v1.4.2_x86_64 ]; then
  wget http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
  tar -zxf finemap_v1.4.2_x86_64.tgz
  rm       finemap_v1.4.2_x86_64.tgz
  mv       finemap_v1.4.2_x86_64 finemap
fi
cd -

# MR-MEGA
mkdir -p ./third_party/mr_mega/
cd       ./third_party/mr_mega/
if ! [ -f MR-MEGA ]; then
  wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip
  unzip MR-MEGA_v0.2.zip
  rm    MR-MEGA_v0.2.zip
  make
fi
cd -

# MsCAVIAR
# If you get errors while compiling, you might be missing these dependencies:
  # FlexiBLAS: https://gitlab.mpi-magdeburg.mpg.de/software/flexiblas-release
  # GSL: https://www.gnu.org/software/gsl/
cd ./third_party/MsCAVIAR/
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
cd ./third_party/PAINTOR_V3.0/
if ! [ -f PAINTOR ]; then
  chmod u+x install.sh
  ./install.sh
fi
cd -

# PolyFun
cd ./third_party/polyfun/
if ! [ -d polyfun_py_env ]; then virtualenv polyfun_py_env; fi
source polyfun_py_env/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy
python -m pip install scipy
python -m pip install sklearn
python -m pip install pandas
python -m pip install tqdm
python -m pip install pyarrow==11.0.0
python -m pip install bitarray
python -m pip install networkx
python -m pip install pandas-plink
python -m pip install rpy2
deactivate
cd -

# SuSiEx
cd ./third_party/SuSiEx/
if ! [ -f bin/SuSiEx ]; then
  cd src
  make all
fi
cd -

# CAVIARBF -- is an R package
# mJAM -- is an R package
