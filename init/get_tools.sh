#!/bin/bash

git submodule update --init --remote # For tools hosted on GitHub, download their source

# Detect if on Compute Canada. If so, use pre-installed Python things. 
PIPFLAGS=""
uname -n | grep calculquebec
if [ $? -eq 0 ]; then PIPFLAGS="--no-index"; fi

# Detect OS & CPU architecture
if   [ "$(uname)" = "Linux"  ] && [ "$(uname -m)" = "x86_64" ]; then OS="Linux_x86"
elif [ "$(uname)" = "Darwin" ] && [ "$(uname -m)" = "x86_64" ]; then OS="Mac_x86"
elif [ "$(uname)" = "Darwin" ] && [ "$(uname -m)" = "arm64"  ]; then OS="Mac_ARM"
else echo "Only Mac x86_64 or ARM64, or Linux x86_64 are supported. Some required third-party programs would fail to run on your machine."; exit 1
fi

# Unfortunately, sed works differently on Linux compared to Mac/BSD. Mac/BSD's sed -i requires an additional extension to create a backup of the edited file.
sed_i() {
  if   [ $OS = "Linux_x86" ]; then sed -i "$@"    #  Linux  sed -i
  else                             sed -i '' "$@" # Mac/BSD sed -i
  fi
}

# Exit if anything fails
set -e 

# fgwas
cd third_party/fgwas/
if ! [ -f src/fgwas-5000_iter ]; then
  autoreconf -f -i # Because https://stackoverflow.com/a/33286344
  ./configure

  # NOTE: cursed editing of fgwas's source to speed up by reducing the iterations before concluding a model doesn't converge.
  # Don't need as many iterations to for the first pass of models with only one annotation, so build separate versions of fgwas.
  sed_i -e 's/iter <5000/iter<100/' -e 's/iter > 4999/iter>99/' src/SNPs.cpp 
  make
  mv src/fgwas src/fgwas-100_iter
  make clean
  sed_i -e 's/iter<1000/iter<1000/' -e 's/iter>999/iter>999/' src/SNPs.cpp
  make
  mv src/fgwas src/fgwas-1000_iter
  make clean
  sed_i -e 's/iter<1000/iter <5000/' -e 's/iter>999/iter > 4999/' src/SNPs.cpp # (back to how it was originally)
  make
  mv src/fgwas src/fgwas-5000_iter
  make clean
fi
cd -

# FINEMAP
cd third_party/
if ! [ -d finemap ]; then
  TMP=""
  if   [ $OS = "Linux_x86" ]; then TMP=finemap_v1.4.2_x86_64
  else                             TMP=finemap_v1.4.2_MacOSX
  fi
  curl -O http://www.christianbenner.com/${TMP}.tgz
  tar -zxf ${TMP}.tgz
  rm       ${TMP}.tgz
  mv       ${TMP} finemap
fi
cd -

# liftOver
mkdir -p third_party/liftover/
cd       third_party/liftover/
if ! [ -f liftOver ]; then
  if   [ $OS = "Linux_x86" ]; then curl -O https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
  elif [ $OS =   "Mac_x86" ]; then curl -O https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/liftOver
  elif [ $OS =   "Mac_ARM" ]; then curl -O https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.arm64/liftOver
  fi
  chmod +x liftOver
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
  if [ $(uname) == "Darwin" ]; then sed_i 's/g++$/g++-13 -I\/opt\/homebrew\/opt\/gsl\/include -L\/opt\/homebrew\/opt\/gsl\/lib/' Makefile; fi # TODO Homebrew gcc & gsl
  set +e
  make
  if [ $? -ne 0 ]; then # Didn't work, try Flexiblas.
    sed_i 's/-llapack -lblas/-lflexiblas/' Makefile # Edit makefile to use flexiblas (lol, cursed)
    make
    sed_i 's/-lflexiblas/-llapack -lblas/' Makefile # Undo, in case this script needs to be run again.
    set -e
  fi
fi
cd -

# PAINTOR
cd third_party/PAINTOR_V3.0/
if ! [ -f PAINTOR ]; then
  if [ $(uname) == "Darwin" ]; then sed_i 's/g++$/g++-13/' Makefile; fi # TODO Homebrew gcc

  # PAINTOR ships with an old version of Eigen. Newer version has better performance!
  rm -r eigen
  git clone https://gitlab.com/libeigen/eigen.git
  sed_i 's/causal_config_bit_vector(causal_index/causal_config_bit_vector.coeffRef(causal_index/' Functions_model.cpp
  sed_i 's/-std=c\+\+11/-std=c\+\+17/' Makefile

  2to3 -wn PAINTOR_Utilities/AnnotateLocus.py
  chmod +x install.sh
  ./install.sh
fi
cd -

# PolyFun
# TODO: major improvements to make for this install:
  # Maybe create new pyvenvs per-job, as apparently this can be faster? (See CCDB Python docs)
cd third_party/polyfun/
if ! [ -d polyfun_py_env ]; then virtualenv polyfun_py_env; fi
source polyfun_py_env/bin/activate
set +e
python -c "import rpy2" # TODO: not a foolproof check but w/e
already_installed=$?
set -e
if [  $already_installed -ne 0 ]; then
  python -m pip install --upgrade pip
  python -m pip install ${PIPFLAGS} numpy==1.25.1
  python -m pip install ${PIPFLAGS} scipy==1.11.1
  python -m pip install ${PIPFLAGS} scikit-learn==1.3.0
  python -m pip install ${PIPFLAGS} pandas==2.0.3
  python -m pip install ${PIPFLAGS} tqdm==4.65.0
  python -m pip install ${PIPFLAGS} pyarrow==11.0.0
  python -m pip install ${PIPFLAGS} bitarray==2.7.6
  python -m pip install ${PIPFLAGS} networkx==3.1
  python -m pip install ${PIPFLAGS} pandas-plink==2.2.9
  python -m pip install ${PIPFLAGS} rpy2==3.5.7

  # TODO: still not working on BU SCC due to rpy2 issues
  set +e
  fix_rpy2.sh # Fix rpy2 bug on BU SCC via a script by the sysadmins :V
  set -e
fi
deactivate
cd -

# SuSiEx
cd third_party/SuSiEx/
if ! [ -f bin/SuSiEx ]; then
  cd src
  if [ $(uname) == "Darwin" ]; then sed_i 's/-fopenmp/-Xclang -fopenmp -I\/opt\/homebrew\/opt\/libomp\/include -L\/opt\/homebrew\/opt\/libomp\/lib -lomp/' makefile; fi # Need extra stuff Mac for clang to get Apple clang to do OpenMP: see https://mac.r-project.org/openmp/ and https://github.com/Rdatatable/data.table/wiki/Installation
  make all
fi
cd -

# CAVIARBF -- is an R package
# mJAM -- is an R package
