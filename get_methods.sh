#!/bin/bash

git submodule update --init

# MR-MEGA
mkdir ./third_party/mr_mega/
cd    ./third_party/mr_mega/
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
  sed -i 's/-llapack -lblas/-lflexiblas/' Makefile # Change makefile to use flexiblas (lol, cursed)
  make
fi
cd -

# PAINTOR
cd ./third_party/PAINTOR_V3.0/
if ! [ -f PAINTOR ]; then
  chmod u+x install.sh
  ./install.sh
fi
cd -

# SuSiEx -- nothing to do
# mJAM -- is an R package
