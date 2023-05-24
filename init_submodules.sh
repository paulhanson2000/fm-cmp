#!/bin/bash

git submodule update --init

# MR-MEGA
mkdir ./third_party/mr_mega/
cd    ./third_party/mr_mega/
wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip
unzip MR-MEGA_v0.2.zip
rm    MR-MEGA_v0.2.zip
make
cd -

# PAINTOR
cd ./third_party/PAINTOR_V3.0/
chmod u+x install.sh
./install.sh
cd -

# SuSiEx -- nothing to do
