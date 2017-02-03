#!/bin/sh

cd preprocess
./preprocess.sh
cd ..
cd src
./ebm.sh
cd ..
cd postprocess
./postprocess.sh
cd ..
rm -rf *~
