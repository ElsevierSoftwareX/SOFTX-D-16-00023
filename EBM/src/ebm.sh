#!/bin/sh
# This script compiles and runs the netcdf version Energy Balance Model
# See compiling detail in Makefile

make
./EBM

rm -rf *~
