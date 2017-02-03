#!/bin/sh

# This step processes the summer temperature in the NH
ferret -script summer.jnl

# This step plot the summer temperature overlay on the basemap
ferret -script figure_demo.jnl
Fprint -o demo.ps demo.plt
ps2pdfwr demo.ps demo.pdf

rm -rf *.plt
rm -rf *~
rm -rf *.ps
