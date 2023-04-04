#!/bin/bash
cd ~/spciaa/iga-ads
cmake --build build -j 9
cd build/examples
cp polution pol/exec
cd pol/data
../exec/polution
cd ..
gnuplot animate.gnuplot
mv animate_tmp.gif gifs
