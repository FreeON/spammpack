#!/usr/bin/gnuplot

set terminal postscript color
set output 'timing.ps'

set xlabel 'core count'
set ylabel 'relative speedup / ideal'

t_1 = 61.366154

plot 'timing.dat' using 1:(t_1/$2/$1) with linespoints title 'freeon-10'
