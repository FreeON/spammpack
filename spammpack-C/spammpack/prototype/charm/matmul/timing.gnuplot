#!/usr/bin/env gnuplot

set terminal postscript color
set output 'timing.ps'

set title 'mustang - matmul 8192x256'

set logscale xy

set xlabel '# nodes (1 core per node)'
set ylabel 'walltime [s]'

set datafile missing '-'

set xrange [1:48]
set xtics ( "1" 1, "2" 2, "4" 4, "8" 8, "16" 16, "32" 32, "48" 48 )

t0 = 22.925205
t1 = 283.885821

plot for [INDEX=0:2] 'mustang.dat' index INDEX using 1:2 with linespoints title columnheader, \
  '' using 1:(t0/$1) with lines linecolor 0 notitle, \
  '' using 1:(t1/$1) with lines linecolor 0 notitle
