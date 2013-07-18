#!/usr/bin/env gnuplot

set terminal postscript color
set output 'timing.ps'

set title 'matmul 2048x64'

set logscale xy

set xlabel '# cores'
set ylabel 'walltime [s]'

set datafile missing '-'

set xrange [1:48]
set xtics ( "1" 1, "2" 2, "4" 4, "8" 8, "16" 16, "32" 32, "48" 48 )

t1 = 22.687797

plot 'timing.dat' using 1:2 with linespoints title 'matmul', \
  '' using 1:(t1/$1) with lines linecolor 0 notitle
