#!/usr/bin/env gnuplot

set terminal postscript eps color lw 3 "Helvetica, 22"
set output 'mustang-8.eps'

set title 'SpAMM 8192x256 on mustang'

set logscale xy

set xlabel '# nodes (1 PE per node)'
set ylabel 'walltime [s]'

set datafile missing '-'

set xrange [1:256]
set xtics ( "1" 1, "2" 2, "4" 4, "8" 8, "16" 16, "32" 32, \
  "64" 64, "128" 128, "256" 256 )

t1 = 844.924630

plot 'mustang-8.dat' using 1:2 with linespoints title 'SpAMM', \
  '' using 1:(t1/$1) with lines linecolor 0 title 'ideal'
