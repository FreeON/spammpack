#!/usr/bin/gnuplot

set terminal postscript color
set output 'timing.ps'

set title 'Futures based multiply'

set logscale xy

set xlabel '# cores'
set ylabel 'walltime [s]'

full_2048_128_t1 = 233.995228
full_1024_64_t1  = 22.848918
full_1024_64_t2  = 22.852871
charm_1024_64_t1 = 0.371398

set datafile missing '-'

set xrange [1:48]
set xtics ( "1" 1, "2" 2, "4" 4, "8" 8, "16" 16, "32" 32, "48" 48 )

plot for [COL=2:4] 'timing.dat' using 1:COL with linespoints title columnheader, \
  '' using 1:(full_2048_128_t1/$1) with lines title 'ideal - 2048x128', \
  '' using 1:(full_1024_64_t1/$1) with lines title 'ideal - 1024x64'

  #'' using 1:(charm_1024_64_t1/$1) with lines title 'charm - 1024x64'
