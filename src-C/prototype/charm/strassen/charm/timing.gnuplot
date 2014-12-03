#!/usr/bin/gnuplot

set terminal postscript color
set output 'timing.ps'

set title 'Futures based multiply'

set logscale xy

set xlabel '# cores'
set ylabel 'walltime [s]'

full_2048_128_t1  = 40.368775
full_1024_64_t1   = 22.848918
full_1024_64_t2   = 3.197707
charm_1024_64_t1  = 0.371398
charm_2048_128_t1 = 1.320883

set datafile missing '-'

set xrange [1:48]
set xtics ( "1" 1, "2" 2, "4" 4, "8" 8, "16" 16, "32" 32, "48" 48 )

plot for [COL=2:3] 'timing.dat' using 1:COL with linespoints title columnheader, \
  '' using 1:($2-$3) with linespoints title 'work', \
  '' using 1:(full_2048_128_t1/$1) with lines linecolor 0 notitle, \
  '' using 1:(charm_2048_128_t1/$1) with lines linecolor 0 notitle

  #'' using 1:(full_1024_64_t1/$1) with lines title 'ideal - 1024x64', \
  #'' using 1:(full_1024_64_t2/$1) with lines title 'ideal cpuaffinity - 1024x64'
  #'' using 1:(charm_1024_64_t1/$1) with lines title 'charm - 1024x64'
