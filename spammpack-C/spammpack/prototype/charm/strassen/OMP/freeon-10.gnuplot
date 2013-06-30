# vim: syntax=gnuplot

set terminal postscript color

set title 'parallel OMP scaling, AMD Opteron 4x12 core'

set xlabel 'cores'

f(x) = x

t0_intel = 72.666956
t0_gcc = 77.539708

set ylabel 'relative speedup to single core'

set output 'freeon-10.scaling.ps'
plot f(x) with lines title 'ideal', \
     'freeon-10.dat' index 0 using 1:(t0_intel/$2) with linespoints title 'Intel', \
     'freeon-10.dat' index 1 using 1:(t0_gcc/$2) with linespoints title 'gcc 4.6.3'

set ylabel 'walltiem [s]'

set logscale xy
set output 'freeon-10.performance.ps'
plot 'freeon-10.dat' index 0 using 1:2 with linespoints title 'Intel', \
     'freeon-10.dat' index 1 using 1:2 with linespoints title 'gcc 4.6.3', \
     '' using 1:(61.502654) with lines title 'Intel - serial', \
     '' using 1:(22.347782) with lines title 'gcc - serial'
