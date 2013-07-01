# vim: syntax=gnuplot

set terminal postscript color

set title 'parallel OMP scaling, AMD Opteron 4x12 core'

set xlabel 'cores'

f(x) = x

t0_intel = 72.666956
t0_gcc_32 = 77.539708
t0_gcc_64 = 604.056920

set ylabel 'relative speedup to single core'

set output 'freeon-10.scaling.ps'
plot f(x) with lines title 'ideal', \
     'freeon-10.dat' index 0 using 1:(t0_intel/$2) with linespoints title 'Intel (2048x32)', \
     'freeon-10.dat' index 1 using 1:(t0_gcc_32/$2) with linespoints title 'gcc 4.6.3 (2048x32)', \
     'freeon-10.dat' index 2 using 1:(t0_gcc_64/$2) with linespoints title 'gcc 4.6.3 (4096x64)'

set ylabel 'walltiem [s]'

set logscale xy
set output 'freeon-10.performance.ps'
plot for [COL=2:4] 'freeon-10.dat' index 0 using 1:COL with linespoints title 'Intel', \
     for [COL=2:4] 'freeon-10.dat' index 1 using 1:COL with linespoints title 'gcc 4.6.3', \
     '' using 1:(61.502654) with lines title 'Intel - serial', \
     '' using 1:(22.347782) with lines title 'gcc - serial'
