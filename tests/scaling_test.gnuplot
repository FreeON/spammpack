# vim: syntax=gnuplot

serial_intel = 19.8928103447
serial_gcc = 1.3729679585

set title '1024x1024 dense random'

set terminal postscript color
set output 'scaling_test.speedup.ps'

set xlabel '# Threads'
set ylabel 'Speedup'

set key left

f(x) = x

plot 'scaling_test.intel.dat' using 1:(serial_intel/$2) title 'SpAMM (Intel)', \
       'scaling_test.gcc.dat' using 1:(serial_gcc/$2) title 'SpAMM (gcc)', \
       f(x) title 'ideal'

set output 'scaling_test.time.ps'

set ylabel 'walltime [s]'

plot 'scaling_test.intel.dat' using 1:2 with linespoints title 'SpAMM (Intel)', \
       'scaling_test.gcc.dat' using 1:2 with linespoints title 'SpAMM (gcc)'
