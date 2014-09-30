# vim: syntax=gnuplot

serial_intel_16_1024 = 19.8928103447
serial_intel_32_2048 = 145.9683935642
serial_intel_64_4096 = 1148.1845097542
serial_gcc_4_4_5 = 1.3729679585
serial_gcc_4_7_0 = 1.4327812195

set title 'dense random'

set terminal postscript color
set output 'scaling_test.speedup.ps'

set xlabel '# Threads'
set ylabel 'Speedup'

set key left

f(x) = x

plot 'scaling_test.intel.16.1024.dat' using 1:(serial_intel_16_1024/$2) title 'SpAMM (Intel, 16, 1024)', \
       'scaling_test.intel.32.2048.dat' using 1:(serial_intel_32_2048/$2) title 'SpAMM (Intel, 32, 2048)', \
       'scaling_test.intel.64.4096.dat' using 1:(serial_intel_64_4096/$2) title 'SpAMM (Intel, 64, 4096)', \
       'scaling_test.gcc-4.4.5.dat' using 1:(serial_gcc_4_4_5/$2) title 'SpAMM (gcc-4.4.5)', \
       'scaling_test.gcc-4.7.0.dat' using 1:(serial_gcc_4_7_0/$2) title 'SpAMM (gcc-4.7.0)', \
       f(x) lc 7 title 'ideal'

set output 'scaling_test.time.ps'

set ylabel 'walltime [s]'

plot 'scaling_test.intel.16.1024.dat' using 1:2 with linespoints title 'SpAMM (Intel, 16, 1024)', \
       'scaling_test.gcc-4.4.5.dat' using 1:2 with linespoints title 'SpAMM (gcc-4.4.5)', \
       'scaling_test.gcc-4.7.0.dat' using 1:2 with linespoints title 'SpAMM (gcc-4.7.0)'
