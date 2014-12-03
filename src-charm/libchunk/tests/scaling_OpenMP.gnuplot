# vim: syntax=gnuplot

set terminal postscript color
set output 'scaling_with_tolerance.ps'

plot '< ./process_scaling_data.py --x-axis tolerance scaling_OpenMP.dat' index 0 using 1:2
