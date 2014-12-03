#!/usr/bin/gnuplot

set terminal postscript color
set output 'diagonal.ps'

set title 'GreedyCommLB performance, 32 nodes, 4096x128 diagonal matrix'

set logscale y

set xlabel 'multiply iteration'
set ylabel 'walltime [s]'

plot for [COL=2:3] 'diagonal.dat' using 1:COL with linespoints title columnheader
