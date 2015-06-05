# vim:syntax=gnuplot

set terminal eps color

set xlabel 'N_basic'
set ylabel 'walltime [s]'

set xtics (16, 32, 64, 128)
set logscale x

set label 1 'L1 - 128 KiB' at 32, 40

set key left
set title 'Single thread performance'
set output 'scan-serial-freeon.eps'
plot 'scan.dat' index 0 using 2:3 with linespoints title 'N_chunk = 2048', \
     '' index 1 using 2:3 with linespoints title 'N_chunk = 1024', \
     '' index 2 using 2:3 with linespoints title 'N_chunk = 512', \
     '' index 3 using 2:3 with linespoints title 'N_chunk = 256', \
     '' index 4 using 2:3 with linespoints title 'N_chunk = 128'

set key center top
set title 'Performance on 48 threads'
set output 'scan-48-freeon.eps'
plot 'scan.dat' index 0 using 2:4 with linespoints title 'N_chunk = 2048', \
     '' index 1 using 2:4 with linespoints title 'N_chunk = 1024', \
     '' index 2 using 2:4 with linespoints title 'N_chunk = 512', \
     '' index 3 using 2:4 with linespoints title 'N_chunk = 256', \
     '' index 4 using 2:4 with linespoints title 'N_chunk = 128'

unset label 1

set key left top
set title 'Speedup on 48 threads'
set ylabel 'speedup'
set output 'scan-speedup-48-freeon.eps'
plot 'scan.dat' index 0 using 2:5 with linespoints title 'N_chunk = 2048', \
     '' index 1 using 2:5 with linespoints title 'N_chunk = 1024', \
     '' index 2 using 2:5 with linespoints title 'N_chunk = 512', \
     '' index 3 using 2:5 with linespoints title 'N_chunk = 256', \
     '' index 4 using 2:5 with linespoints title 'N_chunk = 128'

set key left top
set title 'Speedup on 48 threads'
set ylabel 'parallel efficiency [%]'
set output 'scan-efficiency-48-freeon.eps'
plot 'scan.dat' index 0 using 2:($5/48*100) with linespoints title 'N_chunk = 2048', \
     '' index 1 using 2:($5/48*100) with linespoints title 'N_chunk = 1024', \
     '' index 2 using 2:($5/48*100) with linespoints title 'N_chunk = 512', \
     '' index 3 using 2:($5/48*100) with linespoints title 'N_chunk = 256', \
     '' index 4 using 2:($5/48*100) with linespoints title 'N_chunk = 128'

unset logscale
set xtics autofreq

set xlabel 'tree depth'
set key right
set title 'Single thread performance'
set ylabel 'walltime [s]'
set output 'scan-serial-tiers-freeon.eps'
plot 'scan.dat' index 0 using (log($1/$2)/log(2)):3 with linespoints title 'N_chunk = 2048', \
     '' index 1 using (log($1/$2)/log(2)):3 with linespoints title 'N_chunk = 1024', \
     '' index 2 using (log($1/$2)/log(2)):3 with linespoints title 'N_chunk = 512', \
     '' index 3 using (log($1/$2)/log(2)):3 with linespoints title 'N_chunk = 256', \
     '' index 4 using (log($1/$2)/log(2)):3 with linespoints title 'N_chunk = 128'

set key left bottom
set title 'Performance on 48 threads'
set output 'scan-48-tiers-freeon.eps'
set logscale y
plot 'scan.dat' index 0 using (log($1/$2)/log(2)):4 with linespoints title 'N_chunk = 2048', \
     '' index 1 using (log($1/$2)/log(2)):4 with linespoints title 'N_chunk = 1024', \
     '' index 2 using (log($1/$2)/log(2)):4 with linespoints title 'N_chunk = 512', \
     '' index 3 using (log($1/$2)/log(2)):4 with linespoints title 'N_chunk = 256', \
     '' index 4 using (log($1/$2)/log(2)):4 with linespoints title 'N_chunk = 128'
