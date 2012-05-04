#alias cfgcc 'autoreconf -fis; ./configure --prefix=/home/mchalla/FREEON_BETA_1 F77=gfortran  FC=gfortran CFLAGS="-O3 -ffast-math -ftree-vectorize" FFLAGS="-O3 -ffast-math -ftree-vectorize" FCFLAGS="-O3 -ffast-math -ftree-vectorize"'
#alias cfifc  'autoreconf -i -s; ./configure --disable-random-seed --disable-delete-preprocessed-files --prefix=/home/mchalla/FREEON_BETA_INTEL --enable-static --with-MONDO_SCRATCH=/scratch CC=icc FC=ifort F77=ifort FCFLAGS="-O1 -no-prec-div -static -V -fpe0" FFLAGS="-O1-no-prec-div -static -V -fpe0" '
#alias cfpgi  'autoreconf -i -s; ./configure --build=x86_64-unknown-linux-gnu --disable-random-seed --disable-delete-preprocessed-files --prefix=/home/mchalla/FREEON_BETA_PGI  --enable-phipac --enable-static --with-MONDO_SCRATCH=/scratch CC=pgcc FC=pgf90 F77=pgf90 FCFLAGS="-O2 -fast -fastsse " FFLAGS="-O2 -fast -fastsse " '


FREEON_HOME=/home/mchalla/FREEON_BETA_PGI
FREEON_SRCE=/home/mchalla/freeon_beta_pgi

#FREEON_HOME=/home/mchalla/FREEON_BETA_INTEL
#FREEON_SRCE=/home/mchalla/freeon_beta_intel

.SUFFIXES:	.o .c .f .f90
#
.f90.o:   
#	gfortran -c -O3 -ffast-math -ftree-vectorize -cpp -I$(FREEON_SRCE)/Modules  -I$(FREEON_SRCE)/SCFeqs -I$(FREEON_HOME)/include/ $*.f90
	pgf95 -c -O4 -fast -fastsse -Mpreprocess -I$(FREEON_SRCE)/Modules  -I$(FREEON_SRCE)/SCFeqs -I$(FREEON_HOME)/include/ $*.f90
#	ifort -c  -O3 -fpp -I$(FREEON_SRCE)/Modules  -I$(FREEON_SRCE)/SCFeqs -I$(FREEON_HOME)/include/ $*.f90
#	ifort -c -g -O0 -traceback -fpp -I$(FREEON_SRCE)/Modules  -I$(FREEON_SRCE)/SCFeqs -I$(FREEON_HOME)/include/ $*.f90

.c.o:   
#	gcc -c   -O3 -Iipm-2.0 $*.c
#	icc -c   -O3 -Iipm-2.0 $*.c
	pgcc -c  -O3 -Iipm-2.0 $*.c

Objs = SpAMM_PACKAGE.o SpAMM_PROJECT.o SpAMM_CONVERT.o SpAMM_ALGEBRA.o SpAMM_MNGMENT.o SpAMM_GLOBALS.o SpAMM_TIMER.o SpAMM_DERIVED.o 
#
ipm:	
	cd ipm-2.0; ./configure ; cd ..
	make -C ipm-2.0 

test:	clean ipm $(Objs) SpAMM_TESTONE.o
	gfortran  -o spamm_test.exe SpAMM_TESTONE.o $(Objs) -Lipm-2.0/ -lipm -lm -L/home/mchalla/LIBS/lapack_gcc_4.3.3/ -llapack -lblas

pure:	clean ipm $(Objs) SpAMM_PURIFY2.o
	cp $(FREEON_HOME)/lib/libhdf5.a $(FREEON_HOME)/lib/libfreeonhdf5.a
#	gfortran  -static-libgfortran -o spamm_purify2.exe SpAMM_PURIFY2.o $(Objs) ../SCFeqs/DenMatMethods.o -Lipm-2.0/ -lipm -I../Modules  -I$(FREEON_HOME)/inc/ -L$(FREEON_HOME)/lib/ -L../Modules/.libs/ -L../SCFeqs/  -L../lapack/lapack/.libs/ -L../lapack/blas/.libs/  -L../lapack/install/.libs/ -lfreeonmodules -lfreeonhdf5 -lz -lm -lfreeonlapack -lfreeonblas -lfreeonlapackextra
	pgf95 -Bstatic -o spamm_purify2.exe SpAMM_PURIFY2.o $(Objs) ../SCFeqs/DenMatMethods.o -Lipm-2.0/ -lipm -I../Modules  -I$(FREEON_HOME)/inc/ -L$(FREEON_HOME)/lib/ -L../Modules/.libs/ -L../SCFeqs/  -L../lapack/lapack/.libs/ -L../lapack/blas/.libs/  -L../lapack/install/.libs/ -lfreeonmodules -lhdf5 -lz -lm -lfreeonlapack -lfreeonblas -lfreeonlapackextra
#	ifort  -Bstatic -o spamm_purify2.exe SpAMM_PURIFY2.o $(Objs) ../SCFeqs/DenMatMethods.o -Lipm-2.0/ -lipm -I../Modules  -I$(FREEON_HOME)/inc/ -L$(FREEON_HOME)/lib/ -L../Modules/.libs/ -L../SCFeqs/  -L../lapack/lapack/.libs/ -L../lapack/blas/.libs/  -L../lapack/install/.libs/ -lfreeonmodules -lhdf5 -lz -lm -lfreeonlapack -lfreeonblas -lfreeonlapackextra
#	cp spamm_purify2.exe $(FREEON_HOME)/bin/SP2 

#----------------------------------------------------------------------------
#	openf90 -openmp -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#----------------------------------------------------------------------------
#	ifort -openmp -g -shared-intel -shared-libgcc -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm -L/home/mchalla/LIBS/lapack_gcc_4.3.3/ -llapack -lblas 
#-llal
#	ifort -openmp -g -diag-enable sc3 -o spamm_test.exe $(Objs) ##-Lipm-2.0/ -lipm -lm
#	ifort -openmp -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#	ifort  -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#----------------------------------------------------------------------------
#	pgf95 -mp -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#----------------------------------------------------------------------------
#	gfortran -fopenmp -g -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#	gfortran  -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm -L/home/mchalla/lapack-3.3.1/ -llapack -lblas
#----------------------------------------------------------------------------
ins:	clean ipm spamm
	rm -rf  /home/mchalla/SPAMM/jnk
	mkdir   /home/mchalla/SPAMM/jnk
	inspxe-cl -c mi3 -r /home/mchalla/SPAMM/jnk -- /home/mchalla/SPAMM/spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4 -1182.5623733386851D0
	inspxe-cl -report observations -result-dir /home/mchalla/SPAMM/jnk

dbg:	all
#----------------------------------------------------------------------------
#	idb -args spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4 -1182.5623733386851D0
#----------------------------------------------------------------------------
#	pgdbg spamm_test.exe -c "run"  -program_args /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4  -1182.5623733386851D0
#----------------------------------------------------------------------------
	gdb  --eval-command=run --args spamm_test.exe  /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4  -1182.5623733386851D0
#----------------------------------------------------------------------------
#	pgdbg spamm_test.exe -text -c "run"  -program_args /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 2  -1182.5623733386851D0
#----------------------------------------------------------------------------

#--leak-check=yes 
val:	
	valgrind --tool=helgrind --suppressions=./gfortran_omp.supp ./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 1.0D-3 F  -1182.5623733386851D0
#	valgrind --gen-suppressions=yes --tool=helgrind --suppressions=./ifort_omp.supp ./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 1.0D-3 F  -1182.5623733386851D0
# http://valgrind.org/docs/manual/hg-manual.html

#--gen-suppressions=yes 
#-----------------------------------------------------------------------
clean:	
	make -C ipm-2.0 clean
	rm -f \#*
	rm -f *.mod *.o
	rm -f *~
	rm -f *.exe
	rm -rf jnk
	ls -lS

tarball:clean
	tar -cvf SpAMM_`date '+%B'`_`date '+%d'`_`date +%y`.tar MISC ipm-2.0 \
	SpAMM_DERIVED.f90  SpAMM_GLOBALS.f90  SpAMM_MNGMENT.f90   \
	SpAMM_ALGEBRA.f90  SpAMM.f90 Makefile *.supp *.f90 *.c *.ps *.pdf ;\
	gzip SpAMM_`date '+%B'`_`date '+%d'`_`date +%y`.tar     ;\
#----------------------------------------------------------------------- 
SpAMM_PURIFY2.o:        SpAMM_DERIVED.o  SpAMM_GLOBALS.o  SpAMM_MNGMENT.o  SpAMM_ALGEBRA.o SpAMM_CONVERT.o SpAMM_PROJECT.o SpAMM_PACKAGE.o 
SpAMM_TESTONE.o:        SpAMM_DERIVED.o  SpAMM_GLOBALS.o  SpAMM_MNGMENT.o  SpAMM_ALGEBRA.o SpAMM_CONVERT.o SpAMM_PROJECT.o SpAMM_PACKAGE.o 
SpAMM_PACKAGE.o:        SpAMM_DERIVED.o  SpAMM_GLOBALS.o  SpAMM_MNGMENT.o  SpAMM_ALGEBRA.o SpAMM_CONVERT.o SpAMM_PROJECT.o  
SpAMM_PROJECT.o: 	SpAMM_DERIVED.o  SpAMM_GLOBALS.o  SpAMM_MNGMENT.o  SpAMM_ALGEBRA.o SpAMM_CONVERT.o
SpAMM_CONVERT.o:	SpAMM_DERIVED.o  SpAMM_GLOBALS.o  SpAMM_MNGMENT.o  SpAMM_ALGEBRA.o 
SpAMM_ALGEBRA.o:	SpAMM_DERIVED.o  SpAMM_GLOBALS.o  SpAMM_MNGMENT.o    
SpAMM_MNGMENT.o:	SpAMM_DERIVED.o  SpAMM_GLOBALS.o  
SpAMM_GLOBALS.o:	SpAMM_DERIVED.o


test1:	test
	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4 -1182.5623733386851D0 

scale:	
	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 1  -1182.5623733386851D0 >> spamm_scaling
	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 2  -1182.5623733386851D0 >> spamm_scaling
	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 3  -1182.5623733386851D0 >> spamm_scaling
	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 5  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 6  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 7  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/mchalla/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 8  -1182.5623733386851D0 >> spamm_scaling

tags : *.f90
	ctags --fortran-kinds=+i+L *.f90
