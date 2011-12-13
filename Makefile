.SUFFIXES:	.o .c .f .f90
#
.f90.o:   
#	ifort -c -g -O0 -traceback -fpp -DSPAMM_DOUBLE -DNON_ORTHOGONAL -I/home/matcha/freeon_beta/Modules  -I/home/matcha/freeon_beta/SCFeqs -I/home/matcha/FREEON_HOME_1/include/ $*.f90
	ifort -c  -O3 -fpp -I/home/matcha/freeon_beta/Modules  -I/home/matcha/freeon_beta/SCFeqs -I/home/matcha/FREEON_HOME_1/include/ $*.f90
#	ifort -c -g3 -O3 -fpp -DNON_ORTHOGONAL -I/home/matcha/freeon_beta/Modules  -I/home/matcha/freeon_beta/SCFeqs -I/home/matcha/FREEON_HOME_1/include/ $*.f90
#	ifort -c -g3 -O3 -fpp -DSPAMM_DOUBLE -DNON_ORTHOGONAL -I/home/matcha/freeon_beta/Modules  -I/home/matcha/freeon_beta/SCFeqs -I/home/matcha/FREEON_HOME_1/include/ $*.f90
#	gfortran -c -O3 -cpp -DSPAMM_SINGLE -I/home/matcha/freeon_beta/Modules  -I/home/matcha/freeon_beta/SCFeqs -I/home/matcha/FREEON_HOME_1/include/ $*.f90

.c.o:   
	gcc -c -g -O2 -Iipm-2.0 $*.c
#	icc -c -g3 -O3 -Iipm-2.0 $*.c
Objs = SpAMM_PACKAGE.o SpAMM_PROJECT.o SpAMM_CONVERT.o SpAMM_ALGEBRA.o SpAMM_MNGMENT.o SpAMM_GLOBALS.o SpAMM_TIMER.o SpAMM_DERIVED.o 
#
ipm:	
	cd ipm-2.0; ./configure ; cd ..
	make -C ipm-2.0 

test:	clean ipm $(Objs) SpAMM_TESTONE.o
	gfortran  -o spamm_test.exe SpAMM_TESTONE.o $(Objs) -Lipm-2.0/ -lipm -lm -L/home/matcha/LIBS/lapack_gcc_4.3.3/ -llapack -lblas

#pure:	clean ipm $(Objs) SpAMM_PURIFY2.o
#	gfortran  -o spamm_purify2.exe SpAMM_PURIFY2.o $(Objs) /home/matcha/freeon_beta/SCFeqs/DenMatMethods.o -Lipm-2.0/ -lipm -lm -I/home/matcha/freeon_beta/Modules  -I/home/matcha/FREEON_HOME_1/inc/ -L/home/matcha/LIBS/lapack_gcc_4.5.2//lib -L/home/matcha/LIBS/HDF5//lib -all-static  -llapack -lblas  -L/home/matcha/freeon_beta/Modules/.libs/ -L/home/matcha/freeon_beta/SCFeqs/  -all-static -lfreeonmodules -lz -lhdf5 -lpthread -lm -lm  -llapack -lblas

pure:	clean ipm $(Objs) SpAMM_PURIFY2.o
	ifort -g -Bstatic -o spamm_purify2.exe SpAMM_PURIFY2.o $(Objs) ../SCFeqs/DenMatMethods.o -Lipm-2.0/ -lipm -I../Modules  -I$(MONDO_HOME)/inc/ -L$(MONDO_HOME)/lib/ -L../Modules/.libs/ -L../SCFeqs/  -L../Modules/lapack/lapack/.libs/ -lfreeonmodules -lhdf5 -lz -lm -lfreeonlapack -mkl=sequential 
	cp spamm_purify2.exe $(MONDO_HOME)/bin/SP2 

#-L/home/matcha/LIBS/lapack_ifort

#pure:	clean ipm $(Objs) SpAMM_PURIFY2.o
#	gfortran  -all-static -o spamm_pure.exe SpAMM_PURIFY2.o $(Objs) /home/matcha/freeon_beta/SCFeqs/DenMatMethods.o -Lipm-2.0/ -lipm -lm -I/home/matcha/freeon_beta/Modules  -I/home/matcha/FREEON_HOME_1/inc/ -L/home/matcha/LIBS/LAPACK_GFORT_64  -llapack -lblas  -L/home/matcha/freeon_beta/Modules/.libs/ -lfreeonmodules  -L/home/matcha/FREEON_HOME_1/lib -lhdf5 -L/home/matcha/freeon_beta/SCFeqs/ -lz -lpthread -lm -lm  



#----------------------------------------------------------------------------
#	openf90 -openmp -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#----------------------------------------------------------------------------
#	ifort -openmp -g -shared-intel -shared-libgcc -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm -L/home/matcha/LIBS/lapack_gcc_4.3.3/ -llapack -lblas 
#-llal
#	ifort -openmp -g -diag-enable sc3 -o spamm_test.exe $(Objs) ##-Lipm-2.0/ -lipm -lm
#	ifort -openmp -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#	ifort  -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#----------------------------------------------------------------------------
#	pgf95 -mp -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#----------------------------------------------------------------------------
#	gfortran -fopenmp -g -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm
#	gfortran  -o spamm_test.exe $(Objs) -Lipm-2.0/ -lipm -lm -L/home/matcha/lapack-3.3.1/ -llapack -lblas
#----------------------------------------------------------------------------
ins:	clean ipm spamm
	rm -rf  /home/matcha/SPAMM/jnk
	mkdir   /home/matcha/SPAMM/jnk
	inspxe-cl -c mi3 -r /home/matcha/SPAMM/jnk -- /home/matcha/SPAMM/spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4 -1182.5623733386851D0
	inspxe-cl -report observations -result-dir /home/matcha/SPAMM/jnk

dbg:	all
#----------------------------------------------------------------------------
#	idb -args spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4 -1182.5623733386851D0
#----------------------------------------------------------------------------
#	pgdbg spamm_test.exe -c "run"  -program_args /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4  -1182.5623733386851D0
#----------------------------------------------------------------------------
	gdb  --eval-command=run --args spamm_test.exe  /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4  -1182.5623733386851D0
#----------------------------------------------------------------------------
#	pgdbg spamm_test.exe -text -c "run"  -program_args /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 2  -1182.5623733386851D0
#----------------------------------------------------------------------------

#--leak-check=yes 
val:	
	valgrind --tool=helgrind --suppressions=./gfortran_omp.supp ./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 1.0D-3 F  -1182.5623733386851D0
#	valgrind --gen-suppressions=yes --tool=helgrind --suppressions=./ifort_omp.supp ./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 1.0D-3 F  -1182.5623733386851D0
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
	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4 -1182.5623733386851D0 

scale:	
	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 1  -1182.5623733386851D0 >> spamm_scaling
	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 2  -1182.5623733386851D0 >> spamm_scaling
	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 3  -1182.5623733386851D0 >> spamm_scaling
	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 4  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 5  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 6  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 7  -1182.5623733386851D0 >> spamm_scaling
#	./spamm_test.exe /home/matcha/Desktop/RESEARCH/SCALING/TUBES_3_3/tube_33_132_20432_Geom#1_Base#3_Clone#1  612 732 54 8  -1182.5623733386851D0 >> spamm_scaling

