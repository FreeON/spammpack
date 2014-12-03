#CHARMC = charmc
CHARMC = ~/charm-6.5.1/bin/charmc

openmp : openmp.cc openmp.decl.h
	$(CHARMC) -language charm++ -fopenmp -o $@ $<

%.decl.h : %.ci
	$(CHARMC) $<

run : openmp
	OMP_NUM_THREADS=2 ./openmp
	./charmrun +p2 openmp

clean :
	rm -f charmrun openmp openmp.decl.h openmp.def.h

dist :
	tar cvf - openmp.c{c,i} Makefile | bzip2 > openmp.tar.bz2
