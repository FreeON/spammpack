#EXTRA_CPPFLAGS = -DWITH_DATA

migration : migration.cc migration.decl.h
	charmc -language charm++ -balancer GreedyCommLB $(EXTRA_CPPFLAGS) -o $@ $< -lm

migration.decl.h : migration.ci
	charmc -preprocess $<

run-serial : migration
	./charmrun +p1 migration +LBDebug

run : migration
	./charmrun +p4 migration +LBDebug

dist :
	tar cvf - Makefile migration.c? | bzip2 > migration.tar.bz2
