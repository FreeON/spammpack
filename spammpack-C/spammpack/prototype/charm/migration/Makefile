migration : migration.cc migration.decl.h
	charmc -language charm++ -balancer GreedyCommLB -o $@ $<

migration.decl.h : migration.ci
	charmc $<

run : migration
	./charmrun +p1 migration +LBDebug
	./charmrun +p4 migration +LBDebug

dist :
	tar cvf - Makefile migration.c? | bzip2 > migration.tar.bz2
