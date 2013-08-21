CHARMC = charmc
#CHARMC = /var/tmp/charm-6.5.1/bin/charmc

migration : migration.cc migration.decl.h
	$(CHARMC) -language charm++ -module CommonLBs -balancer GreedyCommLB -o $@ $< -lm

migration.decl.h : migration.ci
	$(CHARMC) $<

run-serial : migration
	./charmrun +p1 migration +LBDebug

run : migration
	./charmrun +p4 migration +LBDebug

dist :
	tar cvf - Makefile migration.c? | bzip2 > migration.tar.bz2

clean :
	rm -f migration.{decl,def}.h charmrun migration
