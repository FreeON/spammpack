#CHARMC = ~/charm/bin/charmc
CHARMC = charmc

reductions : reductions.cc reductions.decl.h
	$(CHARMC) -language charm++ -module CommonLBs -o $@ $< -lrt

%.decl.h : %.ci
	$(CHARMC) $^

run : reductions
	for (( i = 1; i <= 48; i++ )); do ./charmrun +p$$i reductions; done

clean :
	rm -f reductions reductions.decl.h reductions.def.h charmrun
