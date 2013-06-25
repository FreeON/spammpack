CHARMC=/usr/bin/charmc

CXX=$(CHARMC)
CXXFLAGS=$(OPTS)

BINARY=matmul

$(BINARY): matmul.o
	$(CHARMC) -o $@ $< $(OPTS) -lcblas

matmul.o: matmul.C matmul.decl.h
matmul.decl.h: matmul.ci
	$(CHARMC) $<

test:
	./charmrun ++local +p4 ./matmul 50 5

clean:
	rm *.o *.decl.h *.def.h $(BINARY) charmrun
