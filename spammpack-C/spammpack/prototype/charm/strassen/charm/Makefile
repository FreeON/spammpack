.SUFFIXES : .o .cc .ci

CXX = charmc
CXXFLAGS = -O3 -march=native -g -Wall

%.decl.h %.def.h : %.ci
	$(CXX) $<

.cc.o :
	$(CXX) -c -o $@ $(CXXFLAGS) $<

INTERFACEFILES = strassen.ci

GENERATEDHEADERFILES = $(foreach i, $(INTERFACEFILES), $(basename $i).def.h $(basename $i).decl.h)

HEADERFILES = strassen.h

SOURCEFILES = \
	Matrix.cc \
	Node.cc \
	Timer.cc \
	strassen.cc

OBJ = $(foreach i, $(SOURCEFILES), $(basename $i).o)

strassen : $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $^

$(OBJ) : $(HEADERFILES) $(GENERATEDHEADERFILES)

run :
	for P in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		bash -c "OMP_NUM_THREADS=$$P ./strassenOMP -N 2048 -b 32"; \
	done

run-serial :
	./strassenOMP -N 2048 -b 32

clean :
	rm -f $(OBJ)
	rm -f $(GENERATEDHEADERFILES)
	rm -f strassen charmrun
