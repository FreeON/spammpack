.SUFFIXES : .o .cc .ci

CXX = charmc
CXXFLAGS = -O3 -march=native -g -Wall
CTAGS = ctags
LIBS = -lrt

%.decl.h %.def.h : %.ci
	$(CXX) $<

.cc.o :
	$(CXX) -c -o $@ $(CXXFLAGS) $<

INTERFACEFILES = \
	Matrix.ci \
	Messages.ci \
	Node.ci \
	strassen.ci

GENERATEDHEADERFILES = $(foreach i, $(INTERFACEFILES), $(basename $i).def.h $(basename $i).decl.h)

HEADERFILES = \
	Matrix.h \
	Messages.h \
	Node.h \
	Timer.h \
	strassen.h

SOURCEFILES = \
	Matrix.cc \
	Messages.cc \
	Node.cc \
	Timer.cc \
	strassen.cc

OBJ = $(foreach i, $(SOURCEFILES), $(basename $i).o)

strassen : $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $^

ctags : $(INTERFACEFILES) $(SOURCEFILES) $(HEADERFILES)
	$(CTAGS) $^

$(OBJ) : $(HEADERFILES) $(GENERATEDHEADERFILES)

N = 2048
B = 32

run :
	for P in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		bash -c "./strassen +p$$P -N $N -b $B"; \
	done

run-serial :
	./strassen +p1 -N $N -b $B

clean :
	rm -f $(OBJ)
	rm -f $(GENERATEDHEADERFILES)
	rm -f strassen charmrun
