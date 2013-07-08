.SUFFIXES : .o .cc .ci

ALGORITHM = FUTURES
#ALGORITHM = CALLBACKS
#LOGLEVEL = DEBUG
LOGLEVEL = INFO

CXX = charmc
#CXXFLAGS = -O0 -g -Wall
CXXFLAGS = -O3 -g
CTAGS = ctags
LIBS = -lrt
#LDFLAGS = -tracemode projections
CPPFLAGS = -DLOGLEVEL=$(LOGLEVEL) -D$(ALGORITHM)

%.decl.h %.def.h : %.ci
	$(CXX) -preprocess $(CPPFLAGS) $<

.cc.o :
	$(CXX) -c -o $@ $(CPPFLAGS) $(CXXFLAGS) $<

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
	strassen.h \
	Timer.h \
	Utilities.h

SOURCEFILES = \
	Matrix.cc \
	Messages.cc \
	Node.cc \
	strassen.cc \
	Timer.cc \
	Utilities.cc

OBJ = $(foreach i, $(SOURCEFILES), $(basename $i).o)

strassen : $(OBJ)
	$(CXX) -language charm++ -o $@ $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS)

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
