.SUFFIXES : .o .cc .ci

ALGORITHM = FUTURES
#ALGORITHM = CALLBACKS
#LOGLEVEL = DEBUG
LOGLEVEL = INFO
#LOGMETHOD = LOGMETHOD_PRINTF
LOGMETHOD = LOGMETHOD_CKPRINTF

CXX = charmc
#CXX = ~/charm/bin/charmc
#CXXFLAGS = -O0 -g -Wall
CXXFLAGS = -O2 -g
CTAGS = ctags
LIBS = -lrt
#LDFLAGS = -tracemode projections
CPPFLAGS = -DLOGLEVEL=$(LOGLEVEL) -D$(ALGORITHM) -D$(LOGMETHOD)

%.decl.h %.def.h : %.ci
	$(CXX) -preprocess $(CPPFLAGS) $<

.cc.o :
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

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

# Can't go deeper than depth = 4 right now.
#N = 1024
#B = 64
N = 2048
B = 128

run : strassen
	for P in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		bash -c "./strassen +p$$P +setcpuaffinity -N $N -b $B"; \
	done

run-serial : strassen
	./strassen +p1 -N $N -b $B --verify

clean :
	rm -f $(OBJ)
	rm -f $(GENERATEDHEADERFILES)
	rm -f strassen charmrun
