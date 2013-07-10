.SUFFIXES : .o .cc .ci

CXX = charmc
#CXXFLAGS = -O0 -g -Wall
CXXFLAGS = -O2 -g
CTAGS = ctags
#LIBS =
#LDFLAGS = -tracemode projections
#CPPFLAGS =

%.decl.h %.def.h : %.ci
	$(CXX) -preprocess $(CPPFLAGS) $<

.cc.o :
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

INTERFACEFILES = Futures.ci

GENERATEDHEADERFILES = $(foreach i, $(INTERFACEFILES), $(basename $i).def.h $(basename $i).decl.h)

SOURCEFILES = Futures.cc

OBJ = $(foreach i, $(SOURCEFILES), $(basename $i).o)

futures : $(OBJ)
	$(CXX) -language charm++ -o $@ $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS)

ctags : $(INTERFACEFILES) $(SOURCEFILES) $(HEADERFILES)
	$(CTAGS) $^

$(OBJ) : $(HEADERFILES) $(GENERATEDHEADERFILES)

clean :
	rm -f $(OBJ)
	rm -f $(GENERATEDHEADERFILES)
	rm -f futures charmrun
