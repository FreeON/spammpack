.SUFFIXES : .o .cc

#CXX = g++
#CXXFLAGS = -O3 -g -Wall -fopenmp

CXX = icpc
CXXFLAGS = -O3 -xHost -g -Wall -fopenmp

LIBS = -lrt

.cc.o :
	$(CXX) -c -o $@ $(CXXFLAGS) $<

OBJ = strassenOMP.o \
	Matrix.o \
	Node.o \
	Timer.o

strassenOMP : $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $^

clean :
	rm -f $(OBJ)
