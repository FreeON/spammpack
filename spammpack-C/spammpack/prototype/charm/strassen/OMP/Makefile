.SUFFIXES : .o .cc

#CXX = g++
#CXXFLAGS = -O3 -march=native -g -Wall -fopenmp

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

run :
	for P in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		bash -c "OMP_NUM_THREADS=$$P ./strassenOMP -N 2048 -b 32"; \
	done

run-serial :
	./strassenOMP -N 2048 -b 32

clean :
	rm -f $(OBJ)
