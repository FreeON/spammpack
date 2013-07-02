.SUFFIXES : .o .cc

CXX = g++
CXXFLAGS = -O3 -march=native -g -Wall -fopenmp

#CXX = icpc
#CXXFLAGS = -O3 -xHost -g -Wall -fopenmp

LIBS = -lrt

.cc.o :
	$(CXX) -c -o $@ $(CXXFLAGS) $<

OBJ = strassenOMP.o \
	Matrix.o \
	Node.o \
	Timer.o

strassenOMP : $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $^


N = 4096
b = 64

run :
	for P in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		bash -c "OMP_NUM_THREADS=$$P ./strassenOMP -N $N -b $b"; \
	done

run-serial :
	./strassenOMP -N $N -b $b --verify

clean :
	rm -f $(OBJ)
