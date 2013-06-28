.SUFFIXES : .o .cc

.cc.o :
	$(CXX) -c -o $@ $(CXXFLAGS) $<

OBJ = strassenOMP.o \
	Matrix.o \
	Node.o

strassenOMP : $(OBJ)
	$(CXX) -o $@ $(LDFLAGS) $(LIBS) $^

clean :
	rm -f $(OBJ)
