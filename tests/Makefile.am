noinst_PROGRAMS = multiply

TESTS = multiply

AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = -I../src
AM_LDFLAGS = ../src/libspammpack.la

multiply_SOURCES = multiply.F90 load_matrix.F90

TAGS_FILES = ../src/*.F90 ../src/*.c

# Dependencies.
multiply.o : load_matrix.o
