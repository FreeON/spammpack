noinst_PROGRAMS = multiply

TESTS = multiply

AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = -I../src
AM_LDFLAGS = ../src/libspammpack.la

multiply_SOURCES = multiply.F90 test_utilities.F90

TAGS_FILES = ../src/*.F90 ../src/*.c

# Dependencies.
multiply.o : test_utilities.o
