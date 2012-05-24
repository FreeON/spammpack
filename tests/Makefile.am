noinst_PROGRAMS = spamm_multiply

TESTS = spamm_multiply

AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = -I../src
AM_LDFLAGS = ../src/libspammpack.la

spamm_multiply_SOURCES = spamm_multiply.F90 test_utilities.F90

TAGS_FILES = ../src/*.F90 ../src/*.c

# Dependencies.
spamm_multiply.o : test_utilities.o
