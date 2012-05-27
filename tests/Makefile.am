noinst_LTLIBRARIES = libspammtests.la
libspammtests_la_SOURCES = test_utilities.F90 spammtests.F90

noinst_PROGRAMS = \
  spamm_multiply_1 \
  spamm_multiply_2 \
  spamm_multiply_4 \
  spamm_multiply_8

TESTS = \
  spamm_multiply_1 \
  spamm_multiply_2 \
  spamm_multiply_4 \
  spamm_multiply_8

AM_FCFLAGS = $(OPENMP_FCFLAGS)

spamm_multiply_1_SOURCES  = spamm_multiply.F90
spamm_multiply_1_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_1_SOURCES  = spamm_multiply.F90
spamm_multiply_1_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=1
spamm_multiply_2_SOURCES  = spamm_multiply.F90
spamm_multiply_2_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_2_SOURCES  = spamm_multiply.F90
spamm_multiply_2_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=2
spamm_multiply_4_SOURCES  = spamm_multiply.F90
spamm_multiply_4_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_4_SOURCES  = spamm_multiply.F90
spamm_multiply_4_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=4
spamm_multiply_8_SOURCES  = spamm_multiply.F90
spamm_multiply_8_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_8_SOURCES  = spamm_multiply.F90
spamm_multiply_8_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=8

TAGS_FILES = ../src/*.F90 ../src/*.c
CTAGSFLAGS = --fortran-kinds=+L+i

# Dependencies.
spamm_multiply_1.o : test_utilities.o spammtests.o
spamm_multiply_2.o : test_utilities.o spammtests.o
spamm_multiply_4.o : test_utilities.o spammtests.o
spamm_multiply_8.o : test_utilities.o spammtests.o

spammtests.o : test_utilities.o
