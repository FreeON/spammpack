noinst_LTLIBRARIES = libspammtests.la
libspammtests_la_SOURCES = test_utilities.F90 spammtests.F90

noinst_PROGRAMS = \
  spamm_multiply_01 \
  spamm_multiply_02 \
  spamm_multiply_04 \
  spamm_multiply_08 \
  spamm_multiply_16

TESTS = \
  spamm_multiply_01 \
  spamm_multiply_02 \
  spamm_multiply_04 \
  spamm_multiply_08 \
  spamm_multiply_16

AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = -I../src

spamm_multiply_01_SOURCES  = spamm_multiply.F90
spamm_multiply_01_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_01_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=1

spamm_multiply_02_SOURCES  = spamm_multiply.F90
spamm_multiply_02_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_02_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=2

spamm_multiply_04_SOURCES  = spamm_multiply.F90
spamm_multiply_04_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_04_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=4

spamm_multiply_08_SOURCES  = spamm_multiply.F90
spamm_multiply_08_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_08_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=8

spamm_multiply_16_SOURCES  = spamm_multiply.F90
spamm_multiply_16_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_16_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=16

TAGS_FILES = ../src/*.F90 ../src/*.c
CTAGSFLAGS = --fortran-kinds=+L+i

# Dependencies.
spamm_multiply_01.o : test_utilities.o spammtests.o
spamm_multiply_02.o : test_utilities.o spammtests.o
spamm_multiply_04.o : test_utilities.o spammtests.o
spamm_multiply_08.o : test_utilities.o spammtests.o
spamm_multiply_16.o : test_utilities.o spammtests.o

spammtests.o : test_utilities.o
