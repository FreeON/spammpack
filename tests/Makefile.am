noinst_LTLIBRARIES = libspammtests.la
libspammtests_la_SOURCES = test_utilities.F90 spammtests.F90

noinst_SCRIPTS = scaling_test.sh parseTestResults.py

EXTRA_DIST = $(noinst_SCRIPTS)

noinst_PROGRAMS = \
  spamm_multiply_01 \
  spamm_multiply_02 \
  spamm_multiply_04 \
  spamm_multiply_08 \
  spamm_multiply_12 \
  spamm_multiply_16 \
  spamm_multiply_20 \
  spamm_multiply_24 \
  spamm_multiply_28 \
  spamm_multiply_32 \
  spamm_multiply_36 \
  spamm_multiply_40 \
  spamm_multiply_44 \
  spamm_multiply_48

TESTS = scaling_test.sh

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

spamm_multiply_12_SOURCES  = spamm_multiply.F90
spamm_multiply_12_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_12_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=12

spamm_multiply_16_SOURCES  = spamm_multiply.F90
spamm_multiply_16_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_16_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=16

spamm_multiply_20_SOURCES  = spamm_multiply.F90
spamm_multiply_20_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_20_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=20

spamm_multiply_24_SOURCES  = spamm_multiply.F90
spamm_multiply_24_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_24_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=24

spamm_multiply_28_SOURCES  = spamm_multiply.F90
spamm_multiply_28_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_28_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=28

spamm_multiply_32_SOURCES  = spamm_multiply.F90
spamm_multiply_32_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_32_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=32

spamm_multiply_36_SOURCES  = spamm_multiply.F90
spamm_multiply_36_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_36_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=36

spamm_multiply_40_SOURCES  = spamm_multiply.F90
spamm_multiply_40_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_40_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=40

spamm_multiply_44_SOURCES  = spamm_multiply.F90
spamm_multiply_44_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_44_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=44

spamm_multiply_48_SOURCES  = spamm_multiply.F90
spamm_multiply_48_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_48_CPPFLAGS = -I../src -DNUMBER_OF_THREADS=48

TAGS_FILES = ../src/*.F90 ../src/*.c
CTAGSFLAGS = --fortran-kinds=+L+i

# Dependencies.
spamm_multiply_01.o : test_utilities.o spammtests.o
spamm_multiply_02.o : test_utilities.o spammtests.o
spamm_multiply_04.o : test_utilities.o spammtests.o
spamm_multiply_08.o : test_utilities.o spammtests.o
spamm_multiply_12.o : test_utilities.o spammtests.o
spamm_multiply_16.o : test_utilities.o spammtests.o
spamm_multiply_20.o : test_utilities.o spammtests.o
spamm_multiply_24.o : test_utilities.o spammtests.o
spamm_multiply_28.o : test_utilities.o spammtests.o
spamm_multiply_32.o : test_utilities.o spammtests.o
spamm_multiply_36.o : test_utilities.o spammtests.o
spamm_multiply_40.o : test_utilities.o spammtests.o
spamm_multiply_44.o : test_utilities.o spammtests.o
spamm_multiply_48.o : test_utilities.o spammtests.o

spammtests.o : test_utilities.o
