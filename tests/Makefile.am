noinst_LTLIBRARIES = libspammtests.la
libspammtests_la_SOURCES = test_utilities.F90 spammtests.F90

noinst_SCRIPTS = scaling_test.sh parseTestResults.py

EXTRA_DIST = $(noinst_SCRIPTS) \
  freeon-10.Intel.reference \
  freeon-10.gcc.reference \
  scaling_test.reference \
  scaling_test.sh.in \
  set_mode.sh \
  testmatrix_random_1024x1024.coor

if build_test_programs
test_programs = \
  SpAMM_PURIFY2 \
  SpAMM_TESTONE
endif

noinst_PROGRAMS = \
  $(test_programs) \
  spamm_multiply \
  spamm_SP2 \
  convert_to_binary

TESTS = set_mode.sh scaling_test.sh

AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = -I../src

SpAMM_PURIFY2_SOURCES = SpAMM_PURIFY2.F90
SpAMM_PURIFY2_CPPFLAGS = -I$(top_builddir)/.. -I$(top_builddir)/../Modules

SpAMM_TESTONE_SOURCES = SpAMM_TESTONE.F90
SpAMM_TESTONE_CPPFLAGS = -I$(top_builddir)/.. -I$(top_builddir)/../Modules

spamm_multiply_SOURCES  = spamm_multiply.F90
spamm_multiply_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_multiply_CPPFLAGS = -I../src

spamm_SP2_SOURCES  = spamm_SP2.F90
spamm_SP2_LDFLAGS  = ../src/libspammpack.la libspammtests.la
spamm_SP2_CPPFLAGS = -I../src

TAGS_FILES = ../src/*.F90 ../src/*.c
CTAGSFLAGS = --fortran-kinds=+L+i

# Dependencies.
spamm_multiply.o : test_utilities.o spammtests.o
spammtests.o : test_utilities.o

CLEANFILES = *.mod scaling_test.output scaling_test.dat regressionTest.log
