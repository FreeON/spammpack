noinst_LTLIBRARIES = libspammtests.la
libspammtests_la_SOURCES = test_utilities.F90 spammtests.F90

noinst_PROGRAMS = spamm_multiply

TESTS = spamm_multiply

AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = -I../src

spamm_multiply_LDFLAGS = ../src/libspammpack.la libspammtests.la
spamm_multiply_SOURCES = spamm_multiply.F90

TAGS_FILES = ../src/*.F90 ../src/*.c
CTAGSFLAGS = --fortran-kinds=+L+i

# Dependencies.
spamm_multiply.o : test_utilities.o spammtests.o
spammtests.o : test_utilities.o
