AM_CTAGSFLAGS = --Fortran-kinds=+i+L
TAGS_FILES = \
  $(top_srcdir)/src/*c \
  $(top_srcdir)/src/*h

lib_LTLIBRARIES = libspammpack_f90.la

include_HEADERS = spammpack.mod

AM_CPPFLAGS = -I$(top_srcdir)/src

libspammpack_f90_la_SOURCES = \
  spamm_exit.c \
  spamm_interface.c \
  spamm_interface.h \
  spamm_version.F90 \
  spamm_derived.F90 \
  spammpack.F90

CLEANFILES = *.mod
