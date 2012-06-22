AM_CTAGSFLAGS = --Fortran-kinds=+i+L

lib_LTLIBRARIES = libspammpack_f90.la

include_HEADERS = spammpack.mod

libspammpack_f90_la_SOURCES = \
  spamm_convert.F90 \
  spamm_derived.F90 \
  spamm_exit.c \
  spamm_version.F90 \
  spammpack.F90

# Dependencies.
spammpack.lo : spamm_convert.lo
spamm_convert.lo : spamm_derived.lo

CLEANFILES = *.mod
