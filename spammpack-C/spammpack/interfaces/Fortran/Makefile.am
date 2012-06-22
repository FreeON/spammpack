lib_LTLIBRARIES = libspammpack_f90.la

include_HEADERS = spammpack.mod

libspammpack_f90_la_SOURCES = \
  spammpack.F90 \
  spamm_convert.F90 \
  spamm_derived.F90 \
  spamm_version.F90

# Dependencies.
spammpack.lo : spamm_convert.lo
spamm_convert.lo : spamm_derived.lo

CLEANFILES = *.mod
