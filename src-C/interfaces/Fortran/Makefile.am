lib_LTLIBRARIES = libspammpack_fortran.la

libspammpack_fortran_la_LDFLAGS = -version-info 0:3:0
#libspammpack_fortran_la_CPPFLAGS = -I$(top_srcdir)/src -I$(builddir)
#libspammpack_fortran_la_SHORTNAME = sp
AM_CPPFLAGS = -I$(top_srcdir)/src

libspammpack_fortran_la_SOURCES = \
  spamm_fortran_interface.c \
  spammpack_algebra.F90 \
  spammpack_chunk.F90 \
  spammpack_management.F90 \
  spammpack_omp.F90 \
  spammpack_types.F90 \
  spammpack.F90

include_HEADERS = spammpack.mod

BUILT_SOURCES = spammpack.mod

spammpack.mod : libspammpack_fortran.la

# Fortran dependencies.
spammpack.lo : \
	spammpack_algebra.lo \
	spammpack_chunk.lo \
	spammpack_management.lo \
	spammpack_types.lo
spammpack_algebra.lo : \
	spammpack_chunk.lo \
	spammpack_management.lo \
	spammpack_types.lo
spammpack_management.lo : \
	spammpack_types.lo
spammpack_chunk.lo : \
	spammpack_types.lo

CLEANFILES = *.mod
