lib_LTLIBRARIES = libspammpack.la

libspammpack_la_SOURCES = \
  SpAMM_ALGEBRA.F90 \
  SpAMM_CONVERT.F90 \
  SpAMM_DERIVED.F90 \
  SpAMM_GLOBALS.F90 \
  SpAMM_MNGMENT.F90 \
  SpAMM_PROJECT.F90 \
  spamm_exit.c \
  spamm_get_time.c \
  spamm_trap.c \
  spammpack.F90 \
  cmultiply.c

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)

CLEANFILES = *.mod

# Dependencies.
SpAMM_ALGEBRA.lo : SpAMM_DERIVED.lo SpAMM_GLOBALS.lo SpAMM_MNGMENT.lo
SpAMM_CONVERT.lo : SpAMM_DERIVED.lo SpAMM_GLOBALS.lo SpAMM_ALGEBRA.lo SpAMM_MNGMENT.lo
SpAMM_GLOBALS.lo : SpAMM_DERIVED.lo
SpAMM_MNGMENT.lo : SpAMM_DERIVED.lo SpAMM_GLOBALS.lo
spammpack.lo     : SpAMM_DERIVED.lo SpAMM_GLOBALS.lo SpAMM_MNGMENT.lo SpAMM_ALGEBRA.lo SpAMM_PROJECT.lo
SpAMM_PROJECT.lo : SpAMM_DERIVED.lo SpAMM_GLOBALS.lo SpAMM_MNGMENT.lo SpAMM_ALGEBRA.lo

.PHONY : depcheck
depcheck :
	for i in $(libspammpack_la_OBJECTS); do \
		echo "depchecking $$i"; \
		$(MAKE) clean; \
		$(MAKE) $$i || exit 1; \
	done
	$(MAKE) clean

TAGS_FILES = *.F90 *.c
CTAGSFLAGS = --fortran-kinds=+L+i
