SUBDIRS = threaded

lib_LTLIBRARIES = libspammpack-serial.la

libspammpack_serial_la_SOURCES = \
  SpAMM_ALGEBRA.F90 \
  SpAMM_CONVERT.F90 \
  SpAMM_DERIVED.F90 \
  SpAMM_GLOBALS.F90 \
  SpAMM_MNGMENT.F90 \
  SpAMM_PROJECT.F90 \
  spamm_exit.c \
  spamm_get_time.c \
  spamm_trap.c \
  spammpack.F90

nodist_include_HEADERS = \
	spammpack.mod \
	spamm_algebra.mod \
	spamm_convert.mod \
	spamm_derived.mod \
	spamm_globals.mod \
	spamm_mngment.mod \
	spamm_project.mod

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
