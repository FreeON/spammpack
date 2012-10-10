AM_CTAGSFLAGS = --C-kinds=+l+x --Fortran-kinds=+i+L
TAGS_FILES = \
  $(top_srcdir)/src/*.h \
  $(top_srcdir)/src/*.c \
  $(top_srcdir)/src/*.S \
  $(top_srcdir)/interfaces/Fortran/*.F90 \
  $(top_srcdir)/interfaces/Fortran/*.c \
  $(top_srcdir)/interfaces/Fortran/*.h

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/src -I$(top_builddir)/interfaces/Fortran
AM_LDFLAGS = $(top_builddir)/src/libspamm.la
#AM_LDFLAGS = $(top_builddir)/src/libspamm.la -all-static

if HAVE_SSE_IN_CPU
multiply_spamm_SSE = multiply_spamm_kernel_standard_SSE.sh
endif

if HAVE_SSE4_1_IN_CPU
multiply_spamm_SSE4_1 = multiply_spamm_kernel_standard_SSE4_1.sh
endif

multiply_tests = $(multiply_spamm_SSE) $(multiply_spamm_SSE4_1)

check_LTLIBRARIES = libspammtest.la
libspammtest_la_SOURCES = \
  test_utilities.F90 \
  spammtests.F90

noinst_PROGRAMS = \
  convert_to_binary

get_version_SOURCES = get_version.F90

spamm_multiply_SOURCES = spamm_multiply.F90
spamm_multiply_LDFLAGS = $(top_builddir)/interfaces/Fortran/libspammpack_f90.la libspammtest.la
spamm_multiply_DEPENDENCIES = libspammtest.la

SP2_SOURCES = SP2.F90 spamm_project.F90
SP2_LDFLAGS = $(top_builddir)/interfaces/Fortran/libspammpack_f90.la libspammtest.la
SP2_DEPENDENCIES = libspammtest.la

noinst_SCRIPTS = \
  create_spamm_row_major.sh \
  create_spamm_Z_curve.sh \
  convert_dense_row_major.sh \
  convert_dense_Z_curve.sh \
  multiply_spamm_kernel_standard_SSE.sh \
  multiply_spamm_kernel_standard_SSE4_1.sh

EXTRA_DIST = $(noinst_SCRIPTS)

check_PROGRAMS = \
  add_spamm_dense \
  add_spamm_sparse \
  check_spamm \
  convert_dense \
  create_spamm \
  create_spamm_vector \
  get_version \
  hashtable_insert \
  hashtable_remove \
  list \
  list_sort \
  multiply_spamm \
  spamm_multiply \
  stats \
  SP2

TESTS = \
  $(multiply_tests) \
  $(check_PROGRAMS)

CLEANFILES = *.mod
