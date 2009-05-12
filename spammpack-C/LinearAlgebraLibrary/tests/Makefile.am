CLEANFILES = *~

AM_CFLAGS = -I$(top_srcdir)/src

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = \
  matrix_allocation \
  matrix_allocation_f90 \
  matrix_free \
  get_element \
  set_element \
  matrix_equals \
  matrix_transpose \
  matrix_add \
  matrix_multiply

matrix_allocation_SOURCES = matrix_allocation.c
matrix_allocation_LDADD = $(top_builddir)/src/liblal.la

matrix_allocation_f90_SOURCES = matrix_allocation_f90.f90
matrix_allocation_f90_FCFLAGS = -I$(top_builddir)/interface/f90
matrix_allocation_f90_LDADD = $(top_builddir)/src/liblal.la $(top_builddir)/interface/f90/liblalf90.la
matrix_allocation_f90_LDFLAGS = -static

matrix_free_SOURCES = matrix_free.c
matrix_free_LDADD = $(top_builddir)/src/liblal.la
#matrix_free_LDFLAGS = -static

get_element_SOURCES = get_element.c
get_element_LDADD = $(top_builddir)/src/liblal.la
#get_element_LDFLAGS = -static

set_element_SOURCES = set_element.c
set_element_LDADD = $(top_builddir)/src/liblal.la

matrix_equals_SOURCES = matrix_equals.c
matrix_equals_LDADD = $(top_builddir)/src/liblal.la

matrix_transpose_SOURCES = matrix_transpose.c
matrix_transpose_LDADD = $(top_builddir)/src/liblal.la

matrix_add_SOURCES = matrix_add.c
matrix_add_LDADD = $(top_builddir)/src/liblal.la

matrix_multiply_SOURCES = matrix_multiply.c
matrix_multiply_LDADD = $(top_builddir)/src/liblal.la
