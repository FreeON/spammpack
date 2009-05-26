CLEANFILES = *~

AM_CFLAGS = -I$(top_srcdir)/src

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = \
  f90_matrix_allocation \
  f90_set_element

f90_matrix_allocation_SOURCES = f90_matrix_allocation.f90
f90_matrix_allocation_FCFLAGS = -I$(top_builddir)/interface/f90
f90_matrix_allocation_LDADD = $(top_builddir)/interface/f90/liblalf90.la $(top_builddir)/src/liblal.la
#f90_matrix_allocation_LDFLAGS = -static

f90_set_element_SOURCES = f90_set_element.f90
f90_set_element_FCFLAGS = -I$(top_builddir)/interface/f90
f90_set_element_LDADD = $(top_builddir)/interface/f90/liblalf90.la $(top_builddir)/src/liblal.la
