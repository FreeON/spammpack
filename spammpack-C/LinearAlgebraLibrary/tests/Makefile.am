ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = $(top_builddir)/src/libspamm.la

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = create_spamm add_spamm multiply_random_spamm

create_spamm_SOURCES = create_spamm.c
add_spamm_SOURCES = add_spamm.c
multiply_random_spamm_SOURCES = multiply_random_spamm.c
