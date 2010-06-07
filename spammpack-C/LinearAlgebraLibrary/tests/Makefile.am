ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = $(top_builddir)/src/libspamm.la $(top_builddir)/src/libspamm_ll.la

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = \
  ll_add \
  create_spamm \
  nonzero_spamm \
  add_spamm \
  multiply_random_spamm \
  tree_pack

ll_add_SOURCES = ll_add.c
create_spamm_SOURCES = create_spamm.c
nonzero_spamm_SOURCES = nonzero_spamm.c
add_spamm_SOURCES = add_spamm.c
multiply_random_spamm_SOURCES = multiply_random_spamm.c
tree_pack_SOURCES = tree_pack.c

tags :
	ctags ../src/*.{c,h}
