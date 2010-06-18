ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = $(top_builddir)/src/libspamm.la

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = \
  ll_append \
  ll_sort \
  mm_allocate \
  create_spamm \
  nonzero_spamm \
  tree_pack \
  add_spamm \
  add_spamm_linear \
  multiply_random_spamm

ll_append_SOURCES = ll_append.c
ll_sort_SOURCES = ll_sort.c
mm_allocate_SOURCES = mm_allocate.c
create_spamm_SOURCES = create_spamm.c
nonzero_spamm_SOURCES = nonzero_spamm.c
add_spamm_SOURCES = add_spamm.c
add_spamm_linear_SOURCES = add_spamm_linear.c
multiply_random_spamm_SOURCES = multiply_random_spamm.c
tree_pack_SOURCES = tree_pack.c

tags :
	ctags ../src/*.{c,h}
