dist_noinst_SCRIPTS = scaling-mustang.sh

noinst_PROGRAMS = \
	chunk_multiply \
	chunk_tree_multiply

check_PROGRAMS = \
	chunk_add \
	chunk_add_identity \
	chunk_scale \
	chunk_trace \
	chunk_tree_alloc

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = -I $(top_srcdir)/src
AM_LDFLAGS = ../src/libchunk.la

TESTS = \
	$(check_PROGRAMS) \
	$(noinst_PROGRAMS)
