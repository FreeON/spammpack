#SUBDIRS = ll mm add multiply

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = $(top_builddir)/src/libspamm.la
#AM_LDFLAGS = $(top_builddir)/src/libspamm.la -all-static

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = \
  create_spamm \
  nonzero_spamm \
  spamm_norm \
  multiply_spamm \
  multiply_spamm_exponential_decay

#tree_pack

tags :
	ctags ../src/*.{c,h}
