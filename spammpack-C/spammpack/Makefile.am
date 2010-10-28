AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests prototype figures

dist_data_DATA = latex/refman.pdf

EXTRA_DIST = localversion

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/localversion \
  $(top_srcdir)/src/spamm_mainpage.h.in

docs :
	doxygen
	$(MAKE) -C latex
