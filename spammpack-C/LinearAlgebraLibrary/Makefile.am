AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests benchmarks prototype

dist_data_DATA = latex/refman.pdf

EXTRA_DIST = localversion

CONFIG_STATUS_DEPENDENCIES = $(top_srcdir)/localversion

docs :
	doxygen
	$(MAKE) -C latex
