AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests

EXTRA_DIST = localversion

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/localversion

dist: docs

docs:
	doxygen
