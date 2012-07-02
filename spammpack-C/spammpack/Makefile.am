AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src interfaces tests

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/Doxyfile.in \
  $(top_srcdir)/src/spamm_manual_mainpage.h.in

EXTRA_DIST = latex/refman.pdf \
  extract_conf.py

docs :
	doxygen
	$(MAKE) -C latex
