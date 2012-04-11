AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src interfaces tests

CONFIG_STATUS_DEPENDENCIES = $(top_srcdir)/Doxyfile.in

docs :
	doxygen
	$(MAKE) -C latex
