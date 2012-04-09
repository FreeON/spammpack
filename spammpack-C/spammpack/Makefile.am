AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src interfaces tests

#EXTRA_DIST = localversion $(pdf_DATA)

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/Doxyfile.in \
  $(top_srcdir)/src/spamm_manual_mainpage.h.in

#pdfdir = $(docdir)/pdf
#pdf_DATA = latex/refman.pdf

#dist : docs

docs :
	doxygen
	$(MAKE) -C latex
