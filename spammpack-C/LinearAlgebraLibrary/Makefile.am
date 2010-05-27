AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests benchmarks prototype

latex/refman.pdf : docs

dist_data_DATA = latex/refman.pdf

docs :
	doxygen
	$(MAKE) -C latex
