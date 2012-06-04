SUBDIRS = src tests

CONFIG_STATUS_DEPENDENCIES = Doxyfile.in

pdfdir = ${docdir}/pdf
pdf_DATA = latex/refman.pdf

man3dir = ${mandir}/man3
man3_MANS = man/man3/*

EXTRA_DIST = $(pdf_DATA) $(man3_MANS)

docs :
	doxygen
	$(MAKE) -C latex
