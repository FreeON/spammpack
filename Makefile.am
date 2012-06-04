SUBDIRS = src tests

CONFIG_STATUS_DEPENDENCIES = Doxyfile.in

pdfdir = ${docdir}/pdf
pdf_DATA = latex/refman.pdf

EXTRA_DIST = $(pdf_DATA)

docs :
	doxygen
	$(MAKE) -C latex
