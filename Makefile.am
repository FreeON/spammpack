SUBDIRS = src tests

CONFIG_STATUS_DEPENDENCIES = Doxyfile.in

docs :
	doxygen
	$(MAKE) -C latex
