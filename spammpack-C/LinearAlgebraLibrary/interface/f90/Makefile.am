CLEANFILES = *~

lib_LTLIBRARIES = liblalf90.la

liblalf90_la_SOURCES = lal_types.source lal_interface.f90

include_HEADERS = lal.mod

lal.mod : lal_interface.lo
lal_types.mod : lal_types.lo
lal_interface.lo : lal_types.mod

.source.f90 :
	cpp -P -I${top_builddir} $^ > $@
