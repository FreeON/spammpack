bin_PROGRAMS = migration

migration_SOURCES = migration.cc

INTERFACE_FILES = migration.ci

EXTRA_DIST = $(INTERFACE_FILES)

AM_LDFLAGS = -language charm++ -module CommonLBs -balancer GreedyCommLB

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CXX) $(DEFAULT_INCLUDES) -preprocess $<

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun

run-serial : migration$(EXEEXT)
	./charmrun +p1 migration +LBDebug

run : migration$(EXEEXT)
	./charmrun +p4 migration +LBDebug
