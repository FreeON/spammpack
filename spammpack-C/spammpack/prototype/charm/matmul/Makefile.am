TAGS_FILES = config.h

EXTRA_DIST = \
	$(INTERFACE_FILES) \
	matmul.bib \
	configure.LB \
	configure.mustang \
	freeon-submit.sh \
	freeon.job.template \
	moonlight-submit.sh \
	moonlight.sh.template \
	mustang-results.sh \
	mustang-submit.sh \
	mustang.job.template

AM_LDFLAGS = -language charm++ @LOADBALANCER_FLAGS@ @DEFAULT_LB_FLAGS@ @CHARMDEBUG_FLAGS@

bin_PROGRAMS = matmul

INTERFACE_FILES = \
	matmul.ci \
	matrix.ci \
	messages.ci \
	multiply.ci \
	multiplyelement.ci \
	node.ci

matmul_SOURCES = \
	blas_interface.h \
	index.h \
	logger.cc \
	logger.h \
	matmul.cc \
	matmul.h \
	matrix.cc \
	matrix.h \
	messages.cc \
	messages.h \
	multiply.cc \
	multiplyelement.cc \
	multiplyelement.h \
	multiply.h \
	node.cc \
	node.h \
	timer.cc \
	timer.h \
	types.h \
	utilities.cc \
	utilities.h

$(matmul_OBJECTS) : $(BUILT_SOURCES)

# Matrices are blocked into 32z32 at the leaf level.
N = 4096
B = 128

RSH = rsh

run : matmul$(EXEEXT)
	for i in 1 2 4; do \
		./charmrun +p$$i ++remote-shell $(RSH) matmul -N $(N) -b $(B); \
	done

run-local : matmul$(EXEEXT)
	for i in 1 2 4; do \
		./matmul +p$$i -N $(N) -b $(B); \
	done

run-LB : matmul$(EXEEXT)
	./charmrun +p4 matmul -N $(N) -b $(B) --verify +LBDebug

valgrind : matmul$(EXEEXT)
	valgrind --log-file=valgrind.%p.log \
		--track-origins=yes \
		--leak-check=full \
		--show-reachable=yes \
		./$< -N 1

docs :
	$(DOXYGEN)

docs-pdf : docs
	$(MAKE) -C latex

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) $(DEFAULT_INCLUDES) -preprocess $<

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
