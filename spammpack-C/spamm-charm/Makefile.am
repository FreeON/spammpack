SUBDIRS = src tools tests

DISTCHECK_CONFIGURE_FLAGS = --enable-debug

EXTRA_DIST = \
	spamm.bib \
	configure.LB \
	configure.mustang \
	freeon-submit.sh \
	freeon.job.template \
	moonlight-submit.sh \
	moonlight.sh.template \
	mustang-results.sh \
	mustang-submit.sh \
	mustang.job.PP.template \
	mustang.job.decay.template \
	mustang.job.diagonal.template \
	mustang.job.full.template \
	openMP.job

docs :
	$(DOXYGEN)

docs-pdf : docs
	$(MAKE) -C latex
