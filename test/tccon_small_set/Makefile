# Create common data set.
val := $(shell expr $(job_index) + 1)

all: tccon_sounding_$(val)_test/gosat_tccon_sounding_$(val)_test.config

tccon_sounding_$(val)_test/gosat_tccon_sounding_$(val)_test.config: l2_fp $(srcdir)/input/gosat/config/acos_config.lua $(srcdir)/input/gosat/config/gosat_base_config.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/gosat/config/config.lua $(srcdir)/input/gosat/input/l2_gosat_static_input.h5
	mkdir -p tccon_sounding_$(val)_test
	-rm -rf tccon_sounding_$(val)_test/log
	-rm -rf tccon_sounding_$(val)_test/output
	srcdir="${srcdir}" val="${val}" cmd="${cmd}" $(srcdir)/test/tccon_small_set/do_populate.sh
