# Create common data set.
val := $(shell expr $(job_index) + 1)

all: oco2_sounding_$(val)_test/oco_oco2_sounding_$(val)_test.config

oco2_sounding_$(val)_test/oco_oco2_sounding_$(val)_test.config: l2_fp $(srcdir)/input/oco/config/oco_config.lua $(srcdir)/input/oco/config/oco_base_config.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/oco/config/config.lua $(srcdir)/input/oco/input/l2_oco_static_input.h5
	mkdir -p oco2_sounding_$(val)_test
	-rm -rf oco2_sounding_$(val)_test/log
	-rm -rf oco2_sounding_$(val)_test/output
	srcdir="${srcdir}" val="${val}" cmd="${cmd}" $(srcdir)/test/oco2_small_set/do_populate.sh
