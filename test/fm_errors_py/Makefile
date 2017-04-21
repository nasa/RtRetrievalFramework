# Useful to have this as a Makefile so we only run if things have changed

out_dir = fm_errors_py
SHELL := /bin/bash

$(out_dir)/stdout: $(srcdir)/support/utils/l2_fp_fm_errors.py  $(base_dir)/config/config.lua $(srcdir)/input/oco/config/oco_base_config.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/oco/config/oco_config.lua $(prefix)/setup_fp_env.sh
	mkdir -p $(out_dir)
	cp $(base_dir)/out.initial.h5 $(out_dir)/out.h5 # Prime the SV to make test case shorter and to test -u option
	source $(prefix)/setup_fp_env.sh && \
	$(srcdir)/support/utils/l2_fp_fm_errors.py -u -r --num_perturbed 1 -f ils/ABO2 -f spectroscopy/temperature/o2 $(base_dir)/config/config.lua $(out_dir)/out.h5 2>&1 | tee $(out_dir)/stdout

$(out_dir)/stdout_full: $(srcdir)/support/utils/l2_fp_fm_errors.py  $(base_dir)/config/config.lua $(srcdir)/input/oco/config/oco_base_config.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/oco/config/oco_config.lua $(prefix)/setup_fp_env.sh
	mkdir -p $(out_dir)
	source $(prefix)/setup_fp_env.sh && \
	$(srcdir)/support/utils/l2_fp_fm_errors.py -r --num_perturbed 1 -f ils/ABO2 -f spectroscopy/temperature/o2 $(base_dir)/config/config.lua $(out_dir)/out.h5 2>&1 | tee $(out_dir)/stdout_full

compare: $(out_dir)/stdout
    # Compare output with expected
	$(HDF5_BIN)/h5ls -r $(base_dir)/out.expected.h5 | \
              sort > $(out_dir)/out.expected.h5.ls
	$(HDF5_BIN)/h5ls -r $(out_dir)/out.h5 | sort | \
          diff - $(out_dir)/out.expected.h5.ls && \
          echo "No diff" && exit 0; exit 1
	PATH=$(HDF5_BIN):$(PATH) $(RUBY) $(srcdir)/support/utils/compare_hdf \
          $(out_dir)/out.h5 $(base_dir)/out.expected.h5

update_expected_full:
	-rm -r $(out_dir)
	$(MAKE) -f $(base_dir)/Makefile $(out_dir)/stdout_full && \
        cp $(out_dir)/out.h5 $(base_dir)/out.initial.h5 && \
        cp $(out_dir)/out.h5 $(base_dir)/out.expected.h5 && exit 0; exit 1

update_expected:
	-rm -r $(out_dir)
	$(MAKE) -f $(base_dir)/Makefile $(out_dir)/stdout && \
        cp $(out_dir)/out.h5 $(base_dir)/out.expected.h5 && exit 0; exit 1

all: compare
