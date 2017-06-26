# Useful to have this as a Makefile so we only run if things have changed

out_dir = oco2_2stream_fm

$(out_dir)/stdout_lidort: l2_fp $(base_dir)/config/config_lidort.lua $(srcdir)/input/oco/config/oco_base_config.lua $(base_dir)/config/base_test.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/oco/config/oco_config.lua
	mkdir -p $(out_dir)
	(time ./l2_fp $(base_dir)/config/config_lidort.lua $(out_dir)/out_lidort.h5) 2>&1 | tee $(out_dir)/stdout_lidort

$(out_dir)/stdout_2stream: l2_fp $(base_dir)/config/config_2stream.lua $(srcdir)/input/oco/config/oco_base_config.lua $(base_dir)/config/base_test.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/oco/config/oco_config.lua
	mkdir -p $(out_dir)
	(time ./l2_fp $(base_dir)/config/config_2stream.lua $(out_dir)/out_2stream.h5) 2>&1 | tee $(out_dir)/stdout_2stream

compare: $(out_dir)/stdout_lidort $(out_dir)/stdout_2stream
    # Compare LIDORT output with expected
	$(HDF5_BIN)/h5ls -r $(base_dir)/out_lidort.expected.h5 | \
              sort > $(out_dir)/out_lidort.expected.h5.ls
	$(HDF5_BIN)/h5ls -r $(out_dir)/out_lidort.h5 | sort | \
          diff - $(out_dir)/out_lidort.expected.h5.ls && \
          echo "No diff" && exit 0; exit 1
	PATH=$(HDF5_BIN):$(PATH) $(RUBY) $(srcdir)/support/utils/compare_hdf \
          $(out_dir)/out_lidort.h5 $(base_dir)/out_lidort.expected.h5

    # Compare 2stream output with expected
	$(HDF5_BIN)/h5ls -r $(base_dir)/out_2stream.expected.h5 | \
              sort > $(out_dir)/out_2stream.expected.h5.ls
	$(HDF5_BIN)/h5ls -r $(out_dir)/out_2stream.h5 | sort | \
          diff - $(out_dir)/out_2stream.expected.h5.ls && \
          echo "No diff" && exit 0; exit 1
	PATH=$(HDF5_BIN):$(PATH) $(RUBY) $(srcdir)/support/utils/compare_hdf \
          $(out_dir)/out_2stream.h5 $(base_dir)/out_2stream.expected.h5

    # Compare 2stream and LIDORT, there shouldnt be any differences
	PATH=$(HDF5_BIN):$(PATH) $(RUBY) $(srcdir)/support/utils/compare_hdf \
          $(out_dir)/out_2stream.h5 $(out_dir)/out_lidort.h5

update_expected:
	-rm -r $(out_dir)
	$(MAKE) -f $(base_dir)/Makefile $(out_dir)/stdout_lidort && \
        $(MAKE) -f $(base_dir)/Makefile $(out_dir)/stdout_2stream && \
        cp $(out_dir)/out_lidort.h5 $(base_dir)/out_lidort.expected.h5 && \
        cp $(out_dir)/out_2stream.h5 $(base_dir)/out_2stream.expected.h5 && \
        exit 0; exit 1

all: compare
