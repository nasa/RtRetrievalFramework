# Useful to have this as a Makefile so we only run if things have changed

out_dir = fts_test

$(out_dir)/fts_fts_test.config: 
	mkdir -p $(out_dir)
	-rm -rf $(out_dir)/log
	-rm -rf $(out_dir)/output
	srcdir="${srcdir}" val="${val}" cmd="${cmd}" $(srcdir)/test/fts/do_populate.sh	

$(out_dir)/stdout: $(out_dir)/fts_fts_test.config l2_fp $(srcdir)/input/fts/config/config.lua $(srcdir)/input/fts/config/fts_base_config.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/fts/config/fts_config.lua
	$(out_dir)/l2_fp_job.sh -e 0

compare: $(out_dir)/stdout
	$(HDF5_BIN)/h5ls -r $(base_dir)/out.expected.h5 | \
              sort > $(out_dir)/out.expected.h5.ls
	$(HDF5_BIN)/h5ls -r $(out_dir)/output/l2_008.h5 | sort | \
          diff - $(out_dir)/out.expected.h5.ls && \
          echo "No diff" && exit 0; exit 1
	PATH=$(HDF5_BIN):$(PATH) $(RUBY) $(srcdir)/support/utils/compare_hdf \
          $(out_dir)/output/l2_008.h5 $(base_dir)/out.expected.h5

update_expected:
	-rm -r $(out_dir)
	$(MAKE) -f $(base_dir)/Makefile $(out_dir)/stdout cmd=l2_fp && \
        cp $(out_dir)/output/l2_008.h5 $(base_dir)/out.expected.h5 && exit 0; exit 1

all: compare
