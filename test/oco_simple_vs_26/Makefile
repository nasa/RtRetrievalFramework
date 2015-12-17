# Useful to have this as a Makefile so we only run if things have changed

out_dir = oco_simple_run_test

$(out_dir)/stdout: l2_fp $(base_dir)/config/config_ascii_atm.lua $(base_dir)/config/ascii_atm_base_config.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/oco/config/oco_config.lua
	mkdir -p $(out_dir)
	echo $(abscodir)
	echo "# This is the command used to run the full run test, if you want to run it by hand" > $(out_dir)/run_command.txt
	echo "LUA_PATH='$(LUA_PATH)' abscodir=$(abscodir) ./l2_fp $(base_dir)/config/config_ascii_atm.lua $(out_dir)/out.h5" >> $(out_dir)/run_command.txt
	(time ./l2_fp $(base_dir)/config/config_ascii_atm.lua $(out_dir)/out.h5) 2>&1 | tee $(out_dir)/stdout

compare: $(out_dir)/stdout
	$(HDF5_BIN)/h5ls -r $(base_dir)/out.expected.h5 | \
              sort > $(out_dir)/out.expected.h5.ls
	$(HDF5_BIN)/h5ls -r $(out_dir)/out.h5 | sort | \
          diff - $(out_dir)/out.expected.h5.ls && \
          echo "No diff" && exit 0; exit 1
	PATH=$(HDF5_BIN):$(PATH) $(RUBY) $(srcdir)/support/utils/compare_hdf \
          $(out_dir)/out.h5 $(base_dir)/out.expected.h5

update_expected:
	-rm -r $(out_dir)
	$(MAKE) -f $(base_dir)/Makefile $(out_dir)/stdout && \
        cp $(out_dir)/out.h5 $(base_dir)/out.expected.h5 && exit 0; exit 1

all: compare
