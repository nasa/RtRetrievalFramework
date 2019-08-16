# Useful to have this as a Makefile so we only run if things have changed

out_dir = oco2_simulator

SHELL := /bin/bash

$(out_dir)/comparison_plot.pdf: l2_fp $(base_dir)/config/config_orbit_sim_match.lua $(srcdir)/input/oco/config/oco_base_config.lua $(srcdir)/input/common/config/config_common.lua $(srcdir)/input/oco/config/oco_config.lua
	mkdir -p $(out_dir)
	source $(prefix)/setup_fp_env.sh && \
	$(PYTHON) $(base_dir)/compare_simulator.py && \
        $(PYTHON) $(base_dir)/print_values.py

all: $(out_dir)/comparison_plot.pdf
