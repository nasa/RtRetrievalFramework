#! /bin/bash
# We wrap this command in a script so that we know we are using bash rather
# /bin/sh which may or may not be bash.
source ${srcdir}/setup_env.sh && PYTHONPATH=${srcdir}/support:${PYTHONPATH} && \
    cd oco2_sounding_${val}_test && \
    ${srcdir}/support/utils/create_config.py -t oco \
    ${srcdir}/test/oco2_small_set/OCO2_sim_NDa_20120616_168_r69CSUSim02c_001_spliced.h5 \
    ${srcdir}/test/oco2_small_set/OCO2_meteorology_NDa_20120616_168_r69CSUSim02a_001_spliced.h5 \
    ${srcdir}/test/oco2_small_set/OCO2_sim_NDa_20120616_168_r69CSUSim02c_001_spliced.h5
    config_file=${srcdir}/test/oco2_sounding_${val}/config.lua
    if [ -e "$config_file" ]; then
        ${srcdir}/support/utils/populate.py -b `pwd`/../${cmd} oco_oco2_sounding_${val}_test.config --l2_config=${config_file} && exit 0
    else
        ${srcdir}/support/utils/populate.py -b `pwd`/../${cmd} oco_oco2_sounding_${val}_test.config && exit 0
    fi
exit 1


