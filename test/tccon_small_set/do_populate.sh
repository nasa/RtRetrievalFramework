#! /bin/bash
# We wrap this command in a script so that we know we are using bash rather
# /bin/sh which may or may not be bash.
source ${srcdir}/setup_env.sh && PYTHONPATH=${srcdir}/support:${PYTHONPATH} && \
    cd tccon_sounding_${val}_test && \
    ${srcdir}/support/utils/create_config.py -t gosat \
    ${srcdir}/test/tccon_small_set/acos_L1bB2900_tccon_5_good_qual.h5 \
    ${srcdir}/test/tccon_small_set/acos_EcmB2900_tccon_5_good_qual.h5 \
    ${srcdir}/test/tccon_small_set/acos_CldB2900_tccon_5_good_qual.h5
    config_file=${srcdir}/test/tccon_sounding_${val}/config.lua
    if [ -e "$config_file" ]; then
        ${srcdir}/support/utils/populate.py -b `pwd`/../${cmd} gosat_tccon_sounding_${val}_test.config --l2_config=${config_file} && exit 0
    else
        ${srcdir}/support/utils/populate.py -b `pwd`/../${cmd} gosat_tccon_sounding_${val}_test.config && exit 0
    fi
exit 1


