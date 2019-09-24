#! /bin/bash
# We wrap this command in a script so that we know we are using bash rather
# /bin/sh which may or may not be bash.
source ${srcdir}/setup_env.sh && PYTHONPATH=${srcdir}/support:${PYTHONPATH} && \
    cd fts_test && \
    ${srcdir}/support/utils/create_config.py -t fts \
    ${srcdir}/test/fts/input/pa20091103saaaaa_100223160344.008 \
    ${srcdir}/test/fts/input/pa20091103saaaab_100223160344.008 \
    ${srcdir}/test/fts/input/atmosphere_fts.dat \
    ${srcdir}/test/fts/input/tccon_runlog.grl && \
    ${srcdir}/support/utils/populate.py -b `pwd`/../${cmd} fts_fts_test.config && exit 0
exit 1
