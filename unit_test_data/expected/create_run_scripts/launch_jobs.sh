#!/bin/bash
if [[ $# -lt 1 ]]; then
   echo "Usage: $0 <arguments to qsub>"
   echo ""
   echo 'At a minimum, you need to supply the queue to use, e.g., "-q long"'
   exit 1
fi
log_directory="create_run_scripts_test/log"
use_subdir=False
if [ "$use_subdir" = "True" ]; then
    log_directory="create_run_scripts_test/log/qsub"
fi
jid=`qsub $* -S /bin/bash -j oe -o $log_directory -N l2_fp -J 0-4 create_run_scripts_test/l2_fp_job.sh`
do_aggregate=False
if [[ "$do_aggregate" == "True" ]]; then
   # Wait for job to get into torque, otherwise dependency doesn't seem
   # to work right (this appears to be a bug, as of 10/4/2011, but an
   # easy one to work around). Don't think this is an issue for pbs_pro
   # or pleiades, but sleep doesn't hurt them and easier to treat things
   # the same.
   sleep 5                      
   jid2=`qsub $* -S /bin/bash -j oe -o ${log_directory}/qsub_aggregate.log -N aggregate -W depend=afterany:$jid create_run_scripts_test/aggregate.sh`
fi
# The buildbot regression scripts depend on the output from qsub, so go ahead and write that
# out
echo $jid
