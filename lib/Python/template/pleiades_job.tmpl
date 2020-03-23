#!/bin/bash

# These are just sample values, you will probably want to changes these
#PBS -l select=1:ncpus=12:model=has
#PBS -l walltime=60:00

# Stuff that varies
abscodir={0.abscodir}
abscoversion={0.abscoversion}
max_array_index={0.max_array_index}
basedir={0.processing_dir}

# After this, should not need to modify the script
cd $PBS_O_WORKDIR

# Determine number of nodes we have available by finding the number
# of unique names in the list of ssh nodes
nnode=`uniq $PBS_NODEFILE | wc -l`

# Copy absco data to every node local disk
parallel -j $nnode --nonall -u --sshloginfile $PBS_NODEFILE "mkdir -p $TMPDIR/absco && mcp -r $abscodir/$abscoversion $TMPDIR/absco"

# Now run each job. The -q has logging go to log file only, not stdout. 
# The -a points to a different location for the absco data

seq 0 $max_array_index | parallel -j $NCPUS -u --sshloginfile $PBS_NODEFILE "cd $PWD; $basedir/l2_fp_job.sh -q -a $TMPDIR/absco {{}}"

parallel -j $nnode --nonall -u --sshloginfile $PBS_NODEFILE "rm -r -f  $TMPDIR/absco"
