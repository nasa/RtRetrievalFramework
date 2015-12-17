#!/usr/bin/env python

# Do not run directly. Script is appended with run list
# then submitted. This way if the run list is deleted
# the jobs can still run as the script is copied by the 
# queueing system
#
# Indexing is 1 based

import os

run_script = "{run_script}"
run_directories = {run_directories}

if 'PBS_ARRAY_INDEX' in os.environ and os.environ['PBS_ARRAY_INDEX'].isdigit():
    # PBS Pro
    array_index = int(os.environ['PBS_ARRAY_INDEX']) - 1
elif 'PBS_ARRAYID' in os.environ and os.environ['PBS_ARRAYID'].isdigit():
    # Open PBS
    array_index = int(os.environ['PBS_ARRAYID']) - 1
elif 'LSB_JOBINDEX' in os.environ and os.environ['LSB_JOBINDEX'].isdigit():
    # LSF
    array_index = int(os.environ['LSB_JOBINDEX']) - 1
else:
    raise 'array index could not be found from environmental variables'

os.environ['ARRAY_INDEX'] = '%d' % array_index
os.environ['REMOTE_DIR'] = run_directories[array_index]
os.execvp(run_script, [''])
