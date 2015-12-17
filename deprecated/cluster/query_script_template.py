#!/usr/bin/env python

# Do not run directly.
#
# Indexing is 1 based

import sys
import os
import re
import platform

py_minor_ver = int(platform.python_version_tuple()[1])

if py_minor_ver <= 3:
    import popen2
else:
    import subprocess


job_id = "{job_id}"
query_command = "{query_command}"
index_template = "{index_tmpl}"
run_directories = {run_directories}

def query_index(job_index, args=[]):
    run_command = [query_command, index_template % (job_id, job_index)] + args
    
    if py_minor_ver <= 3:
        run_instance = popen2.Popen4(run_command)
        stdout_data = run_instance.fromchild.readlines()
        status_code = run_instance.wait()
        return_code = os.WEXITSTATUS(status_code)
    else:
        p_obj = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
        (stdout_data, stderr_data) = p_obj.communicate()
        return_code = p_obj.returncode

    if return_code == 0:
        return stdout_data
    else:
        return None

    
for run_dir, run_idx in zip(run_directories, range(1,len(run_directories)+1)):
    do_query = False
    addl_args = []
    if len(sys.argv) < 2:
        do_query = True
    elif (sys.argv[1].isdigit() and int(sys.argv[1]) == run_idx) or (not sys.argv[1].isdigit() and re.search(sys.argv[1], run_dir)):
        do_query = True
        addl_args = sys.argv[2:]

    if do_query:
        query_result = query_index(run_idx, addl_args)
    else:
        query_result = None
        
    if query_result != None:
        print '%s:' % run_dir
        sys.stdout.writelines(query_result)
        print ''

    
