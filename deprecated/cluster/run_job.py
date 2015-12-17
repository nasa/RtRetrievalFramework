#!/usr/bin/env python

import os
import sys
import glob
import time
import stat
import copy
import shutil
import signal
import atexit
import dircache
from types import ListType, DictType

if 'CLUSTER_SW_DIR' in os.environ:
    sys.path.append(os.environ['CLUSTER_SW_DIR'])
else:
    sys.path.append(os.path.dirname(sys.argv[0]))

from L2_FP_Util import run_wait_command
from Cluster_Utils import *
from cluster_config import *

cleanup_list = []

def cleanup_job(signum=signal.SIGTERM, frame=None):
    # Cancel signal to avoid infinite loop
    signal.signal(signum, signal.SIG_DFL)

    if VERBOSE_RUN_JOB:
        print 'Executing cleanup'
        
    for file_to_clean in cleanup_list:
        if os.path.exists(file_to_clean):
            if VERBOSE_RUN_JOB:
                print 'Removing: %s' % file_to_clean

            if os.path.isdir(file_to_clean):
                shutil.rmtree(file_to_clean)
            else:
                os.remove(file_to_clean)

    sys.stdout.flush()
    sys.stderr.flush()

def clean_existing_files(file_list, verbose=False):
    for curr_glob in file_list:
        for curr_file in glob.glob(curr_glob):
            if os.path.isdir(curr_file):
                dir_file_list = [ '%s/%s' % (curr_file, dir_file) for dir_file in dircache.listdir(curr_file) ]
                clean_existing_files(dir_file_list, verbose)
            else:
                if verbose:
                    print 'Removing file: %s' % curr_file
                os.remove(curr_file)

def set_absco_path(run_filename, new_absco_path, remote_dir=None):
    run_fo = open(run_filename, 'r')

    inp_filename = None
    for curr_line in run_fo.readlines():
        if curr_line.find('input_file') >= 0:
            inp_filename = curr_line.split('=')[1].strip()

    run_fo.close()

    if inp_filename == None:
        exit_with_error('Could not input file keyword in run file "%s"' % run_filename, REMOTE_DIR)

    if not os.path.exists(inp_filename):
        exit_with_error('Input filename read from run file does not exist "%s"' % inp_filename, REMOTE_DIR)

    bak_inp_filename = inp_filename + '.orig'
    shutil.copyfile(inp_filename, bak_inp_filename)

    inp_fo = open(inp_filename, 'r')
    inpfile_lines = inp_fo.readlines()
    inp_fo.close()

    found_absco = False
    inp_fo = open(inp_filename, 'w')
    for curr_line in inpfile_lines:
        if curr_line.find('absco_path') >= 0:
            curr_absco_path = curr_line.split('=')[1].strip()
            found_absco = True
            
            if os.path.exists(curr_absco_path):
                print 'ABSCO path in %s exists: %s, will not overwrite' % (inp_filename, curr_absco_path)
                inp_fo.write(curr_line)
            else:
                print 'ABSCO path in %s does not exist: %s, setting to: %s' % (inp_filename, curr_absco_path, new_absco_path)
                inp_fo.write('absco_path = %s\n' % new_absco_path)

        else:
            inp_fo.write(curr_line)
    inp_fo.close()

    # Don't make this a fatal error. The code will ultimately complain if it cant find ABSCO
    if not found_absco:
        print >>sys.stderr, 'Could not find absco keyword in run file "%s" or input file "%s"' % (run_filename, inp_filename)
        
def get_run_mode(run_filename):
    run_mode = None
    
    run_fo = open(run_filename, 'r')
    for runfile_line in run_fo.readlines():
        if runfile_line.lower().find(RUN_MODE_KEYWORD.lower()) >= 0:
            line_parts = runfile_line.split('=')

            if line_parts[0].strip().find('#') < 0:
                value_parts = line_parts[1].split('#')
                run_mode = value_parts[0].strip().upper()
    run_fo.close()

    return run_mode

def run_l2_binary(binary_filename, stdout_filename, stderr_filename, parallel_prefix=None, verbose=False):

    l2_command = "%s %s > %s 2> %s" % (TIME_COMMAND, binary_filename, stdout_filename, stderr_filename)
    
    if parallel_prefix != None:
        run_command = ['/bin/sh', '-c', '%s %s' % (parallel_prefix, l2_command)]
    else:
        run_command = ['/bin/sh', '-c', l2_command]

    if verbose:
        print 'Running command: %s' % ' '.join(run_command)

    os.spawnvp(os.P_WAIT, run_command[0], run_command)

    stderr_fo = open(stderr_filename, 'r')

    # Once done check stderr for errors
    stderr_contents = []
    run_error = False
    for stderr_line in stderr_fo.readlines():
        for check_string in ERROR_CHECK_STRINGS:
            if stderr_line.find(check_string) >= 0:
                run_error = True
        stderr_contents.append(stderr_line)
    
    stderr_fo.close()
    
    # Output stderr if verbose and an error was found
    if run_error and verbose:
        print >>sys.stderr, 'ERROR detected in run of %s' % binary_filename
        for stderr_line in stderr_contents:
            print >>sys.stderr, LOG_VERBOSE_PREFIX, stderr_line,

   
def wait_for_sync(queue_file, num_simul, recheck_wait, max_tries, verbose=False):

    print_queue = False

    create_try = 0
    while not os.path.exists(queue_file) and create_try < max_tries:
        queue_fo = open(queue_file, 'w')
        queue_fo.close()
  
        create_try += 1

        time.sleep(5) # Let NFS settle

    if not os.path.exists(queue_file):
        exit_with_error('After %d tries could not create queue file: %s' % (create_try, queue_file))

    if create_try > 1 and verbose:
        print >>sys.stderr, 'Took %d tries to create queue file: %s' % (create_try, queue_file)

    # Get creation time from file
    my_stat = os.stat(queue_file)
        
    sel_queue_files = []
    sel_queue_times = []

    def time_compare(xval, yval):
        xidx = sel_queue_files.index(xval)
        yidx = sel_queue_files.index(yval)
        if sel_queue_times[xidx] == sel_queue_times[yidx]:
            return cmp( xval, yval )
        else:
            return cmp( sel_queue_times[xidx], sel_queue_times[yidx] )

    done_waiting = False
    check_try = 0
  
    while not done_waiting:
        beg_time_win = my_stat[stat.ST_MTIME] - (max_tries-check_try)*recheck_wait
        
        sel_queue_files = []
        sel_queue_times = []
        for curr_queue_basename in dircache.listdir(os.path.dirname(queue_file)):
            curr_queue_fullname = '%s/%s' % (os.path.dirname(queue_file), curr_queue_basename)
            if os.path.exists(curr_queue_fullname):
                # File could disapper on us and cause exception
                try:
                    curr_stat = os.stat(curr_queue_fullname)

                    if curr_stat[stat.ST_MTIME] >= beg_time_win:
                        sel_queue_files.append(curr_queue_basename)
                        sel_queue_times.append(curr_stat[stat.ST_MTIME])
                except OSError:
                    # Ignore disappearing file
                    pass

        # Not sure why this happens but sometimes the code may
        # not find the queue file it wrote earlier
        if not os.path.basename(queue_file) in sel_queue_files:
            sel_queue_files.append(os.path.basename(queue_file))
            sel_queue_times.append(my_stat[stat.ST_MTIME])
            
        sorted_queue_files = copy.copy(sel_queue_files)       
        sorted_queue_files.sort(time_compare)
        sorted_queue_times = copy.copy(sel_queue_times)
        sorted_queue_times.sort()

        queue_pos = sorted_queue_files.index(os.path.basename(queue_file))

        check_try += 1

        if verbose and print_queue:
            print 'Beginning of queue window: %d' % beg_time_win
            print 'Queue:'
            for queue_idx in range(len(sorted_queue_files)):
                if queue_idx == queue_pos:
                    print '-->',
                else:
                    print '-  ',
                print sorted_queue_files[queue_idx], sorted_queue_times[queue_idx]

        if queue_pos < num_simul:
            if verbose:
                print 'Done waiting for rsync, going ahead after attempt %d' % check_try
            done_waiting = True
        else:
            if verbose:
                print 'Time queue position is %d for attempt %d, waiting for %d' % (queue_pos, check_try, recheck_wait)
                
            time.sleep(recheck_wait)


def remove_queue_file(queue_file, verbose=False):
    
    # Catch exception in case accidentally removed or scrubbed by launch_jobs script
    try:
        os.remove(queue_file)
    except OSError:
        if verbose:
            print >>sys.stderr, 'Queue file removed from under us: %s' % queue_file

def determine_node_setup(node_filename):
    node_fo = open(node_filename)
    node_list = []
    node_hash = {}
    for node in node_fo.readlines():
        node = node.strip()
        node_list.append( node )
        node_hash[node] = 1
    node_fo.close()

    num_nodes    = len(node_list)
    num_machines = len(node_hash.keys())

    return (num_machines, num_nodes)

def launch_mpd_daemon(num_machines, node_filename, verbose):
    mpdboot_cmd = 'mpdboot --totalnum=%d -f %s' % (num_machines, node_filename)

    error_code = run_wait_command(mpdboot_cmd, verbose)

    if error_code > 0:
        exit_with_error('Error launching mpd daemon')


def stop_mpd_daemon(verbose):
    mpdkill_cmd = 'mpdallexit'
    
    error_code = run_wait_command(mpdkill_cmd, verbose)

    if error_code > 0 and verbose:
        print >>sys.stderr, 'Error killing MPD daemon'

## Main \/
if __name__ == "__main__":
    signal.signal(signal.SIGTERM, cleanup_job)
    atexit.register(cleanup_job)

    if 'REMOTE_DIR' not in os.environ or len(os.environ['REMOTE_DIR']) == 0:
        exit_with_error('REMOTE_DIR enviromental variable is undefined')

    REMOTE_DIR = os.environ['REMOTE_DIR']
    if REMOTE_DIR[len(REMOTE_DIR)-1] != '/':
        REMOTE_DIR += '/'

    for env_check in ['L2_BINARY', 'CLUSTER_HOSTNAME']:
        if env_check not in os.environ or len(os.environ[env_check]) == 0:
            exit_with_error('%s enviromental variable is undefined' % env_check, REMOTE_DIR)

    L2_BINARY        = os.environ['L2_BINARY']
    CLUSTER_HOSTNAME = os.environ['CLUSTER_HOSTNAME']
    REMOTE_MACHINE   = os.environ['REMOTE_MACHINE']

    if not os.path.exists(L2_BINARY):
        exit_with_error('Specified L2_BINARY does not exist: "%s"' % L2_BINARY, REMOTE_DIR)

    # These already checked for existance in cluster config by launch_jobs
    ABSCO_DIR      = ABSCO_DIRS[CLUSTER_HOSTNAME]
    TIME_QUEUE_DIR = RSYNC_TIME_QUEUE_DIR[CLUSTER_HOSTNAME]

    if os.environ.has_key('PARALLEL_MODE') and len(os.environ['PARALLEL_MODE']) > 0:
        parallel_mode = True
    else:
        parallel_mode = False

    if os.environ.has_key('LOCAL_MODE') and len(os.environ['LOCAL_MODE']) > 0:
        local_mode = True
    else:
        local_mode = False

    # Extract for appropriate mode
    if local_mode:
        if len(REMOTE_MACHINE) > 0:
            exit_with_error('Can not use local mode when REMOTE_MACHINE is defined: %s' % REMOTE_MACHINE, REMOTE_DIR)
        SYNC_BASE_DIRECTORY = None
    elif type(SCRATCH_BASE[CLUSTER_HOSTNAME]) is DictType:
        if parallel_mode:
            SYNC_BASE_DIRECTORY = SCRATCH_BASE[CLUSTER_HOSTNAME]['parallel']
        else:
            SYNC_BASE_DIRECTORY = SCRATCH_BASE[CLUSTER_HOSTNAME]['single']
    else:
        SYNC_BASE_DIRECTORY = SCRATCH_BASE[CLUSTER_HOSTNAME]

    # If no default destination is specified then do not sychronize
    if SYNC_BASE_DIRECTORY == None or len(SYNC_BASE_DIRECTORY) == 0:
        rsync_remote_files = False
    else:
        rsync_remote_files = True

    # If remote machine is not empty then we need to add a : between remote machine and remote dir
    # Otherwise we assume the files is accessible locally
    if REMOTE_MACHINE != None and len(REMOTE_MACHINE) > 0:
        if not rsync_remote_files:
            exit_with_error('rsyncing of remote files disabled when files exist on remote system: %s' % REMOTE_MACHINE, REMOTE_DIR) 
        
        REMOTE_LOCATION = "%s:%s" % (REMOTE_MACHINE, REMOTE_DIR)
        remote_opts = RSYNC_REMOTE_OPTS
    else:
        REMOTE_LOCATION = REMOTE_DIR
        remote_opts = ''


    if rsync_remote_files:
        # Use $ARRAY_INDEX to make run directory as uniq as possible or else possible
        # to have multiple nodes with same name and then overwrite each other
        jobid_string = ""
        for env_name in ['PBS_JOBID', 'LSB_JOBID', 'ARRAY_INDEX']:
            if env_name in os.environ:
                if VERBOSE_RUN_JOB:
                    print '%s = %s' % (env_name, os.environ[env_name])

                if len(jobid_string) > 0:
                    jobid_string += '_'
                jobid_string += os.environ[env_name]

        time_queue_file = '%s/%s/%s' % (TIME_QUEUE_DIR, REMOTE_MACHINE, jobid_string)

        cleanup_list.append(time_queue_file)

        # Directory where remote files are copied, use TMPDIR if cluster software sets it
        if 'TMPDIR' in os.environ and len(os.environ['TMPDIR']) != 0:
            scratch_dir = os.environ['TMPDIR']
        else:
            scratch_dir = '%s/%s.%s.%d' % (SYNC_BASE_DIRECTORY, os.path.basename(REMOTE_DIR.strip('/')), jobid_string, os.getpid())
            scratch_dir = scratch_dir.rstrip('/')

            # Cleanup if we create the directory but let PBS handle TMPDIR itself
            cleanup_list.append(scratch_dir)
    else:
        scratch_dir = REMOTE_DIR

    if VERBOSE_RUN_JOB:
        if 'HOSTNAME' in os.environ:
            print "HOSTNAME = %s" % os.environ['HOSTNAME']

        print "REMOTE_DIR = %s" % REMOTE_DIR
        print "L2_BINARY = %s" % L2_BINARY
        print ""

    if rsync_remote_files and not os.path.exists(SYNC_BASE_DIRECTORY):
        if VERBOSE_RUN_JOB:
            print "Creating SYNC_BASE_DIRECTORY: %s" % SYNC_BASE_DIRECTORY
        try:
            os.makedirs(SYNC_BASE_DIRECTORY)
        except:
            # Dir exists but for some reason check failed
            pass

        if not os.path.exists(SYNC_BASE_DIRECTORY):
            exit_with_error('SYNC_BASE_DIRECTORY not created: %s' % SYNC_BASE_DIRECTORY, REMOTE_DIR)

    if rsync_remote_files:
        if VERBOSE_RUN_JOB:
            print 'rsyncing input files from %s on host %s' % (REMOTE_DIR, REMOTE_MACHINE)

        wait_for_sync(time_queue_file, RSYNC_NUM_SIMUL, QUEUE_RETRY_WAIT, QUEUE_RETRY_COUNT, VERBOSE_RUN_JOB)

        rsync_files(REMOTE_LOCATION, scratch_dir, [RSYNC_INPUT_OPTS[CLUSTER_HOSTNAME], remote_opts], RSYNC_RETRY_WAIT, IN_RETRY_COUNT, verbose=VERBOSE_RUN_JOB)

        remove_queue_file(time_queue_file, verbose=VERBOSE_RUN_JOB)

    if not os.path.exists(scratch_dir):
        exit_with_error('scratch_dir not created: %s' % scratch_dir, REMOTE_DIR)

    if VERBOSE_RUN_JOB:
        print "Changing to scratch dir: %s" % scratch_dir

    os.chdir(scratch_dir)

    if not os.path.exists(RUN_FILENAME):
        exit_with_error('%s does not exist at %s' % (RUN_FILENAME, scratch_dir), REMOTE_DIR)

    if VERBOSE_RUN_JOB:
        print 'Cleaning any existing files'
    clean_existing_files(CLEAN_FILE_LIST, verbose=VERBOSE_RUN_JOB)

    # Fix ABSCO directory
    if VERBOSE_RUN_JOB:
        print 'Setting absco_path in %s in %s' % (ABSCO_DIR, RUN_FILENAME)

    set_absco_path(RUN_FILENAME, ABSCO_DIR, REMOTE_DIR)
    
    # Attempt to set up mpd for parallel mode
    if parallel_mode:
        node_filename = None
        if os.environ.has_key('PBS_NODEFILE') and len(os.environ['PBS_NODEFILE']) > 0:
            node_filename = os.environ['PBS_NODEFILE']
        elif os.environ.has_key('LSB_MCPU_HOSTS') and len(os.environ['LSB_MCPU_HOSTS']) > 0:
            node_filename = 'lsf_machine_file.tmp'
            run_wait_command('machinefile.lsf > %s' % node_filename, verbose=True)

        if node_filename == None:
            exit_with_error('No node filename could be found in enviromental variables', REMOTE_DIR)

        if not os.path.exists(node_filename):
            exit_with_error('Node filename does not exist: "%s"' % node_filename, REMOTE_DIR)

        (num_machines, num_nodes) = determine_node_setup(node_filename)

        print 'Using %d nodes on %d machines' % (num_nodes, num_machines)

        if PARALLEL_LAUNCH_MPD[CLUSTER_HOSTNAME]:
            if VERBOSE_RUN_JOB:
                print 'Launching mpd daemon'

            launch_mpd_daemon(num_machines, node_filename, VERBOSE_RUN_JOB)

        prefix_spec = PARALLEL_RUN_PREFIX[CLUSTER_HOSTNAME]
        if prefix_spec.count('%') == 1:
            parallel_prefix = prefix_spec % (num_nodes)
        else:
            parallel_prefix = prefix_spec
    else:
        parallel_prefix = None

    if VERBOSE_RUN_JOB:
        print 'Running L2 binary %s' % L2_BINARY
        if parallel_mode:
            print 'Using parallel mode for running binary'

    run_l2_binary(L2_BINARY, STDOUT_FILENAME, STDERR_FILENAME, parallel_prefix, verbose=VERBOSE_RUN_JOB)

    if parallel_mode and PARALLEL_LAUNCH_MPD[CLUSTER_HOSTNAME]:
        if VERBOSE_RUN_JOB:
            print 'Stopping mpd daemon'
        stop_mpd_daemon(VERBOSE_RUN_JOB)

    if rsync_remote_files:
        if VERBOSE_RUN_JOB:
            print 'rsyncing output files to %s on host %s' % (REMOTE_DIR, REMOTE_MACHINE)

        run_mode = get_run_mode(RUN_FILENAME)
        if not run_mode in SYNC_OUT_LIST:
            exit_with_error('Could not find run mode of case "%s" in those available for syncing "%s"' % (run_mode, ', '.join(SYNC_OUT_LIST.keys())) , REMOTE_DIR)

        # Build list of files to send back
        root_src_file_list = []
        subs_src_file_list = []
        for sync_glob in SYNC_OUT_LIST[run_mode]:
            for inp_sync_file in glob.glob(sync_glob):
                clean_sync_file = inp_sync_file.strip().rstrip('/')

                if clean_sync_file.find('./') == 0:
                    clean_sync_file = clean_sync_file.replace('./')

                if clean_sync_file.find('/') > 0:
                    subs_src_file_list.append(clean_sync_file)
                else:
                    root_src_file_list.append("%s/%s" % (scratch_dir, clean_sync_file))

        wait_for_sync(time_queue_file, RSYNC_NUM_SIMUL, QUEUE_RETRY_WAIT, QUEUE_RETRY_COUNT, VERBOSE_RUN_JOB)

        if len(root_src_file_list) > 0:
            rsync_files(' '.join(root_src_file_list), REMOTE_LOCATION, [RSYNC_OUTPUT_OPTS[CLUSTER_HOSTNAME], remote_opts], RSYNC_RETRY_WAIT, OUT_RETRY_COUNT, verbose=VERBOSE_RUN_JOB)

        for sync_sub_dir in subs_src_file_list:
            src_filename = '%s/%s' % (scratch_dir, sync_sub_dir)
            dst_filename = '%s/%s' % (REMOTE_LOCATION.rstrip('/'), sync_sub_dir)

            if os.path.isdir(src_filename):
                src_filename += '/'
                dst_filename += '/'

            rsync_files(src_filename, dst_filename, [RSYNC_OUTPUT_OPTS[CLUSTER_HOSTNAME], remote_opts], RSYNC_RETRY_WAIT, OUT_RETRY_COUNT, verbose=VERBOSE_RUN_JOB)

        remove_queue_file(time_queue_file, verbose=VERBOSE_RUN_JOB)

        # Get out of scratch dir before deleting
        if VERBOSE_RUN_JOB:
            print "Getting out of %s" % scratch_dir

        os.chdir(os.path.dirname(scratch_dir))

    # Exit called within
    cleanup_job()
