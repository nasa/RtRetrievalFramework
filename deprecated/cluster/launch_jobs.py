#!/usr/bin/env python

# Load standard modules
import os
import sys
import stat
import time
import socket
import dircache
from optparse import OptionParser
from types import DictType, ListType

from cluster_config import *
from Cluster_Utils import *

from L2_FP_Util import run_wait_command, is_valid_binary

def create_job_array_script(tmpl_filename, script_filename, **template_values):

    with open(tmpl_filename) as tmpl_obj:
        script_lines = [ curr_line.format(**template_values) for curr_line in tmpl_obj.readlines() ]

    with open(script_filename, 'w') as script_obj:
        script_obj.writelines(script_lines)   

def launch_job_array(run_dirs, chunk_count, submit_prefix, job_array_opt, query_script_config, execute=False):
    ja_output_filename = os.path.abspath(JOB_ARRAY_SCRIPT_OUTPUT % ("%d_%d" % (os.getpid(), chunk_count)))

    print 'Writing job array script: %s' % ja_output_filename
    create_job_array_script(JOB_ARRAY_SCRIPT_TEMPLATE, ja_output_filename, run_directories=run_dirs, run_script=RUN_SCRIPT)
    
    ja_opts = job_array_opt % (1, len(run_dirs))

    sub_exec = "%s %s %s" % (submit_prefix, ja_opts, ja_output_filename)
    job_name = "%d index job array" % len(run_dirs)
    if execute:
        print 'Launching %s' % job_name
        os.chmod(ja_output_filename, stat.S_IRWXU)

        with os.popen(sub_exec) as sub_o:
            submit_lines = sub_o.readlines()
            sys.stderr.writelines(submit_lines)

        job_id, job_host = submit_lines[0].split('.', 1)
        job_id = job_id.replace('[]', '') # Remove PBS pro job arrays nomenclature
        print 'Job ID: %s' % job_id

        # Remove run script since the scheduler should have copied it from us
        os.remove(ja_output_filename)

        query_script = QUERY_SCRIPT_OUTPUT % job_id
        print 'Creating query script: %s' % query_script
        create_job_array_script(QUERY_SCRIPT_TEMPLATE, query_script, job_id=job_id, run_directories=run_dirs, **query_script_config)
        
    else:
        print 'Launch command for: %s' % job_name
        print sub_exec

def launch_individual(run_dirs, submit_prefix, single_opt, parallel_opt=None, execute=False):
    for curr_dir in run_dirs:
        os.environ['REMOTE_DIR'] = curr_dir

        if parallel_opt != None:
            sub_exec = "%s %s %s" % (submit_prefix, parallel_opt, RUN_SCRIPT)
        else:
            sub_exec = "%s %s %s" % (submit_prefix, single_opt, RUN_SCRIPT)
        
        job_name = curr_dir[-60:]
        if execute:
            print 'Launching ...%s' % job_name
            
            sub_o = os.popen(sub_exec)
            for line in sub_o.readlines():
                print line,
            sub_o.close()
        else:
            print 'Launch command for: ...%s' % job_name
            print sub_exec

def check_for_remote_dir(remote_machine, remote_dir):
    if remote_machine != None and len(remote_machine) > 0:
        test_command = '%s %s test -d "%s"' % (SSH_EXEC, remote_machine, remote_dir)
        print "Sanity checking on %s for existance of: %s" % (remote_machine, remote_dir)
    else:
        test_command = 'test -d "%s"' % (remote_dir)
        print "Sanity checking for existance of: %s" % (remote_dir)

    error_code = run_wait_command(test_command, True)
    
    if error_code == 0:
        does_exist = True
    else:
        does_exist = False
        
    return does_exist

def process_run_spec(run_spec_file, run_dirs):
    rs_binary = None
    rs_work_path = None

    if run_spec_file != None:
        if not os.path.exists(run_spec_file):
            parser.error("Run directory file '%s' does not exist" % run_spec_file)

        run_spec_obj = open(run_spec_file)

        in_header = False
        for spec_line in run_spec_obj.readlines():
            spec_line = spec_line.strip()

            if len(spec_line) == 0 or spec_line.find('#') == 0:
                continue

            if spec_line.lower().find('begin HEADER'.lower()) == 0:
                in_header = True
            elif spec_line.lower().find('end HEADER'.lower()) == 0:
                in_header = False
            elif in_header:
                head_parts = [ pstr.strip() for pstr in spec_line.split('=') ]

                if len(head_parts) != 2:
                    raise ValueError('Invalid keyword value pair string "%s" in run spec file %s' % (spec_line, run_spec_file))

                (key_name, key_val) = head_parts

                if key_name.lower().find('binary_filename') == 0:
                    rs_binary = key_val
                elif key_name.lower().find('work_path') == 0:
                    rs_work_path = key_val
            else:
                run_dirs.append(spec_line)
                
    return (rs_binary, rs_work_path)

def check_cluster_config(cluster_name, use_parallel):
    
    hash_check_vars = ['SUBMIT_EXEC', 'ABSCO_DIRS', 'WORK_REMOTE_MACHINE', 'RUNS_REMOTE_MACHINE', 'SCRATCH_BASE', 'RSYNC_TIME_QUEUE_DIR', 'SINGLE_MODE_OPT', 'QUERY_SCRIPT_CONFIG']
    if use_parallel:
        hash_check_vars.append('PARALLEL_SUBMIT_OPT')
        hash_check_vars.append('PARALLEL_RUN_PREFIX')
        hash_check_vars.append('PARALLEL_PROC_PER_NODE')
        hash_check_vars.append('PARALLEL_LAUNCH_MPD')
        
    for hash_name in hash_check_vars:
        if not cluster_name in eval("%s.keys()" % hash_name):\
           raise RuntimeError('%s not found for cluster "%s"' % (hash_name, cluster_name))

    if not cluster_name in JOB_ARRAY_OPT.keys() and options.job_array:
        raise RuntimeError('JOB_ARRAY_OPT not found or supported for cluster "%s"' % cluster_name)

    print 'Recognized cluster as: %s' % cluster_name

def get_cluster_submit_prefix(cluster_name):

    submit_prefix = None
    if type(SUBMIT_EXEC[cluster_name]) is DictType and options.queue == None:
        print 'queue must be specified on cluster %s, valid options: %s' % (cluster_name, SUBMIT_EXEC[cluster_name].keys())
        sys.exit()
    elif type(SUBMIT_EXEC[cluster_name]) is DictType and not options.queue in SUBMIT_EXEC[cluster_name]:
        print 'queue name "%s" not found for cluster %s, valid options: %s' % (options.queue, cluster_name, SUBMIT_EXEC[cluster_name].keys())
        sys.exit()
    elif type(SUBMIT_EXEC[cluster_name]) is DictType and options.queue in SUBMIT_EXEC[cluster_name]:
        print 'Using queue: %s' % options.queue
        submit_prefix = SUBMIT_EXEC[cluster_name][options.queue]
    else:
        submit_prefix = SUBMIT_EXEC[cluster_name]

    return submit_prefix

def get_remote_machine_ip(cmd_line_name, config_remote_machine):

    if cmd_line_name != None:
        remote_name = cmd_line_name
        print 'Using command line specified remote machine "%s"' % remote_name
    else:
        remote_name = config_remote_machine
        print 'Using default remote machine "%s"' % remote_name

    if len(remote_name) > 0:
        if remote_name.find('@') >= 0:
	    username, remote_hostname = remote_name.split('@', 1)
        else:
	    username = None
            remote_hostname = remote_name

        remote_ip = socket.gethostbyname(remote_hostname)

        if username != None:
            remote_ip = '%s@%s' % (username, remote_ip)

        print 'Resolved remote machine %s as %s' % (remote_name, remote_ip)
    else:
        # Local, no remote machine
        print 'Using local filesystem'
        remote_ip = ''

    return remote_ip

def init_remote_queue_dir(remote_machine):
    remote_queue_dir = '%s/%s' % (RSYNC_TIME_QUEUE_DIR[hostname], remote_machine)
    if not os.path.exists(remote_queue_dir):
        print 'Creating time queue dir: %s' % remote_queue_dir
        os.makedirs(remote_queue_dir)
    else:
        # Clean old queue files here
        remove_count = 0
        for curr_queue_file in dircache.listdir(remote_queue_dir):
            full_queue_filename = '%s/%s' % (remote_queue_dir, curr_queue_file)
            curr_stat = os.stat(full_queue_filename)
            if curr_stat[stat.ST_MTIME] < (time.time() - (24*60)):
                remove_count += 1
                if os.path.isfile(full_queue_filename):
                    os.remove(full_queue_filename)
                    remove_count += 1
        if remove_count > 0:
            print 'Removed %d old time queue files' % remove_count

def resolve_work_path(opt_work_path, rs_work_path):

    if 'L2_WORK_PATH' in os.environ:
        env_work_path = os.environ['L2_WORK_PATH']
    else:
        env_work_path = None

    check_names = ['environment', 'command line', 'run spec file']
    check_values = [ env_work_path, opt_work_path, rs_work_path ]

    work_path = None
    prev_name = None
    for (curr_src_name, curr_path_val) in zip(check_names, check_values):
        if curr_path_val != None:
            if work_path != None and work_path.rstrip('/') != curr_path_val.rstrip('/'):
                raise ValueError('Work path from %s with value "%s" does not match path from %s with value "%s"' % (prev_name, work_path, curr_src_name, curr_path_val))
            else:
                work_path = curr_path_val
                prev_name = curr_src_name

    return work_path

def resolve_binary_name(opt_binary, rs_binary):

    binary = None

    if opt_binary != None:
            binary = opt_binary
    elif rs_binary != None:
        if os.path.exists(rs_binary):
            binary = rs_binary
        else:
            raise IOError('Can not use non-existant run spec binary %s, specify binary on command line' % rs_binary)

    if binary == None:
        raise ValueError('No binary specified through run spec or on command line')

    if not os.path.exists(binary):
        raise IOError('specified binary does not exist: "%s"' % binary)

    return binary

def sync_work_path(work_path, remote_machine, hostname):
    # At this point work_path is the same on remote and dest system

    if len(remote_machine) == 0:
        print >>sys.stderr, 'Work path is local not syncing'
        return

    if not check_for_remote_dir(remote_machine, work_path):
        raise IOError('work_path %s does not exist on remote machine %s' % (work_path, remote_machine))

    if not os.path.exists(work_path):
        os.mkdir(work_path)

    remote_loc = '%s:%s/' % (remote_machine, work_path)

    rsync_files(remote_loc, work_path, [RSYNC_INPUT_OPTS[hostname], RSYNC_REMOTE_OPTS, '-v'], wait_time=60,  max_retry_count=1, verbose=True)

if __name__ == "__main__":
    # Load command line options
    parser = OptionParser(usage="""usage: %prog [options] [<run_dir> <run_dir> ...]

Run directories can be specified optionally on the command line in addition or
instead of using a file with a list of run directories. Run directories must
always be fully qualified paths.""")

    parser.add_option( "-r", "--run_spec", dest="run_spec",
                       metavar="FILE",
                       help="file containing list of remote run directories to launch. may optional header information")

    parser.add_option( "-w", "--work_path", dest="work_path",
                       metavar="FILE",
                       help="work path that matches remote system if not specfied in run_spec")

    parser.add_option( "-q", "--queue", dest="queue",
                       help="which queue to use on target system if supported"
                       )

    parser.add_option( "-b", "--binary", dest="binary",
                       help="location of L2 FP binary to use for jobs"
                       )

    parser.add_option( "-m", "--remote_machine", dest="remote_machine",
                       help="computer system to sync run directories with other than default"
                       )

    parser.add_option( "-i", "--individual", dest="job_array",
                       default=True,
                       action="store_false",
                       help="Submit directories individually instead of as job array"
                       )

    parser.add_option( "-l", "--local", dest="local",
                       default=False,
                       action="store_true",
                       help="run job in location where it lives instead of copying to local scratch"
                       )

    parser.add_option( "-p", "--parallel", dest="parallel",
                       metavar="NUM",
                       type="int",
                       default=0,
                       help="run parallel mode binary with NUM instances"
                       )

    parser.add_option( "-c", "--chunk", dest="chunk_size",
                       metavar="NUM",
                       type="int",
                       help="split run dirs into chunks for job array submission"
                       )

    parser.add_option( "-d", "--dry_run", dest="execute",
                       default=True,
                       action="store_false",
                       help="Do not really submit jobs"
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    # Gather list of remote run directories to launch
    run_dirs = []

    if (len(args) > 0):
        for arg_dir in args:
            # Expand path, if this does not exist locally then it must
            # be a remote path that will be checked later
            expanded_path = os.path.realpath(os.path.expanduser(arg_dir))
            if os.path.exists(expanded_path):
                run_dirs.append(expanded_path)
            else:
                run_dirs.append(arg_dir)

    # Load run spec and possibly the binary specified there and work path
    (rs_binary, rs_work_path) = process_run_spec(options.run_spec, run_dirs)

    # If no directories are specified on command line and or run spec
    # can not continue
    if len(run_dirs) == 0:
        parser.error('at least one run directory must be specified')

    # Cluser name is hostname without full domain
    hostname = os.environ['HOSTNAME'].split('.')[0]

    # Whether or not to use parallel mode when running
    use_parallel = options.parallel > 0
    
    if use_parallel:
        options.job_array = False

    # Make sure 
    check_cluster_config(hostname, use_parallel)

    # Get prefix for command line to call cluster queuing system
    submit_prefix = get_cluster_submit_prefix(hostname)

    # Find where we are rsyncing files from
    runs_remote_machine = get_remote_machine_ip(options.remote_machine, RUNS_REMOTE_MACHINE[hostname])

    # Make rsync time queue directory if it does not exists
    init_remote_queue_dir(runs_remote_machine)

    # Resolve binary and work path to use
    try:
        work_path = resolve_work_path(options.work_path, rs_work_path)
        binary = resolve_binary_name(options.binary, rs_binary)
    except:
        exctype, value = sys.exc_info()[:2]
        if value == None:
            parser.error(exctype)
        else:
            parser.error(value)

    # Make sure the filename passed to us really is a L2 FP binary
# Name of binary has changed, so don't bother checking this.
#    if not is_valid_binary(binary, use_parallel):
#        parser.error('%s is not a valid binary, parallel mode = %s' % (binary, use_parallel))

    # Check that at least first run dir can be retrieved, saves time debugging
    # why jobs are failing due to bad run_spec setup
    if not check_for_remote_dir(runs_remote_machine, run_dirs[0]):
        if runs_remote_machine != None and len(runs_remote_machine) > 0:
            raise IOError('Remote directory %s does not exist on %s' % (run_dirs[0], runs_remote_machine))
        else:
            raise IOError('Local directory %s does not exist' % (run_dirs[0]))

    # Synchronize the work path if it is set
    if work_path != None:
        sync_work_path(work_path, WORK_REMOTE_MACHINE[hostname], hostname)

    # Set invariant environmental variables
    os.environ['L2_BINARY']        = binary
    os.environ['CLUSTER_HOSTNAME'] = hostname
    os.environ['CLUSTER_SW_DIR']   = os.path.dirname(sys.argv[0])
    os.environ['REMOTE_MACHINE']   = runs_remote_machine

    if use_parallel:
        os.environ['PARALLEL_MODE'] = 'TRUE'

        if (options.parallel % PARALLEL_PROC_PER_NODE[hostname]) > 0:
            parser.error('Number of parallel instances must be a multiple of %d' % PARALLEL_PROC_PER_NODE[hostname])
        num_nodes = (options.parallel / PARALLEL_PROC_PER_NODE[hostname])

        submit_opt = PARALLEL_SUBMIT_OPT[hostname]

        if submit_opt.count('%') == 2:
            parallel_opt = submit_opt % (num_nodes, PARALLEL_PROC_PER_NODE[hostname])
        elif submit_opt.count('%') == 1:
            parallel_opt = submit_opt % (num_nodes)
        else:
            raise ValueError('Improper submit options string: "%s" for hostname: %s' % (submit_opt, hostname))
    else:
        parallel_opt = None

    if options.local:
        print 'Running jobs locally at run directories'
        os.environ['LOCAL_MODE'] = 'TRUE'

    # Launch jobs possibly in chunks
    if options.chunk_size == None:
        options.chunk_size = len(run_dirs)

    chunk_count = 1
    while len(run_dirs) > 0:
        chunk_amount = min(len(run_dirs),options.chunk_size)
        chunk_dirs = run_dirs[0:chunk_amount]

        if len(run_dirs) >= chunk_amount:
            run_dirs = run_dirs[chunk_amount:]

	# Removed not passed to node and it can be set by
	# the node. This is necessary since we pass all 
	# enviromental variables from the head node to
	# the processing nodes
	del os.environ['HOSTNAME']

        if options.job_array:
            launch_job_array(chunk_dirs, chunk_count, submit_prefix, JOB_ARRAY_OPT[hostname], QUERY_SCRIPT_CONFIG[hostname], options.execute)
        else:
            launch_individual(chunk_dirs, submit_prefix, SINGLE_MODE_OPT[hostname], parallel_opt, options.execute)

        chunk_count += 1
