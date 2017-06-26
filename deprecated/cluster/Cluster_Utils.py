import os
import time
from types import DictType, ListType

from cluster_config import *

from L2_FP_Util import run_wait_command, BOGUS_ERROR_CODE

def exit_with_error(error_msgs, remote_dir=None):
    if not type(error_msgs) is ListType:
        error_msgs = [error_msgs]

    for curr_msg in error_msgs:
        print >>sys.stderr, curr_msg
        
    if remote_dir != None:
        print >>sys.stderr, 'ERROR for REMOTE_DIR = %s' % remote_dir

    sys.exit(1)

def rsync_files(src_location, dest_location, rsync_opts, wait_time, max_retry_count, verbose=False):

    if type(rsync_opts) is ListType:
        opts_string = ' '.join(rsync_opts)
    else:
        opts_string = rsync_opts
        
    retry = 0
    error_code = BOGUS_ERROR_CODE # bogus code to get loop started
    while retry < max_retry_count:
        if retry > 0:
            time.sleep(wait_time)

        rsync_command = RSYNC_BIN + " " + opts_string + " " + src_location + " " + dest_location
        if verbose:
            print "Running rsync command: %s" % rsync_command

        error_code = run_wait_command(rsync_command, verbose)

        if verbose:
            # In case we need to debug which error codes are not
            # being used as errors
            print "rsync error code: %d" % error_code
            sys.stdout.flush()

        if error_code not in RSYNC_RETRY_EXIT_CODES:
            break
        
        retry += 1

    if error_code > 0 or retry > 1:
        error_strings = []

        if error_code not in RSYNC_RETRY_EXIT_CODES:
            error_strings.append( "rsync failed after %d tries" % retry )
            error_strings.append( "rsync exit code: %d" % error_code )
        else:
            error_strings.append( "rsync succeeded after %d tries" % retry )
            
        error_strings.append( "rsync source: %s" % src_location )
        error_strings.append( "rsync destination: %s" % dest_location )

        if error_code not in RSYNC_RETRY_EXIT_CODES:
            exit_with_error(error_strings)
        else:
            for err_line in error_strings:
                print >>sys.stderr, err_line
