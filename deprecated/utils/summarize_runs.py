#!/usr/bin/env python

import os
import re
import sys
import subprocess
import time
import glob
import logging
from contextlib import contextmanager
from optparse import OptionParser
import logging
import traceback

import L2_Log_Util
from OCO_Matrix import OCO_Matrix

import stats
from testset_summary import testset_summary
from compute_runtimes import compute_runtimes

LOGGER_NAME = 'summarize_runs'
LOG_FILE_BASE = 'summarize_runs.log'

XVFB_BIN = "/usr/bin/Xvfb"
XVFB_CHECK_CMD = '/bin/ps w -u %s | grep Xvfb' % os.environ['USER']
XVFB_COMMAND_FMT = XVFB_BIN + " :%d -screen 0 1600x1200x16 -fbdir /tmp -pixdepths 8 16 32"
XVFB_DISPLAY_NUM_START = 5
XVFB_LAUNCH_MAX_ATTEMPTS = 5
RESULTS_DIR_BASENAME = "results"

NEVER_RAN_LIST_BASENAME = 'never_ran.list'
ALL_LIST_BASENAME = 'all.list'

RUN_SPEC_BASENAME = 'run_spec.dat'
SET_DIR_SEARCH_FILES = [ RUN_SPEC_BASENAME, 'input_config.dat', 'orbitsim_.*.config', 'oco_.*.config' ]

SUMM_FILE_FMT = "%s/summary_%s.dat"
STAT_FILE_FMT = "%s/stats_%s.dat"
PLOT_FILE_FMT = "%s/summary_plots_%s.ps"
TIME_FILE_FMT = "%s/timing_%s.dat"
PLOTS_DIR_FMT = "%s/plots_%s/"

SUMMARY_TYPES = ['converged', 'exceeded_max', 'failed', 'forward_model', 'jacobian_only']

HIST_SUMM_TYPES = { 'converged':    ['converged'],
                    'exceeded_max': ['max_iter', 'max_div'],
                    'failed'      : ['calc_error'],
                    }

SZA_HIST_FILE_TMPL = '%s/sza_histogram.ps'
SUMM_FILE_TMPL = '%s/summary_%s.dat'
LIST_FILE_TMPL = '%s/%s.list'
BRDF_MAP_FILE_TMPL = '%s/setup/sounding_id_brdf.map'

MAX_RECURSE_SEARCH_DEPTH = 5

DO_RECOMPILE = False

@contextmanager
def error_trapping(fail_on_error=True, processing_msg=None):
    logger = logging.getLogger('error_trapping')
    try:
        yield None
    except:
        if fail_on_error:
            raise
        else:
            if processing_msg != None:
                logger.error('ERROR: Failed during: %s' % processing_msg)
            for error_line in traceback.format_exception(*sys.exc_info()):
                if error_line != None:
                    logger.error(error_line.strip())

def launch_xvfb():
    logger = logging.getLogger('launch_xvfb')
    
    if not os.path.exists(XVFB_BIN):
        logger.error('Could not find Xvfb program, not changing DISPLAY environment')
        return None
    
    launch_try = 0
    display_index = None
    while display_index == None and launch_try < XVFB_LAUNCH_MAX_ATTEMPTS:
        logger.debug('Xvfb launch attempt %d' % (launch_try+1))
        
        check_obj = os.popen(XVFB_CHECK_CMD)
        check_results = check_obj.readlines()
        check_obj.close()

        xvfb_line_str = None
        for check_line in check_results:
            if check_line.find(XVFB_BIN) >= 0:
                xvfb_line_str = check_line
                break

        if xvfb_line_str != None:
            logger.debug('ps check line: %s' % xvfb_line_str.strip())
            disp_loc_beg = xvfb_line_str.find(XVFB_BIN) + len(XVFB_BIN)+1
            disp_loc_end = xvfb_line_str.find(' ', disp_loc_beg)

            check_str = xvfb_line_str[disp_loc_beg:disp_loc_end].replace(':','')
            logger.debug('Trying to use check string: %s' % check_str)
            
            try:
                display_index = int(check_str)
            except:
                check_results = None
            
            if check_results != None:
                logger.debug('Xvfb found running on display %d' % display_index)
                
        else:
        
            try_index = XVFB_DISPLAY_NUM_START + launch_try * 2

            logger.debug('Xvfb not running, attempting to load new Xvfb session for display %d' % try_index)
                
            xvfb_command = XVFB_COMMAND_FMT % try_index

            run_instance = subprocess.Popen(xvfb_command.split(), stdin=subprocess.PIPE)

            # Give Xvfb time to load
            time.sleep(2)

        # Keep track of how many times tried
        launch_try += 1
           
    return display_index

def summarize_recursively(search_base_dir, results_dir, summ_obj=None, tstplot_obj=None, vis_opts=None, hist_obj=None, overwrite=False, fail_on_error=False):
    logger = logging.getLogger(LOGGER_NAME)

    # Locate each base dir for a run set based on certain key files
    logger.info('Summarizing recursively from %s' % search_base_dir)
    for root, dirs, files in os.walk(search_base_dir, topdown=True):
        if len(root.split('/')) > len(search_base_dir.split('/')) + MAX_RECURSE_SEARCH_DEPTH:
            continue

        is_set_dir = False
        for curr_file in files:
            for search_re in SET_DIR_SEARCH_FILES:
                if re.search(search_re, curr_file):
                    is_set_dir = True
                    break
            if is_set_dir:
                break

        if is_set_dir:
            # Remove subdirs of root so they will not be followed
            while len(dirs) > 0: dirs.pop(0)

            if results_dir.find('%') >= 0:
                curr_result_dir = results_dir % root
            else:
                curr_result_dir = results_dir

            all_list_file = os.path.join(curr_result_dir, ALL_LIST_BASENAME)
            never_ran_file = os.path.join(curr_result_dir, NEVER_RAN_LIST_BASENAME)
            
            if overwrite or not os.path.exists(curr_result_dir) or os.path.exists(never_ran_file) or not os.path.exists(all_list_file):
                logger.info('Summarizing for set dir: %s' % root)
                summarize_set_dir(root, curr_result_dir, summ_obj=summ_obj, tstplot_obj=tstplot_obj, vis_opts=vis_opts, hist_obj=hist_obj, overwrite=overwrite, fail_on_error=fail_on_error)
            else:
                logger.info('Skipping completed set dir: %s' % root)


class sza_histogram():
    idl_routines = ['matrix_file__define', 'hist_l2_run_stats']

    def __init__(self, recompile=False, **addl_args):
        from IDL_Util import IDL_Util
        
        self.idl_util = IDL_Util()
        self.idl_util.idl_routines = self.idl_routines
        self.idl_util.load_idl(recompile)

    def plot(self, set_dir, results_dir, verbose=False):
        logger = logging.getLogger('sza_histogram')
        
        plot_file = SZA_HIST_FILE_TMPL % results_dir

        summ_files = []
        hist_list_types = []
        for curr_type, curr_lists in HIST_SUMM_TYPES.items():
            found_files = glob.glob(SUMM_FILE_TMPL % (results_dir, curr_type))
            if len(found_files) > 0:
                hist_list_types += curr_lists 
            for curr_filename in found_files:
                summ_files.append(curr_filename)

        if len(summ_files) == 0:
            logger.error('Could not find any summary files for histogram plot for set dir: %s' % set_dir)
            return

        list_files = []
        for curr_type in hist_list_types:
            for found_file in glob.glob(LIST_FILE_TMPL % (results_dir, curr_type)):
                list_files.append(found_file)

        # This check causes unnecessary problems when both max_div and max_iter are not present
        # for exceeded_max
        #if len(list_files) != len(hist_list_types):
        #    raise IOError('Could not find all list file types: %s in results dir: %s, found: %s' % (hist_list_types, results_dir, list_files))

        brdf_map_file = BRDF_MAP_FILE_TMPL % set_dir

        if not os.path.exists(brdf_map_file):
            raise IOError('Could not find brdf map file: %s' % brdf_map_file)

        for put_name in ['list_files', 'summ_files', 'brdf_map_file', 'set_dir', 'plot_file']:
            eval(' self.idl_util.idl.put("{0}", {0}) '.format(put_name))

        logger.info('Building SZA histogram plot')
        logger.info('summary files: %s' % summ_files)
        logger.info('list files: %s' % list_files)
        
        self.idl_util.idl.eval('hist_l2_run_stats, list_files, summ_files, brdf_map_file, set_dir, plot_file')

def summarize_set_dir(set_dir, results_dir, summ_obj=None, tstplot_obj=None, vis_opts=None, hist_obj=None, overwrite=False, fail_on_error=False,  verbose=False):
    logger = logging.getLogger(LOGGER_NAME)

    # Remove trailing slashes
    set_dir = set_dir.rstrip('/')
    results_dir = results_dir.rstrip('/')

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Open log file for set dir
    log_handler = L2_Log_Util.open_log_file( os.path.join(results_dir, LOG_FILE_BASE) )

    # Create necessary objects if not passed in
    if summ_obj == None:
        summ_obj = testset_summary()

    # Create status files
    run_lists = stats.gather_status(set_dir, results_dir)

    for (run_type, run_directories) in run_lists.items():
        if run_type not in SUMMARY_TYPES:
            continue

        input_files_only = run_type == 'forward_model' or run_type == 'jacobian_only'
        
        summ_file = SUMM_FILE_FMT % (results_dir, run_type)
        stat_file = STAT_FILE_FMT % (results_dir, run_type)
        plot_file = PLOT_FILE_FMT % (results_dir, run_type)
        time_file = TIME_FILE_FMT % (results_dir, run_type)

        plots_dir = PLOTS_DIR_FMT % (results_dir, run_type)

        processing_msg = 'generating summary files for run type: %s for set dir: %s' % (run_type, set_dir)
        with error_trapping(fail_on_error, processing_msg):
            summ_obj.make_summary_output(run_directories, stat_filename=stat_file, summ_filename=summ_file, overwrite=overwrite)

        ####

        if tstplot_obj != None:
            processing_msg = 'generating summary plots for run type: %s for set dir: %s' % (run_type, set_dir)
            with error_trapping(fail_on_error, processing_msg):
                tstplot_obj.make_summary_output(stat_filename=stat_file, summ_filename=summ_file, plot_filename=plot_file, overwrite=overwrite)

        ####

        if vis_opts != None and not input_files_only:
            import l2_vis
            if not os.path.exists(plots_dir):
                logger.info('Creating plots dir:', plots_dir)
                os.makedirs(plots_dir)

            processing_msg = 'creating per run directory plots for run type: %s for set dir: %s' % (run_type, set_dir)
            with error_trapping(fail_on_error, processing_msg):
                l2_vis.make_run_plots(run_directories, plots_dir, overwrite=overwrite, **vis_opts)

        ####

# Commented out as the new L2 code does not have the output that this routine reads anymore
#         processing_msg = 'generating summary timings info file for run type: %s for set dir: %s' % (run_type, set_dir)
#         with error_trapping(fail_on_error, processing_msg):
#             make_time_file = True
#             if not overwrite and os.path.exists(time_file):
#                 time_obj = OCO_Matrix()
#                 time_obj.read(time_file, read_data=False)

#                 if time_obj.dims[0] == len(run_directories):
#                     if verbose:
#                         logger.info('Skipping recreating timings file: %s' % time_file)
#                     make_time_file = False

#             if make_time_file:
#                 compute_runtimes(run_directories, time_file)
    ####

    if hist_obj != None:
        processing_msg = 'creating solar zenith angle histogram plot for set dir: %s' % set_dir
        with error_trapping(fail_on_error, processing_msg):
            hist_obj.plot(set_dir, results_dir, verbose)

    L2_Log_Util.close_log_file(log_handler)

def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] [run_dir_base]...")

    parser.add_option( "-o", "--overwrite", dest="overwrite",
                       default=False,
                       action="store_true",
                       help="overwrite any existing files instead of skipping")      

    parser.add_option( "-d", "--output_dir", dest="output_dir",
                       metavar="DIR",
                       type="string",
                       help="directory for outputting plots and summary files",
                       )

    parser.add_option( "-r", "--recurse", dest="recurse",
                       default=False,
                       action="store_true",
                       help="recurse from base dir looking for set directories")      

    parser.add_option( "-f", "--fail_on_error", dest="fail_on_error",
                       default=False,
                       action="store_true",
                       help="fail if any sub-programs fail")

    parser.add_option( "-p", "--make_plots", dest="make_plots",
                       default=False,
                       action="store_true",
                       help="make plots for run directories")

    parser.add_option( "-X", "--use_xvfb", dest="use_xvfb",
                       default=False,
                       action="store_true",
                       help="use X frame buffer in case no X server available")

    parser.add_option( "--use_idl", dest="use_idl",
                       default=False,
                       action="store_true",
                       help="use idl based routines, requires pyIDL Python module")

    parser.add_option( "--opt", dest="pass_option",
                       metavar="KEY=VALUE",
                       type="string",
                       action="append",
                       help="additional arguments to pass along to modules",
                       )

    # Initialize logging
    L2_Log_Util.init_logging()

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) >= 1:
        base_dir_list = args
    else:
        base_dir_list = [ './' ]

    addl_args = {}
    if options.pass_option != None:
        for arg_string in options.pass_option:
            if arg_string.find('=') < 0:
                parser.error('pass_option must be formated KEY=VALUE')
            (key, value) = arg_string.split('=')
            addl_args[key.strip()] = eval(value.strip())

    # Need a place to draw for IDL commands
    if options.use_xvfb:
        display_index = launch_xvfb()
    
        if display_index == None:
            raise Exception('Was unable to initialize or use xvfb')
        else:
            os.environ['DISPLAY'] = ':%d' % display_index

    summ_obj    = testset_summary(**addl_args)

    if options.make_plots:
        import l2_vis
        
        vis_opts = dict(idl_recompile=DO_RECOMPILE, do_sv_profiles=options.use_idl, **addl_args)
        
    else:
        vis_opts = None

    if options.make_plots and options.use_idl:
        from standard_testplots import standard_testplots
        tstplot_obj = standard_testplots(recompile=DO_RECOMPILE, **addl_args)

        hist_obj    = None
        # Disabled for now
        #hist_obj = sza_histogram(recompile=DO_RECOMPILE, **addl_args)
    else:       
        tstplot_obj = None
        hist_obj    = None
    
    for base_dir in base_dir_list:
        if options.output_dir == None:
            if options.recurse:
                output_dir = '%s/' + RESULTS_DIR_BASENAME
            else:
                output_dir = '%s/%s' % (base_dir.rstrip('/'), RESULTS_DIR_BASENAME)
        else:
            output_dir = options.output_dir
        
        if options.recurse:
            summarize_recursively(base_dir, output_dir, summ_obj=summ_obj, tstplot_obj=tstplot_obj, vist_opts=vis_opts, hist_obj=hist_obj, overwrite=options.overwrite, fail_on_error=options.fail_on_error)
        else:
            summarize_set_dir(base_dir, output_dir, summ_obj=summ_obj, tstplot_obj=tstplot_obj, vis_opts=vis_opts, hist_obj=hist_obj, overwrite=options.overwrite, fail_on_error=options.fail_on_error)


if __name__ == "__main__":
    standalone_main()

