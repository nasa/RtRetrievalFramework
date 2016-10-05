#!/usr/bin/env python

import os
import sys
import re
import logging
from time import gmtime, strftime
from optparse import OptionParser

import numpy

import L2_Log_Util
from OCO_TextUtils import *
from OCO_Matrix import OCO_Matrix

LOGGING_NAME = 'testset_summary'

STDOUT_BASE = "stdout"
RUNFILE_BASE = "oco_l2.run"
SUMMARY_BASE = "out/summary.dat"

# Define the different types of labels
DEFAULT_LABELS = ["Testcase_Name", "Total_Runtime", "Num_Iter"]
PARALLEL_LABELS = ["Num_Nodes", "Avg_Node_Runtime", "Variance", "Min_Runtime", "Max_Runtime"]
ADDITIONAL_LABELS = ["Single_Scatter", "Polarization", "Num_Moments", "SV_Size"]


def compute_runtimes(run_dirs, output_data_file, additional_cols=False, parallel_cols=False):
    logger = logging.getLogger(LOGGING_NAME)
    
    # Names of runs for first column of timing file
    run_names = extract_run_names(run_dirs)

    default_data = numpy.zeros((len(run_dirs), len(DEFAULT_LABELS)), dtype=numpy.chararray)
    parallel_data = numpy.zeros((len(run_dirs), len(PARALLEL_LABELS)), dtype=numpy.chararray)
    additional_data = numpy.zeros((len(run_dirs), len(ADDITIONAL_LABELS)), dtype=numpy.chararray)
    row_idx = 0
    for (test_data_loc, tc_name) in zip(run_dirs, run_names):
        test_data_loc = test_data_loc.strip()
        tc_name = tc_name.strip()

        logger.debug("Gathering runtimes for %s" % test_data_loc)

        # Check if supplied location is the stdout filename itself or a path to
        # a test case directory underwhich it is expect to have a file called
        # stdout
        if test_data_loc.find(STDOUT_BASE) >= 0:
            stdout_file = test_data_loc
        else:
            stdout_file = test_data_loc + "/" + STDOUT_BASE

        # Do nothing if a stdout file does not exist
        if os.path.exists(stdout_file):
            # Collect stdout info
            stdout_f = open(stdout_file, 'r')
            so_lines = stdout_f.readlines()
            stdout_f.close()

            total_runtime = 0
            iter_count    = 0

            sum_sq_runtime = 0
            min_runtime    = sys.maxint
            max_runtime    = -sys.maxint - 1
            rt_line_count  = 0
            num_pf_mom     = 0

            retr_out_second_line = False
            for line in so_lines:
                if( line.find("RUNTIME RETRIEVAL") >= 0 ):
                    num_re = re.findall('([\d+.]+)', line)
                    curr_time = float(num_re[0])

                    total_runtime  = total_runtime + curr_time
                    sum_sq_runtime = sum_sq_runtime + curr_time*curr_time
                    min_runtime    = min(min_runtime, curr_time)
                    max_runtime    = max(max_runtime, curr_time)
                    rt_line_count  = rt_line_count + 1
                    continue

                if retr_out_second_line:
                    # Try and parse the next line as retr_out information, due to
                    # ifort seperating lines
                    retr_out_second_line = False
                    line_parts = line.split()
                    try:
                        iter_count = int(line_parts[1])
                    except IndexError:
                        iter_count = None
                    continue
                
                if( line.find("Retr_out_para") >= 0 ):
                    line_parts = line.split()
                    try:
                        iter_count = int(line_parts[8])
                    except IndexError:
                        iter_count = None
                        retr_out_second_line = True
                    continue

                if( line.find("max_pf_mom") >= 0 ):
                    line_parts = line.split()
                    num_pf_mom = int(line_parts[3])
                    continue

            if rt_line_count > 0:
                node_avg = total_runtime / rt_line_count
            else:
                node_avg = 0.0

            if(rt_line_count > 1):
                n_count = float(rt_line_count)
                node_var = (1 / (n_count-1)) * sum_sq_runtime - (n_count / (n_count-1)) * node_avg*node_avg
            else:
                node_var = 0.0

            # Collect additonal info
            single_scat = 'false'
            polarization = 'false'
            sv_size = 0
            if additional_cols:
                # Collect oco_l2.run info
                run_file = test_data_loc + "/" + RUNFILE_BASE
                if os.path.exists(run_file):
                    run_f = open(run_file, 'r')
                    run_lines = run_f.readlines()
                    run_f.close()

                    for line in run_lines:
                        if( line.find("single_scatter_correction") >= 0 ):
                            line_parts = line.split()
                            single_scat = line_parts[2].lower()

                        if( line.find("polarization") >= 0 ):
                            line_parts = line.split()
                            polarization = line_parts[2].lower()

                # Collect info from summary.dat file
                summary_file = test_data_loc + "/" + SUMMARY_BASE
                if os.path.exists(summary_file):
                    summ_f = open(summary_file, 'r')
                    summ_lines = summ_f.readlines()
                    summ_f.close()

                    for line in summ_lines:
                        if( line.find("Total") >= 0 ):
                            line_parts = line.split()
                            sv_size = int(line_parts[2])


            # Save data now in individual variables into matricies
            default_row = [tc_name, total_runtime, iter_count]
            for col_idx in range(len(default_row)):
                default_data[row_idx][col_idx] = default_row[col_idx]

            parallel_row = [rt_line_count, node_avg, node_var, min_runtime, max_runtime]
            for col_idx in range(len(parallel_row)):
                parallel_data[row_idx][col_idx] = parallel_row[col_idx]

            addl_row = [single_scat, polarization, num_pf_mom, sv_size]
            for col_idx in range(len(addl_row)):
                additional_data[row_idx][col_idx] = addl_row[col_idx]

            row_idx += 1


    # Put together the final output label list
    file_labels = copy.copy(DEFAULT_LABELS)
    if parallel_cols:
        file_labels += PARALLEL_LABELS
    if additional_cols:
        file_labels += ADDITIONAL_LABELS

    # Create a new data matrix for concatenated data
    file_data = numpy.zeros((len(run_dirs), len(file_labels)), dtype=numpy.chararray)

    # Concatenate various types of data
    for row_idx in range(len(run_dirs)):
        dflt_beg = 0
        dflt_end = len(DEFAULT_LABELS)
        file_data[row_idx][dflt_beg:dflt_end] = default_data[row_idx][:]

        par_end = dflt_end
        if parallel_cols:
            par_beg = dflt_end
            par_end = par_beg + len(PARALLEL_LABELS)
            file_data[row_idx][par_beg:par_end] = parallel_data[row_idx][:]

        if additional_cols:
            addl_beg = par_end
            addl_end = addl_beg + len(ADDITIONAL_LABELS)
            file_data[row_idx][addl_beg:addl_end] = additional_data[row_idx][:]

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = "Testset Timing Results"
    out_mat_obj.labels = file_labels
    out_mat_obj.format_col = [False]
    out_mat_obj.data = file_data

    out_mat_obj.write(output_data_file, auto_size_cols=True)


def standalone_main():
    # Load command line options
    parser = OptionParser(usage="usage: %prog [options] [run_dir_file]")

    parser.add_option( "-t", "--timings_file", dest="timings_file",
                       metavar="FILE",
                       default="timings.dat",
                       help="filename to output timing information other than default")

    parser.add_option( "-a", "--additional_info", dest="additional_info",
                       default=False,
                       action="store_true",
                       help="output additional run information columns")

    parser.add_option( "-p", "--parallel_info", dest="parallel_info",
                       default=False,
                       action="store_true",
                       help="output parallel statistics information columns"
                       )
    
    parser.add_option( "-v", "--verbose", dest="verbose",
                       default=False,
                       action="store_true",
                       help="Output more information to screen on processing"
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    # Initialize logging
    if options.verbose:
        L2_Log_Util.init_logging(logging.DEBUG)
    else:
        L2_Log_Util.init_logging(logging.INFO)

    if (len(args) < 1):
        parser.error('A file containing a list of run directories is needed')

    run_dirs_file = args[0]
    output_data_file = options.timings_file

    test_data_locs_f = open(run_dirs_file, 'r')
    run_dirs = test_data_locs_f.readlines()
    test_data_locs_f.close()

    compute_runtimes(run_dirs, output_data_file, additional_cols=options.additional_info, parallel_cols=options.parallel_info)
    

if __name__ == "__main__":
    standalone_main()
