#!/usr/bin/env python

# Load standard modules
import os
import re
import sys
import glob
import logging
import traceback
from optparse import OptionParser

# Load additional modules
import numpy

# Load L2 modules
from L2_Plot_Util import *
import L2_Log_Util
import OCO_TextUtils

# Map plotting routines to globs for run directories to match files for routine
# This definition must be declared after the routines are declared
# These routines are generally placed in L2_Plot_Util and have the basic interface:
# (filenames, plot_title, legend_labels)
PLOT_ROUTINE_GLOBS = { 'out/control1/*_pd*.dat': plot_jacobians,
                       'out/rad_conv.dat':       plot_spec_residuals }

#################################

def get_path_groups(paths_to_grp):
    """Groups a set of paths based on a common basename"""

    # Can not group 1 or less items
    if len(paths_to_grp) <= 1:
        return []

    # Get list of tuples with basename for each input path along with index into original list
    base_dir_names = [ (curr_idx, os.path.basename(curr_dir)) for curr_idx, curr_dir in enumerate(paths_to_grp) ]

    found_groups = []
    while len(base_dir_names) > 0:
        # Pop off first item to match to other times in list
        srch_dir_info = base_dir_names.pop()

        # If we ever pop off the last item then we have a straggler that did not matching anything previously
        if len(base_dir_names) == 0:
            raise Exception('Exhausted path list without finding match for last item')

        # Results in a list with the size of the matching basename between the previously popped basename
        # and all other basenames
        matching_count = [ len(os.path.commonprefix( (srch_dir_info[1], curr_srch[1]) )) for curr_srch in base_dir_names ]

        # Find the largest matching length and where else that match occurs
        largest_matches = numpy.where(numpy.array(matching_count) == max(matching_count))
        matched_items = [ base_dir_names[curr_idx] for curr_idx in  largest_matches[0] ]

        # Now group together the items that have the same longest matching basename
        new_group = [ paths_to_grp[ srch_dir_info[0] ] ]
        for curr_item in matched_items:
            new_group.append( paths_to_grp[ curr_item[0] ] )
            base_dir_names.remove(curr_item)

        found_groups.append( new_group )

    return found_groups


def plot_group(run_directories, output_directory=None, filename_filters=None, overwrite=False, debugging_logs=False, **kwargs):
    """Plot files for a group of supposedly similar run directories"""

    logger = logging.getLogger()

    # Find common name among basenames of run directories
    dir_group_name = os.path.commonprefix( [os.path.basename(curr_dir.strip('/')) for curr_dir in run_directories] )

    logger.info('Ploting run group: ' + dir_group_name)

    for curr_glob, curr_routine in PLOT_ROUTINE_GLOBS.items():

        # For each type of plotting glob find matching files
        plot_files = []
        for curr_dir in run_directories:
            # glob may return multiple
            glob_results = glob.glob( os.path.join(curr_dir, curr_glob) )

            if filename_filters == None:
                plot_files += glob_results
            else:
                for curr_filter in filename_filters:
                    for curr_result in glob_results:
                        if re.search(curr_filter, curr_result):
                            plot_files.append(curr_result)

        # Group together globbed files according to basename
        # Safe way to only plot files that appear in all grouped dirs
        file_groups = get_path_groups(plot_files)

        for curr_group_files in file_groups:
            group_basename = os.path.splitext(os.path.basename(curr_group_files[0]))[0]
            group_title = '%s %s' % (dir_group_name, group_basename)

            output_filename = '%s_%s.%s' % (dir_group_name, group_basename, OUTPUT_FORMAT)
            if output_directory != None:
                if not os.path.exists(output_directory):
                    os.makedirs(output_directory)
                
                output_filename = os.path.join(output_directory, output_filename)

            if debugging_logs:
                log_filename = '%s.%s' % (os.path.splitext(output_filename)[0], 'log')
                logger.info('Writing log file: %s' % log_filename)
                log_handler = L2_Log_Util.open_log_file(log_filename)

            logger.info('Plotting file group: ' + group_basename)
            logger.debug('Group files:')
            for curr_file in curr_group_files:
                logger.debug(os.path.realpath(curr_file))
            
            # Remove common parts from end of filenames
            group_file_labels = [ os.path.dirname(re.sub(dir_group_name + '.*', '', curr_file)) for curr_file in curr_group_files ]

            # Remove common directory part from beginning of labels
            common_dir_prefix = os.path.commonprefix([ os.path.dirname(curr_lbl) for curr_lbl in group_file_labels])
            if len(common_dir_prefix) > 0:
                group_file_labels = [ curr_lbl.replace(common_dir_prefix, '') for curr_lbl in group_file_labels ]

            if os.path.exists(output_filename):
                logger.info('Skipping existing plot file: %s' % output_filename)
            else:
                logger.info('Creating plot file: %s' % output_filename)

                # Plot files using plotting routine associated with file glob
                try:
                    plotted_filenames = curr_routine(curr_group_files, group_title, group_file_labels, **kwargs)
                    combine_pdfs(plotted_filenames, output_filename, remove_source_files=True)

                    logger.debug('Combined %d plot files' % len(plotted_filenames))
                except:
                    logger.error('Failed to create plot file: %s' % output_filename)
                    logger.error('Error running routine %s for files %s:' % (curr_routine, curr_group_files))
                    logger.error(''.join(traceback.format_exception(*sys.exc_info(), limit=2)))

            if debugging_logs:
                L2_Log_Util.close_log_file(log_handler)

def plot_run_comparisons(run_directories, output_directory=None, **kwargs):
    """Plot comparisons for run directories that should be grouped together"""

    for curr_group in get_path_groups(run_directories):
        plot_group(curr_group, output_directory, **kwargs)

def standalone_main():
    # Load command line options
    parser = OptionParser(usage="usage: %prog [options] [summary_files]...")

    parser.add_option( "-r", "--run_dir_file", dest="run_dir_file",
                       metavar="FILE",
                       help="file to read list of run directories from")

    parser.add_option( "-o", "--output_dir", dest="output_dir",
                       metavar="DIR",
                       help="directory where to write plot files")

    parser.add_option( "-f", "--filter", dest="filter_names",
                       metavar="REGEXP",
                       type="string",
                       action="append",
                       help="filter matching filenames according to the supplied regular expression",
                       )

    parser.add_option( "-i", "--indexes", dest="indexes",
                       metavar="NUM",
                       type="string",
                       action="append",
                       help="which indexes that would generate multiple pages of plots to use",
                       )

    parser.add_option( "--log", dest="make_logs",
                       default=False,
                       action="store_true",
                       help="make debugging log files per plot file")


    # Parse command line arguments
    (options, args) = parser.parse_args()

    # Initialize logging
    L2_Log_Util.init_logging()

    # Gather list of run directories to gather information from
    run_dirs = []

    for arg_dir in args:
        run_dirs.append(arg_dir)

    if options.indexes != None:
        indexes = []
        for idx_arg in options.indexes:
            indexes += OCO_TextUtils.index_range_list(idx_arg)
    else:
        indexes = None

    if options.run_dir_file != None:
        if not os.path.exists(options.run_dir_file):
            parser.error("Run directory file '%s' does not exist" % options.run_dir_file)

        with open(options.run_dir_file, 'r') as run_dir_fh:
            run_dirs += [ file_dir.strip() for file_dir in run_dir_fh.readlines() ]

    plot_run_comparisons(run_dirs, output_directory=options.output_dir, filename_filters=options.filter_names, page_indexes=indexes, debugging_logs=options.make_logs)

if __name__ == "__main__":
    standalone_main()


