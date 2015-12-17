#!/usr/bin/env python

# Load standard modules
import os
import sys
import copy
import glob
from optparse import OptionParser

# Load L2 modules
from OCO_TextUtils import *
from IDL_Util import IDL_Util

## For standard test plots
TITLE_FILTER_STRS = ['scratch', 'groups', 'algorithm', os.environ['USER'], 'runs', 'l2_fp']

MAX_LOG_SEARCH_DEPTH = 2
TITLE_SEP_STR = ' : '
MAX_TITLE_LEN = 80

class standard_testplots(IDL_Util):
    idl_routines = ['matrix_file__define', 'compile_dummies', 'standard_testplots']

    def __init__(self, recompile=False, **addl_args):
        IDL_Util.__init__(self)
        
        # Load IDL interpreter
        self.setup_idl()

        if not recompile and os.path.exists(self.sav_filename):
            self.load_sav()
        else:
            self.compile_sav()

    def is_orbit_sim_log(self, log_filename):
        log_obj = open(log_filename)

        found_log = False
        line_number = 1
        for line in log_obj.readlines():
            if line.find('# FILES MENU') >= 0:
                found_log = True
                break
            elif line_number > 20:
                break
            line_number += 1
               
        log_obj.close()

        return found_log           

    def get_common_path(self, run_dirs):
        common_path = os.path.commonprefix([ os.path.realpath(rd) for rd in run_dirs ])
        
        # Remove last common part in case run dirs have common part that is not a full
        # real path
        if not os.path.exists(common_path):
            common_path = os.path.dirname(common_path)
            
        return common_path

    def find_log_file(self, common_path):
        #full_run_dirs = [ os.path.realpath(rdir) for rdir in run_dirs ]
        
        for root, dirs, files in os.walk(common_path, topdown=True):
            if len(root.split('/')) > len(common_path.split('/')) + MAX_LOG_SEARCH_DEPTH:
                continue
                                         
            #for curr_dir in copy.copy(dirs):
                #if os.path.join(root, curr_dir) in full_run_dirs:
                #    dirs.remove(curr_dir)

            for curr_file in files:
                curr_file = os.path.join(root, curr_file)
                if curr_file[-4:] == '.log' and self.is_orbit_sim_log(curr_file):
                    return os.path.realpath(curr_file)
                elif curr_file[-4:] == '.hdf' and os.path.islink(curr_file):
                    link_loc = os.path.dirname(os.readlink(curr_file))
                    for hdf_dir_file in os.listdir(link_loc):
                        hdf_dir_file = os.path.join(link_loc, hdf_dir_file)
                        if hdf_dir_file[-4:] == '.log' and self.is_orbit_sim_log(hdf_dir_file):
                            return hdf_dir_file
        return ''

    def determine_plot_title(self, common_path):
        path_parts = common_path.split('/')
        path_parts.reverse()
        
        total_len = 0
        title_parts = []
        for common_part in path_parts:
            if len(common_part) > 0 and common_part.lower() not in TITLE_FILTER_STRS:
                total_len += len(common_part) + len(TITLE_SEP_STR)
                if total_len > MAX_TITLE_LEN:
                    break
                else:
                    title_parts.insert(0, common_part)

        title = TITLE_SEP_STR.join(title_parts)

        return title

    def make_summary_output(self, summ_filename, stat_filename, plot_filename, overwrite=False):

        # Put variables directly into IDL space or else direct to IDL
        # code with too many items in run_dirs will cause an error
        self.idl.put('stat_filename', stat_filename)
        self.idl.put('summ_filename', summ_filename)

        if overwrite:
            self.idl.put('overwrite', 1)
        else:
            self.idl.put('overwrite', 0)

        log_srch_base = os.path.split( os.path.dirname(summ_filename) )[0]
            
        self.idl.put('plot_filename', plot_filename)
        self.idl.put('log_file',      self.find_log_file(log_srch_base))
        self.idl.put('plot_title',    self.determine_plot_title(log_srch_base))
        
        self.idl.eval("standard_testplots, summ_filename, stat_filename, plot_filename, log_file, plot_title=plot_title")
        
               
def standalone_main():

    # Load command line options
    parser = OptionParser(usage="usage: %prog <summ_filename> <stat_filename>...")

    parser.add_option( "--plot_file", dest="plot_filename",
                       metavar="FILE", default="./summary.ps",
                       help="name of summary plot file other than default")

    parser.add_option( "--overwrite", dest="overwrite",
                       default=False,
                       action="store_true",
                       help="Overwrite existing summary and stats file instead of adding to existing contents"
                       )

    parser.add_option( "--recompile", dest="recompile",
                       default=False,
                       action="store_true",
                       help="recompile IDL routines regardless if .sav file already exists"
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    # Gather list of run directories to gather information from
    run_dirs = []

    if (len(args) < 2):
        parser.error('Need both summary and stat filenames')

    summ_filename = args[0]
    stat_filename = args[1]

    summ_plot = standard_testplots(recompile=options.recompile)
    summ_plot.make_summary_output(summ_filename, stat_filename, options.plot_filename, overwrite=options.overwrite)

if __name__ == "__main__":
    standalone_main()

