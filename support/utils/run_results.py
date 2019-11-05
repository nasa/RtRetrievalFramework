#!/usr/bin/env python

from __future__ import print_function
from past.builtins import cmp
from builtins import str
from builtins import object
import sys
import os
import re
import inspect
import glob
from optparse import OptionParser

if int(sys.version[0]) == 3:
    DictType = dict
    StringType = bytes
    import functools
else:
    from types import DictType, StringType

import h5py

from full_physics.l2_input import L2InputFile

from full_physics.text_util import BackwardsReader

# Regular expression used for searching from current path
# for run ids. The first group should parse out the run id
RUN_ID_SEARCH = "sounding_id.list"

CONFIG_RUN_ID_SEC_RE = "input->.*FullPhysics->LIST->VALUES"

# Place limits on the file find to speed things up
SRCH_IGNORE_DIRS = ["config", "input", "populate"]
DEFAULT_MAX_SRCH_DEPTH = 2

# These are globs as well as will be formated with {run_id} and {base_path} located under base path
HDF_OUTPUT_FILENAME_SRCH = "output/l2_{run_id}.h5"
HDF_ERROR_FILENAME_SRCH = "output/l2_{run_id}.h5.error"
LOG_OUTPUT_FILENAME_SRCH = "log/l2_{run_id}.log"
LOG_RUNNING_FILENAME_SRCH = "log/l2_{run_id}.log.running"

# Look into output HDF file for state of execution
OUTCOME_FLAG_DATASET = "RetrievalResults/outcome_flag"

# String output into L2 code meaning that any error
# causing premature execution was handled properly
HANDLED_EXCEPTION_STRS = [r"Exception thrown by Full Physics code" ]
FILE_EXCEPTION_STRS = [r"While trying to open file '[^']*' a HDF 5 Exception thrown:"]

# Use state vector names to add some debugging counts of configuration type
STATEVECTOR_NAMES_DATASET = "RetrievedStateVector/state_vector_names"

# The search strings are not regular expressions, just used in a call to "".find() 
SV_CONFIG_TYPE_STRINGS = { "lambertian": "Ground Lambertian",
                           "coxmunk"   : "Ground Coxmunk",
                           }

# Only look at most 20 lines from end, do not search whole file
LOG_EXCEPTION_SRCH_LIMIT = 35

class results_count(object):
    def __init__(self, verbose=False):
        self.completed_counts = { 
                                  "outcome_1"   : {"num": 0, "prefix": "with", "order": 3, "indent": 4, "add_to_total": False},
                                  "outcome_2"    : {"num": 0, "prefix": "with", "order": 4, "indent": 4, "add_to_total": False},

                                  "converged"              : {"num": 0, "order": 0},
                                  "max_iter"               : {"num": 0, "order": 1, "alt": "exceeded maximum iterations"},
                                  "max_div"                : {"num": 0, "order": 2, "alt": "exceeded maximum divergence"},
                                  }

        self.cfg_type_counts = { "lambertian"     : {"num": 0, "indent": 4, "add_to_total": False, "alt": "completed lambertian type"},
                                 "coxmunk"        : {"num": 0, "indent": 4, "add_to_total": False, "alt": "completed coxmunk type"},
                                 }

        self.error_counts   = { "unknown_quality" : {"num": 0, "prefix": "with"},
                                "handled_error"   : {"num": 0, "alt": "handled errors"},
                                "file_open_error"   : {"num": 0, "alt": "opening file errors"},
                                "exec_error"      : {"num": 0, "alt": "execution errors"},
                                }

        self.other_counts    = { "unknown"       : {"num": 0, "alt": "unknown result"},
                                 "running"       : {"num": 0, "alt": "running"},
                                 "never_ran"     : 0,
                                 }

        # For the special case of producing an all list
        self.all_run_dirs = []

        self.type_output_files = {}

        self.count_dicts = (self.completed_counts, self.cfg_type_counts, self.error_counts, self.other_counts)

        self.verbose = verbose

    def __del__(self):
        for output in list(self.type_output_files.values()):
            output["object"].close()

    def register_output(self, type_name, output):
        if not type_name in self.type_names(unique=False):
            raise ValueError("Type name %s not valid" % (type_name))

        if type(output) is StringType:
            found_existing_file = False
            for out_file in list(self.type_output_files.values()):
                if out_file["filename"] == output:
                    self.type_output_files[type_name] = out_file
                    found_existing_file = True
                    break

            if not found_existing_file:
                self.type_output_files[type_name] = {"filename": output, "object": open(output, "w")}
        else:
            self.type_output_files[type_name] = {"filename": None, "object": output}

    def type_names(self, unique=True):
        names = []
        for count_dict in self.count_dicts:
            for count_name, count_desc in list(count_dict.items()):
                if unique == False or type(count_desc) is not DictType or "add_to_total" not in count_desc or count_desc["add_to_total"] == True:                    
                    names.append(count_name)
        if unique == False:
            names.append("all")
        return names

    def add_count(self, result_obj, count_duplicates=False):

        for result_type in result_obj.result_types:
            was_counted = False
            for count_dict in self.count_dicts:
                if result_type in count_dict:
                    if type(count_dict[result_type]) is DictType:
                        count_dict[result_type]["num"] += 1
                    else:
                        count_dict[result_type] += 1
                    was_counted = True
            
            if not was_counted:
                raise ValueError('Unknown result type "%s" for run id "%s"' % (result_type, result_obj.run_id))

            if result_type in self.type_output_files:
                print(result_obj.run_id, file=self.type_output_files[result_type]["object"])

        if result_obj.run_id in self.all_run_dirs and self.verbose:
            print("Duplicate sounding id: %s" % result_obj.run_id)                

        if count_duplicates or (not result_obj.run_id in self.all_run_dirs):
            self.all_run_dirs.append(result_obj.run_id)
            if "all" in self.type_output_files:
                print(result_obj.run_id, file=self.type_output_files["all"]["object"])
                    
    def print_overall_stats(self, out_obj=sys.stdout):
        stats_format = "%4s %s\n"

        total_count = 0

        count_sect_strs = []
        for count_dict in self.count_dicts:
            count_items = list(count_dict.items())
            count_items.reverse()

            def count_compare(x, y):
                (x_id, x_desc) = x
                (y_id, y_desc) = y

                if type(x_desc) is DictType and "order" in x_desc:
                    x_order = x_desc["order"]
                else:
                    x_order = -1

                if type(y_desc) is DictType and "order" in y_desc:
                    y_order = y_desc["order"]
                else:
                    y_order = -1

                if x_order >= 0 and y_order >= 0:
                    return cmp(x_order, y_order)
                else:
                    return cmp(x_id, y_id)
                
            if int(sys.version[0]) == 3:
                count_items.sort(key=functools.cmp_to_key(count_compare))
            else:
                count_items.sort(count_compare)

            sect_string = ""
            for (count_id, count_desc) in count_items:
                if type(count_desc) is DictType:
                    count_num = count_desc["num"]
                else:
                    count_num = count_desc

                if count_num > 0:
                    if type(count_desc) is DictType and "add_to_total" in count_desc:
                        if count_desc["add_to_total"]:
                            total_count += count_num
                    else:
                        total_count += count_num
                
                    if type(count_desc) is DictType and "alt" in count_desc:
                        count_name = count_desc["alt"]
                    else:
                        count_name = count_id.replace("_", " ")

                    if type(count_desc) is DictType and "prefix" in count_desc:
                        count_name = "%s %s" % (count_desc["prefix"], count_name)

                    if type(count_desc) is DictType and "indent" in count_desc:
                        indent = " " * count_desc["indent"]
                    else:
                        indent = ""

                    sect_string += indent
                    sect_string += stats_format % (count_num, count_name)

            if len(sect_string) > 0:
                count_sect_strs.append(sect_string)

        if total_count != len(self.all_run_dirs):
            raise Exception("Length of all run ids processed: %d does not match count from all count types: %d. Possibly because counting of duplicates was disabled." % (len(self.all_run_dirs), total_count))
        print("----\n".join(count_sect_strs), end=' ', file=out_obj)
        print("====", file=out_obj)
        print(stats_format % (total_count, "total runs"), file=out_obj)

class run_id_result(object):
    
    def __init__(self, base_path, run_id, aggregate_file=None, verbose=False):
        self.result_types  = []
        self.iterations    = 0

        self.base_path = base_path
        self.run_id = run_id
        self.aggregate_file = aggregate_file
        self.verbose = verbose

        class_members = inspect.getmembers(self, inspect.ismethod)
        class_members.sort()

        if self.verbose:
            print("run_id: ", run_id)

        for member in class_members:
            if member[0].find("check_") == 0:
                check_result = member[1]()
                if check_result != None and check_result != False:
                    if self.verbose:
                        print("Result types: %s" % (", ".join(self.result_types)))
                    break
                
    def find_check_file(self, file_glob):
        file_search = (self.base_path + "/" + file_glob).format(base_path=self.base_path, run_id=self.run_id)
        file_result = glob.glob( file_search )
        if(len(file_result) == 0):
            file_search2 = (self.base_path + "/*/" + file_glob).format(base_path=self.base_path, run_id=self.run_id)
            file_result = glob.glob( file_search2 )

        if self.verbose: print("Checking for file %s" % file_search)
        if len(file_result) > 0 and file_result[0] != "":
            if self.verbose: print("Found file %s" % file_result[0])
            return file_result[0]
        else:
            return None

    def check_hdf_results(self, hdf_obj, run_index):
        # Try and detect brdf type
        try:
            sv_names = hdf_obj[STATEVECTOR_NAMES_DATASET][run_index, :]
                
            for cfg_type_name, cfg_type_srch in SV_CONFIG_TYPE_STRINGS.items():
                fnd_of_type = [sv_elem for sv_elem in sv_names if str(sv_elem).find(cfg_type_srch) >= 0]
                if len(fnd_of_type) > 0:
                    self.result_types.append(cfg_type_name)
                    break
        except KeyError:
            pass

        try:
            outcome = hdf_obj[OUTCOME_FLAG_DATASET][run_index]
        except KeyError:
            outcome = None

        # Check outcome and use master qual value if available.
        # Do not use master qual for max iter/div types or else
        # total count would come out wrong
        if outcome == 1:
            self.result_types.append( "converged" )
            return True
        elif outcome == 2:
            self.result_types.append( "converged" )
            return True
        elif outcome == 3:
            self.result_types.append( "max_iter" )
            return True
        elif outcome == 4:
            self.result_types.append( "max_div" )
            return True
        else:
            return False

    def check_01_aggregate(self):
        if self.verbose: print("Checking aggregate file")

        if self.aggregate_file != None:
            try:
                index = self.aggregate_file.get_sounding_indexes(self.run_id)
            except KeyError:
                return False

            return self.check_hdf_results(self.aggregate_file, index)
        else:
            return False

    def check_02_if_running(self):
        if self.verbose: print("Checking if running")
        running_find = self.find_check_file(LOG_RUNNING_FILENAME_SRCH)

        if running_find != None:
            # If found HDF file with .error on end then some exception
            # has handled correctly
            self.result_types.append( "running" )
            return True

        return False

    def check_03_if_ran(self):
        if self.verbose: print("Checking if ran")
        self.hdf_find = self.find_check_file(HDF_OUTPUT_FILENAME_SRCH)
        self.log_find = self.find_check_file(LOG_OUTPUT_FILENAME_SRCH)

        # No log file or hdf file means we were passed
        # a run id which was never processed
        if self.hdf_find == None and self.log_find == None:
            self.result_types.append( "never_ran" )
            return True

        return False

    def check_04_hdf_outcome(self):
        if self.verbose: print("Checking hdf outcome")

        if self.hdf_find != None:
            with h5py.File(self.hdf_find, "r") as hdf_obj:
                return self.check_hdf_results(hdf_obj, 0)

        # Could not determine anything 
        return False

    def check_05_error_check(self):
        if self.verbose: print("Checking for an error")
        err_find = self.find_check_file(HDF_ERROR_FILENAME_SRCH)

        if err_find != None:
            # If found HDF file with .error on end then some exception
            # has handled correctly
            self.result_types.append( "handled_error" )
            return True

        # Now check the log, maybe exception was caught, but
        # output files not named appropriately
        if self.log_find != None:
            if self.verbose: print("Searching %s for errors" % self.log_find)
            # Search backwards from end of file (for speed)
            # to find exception handling message
            with open(self.log_find) as log_fo:
                line_count = 0
                for log_line in BackwardsReader(log_fo):
                    if line_count > LOG_EXCEPTION_SRCH_LIMIT:
                        break
                    line_count += 1

                    for error_str in FILE_EXCEPTION_STRS:
                        if re.search(error_str, log_line):
                            if self.verbose: "Found error message:", log_line
                            self.result_types.append( "file_open_error" )
                            return True
                    for error_str in HANDLED_EXCEPTION_STRS:
                        if re.search(error_str, log_line):
                            if self.verbose: "Found error message:", log_line
                            self.result_types.append( "handled_error" )
                            return True

        # A log file w/ no HDF file means it started running
        # and where no exception handling message is printed
        # but must have seg faulted without creating hdf file
        #
        # And yes this is same condition check as above, if
        # log file above is checked and nothing found we check here too
        if self.log_find != None and self.hdf_find == None:
            self.result_types.append( "exec_error" )
            return True

        return False

    def check_99_unknown_fall_through(self):
        if self.verbose: print("No other sufficient result recorded. Marking run unknown")
        self.result_types.append( "unknown" )
        return True


def find_run_ids(base_dir, run_dir_ids, search_depth=DEFAULT_MAX_SRCH_DEPTH, filter=None, verbose=False):
    
    for root, dirs, files in os.walk(base_dir, topdown=True):
        if len(root.split("/")) > len(base_dir.split("/")) + DEFAULT_MAX_SRCH_DEPTH:
            continue
        elif os.path.basename(root) in SRCH_IGNORE_DIRS:
            while len(dirs) > 0: dirs.pop(0)
            continue

        for curr_file in files:
            if curr_file == RUN_ID_SEARCH:
                id_list = open(os.path.join(root, curr_file)).read().split()
                for id_name in id_list:
                    run_dir_ids.append( (root, id_name) )

def read_run_id_file(run_id_file):
    # Use L2_Input to parse XML / heritage format file
    # as well as a regular ascii file
    parsed_file = L2InputFile(run_id_file)
    section_names = parsed_file.get_all_section_names()

    # If there are section names then this is an XML or heritage format
    # used to drive the populator
    if len(section_names) > 0:
        # Try to find section where run ids live
        list_sections = [sec_name for sec_name in section_names if re.search(CONFIG_RUN_ID_SEC_RE, sec_name)]

        id_lines = []
        for sec_name in list_sections:
            found_sec_objs = parsed_file.get_section(sec_name)
            if len(found_sec_objs) > 0:
                id_lines += found_sec_objs[0].get_matrix_data()
    else:
        # Otherwise this is just a plain ascii file we
        # have parsed into a list of lists
        id_lines = parsed_file.get_matrix_data()

    # Pull run ids out of list of lists parsed above
    run_id_list = []
    for id_line in id_lines:
        # Parse out first part of space seperated line string or first item from list
        if isinstance(id_line, str):
            run_id = id_line.split()[0].strip()
        elif hasattr(id_line, "__iter__"):
            run_id = id_line[0]

        run_id_list.append(run_id)

    return run_id_list

def standalone_main():
    # Set up command line arguments
    parser = OptionParser(usage="usage: %prog [options] [run_dir_base] [run_dir_base]...")

    parser.add_option( "-t", "--list_all_types", dest="list_all_types",
                       action="store_true",
                       help="list possible types for list_type option",
                       default=False)                    

    parser.add_option( "-d", "--type_names", dest="type_names",
                       metavar="TYPE_NAME",
                       type="string",
                       action="append",
                       help="outputs run ids of a specific type",
                       )
    
    parser.add_option( "-o", "--output_file", dest="output_file",
                       metavar="FILE",
                       type="string",
                       action="append",
                       help="filename for each associated type_names argument",
                       )

    parser.add_option( "-s", "--stats_file", dest="stats_file",
                       metavar="FILE",
                       type="string",
                       help="filename for overall stats instead of STDOUT",
                       )

    parser.add_option( "-v", "--verbose", dest="verbose",
                       action="store_true",
                       help="output information for each test case seen",
                       default=False)                    

    parser.add_option( "-f", "--filter", dest="filter",
                       metavar="REGEXP",
                       action="append",
                       help="filter out tests that match the specified regular expression")

    parser.add_option( "-r", "--run_id_file", dest="run_id_file",
                       metavar="FILE",
                       help="file to read with list of run ids instead of relying on file searching, assumed ids located under current directory")

    parser.add_option( "-m", "--max_search_depth", dest="max_search_depth",
                       metavar="NUM",
                       type="int",
                       help="maximum directory depth to search for output files",
                       default=DEFAULT_MAX_SRCH_DEPTH
                       )

    parser.add_option( "--dups", "--allow_duplicates", dest="allow_duplicates",
                       action="store_true",
                       help="Do not allow duplicate sounding ids to be counted",
                       default=False)                    


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if (len(args) > 0):
        search_base_dirs = [ curr_dir.rstrip("/") for curr_dir in args ]
    elif options.run_id_file == None:
        search_base_dirs = [ "." ]
    else:
        # Ignore search dirs if run id file specified
        search_base_dirs = []

    r_count = results_count(verbose=options.verbose)
    
    if options.list_all_types:
        for type_name in r_count.type_names(unique=False):
            print(type_name)
        sys.exit(0)

    if options.output_file != None and len(options.output_file) > 0:
        last_filename = None
        for type_name in options.type_names:
            if len(options.output_file) > 0:
                out_filename = options.output_file.pop(0)
                r_count.register_output( type_name, out_filename )
                last_filename = out_filename
            else:
                r_count.register_output( type_name, last_filename )

    elif options.type_names != None and len(options.type_names) > 0:
        for type_name in options.type_names:
            r_count.register_output(type_name, sys.stdout)

    run_dir_ids = []
    # Add any run ids specified in a file
    # Handle parsing either the SDOS XML config or an ascii
    # file where the run id is the first part of the line
    if options.run_id_file != None:
        if not os.path.exists(options.run_id_file):
            parser.error("Run id file '%s' does not exist" % options.run_id_file)
        if options.verbose:
            print("Reading run ids from %s" % options.run_id_file)
        curr_dir = os.curdir
        for run_id in read_run_id_file(options.run_id_file):
            run_dir_ids.append( (curr_dir, run_id) )

    # Append any directories found by searching from command line specified base paths

    for curr_base_dir in search_base_dirs:
        if options.verbose:
            print("Finding run id starting at %s with depth %d" % (curr_base_dir, options.max_search_depth))

        find_run_ids(curr_base_dir, run_dir_ids, options.max_search_depth, options.filter, options.verbose)

    # Go through all found run ids and find the result of the run
    count = 1
    for id_dir, run_id in run_dir_ids:
        if options.verbose:
            print('Processing run id #%04d %s located under directory: "%s"' % (count, run_id, id_dir))
            count += 1

        d_result = run_id_result(id_dir, run_id, verbose=options.verbose)
        r_count.add_count( d_result, count_duplicates=options.allow_duplicates )

    if options.verbose:
        print("")

    if options.type_names == None or len(options.type_names) <= 0 or options.stats_file != None:
        stats_out_obj = sys.stdout
        if options.stats_file != None:
            stats_out_obj = open(options.stats_file, "w")
            
        r_count.print_overall_stats(stats_out_obj)

        if options.stats_file != None:
            stats_out_obj.close()

    if options.output_file != None and len(options.output_file) > 0:
        for out_obj in type_dir_out_objs:
            out_obj.close()
    
if __name__ == "__main__":
    standalone_main()
