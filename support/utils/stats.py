#!/usr/bin/env python

from future import standard_library
standard_library.install_aliases()

import os
from glob import glob
from optparse import OptionParser
from io import StringIO

from full_physics import acos_file, docopt_simple

import run_results

version = "July 18, 2017"
usage = '''Usage:
  stats.py [options] [<run_dir_base>]
  stats.py -h | --help
  stats.py -v | --version

Return stats on the jobs run in the given directory (default is the current
director). 

Options:
  -h --help
     Print this message

  -o --output-dir=s
     Directory for outputting status files. [default: ./]

  -r --run-id-file=s
     File to read with the list of run ids instead of relying on file 
     searching, assumed ids located under current directory.

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)

SOUNDING_ID_FILENAME_GLOBS = ['rmgr_sounding_ids.list', '*.config']

status_count_filename = 'run_results.txt'
aggregate_types = { 'exceeded_max': [ 'max_iter', 'max_div' ],
                    'failed' : ['handled_error', 'file_open_error', 'exec_error', 'other', 'unknown_quality']
                    }

def gather_status(base_dir, output_dir, run_id_file=None,
                  aggregate_file=None, verbose=False):

    # Strip trailing slash
    output_dir = output_dir.rstrip('/')

    # Make output location if not already existing
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load aggregate file
    if aggregate_file != None:
        aggregate_file = acos_file.L2(aggregate_file)

    r_count = run_results.results_count()

    out_filenames = {}
    out_objects = {}
    for type_name in r_count.type_names(unique=False):
        filename = '%s/%s.list' % (output_dir, type_name)
        out_obj  = StringIO()
        out_filenames[type_name] = filename
        out_objects[type_name] = out_obj
        r_count.register_output(type_name, out_obj)

    run_id_list = []
    if (run_id_file != None):
        for run_id in run_results.read_run_id_file(run_id_file):
            run_id_list.append( (base_dir, run_id) )
    else:
        run_results.find_run_ids(base_dir, run_id_list)

    for id_dir, run_id in run_id_list:
        d_result = run_results.run_id_result(id_dir, run_id, aggregate_file=aggregate_file)
        r_count.add_count(d_result)

    # Output stats to file
    stats_out_obj = open('%s/%s' % (output_dir, status_count_filename), 'w')
    r_count.print_overall_stats(stats_out_obj)
    stats_out_obj.close()

    # Output stats to screen
    if verbose:
        r_count.print_overall_stats()

    # Output non-empty files
    run_lists = {}
    aggregated_dirs = {}
    for key in list(aggregate_types.keys()): aggregated_dirs[key] = ''
    
    for type_name, filename in list(out_filenames.items()):
        out_str_obj = out_objects[type_name]
        out_strings = out_str_obj.getvalue()

        if len(out_strings) > 0:
            run_lists[type_name] = out_strings.split('\n')[0:-1] # Loose last empty one due to \n
            
            out_file_obj = open(filename, 'w')
            out_file_obj.write(out_strings)
            out_file_obj.close()

            for agg_name, agg_types in list(aggregate_types.items()):
                if type_name in agg_types:
                    aggregated_dirs[agg_name] += out_strings
        elif os.path.exists(filename):
            os.remove(filename)
            
    # Output non empty aggregated files
    for agg_name, agg_strings in list(aggregated_dirs.items()):
        agg_filename = '%s/%s.list' % (output_dir, agg_name)

        if len(agg_strings) > 0:
            run_lists[agg_name] = agg_strings.split('\n')[0:-1] # Loose last empty one due to \n
            
            out_file_obj = open(agg_filename, 'w')
            out_file_obj.write(agg_strings)
            out_file_obj.close()
        elif os.path.exists(agg_filename):
            os.remove(agg_filename)
            
    return run_lists


if(args.run_dir_base):
    base_dir = args.run_dir_base
else:
    base_dir = "./"

# Try and find a sounding id file to use
run_id_file = args.run_id_file
if run_id_file is None:
    for srch_glob in SOUNDING_ID_FILENAME_GLOBS:
        res = glob(os.path.join(base_dir, srch_glob))
        if len(res) > 0:
            run_id_file = res[0]
            break

aggregate_srch = glob(os.path.join(base_dir, "l2_aggregate.h5"))
if len(aggregate_srch) > 0:
    aggregate_file = aggregate_srch[0]
else:
    aggregate_file = None

gather_status(base_dir, args.output_dir, run_id_file=run_id_file,
              aggregate_file=aggregate_file, verbose=True)

