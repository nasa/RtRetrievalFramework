#!/usr/bin/env python

import os
from optparse import OptionParser
from StringIO import StringIO

import run_results

status_count_filename = 'run_results.txt'
aggregate_types = { 'exceeded_max': [ 'max_iter', 'max_div' ],
                    'failed' : ['handled_error', 'exec_error', 'other', 'unknown_quality']
                    }

def gather_status(base_dir, output_dir, run_id_file=None, verbose=False):

    # Strip trailing slash
    output_dir = output_dir.rstrip('/')

    # Make output location if not already existing
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

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
    if(run_id_file != None):
        for run_id in run_results.read_run_id_file(run_id_file):
            run_dir_ids.append( (base_dir, run_id) )
    else:
        run_results.find_run_ids(base_dir, run_id_list)

    for id_dir, run_id in run_id_list:
        d_result = run_results.run_id_result(id_dir, run_id)
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
    for key in aggregate_types.keys(): aggregated_dirs[key] = ''
    
    for type_name, filename in out_filenames.items():
        out_str_obj = out_objects[type_name]
        out_strings = out_str_obj.getvalue()

        if len(out_strings) > 0:
            run_lists[type_name] = out_strings.split('\n')[0:-1] # Loose last empty one due to \n
            
            out_file_obj = open(filename, 'w')
            out_file_obj.write(out_strings)
            out_file_obj.close()

            for agg_name, agg_types in aggregate_types.items():
                if type_name in agg_types:
                    aggregated_dirs[agg_name] += out_strings
        elif os.path.exists(filename):
            os.remove(filename)
            
    # Output non empty aggregated files
    for agg_name, agg_strings in aggregated_dirs.items():
        agg_filename = '%s/%s.list' % (output_dir, agg_name)

        if len(agg_strings) > 0:
            run_lists[agg_name] = agg_strings.split('\n')[0:-1] # Loose last empty one due to \n
            
            out_file_obj = open(agg_filename, 'w')
            out_file_obj.write(agg_strings)
            out_file_obj.close()
        elif os.path.exists(agg_filename):
            os.remove(agg_filename)
            
    return run_lists


def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] [run_dir_base]")

    parser.add_option( "-o", "--output_dir", dest="output_dir",
                       metavar="DIR",
                       type="string",
                       default='./',
                       help="directory for outputting status files",
                       )

    parser.add_option( "-r", "--run_id_file", dest="run_id_file",
                       metavar="FILE",
                       help="file to read with list of run ids instead of relying on file searching, assumed ids located under current directory")

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) > 1:
        raise ValueError('Only one base directory should be specified')
    elif len(args) == 1:
        base_dir = args[0]
    else:
        base_dir = './'

    gather_status(base_dir, options.output_dir, run_id_file=options.run_id_file, verbose=True)

if __name__ == "__main__":
    standalone_main()
