#!/usr/bin/env python

import os
import sys
import re
from time import gmtime, strftime
from optparse import OptionParser
from sets import Set
import glob
from types import ListType

from OCO_TextUtils import *
from OCO_Matrix import OCO_Matrix
import numpy

fill_value = -999.0
default_out_filename = 'aggregate_data.dat'

##
# Return format used for naming columns that have a digit component

def get_column_format(num_cols):
    if num_cols >= 1000:
        return "%s_%04d"
    elif num_cols >= 100:
        return "%s_%03d"
    elif num_cols >= 10:
        return "%s_%02d"
    else:
        return "%s_%0d"

##
# Tries to load a data object using OCO_Matrix format. Failing that the file
# is read a tabled data where each space indicates a new column

def get_data_object(data_filename):

    # Try to load data using OCO_Matrix class
    try:
        data_obj = OCO_Matrix(data_filename)
        return data_obj
    except:
        pass

    # Now load file as tabled data
    table_file_obj = open(data_filename, 'r')
    file_lines = table_file_obj.readlines()
    table_file_obj.close()

    # Seperate each line by spaces. Keep count of maximum
    # number of columns seen for when file is added so we can
    # know how to size the resultng matrix
    max_cols = 0
    file_rows = []
    for line in file_lines:
        if line.find('#') < 0 and len(line.strip()) != 0:
            line_cols = line.strip().split()
            file_rows.append(line_cols)
            max_cols = max(max_cols, len(line_cols))

#    data_mat = numpy.zeros((len(file_rows), max_cols), dtype=float)
    data_mat = numpy.zeros((len(file_rows), max_cols), dtype=numpy.chararray)

    for row_idx in range(len(file_rows)):
        num_cols = len(file_rows[row_idx])
        for col_idx in range(num_cols):
            col_value = file_rows[row_idx][col_idx]
            data_mat[row_idx][col_idx] = col_value
#            try:
#                data_mat[row_idx][col_idx] = float(col_value)
#            except:
#                data_mat[row_idx][col_idx] = fill_value

    # Create label names based on filename and index or else can
    # not select specific columns
    label_base = os.path.basename(data_filename)
    label_base = label_base[0:label_base.rfind('.')] # Remove extension

    data_labels = []    
    for col_idx in range(max_cols):
        data_labels.append( get_column_format(max_cols) % (label_base, col_idx) )
    
    # Save data into OCO Matrix object
    data_obj = OCO_Matrix()
    data_obj.dims = [len(file_rows), max_cols]
    data_obj.labels = data_labels
    data_obj.data = data_mat
    
    return data_obj
    
##
# main

def standalone_main(arg_list):
    parser = OptionParser(usage="usage: %prog [options] [run_dir] [run_dir]...")
    parser.add_option( "-d", "--data_file", dest="data_files",
                       metavar="FILE",
                       action="append",
                       help="relative filename to extract data from",
                       )

    parser.add_option( "-c", "--column", dest="columns",
                       metavar="NAME",
                       action="append",
                       help="name of a column to extract from files")

    parser.add_option( "-r", "--run_dir_file", dest="run_dir_file",
                       metavar="FILE",
                       help="file to read list of run directories from")

    parser.add_option( "-o", "--output_data_file", dest="output_data_file",
                       metavar="FILE", default=default_out_filename,
                       help="file where output table data will be saved, default: %s" % default_out_filename)

    parser.add_option( "--multi_row", dest="multi_row",
                       default=False,
                       action="store_true",
                       help="Allow multiple rows per run directory for a certain column"
                       )

    parser.add_option( "--no_case_name", dest="add_case_name",
                       default=True,
                       action="store_false",
                       help="do not add case name column"
                       )


    ##
    ## Parse command line arguments
    (options, args) = parser.parse_args(args=arg_list)

    # Make sure at least one file type is specified
    if options.data_files == None or len(options.data_files) <= 0:
        parser.error("At least one data file must be specified by the --data_file option")

    # Make sure we have at least one column
    if options.columns == None or len(options.columns) == 0:
        parser.error("At least one column must be defined")

    ##
    ## Gather list of run directories to gather information from either
    ## from non option arguments or from specified run dir file
    run_dir_list = []

    if options.run_dir_file != None:
        if not os.path.exists(options.run_dir_file):
            parse.error("Run directory file '%s' does not exist" % options.run_dir_file)

        run_dir_fh = open(options.run_dir_file, 'r')
        for file_dir in run_dir_fh.readlines():
            run_dir_list.append(file_dir.strip())
        run_dir_fh.close()

    if (len(args) > 0):
        for arg_dir in args:
            run_dir_list.append(arg_dir)
    else:
        if len(run_dir_list) == 0:
            run_dir_list.append('./')

    ##
    ## Make sure there are some run files available
    if run_dir_list == None or len(run_dir_list) <= 0:
        parser.error("No run directories have been specified")

    runs_data = {}
    column_names = {}
    total_num_cols = 0
    for run_dir in run_dir_list:   
        run_name = os.path.dirname(run_dir)

        for data_name in Set(options.data_files):
            if data_name.find('/') == 0:
                data_file_search = glob.glob("%s" % (data_name))
            else:
                data_file_search = glob.glob("%s/%s" % (run_dir, data_name))

            # Skip this data file for this run dir
            if len(data_file_search) == 0:
                continue

            data_full_path = data_file_search[0]

            # Parse file and load into memory
            print 'Opening: %s' % data_full_path
            data_obj = get_data_object(data_full_path)

            # If there was no way to parse the file then skip this data file
            if data_obj == None:
                continue        

            # Get list of columns to take data from
            used_columns = []
            row_specifiers = []
            renamed_columns = []
            if options.columns == None or len(options.columns) <= 0:
                used_columns = data_obj.labels
            else:
                for col_option in options.columns:
                    col_parts = col_option.split('#')
                    col_name = col_parts[0]

                    if len(col_parts) >= 2:
                        row_spec = col_parts[1]
                        if row_spec.isdigit():
                            row_spec = '%d' % int(row_spec)
                    else:
                        row_spec = ':'

                    if len(col_parts) >= 3:
                        col_new_name = col_parts[2]
                    else:
                        col_new_name = col_name

                    if col_name.isdigit():
                        # If column name is a integer then try and look it up
                        # in the labels, failing that use the index if it is not
                        # larger than the number of columns
                        if int(col_name) >= 0 and int(col_name) < len(data_obj.labels):
                            used_columns.append(data_obj.labels[int(col_name)])
                            if col_new_name == col_name:
                                col_new_name = data_obj.labels[int(col_name)]
                        elif int(col_name) < data_obj.dims[1]:
                            used_columns.append(col_name)
                        row_specifiers.append(row_spec)

                        renamed_columns.append(col_new_name)
                    elif col_name in data_obj.labels:
                        # Use the column name as is since it appears in the file's label list
                        used_columns.append(col_name)
                        row_specifiers.append(row_spec)

                        renamed_columns.append(col_new_name)

            # Get data for each used column
            for (col_orig_name, row_spec, col_new_name) in zip(used_columns, row_specifiers, renamed_columns):
                # Find the index for the column so we know how to extract it
                if col_orig_name.isdigit():
                    col_index = int(col_orig_name)
                else:
                    col_index = data_obj.labels.index(col_orig_name)

                col_data = data_obj.data[:, col_index]

                # Try the row_spec as a range for an array otherwise use as a filter
                all_data_range = range(0, data_obj.dims[0])
                try:
                    used_data_range = eval('all_data_range[' + row_spec + ']')
                except:
                    used_data_range = []
                    for row_index in all_data_range:
                        for col_index in range(data_obj.dims[1]):
                            if re.search(row_spec, str(data_obj.data[row_index, col_index])):
                                used_data_range.append(row_index)
                                break

                if not type(used_data_range) is ListType:
                    used_data_range = [ used_data_range ]

                data_count = 1
                for data_idx in used_data_range:
                    # If this data has more than one row, then add same column name
                    # multiple times with an index appended to the name for each row
                    # so long as multi_row is not enabled
                    if len(used_data_range) > 1 and not options.multi_row:
                        data_col_name = get_column_format(data_obj.dims[0]) % (col_new_name, data_count)
                    else:
                        data_col_name = col_new_name

                    # If multi_row is enabled then add new rows for multiple items in a file
                    if len(used_data_range) > 1 and options.multi_row:
                        data_run_name = get_column_format(data_obj.dims[0]) % (run_name, data_count)
                    elif len(used_data_range) == 1 and options.multi_row:
                        if column_names.has_key(data_col_name) and data_name in column_names[data_col_name]:
                            row_count = len(column_names[data_col_name]) + 1
                        else:
                            row_count = 0
                        data_run_name = get_column_format(data_obj.dims[1]) % (run_name, row_count)
                    else:
                        data_run_name = run_name

                    # Store name of file associated with column in case we need
                    # to go back and rename a column to append its filename
                    # if we encountered the column name twice
                    if not column_names.has_key(data_col_name):
                        column_names[data_col_name] = []
                    if not data_name in column_names[data_col_name]:
                        column_names[data_col_name].append(data_name)
                        total_num_cols = total_num_cols + 1

                    # Store column data indexed by run dir and column name
                    if not runs_data.has_key(data_run_name):
                        runs_data[data_run_name] = {}
                    if not runs_data[data_run_name].has_key(data_col_name):
                        runs_data[data_run_name][data_col_name] = {}

                    runs_data[data_run_name][data_col_name][data_name] = col_data[data_idx]

                    data_count += 1

    ##
    # Now that data from all the run dirs has been collected we can loop over the collected
    # data and place it into a new matrix
    print '%d directories found' % len(run_dir_list)
    print '%s rows used' % len(runs_data)
    print '%s columns found' % total_num_cols

    if options.add_case_name:
        total_num_cols = total_num_cols+1

    output_data_matrix = numpy.zeros((len(runs_data), total_num_cols), dtype=numpy.chararray)

    # Get column name list that will be used for looping over data and placing into file
    # Sort based upon the ordet the columns were input from the command line
    supplied_col_order = []
    for col_opt in options.columns:
        col_parts = col_opt.split('#')

        if len(col_parts) >= 3:
            col_name = col_parts[2]
        else:
            col_name = col_parts[0]

        supplied_col_order.append(col_name)

    def col_compare(x, y):
        if supplied_col_order == None:
            return cmp(x, y)
        if x in supplied_col_order and not y in supplied_col_order:
            return -1
        elif not x in supplied_col_order and y in supplied_col_order:
            return 1
        elif not x in supplied_col_order and y not in supplied_col_order:
            return cmp(x, y)
        else:
            return cmp(supplied_col_order.index(x), supplied_col_order.index(y))

    col_name_list = column_names.keys()
    col_name_list.sort(col_compare)

    used_run_names = runs_data.keys()
    used_run_names.sort()
    testcase_names = extract_run_names(used_run_names)

    row_index = 0
    for (run_name, tc_name) in zip(used_run_names, testcase_names):   

        if options.add_case_name:
            output_data_matrix[row_index, 0] = tc_name
            col_index = 1
        else:
            col_index = 0

        for col_name in col_name_list:
            data_name_list = column_names[col_name]

            for data_name in data_name_list:
                if runs_data.has_key(run_name) and runs_data[run_name].has_key(col_name) and \
                   runs_data[run_name][col_name].has_key(data_name):
                    output_data_matrix[row_index, col_index] = runs_data[run_name][col_name][data_name]
                else:
                    output_data_matrix[row_index, col_index] = fill_value

                col_index = col_index + 1

        row_index = row_index + 1

    ##
    # Flatten out column names for output file label
    if options.add_case_name:
        flat_column_names = ['case_name']
    else:
        flat_column_names = []

    for col_name in col_name_list:
        data_name_list = column_names[col_name]

        if len(data_name_list) > 1:
            for data_name in data_name_list:
                flat_column_names.append( col_name + "-" + os.path.basename(data_name) )
        else:
            flat_column_names.append(col_name)

    print 'Writing: %s' % options.output_data_file

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = 'Data aggregated from: %s' % ', '.join(options.data_files)
    out_mat_obj.dims = [len(runs_data), total_num_cols]
    out_mat_obj.labels = flat_column_names
    out_mat_obj.data = output_data_matrix

    out_mat_obj.write(options.output_data_file)

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):
    if type(scriptOptions) is not ListType:
        scriptOptions = scriptOptions.split()
    
    scriptOptions.append('-d')
    scriptOptions.append(fileObj.filename)
    scriptOptions.append('--no_case_name')
    standalone_main(scriptOptions)

    print '%s: %s' % ('gather_data', ' '.join(scriptOptions))

if __name__ == "__main__":
    standalone_main(sys.argv[1:])
