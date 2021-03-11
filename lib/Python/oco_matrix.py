#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
import bisect
import datetime
import inspect
import math
import sys
import re
import six
my_list_type = None
try:
    # This is for python 2
    from types import ListType
    my_list_type = ListType
except ImportError:
    # This is for python 3
    my_list_type = list

import numpy

import copy as copy_module

from .text_util import TeX_string

# IEEE 754 double-precision
EPS = 1.11022302463e-16

FILL_VALUE = 'N/A'

class OcoMatrix(object):
    """Read a matrix file output by the L2 code
    """

    def __init__(self, file_input=None, **read_args):

        # Encapsulate old attributes as properies
        self._dims = numpy.zeros(2, dtype=int)

        self.labels = []
        self.units = []

        self.header = {}
        self.pixels = []

        self.file_id = []
        self.file_creation = ""

        self.verbose_pre = "OcoMatrix:"

        self.filename = None
        if file_input != None:
            self.read(file_input, **read_args)

    @property
    def dims(self):
        if hasattr(self, 'data') and self.data != None and hasattr(self, 'shape'):
            return self.data.shape
        else:
            return self._dims

    @dims.setter
    def dims(self, new_value):
        self._dims[:] = new_value[:]

    def pixel_ranges(self):
        beg_idxs = self.pixels
        end_idxs = self.pixels[1:] + [self.dims[0]]

        return list(zip(beg_idxs, end_idxs))

    @property
    def labels_lower(self):
        return [ curr_lbl.lower() for curr_lbl in self.labels ]

    def get_header_value(self, keyword_name):

        matched_values = []
        for curr_name in list(self.header.keys()):
            if curr_name.lower().strip() == keyword_name.lower().strip():
                matched_values.append(self.header[curr_name])

        if len(matched_values) == 0:
            return None
        elif len(matched_values) == 1:
            return matched_values[0]
        else:
            return matched_values

    def __getitem__(self, column_ident):

        if isinstance(column_ident, six.integer_types) or column_ident.isdigit():
            column_indexes = int(column_ident)
        else:
            column_indexes = self.find_labels(column_ident, match_case=True, indexes=True)
            if len(column_indexes) == 0:
                column_indexes = self.find_labels(column_ident, match_case=False, indexes=True)
            if len(column_indexes) == 0:
                raise LookupError('Could not find column(s) named like "%s" among labels: %s' % (column_ident, self.labels))

        return self.data[:, column_indexes]

    def find_labels(self, label_re, match_case=True, indexes=False):

        if not match_case:
            labels_srch = [ curr_lbl.lower() for curr_lbl in self.labels ]
            label_re = label_re.lower()
        else:
            labels_srch = self.labels

        found_labels = []
        for curr_lbl in labels_srch:
            if re.search(label_re, curr_lbl):
                if indexes:
                    found_labels.append( labels_srch.index(curr_lbl) )
                else:
                    found_labels.append( self.labels[ labels_srch.index(curr_lbl) ] )

        return found_labels

    def as_dict(self, key_label=0, existing_dict=None):

        if existing_dict == None:
            existing_dict = {}

        if isinstance(key_label, six.integer_types) or key_label.isdigit():
            key_index = int(key_label)
        else:
            key_index = self.find_labels(key_label, match_case=False, indexes=True)
            if len(key_index) == 0:
                raise LookupError('Could not find label: %s' % key_label)
            elif len(key_index) > 1: 
                raise LookupError('Found too many matches for label: %s => %s' % (key_label, self.labels[key_index]))
            else:
                key_index = key_index[0]

        dict_labels = []
        for col_idx in range(self.data.shape[1]):
            if col_idx < len(self.labels):
                dict_labels.append( self.labels[col_idx] )
            else:
                dict_labels.append( 'column_%04d' % col_idx )

        for row_idx in range(self.data.shape[0]):
            row_key_value = self.data[row_idx, key_index]

            if row_key_value not in existing_dict:
                existing_dict[row_key_value] = {}
            
            for col_idx in range(self.data.shape[1]):
                existing_dict[row_key_value][dict_labels[col_idx]] = self.data[row_idx, col_idx]

        return existing_dict

    def add_column(self, column_name):
        self.labels.append(column_name)
        self.units.append("")
        new_data = numpy.zeros((self.data.shape[0], self.data.shape[1]+1), dtype=self.data.dtype)
        new_data[:, 0:self.data.shape[1]] = self.data[:,:]
        self.data = new_data
        self._dims[:] = self.data.shape[:]

    def read(self, file_input, read_data=True, ignore_conv_err=False, as_strings=False, strip_nonnumbers=False):

        if isinstance(file_input, six.string_types):
            # Given a filename so try to open it
            self.filename = file_input
            try:
                file_obj = open(file_input, "r")
            except IOError:
                raise IOError("Error opening file: %s" % file_input)
        elif hasattr(file_input, 'readline'):
            # Given a file like object
            file_obj = file_input
        else:
            raise IOError('Unknown object passed for reading: %s' % file_input)
          
        # Read the header
        in_header = False
        header_less = True
        file_lines = []
        while(True):
            this_line = file_obj.readline()
            if (len(this_line) == 0): break

            # Skip commented lines
            hash = this_line.find("#")
            if (hash != -1):
                this_line = this_line[:hash]
            if (len(this_line) == 0): continue

            line = this_line.lower()
            if (line.find("begin header") != -1):
                in_header = True
                header_less = False
                continue

            if (line.find("end header") != -1):
                break
            tokens = line.split()
            if (len(tokens) == 0): continue
            if (tokens[0] == "num_rows"):
                self._dims[0] = int(tokens[2])
            elif (tokens[0] == "num_columns"):
                self._dims[1] = int(tokens[2])
            elif (tokens[0] == "start_pixels"):
                tokens = this_line.split()
                for j in range(2, len(tokens)):
                    self.pixels.append(int(tokens[j]) - 1)
            elif (tokens[0] == "labels"):
                tokens = this_line.split('"')
                for j in range(1, len(tokens)):
                    if (tokens[j].isspace() == False):
                        self.labels.append(tokens[j])
            elif (tokens[0] == "units"):
                tokens = this_line.split('"')
                for j in range(1, len(tokens)):
                    if (tokens[j].isspace() == False):
                        self.units.append(tokens[j])
            elif (tokens[0] == "file_id"):
                self.file_id = this_line.split('=')[1].strip().strip('"')
            elif (tokens[0] == "file_creation"):
                self.file_creation = this_line.split('=')[1].strip().strip('"')
            elif in_header:
                (head_key, head_val) = [item.strip().strip('"').strip() for item in this_line.split('=')]
                self.header[head_key] = head_val
            elif not in_header:
                file_lines.append( this_line )

        # Read any additional data not read above
        while (True):
            this_line = file_obj.readline()
            if this_line == '':
                break
            file_lines.append( this_line )

        # Only close if we were given a filename we opened ourselves
        if isinstance(file_input, six.string_types):
            file_obj.close()

        # Get dimensions from read data if not found in header
        if (self._dims[0] == 0):
            self._dims[0] = len(file_lines)

        if (self._dims[1] == 0):
            for this_line in file_lines:
                line_parts = this_line.strip().split()
                self._dims[1] = max(self._dims[1], len(line_parts))


        if as_strings:
            data_type = numpy.chararray
        else:
            data_type = float

        if header_less:
            self.data = []
        else:
            self.data = numpy.zeros((self._dims[0], self._dims[1]), dtype=data_type, order='F')

        dst_row = 0
        dst_col = 0
        for this_line in file_lines:
            if(self._dims[0] < dst_row):
                raise ValueError("more rows read than specified in file for file: %s" % filename)

            # Skip commented lines
            hash = this_line.find("#")
            if (hash != -1):
                this_line = this_line[:hash]
            if (len(this_line) == 0): continue

            line = this_line.split()

            if header_less:
                self.data.append( [] )

            for src_col in range (0, len(line)):
                try:
                    if as_strings:
                        new_value = line[src_col]
                    else:
                        new_value = float(line[src_col])
                except ValueError:
                    # Fix malformed exponentials that ran out of string space
                    if re.search('\d+[.]\d*[dD]([-+]?\d+)$', line[src_col]):
                        exp_start = re.search('[dD]([-+]?\d+)$', line[src_col]).start()
                        new_value = line[src_col].lower()
                        new_value = new_value.replace('d', '')
                        new_value = float( new_value.replace(new_value[exp_start:], 'E' + new_value[exp_start:]) )
                    elif re.search('[A-Za-z_/.]+', line[src_col]) and strip_nonnumbers:
                        cleaned_val = re.sub('[A-Za-z_/.]+', '', line[src_col])
                        if len(cleaned_val) > 0:
                            new_value = float(cleaned_val)
                        else:
                            new_value = None
                    else:
                        if not ignore_conv_err:
                            print('invalid literal for float(): "%s" at row: %d, column: %d' % (line[src_col], dst_row, src_col), file=sys.stderr)
                        new_value = None

                if header_less:
                    self.data[dst_row].append(new_value)
                else:
                    self.data[dst_row, dst_col] = new_value
                dst_col += 1

            # Only advance destination row if completely filled up columns
            if (dst_col >= self._dims[1] or header_less):
                dst_row += 1
                dst_col = 0

        # Create the y axis labels
        self.ytitle = []
        if (len(self.labels) == 0):
            for i in range(0, self._dims[1]):
                self.ytitle.append("")
        else:
            for i in range(0, len(self.labels)):
                if (i >= len(self.units)):
                    self.ytitle.append("%s" % self.labels[i])
                else:
                    self.ytitle.append("%s (%s)" % (self.labels[i], 
                                                    TeX_string(self.units[i])))
            if (len(self.labels) == 1):
                for i in range(1, self._dims[1]):
                    self.ytitle.append(self.ytitle[0])

    def get_column_format_info(self, column_index):

        if type(self.data) is my_list_type:
            column_data = self.data
        else:
            column_data = self.data[:, column_index]

        is_num = []
        for item in column_data:
            if isinstance(item, six.string_types):
                cast_failed = False
                try:
                    float(item)
                except ValueError:
                    cast_failed = True
                    
                if not cast_failed:
                    is_num.append(True)
                else:
                    is_num.append(False)
            else:
                is_num.append(True)
    
        if is_num.count(True) == len(column_data):
            # Ensure that column data is values and not strings
            float_data = []
            for item in column_data:
                if item == None:
                    float_data.append( 0.0 )
                else:
                    float_data.append( float(item) )

            value_small = [ elem < EPS for elem in float_data ].count(True) == len(float_data)
            
            holds_precision = True
            best_precision = 10
            curr_precision = best_precision
            while not value_small and holds_precision and curr_precision >= 0:
               # Convert column to string w/ current precision then back to
               # a double in order to compare against the original column data and
               # machine precision
               column_compare = [ float(('%0.'+str(curr_precision)+'e') % elem) for elem in float_data ]

               # Make sure all are less than machine noise in differents and
               # mark the current precision as best seen
               compare_diff = [ abs(a-b) for a, b in zip(float_data, column_compare) ]

               diff_small = [ elem < EPS for elem in compare_diff ].count(True) == len(float_data)

               if diff_small:
                   best_precision = curr_precision
               else:
                   holds_precision = False

               curr_precision -= 1

            # Make sure that item is really an int by comparing the conversion
            # of the column data to ints with the real numbers
            if len(float_data) > 0:
                int_column_data    = []
                for elem in float_data:
                    if str(elem) == "nan":
                        int_column_data.append(0)
                    else:
                        int_column_data.append(int(elem))
                        
                int_diff = [ abs(a-b) for a, b in zip(int_column_data, float_data) ]

                col_all_ints = ([ elem < EPS for elem in int_diff ].count(True) == len(float_data))# and (max(int_column_data) < 1.0E6)
            else:
                col_all_ints = []
            
            column_width = 0
            if not value_small and col_all_ints:
                for col_elem in float_data:
                    column_width = max(column_width, len('%d' % col_elem))
                column_width += 1
                prec_type_code = 'd'
            else:
                column_width = curr_precision + 8
                prec_type_code = 'e'
        else:
            lens = [ len(str(item)) for item in column_data ]
            best_precision = max(lens)
            column_width = best_precision
            prec_type_code = 's'
       
        return (column_width, best_precision, prec_type_code)

    def write(self, file_output, use_set_dims=False, auto_size_cols=False, default_precision=8, verbose=False):

        if isinstance(file_output, six.string_types):
            self.filename = file_output
            file_obj = open(file_output, "w")
        elif hasattr(file_output, 'write'):
            file_obj = file_output
        else:
            raise IOError('Unknown object passed for writing: %s' % file_output)

        if not use_set_dims and hasattr(self, 'data'):
            if type(self.data) is my_list_type:
                self._dims = [len(self.data), 1]
            else:
                self._dims[:] = self.data.shape[:]

        if verbose:
            print(self.verbose_pre, 'Dims: ', self._dims)

        header_wrt = copy_module.copy(self.header)
        header_wrt['File_ID']       = '"%s"'   % (self.file_id)
        header_wrt['File_Creation'] = '"%s"' % (datetime.datetime.utcnow().strftime("%FT%T%Z"))

        # If there is no data attribute then this is probably
        # a header only file
        if hasattr(self, 'data'):
            header_wrt['Num_Rows']      = '%d'     % (self._dims[0])
            header_wrt['Num_Columns']   = '%d'     % (self._dims[1])
            
        if (len(self.labels) > 0):
            header_wrt['Labels'] = ' '.join([ '"%s"' % val for val in self.labels ])

        if (len(self.units) > 0):
            header_wrt['Units'] = ' '.join([ '"%s"' % val for val in self.units ])

        if (len(self.pixels) > 0):
            header_wrt['Start_Pixels'] = ' '.join([ '%d' % (pix+1) for pix in self.pixels ])

        header_names = list(header_wrt.keys())
        header_names.sort()

        max_name_len = max([ len(name) for name in header_names ])

        file_obj.write("begin HEADER\n")
        for curr_name in header_names:
            num_space = max_name_len - len(curr_name)
            if hasattr(header_wrt[curr_name], 'shape'):
                value_str = ' '.join([str(curr_val) for curr_val in numpy.ravel(header_wrt[curr_name])])
            else:
                value_str = header_wrt[curr_name]
            file_obj.write("  %s%s = %s\n" % (curr_name, ' ' * num_space, value_str))
        file_obj.write("end HEADER\n\n")

        # Write data portion of file
        if hasattr(self, 'data'):
            # Get column formats that take in consideration element precision
            # as well as label widths
            data_col_formats = []
            data_col_types = []
            lbl_col_formats = []

            if auto_size_cols:

                if verbose:
                    print(self.verbose_pre, 'Auto sizing columns')

                for j in range(self._dims[1]):
                    (col_width, prec, prec_type) = self.get_column_format_info(j)

                    # Find which is bigger in width, the data or the label and
                    # add some white space padding for good measure
                    if (len(self.labels) > j):
                        col_width = max(col_width, len(self.labels[j]) + 1)

                    # Add spacing between items
                    col_width += 2

                    lbl_col_formats.append('%%%ds' % (col_width))

                    # Add one extra space for first data column due to comment mark
                    # of labels. Now data column lines up with data label
                    if j == 0:
                        col_width += 1

                    if prec_type == 'e':
                        data_col_formats.append('%%%d.%de' % (col_width, prec))
                    else:
                        data_col_formats.append('%%%d%s' % (col_width, prec_type))

                    data_col_types.append(prec_type)

            else:
                if verbose:
                    print(self.verbose_pre, 'Using standard column sizes')

                if type(self.data) is my_list_type:
                    data_col_formats   = [ '%s' ]
                    data_col_types     = [ 's' ]
                    lbl_col_formats    = [ '%s' ]
                    lbl_col_formats[0] = '%s'
                else:
                    default_col_len = default_precision + 9
                    col_lens = [ max(default_col_len, len(self.labels[col_count])+1) for col_count in range(len(self.labels)) ]
                    for x in range(self._dims[1]-len(self.labels)):
                        col_lens.append(default_col_len)

                    if ( len(self.data.shape) == 2 and self.data.shape[1] != 0 and isinstance(self.data[0,0], six.string_types)) or ( len(self.data.shape) == 1 and isinstance(self.data[0], six.string_types)):
                        data_col_formats   = [ '%' + '%s' % (col_lens[col_count]) + 's' for col_count in range(self._dims[1]) ]
                        data_col_types     = [ 's' for col_count in range(self._dims[1]) ]
                    else:
                        data_col_formats   = [ '%' + '%d' % (col_lens[col_count]) + '.%de' % default_precision for col_count in range(self._dims[1]) ]
                        data_col_types     = [ 'e' for col_count in range(self._dims[1]) ]
                    lbl_col_formats    = [ '%' + '%d' % (col_lens[col_count]) + 's' for col_count in range(self._dims[1]) ]
                    if len(lbl_col_formats) > 0:
                        lbl_col_formats[0] = '%' + '%d' % (col_lens[col_count]-1) + 's' # account for comment char

            # Write labels column comment
            if len(self.labels) > 0:
                file_obj.write("#")
                for i in range(len(self.labels)):
                    file_obj.write(lbl_col_formats[i] % self.labels[i])
                file_obj.write("\n")

            # Write units column comment
            if len(self.units) > 0:
                file_obj.write("#")
                for i in range(self._dims[1]):
                    if i < len(self.units):
                        file_obj.write(lbl_col_formats[i] % self.units[i])
                    else:
                        file_obj.write(lbl_col_formats[i] % 'unknown')
                file_obj.write("\n")

            # Write data values
            if verbose:
                print(self.verbose_pre, 'Writing data values')

            for i in range(self._dims[0]):
                for j in range(self._dims[1]):
                    if type(self.data) is my_list_type:
                        data_value = self.data[i]
                    else:
                        data_value = self.data[i,j]

                    if data_value == None or (isinstance(data_value, six.string_types) and len(data_value) == 0):
                        file_obj.write(lbl_col_formats[j] % FILL_VALUE)
                    elif data_col_types[j] == 'e' or data_col_types[j] == 'd':
                        file_obj.write(data_col_formats[j] % float(data_value))
                    else:
                        file_obj.write(data_col_formats[j] % str(data_value))
                file_obj.write("\n")

        # Close file only if we opened it
        if isinstance(file_output, six.string_types):
            file_obj.close()
