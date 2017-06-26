#!/usr/bin/env python

import os
import re
import sys
from optparse import OptionParser

from OCO_Matrix import OCO_Matrix

import numpy

# For reading excel files
import xlrd


ERROR_CODE_FROM_TEXT = dict([ (b, a) for a, b in xlrd.error_text_from_code.items() ])
FILL_CODE = ERROR_CODE_FROM_TEXT['#N/A!']

FILE_ID           = 'GOSAT TANSO FTS Conversion Factors'
FILE_SUFFIX = [ 'B1', 'B2', 'B3' ]

WAVENUMBER_COLUMN_NAME = 'Wavenumber'
DATA_COLUMN_NAMES_TMPL = 'CNV({band}{pol}_gain{gain})'
                                                                                 
NIES_DATA_SHEETS = ['CNV(B1)', 'CNV(B2)', 'CNV(B3)']

KUZE_BAND_ROW_MATCH = 'Band \d'
KUZE_GAIN_ROW_BEG = 'Gain'
KUZE_COLUMN_EXTRACT_MATCH = '(corrected)'
KUZE_LASER_ROW_BEG = 'Primary Sampling laser'
KUZE_LASER_ROW_COL_OFFSET = 5

def get_cell_data(sheet, row_idx, col_idx):
    cell_value = sheet.cell_value(row_idx, col_idx)
    cell_type  = sheet.cell_type(row_idx, col_idx)

    if cell_type == xlrd.XL_CELL_ERROR and cell_value == FILL_CODE:
        cell_data = None
    else:
        cell_data = cell_value

    return cell_data

def write_radcnv_file(output_filename, column_names, sheet, input_rows, input_cols, scaling=1.0):
    file_obj = OCO_Matrix()
    file_obj.labels = column_names
    file_obj.file_id = FILE_ID

    file_obj.data = numpy.zeros((len(input_rows), len(input_cols)), dtype=numpy.chararray)
    out_row_idx = 0
    for in_row_idx, out_row_idx in zip(input_rows, range(len(input_rows)) ):
        for in_col_idx, out_col_idx in zip(input_cols, range(len(input_cols)) ):
            cell_data = get_cell_data(sheet, in_row_idx, in_col_idx)
            if cell_data != None and column_names[out_col_idx] != WAVENUMBER_COLUMN_NAME:
                try:
                    cell_data *= scaling
                except TypeError as e:
                    raise TypeError('%s: cell_data = "%s", scaling = "%s" at row: %d column %d' % (e, cell_data, scaling, in_row_idx, in_col_idx))
            file_obj.data[out_row_idx, out_col_idx] = cell_data

    print 'Writing output filename: %s' % output_filename
    file_obj.write(output_filename)

def convert_radcnv_nies(xls_book, out_file_base):
   
    for sheet_name, suffix in zip(NIES_DATA_SHEETS, FILE_SUFFIX):
        curr_sheet = xls_book.sheet_by_name(sheet_name)

        column_names = [ WAVENUMBER_COLUMN_NAME ]
        for curr_col in curr_sheet.row(0)[1:]:
            curr_name = str(curr_col.value)
            column_names.append(curr_name)

        output_filename = '%s_%s.dat' % (out_file_base, suffix)

        input_rows = range(1, curr_sheet.nrows)
        input_cols = range(curr_sheet.ncols)
        write_radcnv_file(output_filename, column_names, curr_sheet, input_rows, input_cols)

def find_cell_in_sheet(sheet, cell_re, rows, columns, stop_on_first=True):

    locations = []
    for row_idx in rows:
        for col_idx in columns:
            cell_data = get_cell_data(sheet, row_idx, col_idx)
            if re.search(cell_re, str(cell_data)):
                locations.append( (row_idx, col_idx) )
                if stop_on_first:
                    break

    return locations

def convert_radcnv_kuze(xls_book, out_file_base):

    curr_sheet = xls_book.sheet_by_index(0)

    laser_row_loc = find_cell_in_sheet(curr_sheet, KUZE_LASER_ROW_BEG, range(5), range(5))
    if len(laser_row_loc) == 0:
        raise Exception('Searched worksheet and could not find where laser frequency is identified')
    laser_value = get_cell_data(curr_sheet, laser_row_loc[0][0], laser_row_loc[0][1] + KUZE_LASER_ROW_COL_OFFSET)

    gain_row_loc = find_cell_in_sheet(curr_sheet, KUZE_GAIN_ROW_BEG, range(20), range(5))

    if len(gain_row_loc) == 0:
        raise Exception('Searched worksheet and could not find where gains are identified')
    

    # Find where band data columns begin
    band_name_loc = find_cell_in_sheet(curr_sheet, KUZE_BAND_ROW_MATCH, range(1, gain_row_loc[0][0]), range(curr_sheet.ncols), stop_on_first=False)

    if len(band_name_loc) == 0:
        raise Exception('Searched worksheet and could not find where bands are identified')

    band_col_extents = []
    for band_idx in range(len(band_name_loc)):
        col_beg = band_name_loc[band_idx][1]
        if band_idx < len(band_name_loc)-1:
            col_end = band_name_loc[band_idx+1][1] - 1
        else:
            col_end = curr_sheet.ncols-1
        band_col_extents.append( (col_beg, col_end) )

    for band_cols, band_suffix in zip(band_col_extents, FILE_SUFFIX):
        beg_col = band_cols[0]
        end_col = band_cols[1]
        
        data_row_beg = gain_row_loc[0][0] + 2
        data_row_end = curr_sheet.nrows-1
        for row_idx in range(data_row_beg, curr_sheet.nrows):
            cell_data = get_cell_data(curr_sheet, row_idx, beg_col)
            if cell_data == None or str(cell_data).strip() == '':
                break
            else:
                data_row_end = row_idx

        gain_names = []
        for gain_col_idx in range(beg_col, end_col+1):
            gain_val = get_cell_data(curr_sheet, gain_row_loc[0][0], gain_col_idx)
            if gain_val == None or gain_val == '' and len(gain_names) > 0:
                gain_val = gain_names[-1]
            gain_names.append( gain_val )

        extract_find = find_cell_in_sheet(curr_sheet, KUZE_COLUMN_EXTRACT_MATCH, [data_row_beg-1], range(beg_col, end_col+1), stop_on_first=False)

        desired_cols = [ beg_col ]
        column_names = [ WAVENUMBER_COLUMN_NAME ]
        for extract_col in extract_find:
            row_idx = extract_col[0]
            col_idx = extract_col[1]
            
            gain_code = gain_names[ col_idx - beg_col ]
            pol_code  = get_cell_data(curr_sheet, row_idx, col_idx)[0:1]

            column_names.append( DATA_COLUMN_NAMES_TMPL.format(band=band_suffix, gain=gain_code, pol=pol_code) )
            desired_cols.append( col_idx )

        output_filename = '%s_%s.dat' % (out_file_base, band_suffix)

        input_rows = range(data_row_beg, data_row_end+1)
        write_radcnv_file(output_filename, column_names, curr_sheet, input_rows, desired_cols, scaling=laser_value*1e-4 / 2.0)
        

def convert_radcnv(xls_file, output_dir):

    out_file_base = os.path.join(output_dir, re.sub('\..*$', '', os.path.basename(xls_file)))

    print 'Reading input file: %s' % xls_file
    xls_book = xlrd.open_workbook(xls_file)

    first_sheet = xls_book.sheet_by_index(0)        
    
    if first_sheet.name == NIES_DATA_SHEETS[0]:
        print 'Using NIES format'
        convert_radcnv_nies(xls_book, out_file_base)
    else:
        print 'Using Kuze format'
        convert_radcnv_kuze(xls_book, out_file_base)

def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] rad_cnv_xls [output_dir]")

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if (len(args) < 1):
        parser.error("Need input_flename")

    xls_file  = args[0]

    if not os.path.exists(xls_file):
        parser.error("%s does not exist" % xls_file)

    if len(args) >= 2:
        output_dir = args[1]
    else:
        output_dir = '.'

    convert_radcnv(xls_file, output_dir)

if __name__ == "__main__":
    standalone_main()

