#!/usr/bin/env python

import os
import shutil
import itertools
from optparse import OptionParser

import numpy

from OCO_Matrix import OCO_Matrix
import L2_Input
from OCO_TextUtils import extract_run_names

SOUNDINGINFO_DEST_NAME = "soundinginfo_%s.dat"
ATMOSPHERE_DEST_NAME = "atmosphere_%s.dat"
PSURF_DEST_NAME = "psurf_%s.dat"

AEROSOL_DEST_NAME = "aerosol_%s.dat"
ALBEDO_DEST_NAME = "albedo_%s.dat"

def get_inp_file_path(inp_fileobj, in_filename):
    if not os.path.exists(in_filename):
        out_filename = os.path.realpath( os.path.join(os.path.dirname(inp_fileobj.filename), in_filename) )
    else:
        out_filename = in_filename

    if not os.path.exists(out_filename):
        raise IOError("Could not find %s either as absolute or relative from current path or relative to input file location: %s" % (in_filename, os.path.dirname(inp_fileobj.filename)))

    return out_filename

def write_source_into_header(output_filename, src_filename):

    file_obj = OCO_Matrix(output_filename, as_strings=True)
    file_obj.header["source_filename"] = os.path.realpath(src_filename)
    file_obj.write(output_filename)

def aggregate_matrix_files(input_files, output_filename, src_column_names=None, dst_column_names=None):

    # Create combined atmosphere file from values in all GAS sections
    # Will overwrite columns of the same name in sequential order

    # First gather contents of various files
    uniq_file_ids = []
    uniq_labels = []
    uniq_units = []
    uniq_filenames = []
    file_objs = []
    max_num_rows = -1
    for curr_input_filename in input_files:
        if curr_input_filename not in uniq_filenames:
            print "%s -> %s" % (curr_input_filename, output_filename)
            uniq_filenames.append(curr_input_filename)
            curr_matobj = OCO_Matrix(curr_input_filename)
            file_objs.append(curr_matobj)

            max_num_rows = max(max_num_rows, curr_matobj.dims[0])

            for (curr_lbl, curr_unit) in itertools.izip_longest(curr_matobj.labels, curr_matobj.units, fillvalue=""):
                if not curr_lbl in uniq_labels:
                    uniq_labels.append(curr_lbl)
                    uniq_units.append(curr_unit)
                    
            if not curr_matobj.file_id in uniq_file_ids:
                uniq_file_ids.append(curr_matobj.file_id)
            
    out_matobj = OCO_Matrix()
    if len(uniq_file_ids) == 1:
        out_matobj.file_id = uniq_file_ids[0]
    else:
        out_matobj.file_id = uniq_file_ids

    if src_column_names == None:
        src_column_names = uniq_labels

    if dst_column_names == None:
        dst_column_names = uniq_labels

    src_column_names = [ curr_col.upper() for curr_col in src_column_names ]
    dst_column_names = [ curr_col.upper() for curr_col in dst_column_names ]

    out_matobj.labels = dst_column_names
    out_matobj.units = uniq_units

    out_matobj.data = numpy.zeros((max_num_rows, len(dst_column_names)), dtype=float)

    # Now aggregate data
    for curr_matobj in file_objs:
        out_matobj.header.update(curr_matobj.header)
        n_src_rows = curr_matobj.dims[0]
        
        for (src_col_idx, col_name) in enumerate(curr_matobj.labels):
            if col_name.upper() in src_column_names:
                dst_col_idx = src_column_names.index(col_name.upper())

                out_matobj.data[:,dst_col_idx] = 0 # make sure not to overlay mismatched row sizes
                out_matobj.data[:n_src_rows, dst_col_idx] = curr_matobj.data[:, src_col_idx]

    out_matobj.write(output_filename)

    if len(uniq_filenames) == 1:
        src_filename = uniq_filenames[0]
    else:
        src_filename = ", ".join(uniq_filenames)

    write_source_into_header(output_filename, src_filename)

def extract_soundinginfo_file(inp_fileobj, output_filename):
    soundinginfo_sect = inp_fileobj.Get_Section("SOUNDING_INFO")
    soundinginfo_file = get_inp_file_path(inp_fileobj, soundinginfo_sect[0].Get_Keyword_Value("soundinginfo_file"))

    print "%s -> %s" % (soundinginfo_file, output_filename)

    # Handle converting sounding info files without a HEADER block
    si_fileobj = L2_Input.Input_File(soundinginfo_file)

    header_section = si_fileobj.Get_Section("HEADER")
    if len(header_section) > 0:
        key_dict = header_section[0].Get_Keywords_Dict()
    else:
        key_dict = si_fileobj.Get_Keywords_Dict()

    output_matobj = OCO_Matrix()
    output_matobj.header.update(key_dict)
    output_matobj.file_id = "Sounding Information"
    output_matobj.write(output_filename)

    write_source_into_header(output_filename, soundinginfo_file)
    
def extract_atmosphere_file(inp_fileobj, output_filename):
    atm_sections = inp_fileobj.Get_Section("PARAMETER_DEFINITION->GAS")
    atm_sections += inp_fileobj.Get_Section("PARAMETER_DEFINITION->TEMPERATURE")

    a_priori_filenames = []
    for curr_sect in atm_sections:
        a_priori_filenames.append( get_inp_file_path(inp_fileobj, curr_sect.Get_Keyword_Value("a_priori")) )

    aggregate_matrix_files(a_priori_filenames, output_filename)

def extract_aerosol_file(inp_fileobj, output_filename):
    aer_sections = inp_fileobj.Get_Section("PARAMETER_DEFINITION->AEROSOL")

    a_priori_filenames = []
    src_column_names = ['PRESSURE']
    dst_column_names = ['PRESSURE']
    for curr_sect in aer_sections:
        a_priori_filenames.append( get_inp_file_path(inp_fileobj, curr_sect.Get_Keyword_Value("a_priori")) )
        
        src_column_names.append( curr_sect.Get_Keyword_Value("name") )
        mie_file = curr_sect.Get_Keyword_Value("mie_file")
        dst_column_names.append( os.path.basename(mie_file).replace('.mie','').upper() )

    aggregate_matrix_files(a_priori_filenames, output_filename, src_column_names, dst_column_names)

def extract_psurf_file(inp_fileobj, output_filename):
    psurf_section = inp_fileobj.Get_Section("PARAMETER_DEFINITION->SURFACE_PRESSURE")
    psurf_filename = get_inp_file_path(inp_fileobj, psurf_section[0].Get_Keyword_Value("a_priori"))

    print "%s -> %s" % (psurf_filename, output_filename)
    shutil.copyfile(psurf_filename, output_filename)

    write_source_into_header(output_filename, psurf_filename)

def extract_albedo_file(inp_fileobj, output_filename):
    albedo_section = inp_fileobj.Get_Section("PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT")
    albedo_filename = get_inp_file_path(inp_fileobj, albedo_section[0].Get_Keyword_Value("a_priori"))

    print "%s -> %s" % (albedo_filename, output_filename)
    shutil.copyfile(albedo_filename, output_filename)

    write_source_into_header(output_filename, albedo_filename)

def create_dataset_from_inp_files(inp_files, output_dir):

    output_names = []
    for curr_name in extract_run_names(inp_files):
        curr_name = curr_name.replace("oco_l2_","").replace(".inp","")
        output_names.append(curr_name)
   
    for curr_inp_filename, curr_out_name in zip(inp_files, output_names):
        # Directory where to place each sounding's files
        curr_out_dir = os.path.join(output_dir, curr_out_name)
        if not os.path.exists(curr_out_dir):
            os.makedirs(curr_out_dir)

        print "Creating %s from %s" % (curr_out_dir, curr_inp_filename)

        inp_fileobj = L2_Input.Input_File(curr_inp_filename)
        
        extract_soundinginfo_file(inp_fileobj, os.path.join(curr_out_dir, SOUNDINGINFO_DEST_NAME % curr_out_name))
        extract_atmosphere_file(inp_fileobj, os.path.join(curr_out_dir, ATMOSPHERE_DEST_NAME % curr_out_name))
        extract_psurf_file(inp_fileobj, os.path.join(curr_out_dir, PSURF_DEST_NAME % curr_out_name))
        extract_aerosol_file(inp_fileobj, os.path.join(curr_out_dir, AEROSOL_DEST_NAME % curr_out_name))
        extract_albedo_file(inp_fileobj, os.path.join(curr_out_dir, ALBEDO_DEST_NAME % curr_out_name))

        print ""

def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] <inp_file> [ <inp_file> ... ]")

    parser.add_option( "-o", "--output_dir", dest="output_dir",
                       metavar="DIR", default="./",
                       help="location where output files are written")
        
    # Parse command line arguments
    (options, args) = parser.parse_args()

    create_dataset_from_inp_files(args, options.output_dir)

if __name__ == "__main__":
    standalone_main()


