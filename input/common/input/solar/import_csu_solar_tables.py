#!/usr/bin/env python
#
# Creates the solar continuum and absorption datasets for inclusion in a l2 static input file
# It uses a combination of the old solar model and tables from CSU to create the tables

import os
import sys
import re
from optparse import OptionParser
from collections import namedtuple

import numpy
import h5py
from scipy.constants import speed_of_light, h as planck

from full_physics import *

OLD_SOLAR_MODEL_INPUT_FILE = os.path.join(os.path.dirname(sys.argv[0]), 'old_solar_model_input.h5')

BAND_RANGES = [ numpy.arange(12700, 13300, 0.01),
                numpy.arange(5800, 6500, 0.01),
                numpy.arange(4700, 5000, 0.01),
              ]

ABS_DATASET_WN = "/Solar/Absorption/Absorption_%d/wavenumber"
ABS_DATASET_SPEC = "/Solar/Absorption/Absorption_%d/spectrum"
ABS_SPEC_UNIT = "dimensionless"

CONT_DATASET_WN = "/Solar/Continuum/Continuum_%d/wavenumber"
CONT_DATASET_SPEC = "/Solar/Continuum/Continuum_%d/spectrum"
CONT_SPEC_UNIT = "ph / s / m^2 / micron"

FRAME_WN_ATTR = "Solar rest frame"
WN_UNIT = "cm^-1"

class FpSolarModelCompute(object):
    def __init__(self, static_inp_file):
        self.si_hdf = HdfFile(static_inp_file, HdfFile.READ_WRITE)
    
    def continuum(self, wn):
        param = self.si_hdf.read_double_with_unit_1d("/Solar/a_priori")
        spoly = SolarContinuumPolynomial(param, False)
        scont = spoly.solar_continuum_spectrum(SpectralDomain(wn)).value

        return scont

    def absorption(self, wn):
        sfil = SolarAbsorptionOcoFile(self.si_hdf, "Solar", 1.0)
        sabs = sfil.solar_absorption_spectrum(SpectralDomain(wn)).value
        
        return sabs

def read_csu_ascii_table(filename, header_size=5):
    # Read header values
    doppler_vel = None
    solar_dist = None
    column_names = None
    cont_units = None
    with open(filename) as table_fn_obj:
        line_no = 0
        while line_no < header_size:
            line_no += 1
            line_txt = table_fn_obj.readline()
            if line_txt.find("Solar Doppler Velocity") >= 0:
                part1, part2 = line_txt.split(":", 1)
                doppler_vel = float(re.sub('#.*$', '', part2))
            elif line_txt.find("Distance to Sun in AU") >= 0:
                part1, part2 = line_txt.split(":", 1)
                solar_dist = float(part2)
            elif line_txt.find("Wavenumber") >= 0:
                column_names = filter(lambda s: len(s) > 0, [s.strip() for s in line_txt.strip().split('  ')])
            elif line_txt.find("Solar Intensity Units") >= 0:
                heading, cont_units = line_txt.split(":", 1)
                cont_units = cont_units.strip()

        if not doppler_vel:
            raise ValueError("Failed to parse doppler velocity")
        elif not solar_dist:
            raise ValueError("Failed to parse solar distance")
        elif not column_names:
            raise ValueError("Failed to parse column names")
        elif not cont_units: 
            raise ValueError("Failed to parse continuum units")

    table_data = numpy.loadtxt(filename, skip_header=header_size)

    SolarTableData = namedtuple('SolarTableData', 'doppler_velocity solar_distance table_columns table_data continuum_units')
    return SolarTableData(doppler_vel, solar_dist, column_names, table_data, cont_units)

def shift_wn(wn, doppler_vel):
    return wn / (1 + doppler_vel/speed_of_light)

def merge_values(inp_wn, inp_vals, dst_wn, dst_vals):
    # Replace a whole section of the fp_cont with the table data
    where_repl_1 = numpy.where(dst_wn >= min(inp_wn))
    where_repl_2 = numpy.where(dst_wn <= max(inp_wn))
    where_repl = list(set(where_repl_1[0]).intersection(set(where_repl_2[0])))

    dst_wn[where_repl] = inp_wn
    dst_vals[where_repl] = inp_vals

    return dst_wn, dst_vals
   
def compute_merged_continuum(input_tables, fp_model):
    wn = []
    cont = []

    for band_idx, band_range in enumerate(BAND_RANGES):
        band_table = input_tables[band_idx]

        # Shift wn to match the input table
        dst_wn = shift_wn(band_range, band_table.doppler_velocity)

        # Convert range to same as used by input table
        fp_cont = numpy.array(fp_model.continuum(dst_wn))

        # Read data and sort to increasing wn
        sort_inds = band_table.table_data[:,1].argsort()
        inp_wn, inp_abs, inp_cont = band_table.table_data[:,1][sort_inds], band_table.table_data[:,3][sort_inds], band_table.table_data[:,4][sort_inds]

        # Convert units
        cont_awu = ArrayWithUnit_double_1()
        cont_awu.data = inp_cont
        cont_awu.units = Unit(band_table.continuum_units)
        if not cont_awu.units.is_commensurate(Unit(CONT_SPEC_UNIT)):
            output_units = Unit("(W / cm^-1) / (ph / s / micron)")
            wavenumber_unit = Unit("cm^-1")
            alpha = DoubleWithUnit(speed_of_light * planck, Unit("m s^-1 J s"))
            convert_awu = ArrayWithUnit_double_1()
            convert_awu.data = alpha.value() / inp_wn 
            convert_awu.units = alpha.units() / wavenumber_unit / Unit("ph")
            convert_awu = convert_awu.convert(output_units)
             
            cont_awu.data = cont_awu.data / convert_awu.data
            cont_awu.units = cont_awu.units / convert_awu.units

        inp_cont = cont_awu.convert_wave(Unit(CONT_SPEC_UNIT)).data

        # Convert input continuum back to 1 AU and remove absorption lines
        inp_cont = inp_cont * (band_table.solar_distance**2) / inp_abs

        out_wn, out_cont = merge_values(inp_wn, inp_cont, dst_wn, fp_cont) 

        wn.append(out_wn)
        cont.append(out_cont)

    return wn, cont

def compute_merged_absorption(input_tables, fp_model):
    wn = []
    absor = []

    for band_idx, band_range in enumerate(BAND_RANGES):
        band_table = input_tables[band_idx]

        # Shift wn to match the input table
        dst_wn = shift_wn(band_range, band_table.doppler_velocity)

        # Convert range to same as used by input table
        fp_absor = numpy.array(fp_model.absorption(dst_wn))

        # Read data and sort to increasing wn
        sort_inds = band_table.table_data[:,1].argsort()
        inp_wn, inp_absor = band_table.table_data[:,1][sort_inds], band_table.table_data[:,3][sort_inds]

        out_wn, out_absor = merge_values(inp_wn, inp_absor, dst_wn, fp_absor) 

        wn.append(out_wn)
        absor.append(out_absor)

    return wn,absor
 
def write_continuum(out_hdf_obj, wn, cont):
    for band_idx, (band_wn, band_cont) in enumerate(zip(wn,cont)):
        wn_ds = out_hdf_obj.require_dataset(CONT_DATASET_WN % (band_idx+1), shape=band_wn.shape, dtype=band_wn.dtype)
        wn_ds[:] = band_wn
        wn_ds.attrs['Units'] = [numpy.array(WN_UNIT+'\0')]
        wn_ds.attrs['Frame'] = [numpy.array(FRAME_WN_ATTR+'\0')]

        spec_ds = out_hdf_obj.require_dataset(CONT_DATASET_SPEC % (band_idx+1), shape=band_cont.shape, dtype=band_cont.dtype)
        spec_ds[:] = band_cont
        spec_ds.attrs['Units'] = [numpy.array(CONT_SPEC_UNIT+'\0')]

def write_absorption(out_hdf_obj, wn, absor):
     for band_idx, (band_wn, band_absor) in enumerate(zip(wn,absor)):
        wn_ds = out_hdf_obj.require_dataset(ABS_DATASET_WN % (band_idx+1), shape=band_wn.shape, dtype=band_wn.dtype)
        wn_ds[:] = band_wn
        wn_ds.attrs['Units'] = [numpy.array(WN_UNIT+'\0')]
        wn_ds.attrs['Frame'] = [numpy.array(FRAME_WN_ATTR+'\0')]

        spec_ds = out_hdf_obj.require_dataset(ABS_DATASET_SPEC % (band_idx+1), shape=band_absor.shape, dtype=band_absor.dtype)
        spec_ds[:] = band_absor
        spec_ds.attrs['Units'] = [numpy.array(ABS_SPEC_UNIT+'\0')]

def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] <solar_input_filenames>")
    
    parser.add_option( "-i", "--old_model_input_file", dest="old_model_input_file",
                       metavar="FILE",
                       default=OLD_SOLAR_MODEL_INPUT_FILE,
                       help="file name with inputs for old solar model to fill in the gaps")

    parser.add_option( "-o", "--output_file", dest="output_file",
                       metavar="FILE",
                       default="solar_tables.h5",
                       help="name for main output filename")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("Need to specify satic input filename")

    # Read CSU data
    print "Reading input tables"
    input_tables = []
    for band_idx, filename in enumerate(args):
        if not filename.find(str(band_idx+1)) >= 0:
            raise ValueError("Could not find index %d in filename: %s. Possibly filenames supplied out of order" % (band_idx+1, filename))
        input_tables.append( read_csu_ascii_table(filename) )

    print "Reading old model data"
    fp_model = FpSolarModelCompute(options.old_model_input_file)

    # Create the new hdf file if it doesnt exist
    print "Creating output file: %s" % options.output_file
    out_hdf_obj = None 
    if not os.path.exists(options.output_file):
        out_hdf_obj = h5py.File(options.output_file, "w")
    else:
        out_hdf_obj = h5py.File(options.output_file, "r+")

    print "Writing continuum data"
    wn_cont, data_cont = compute_merged_continuum(input_tables, fp_model)
    write_continuum(out_hdf_obj, wn_cont, data_cont)

    print "Writing absorption data"
    wn_abs, data_abs = compute_merged_absorption(input_tables, fp_model)
    write_absorption(out_hdf_obj, wn_abs, data_abs)

    print "Done. Closing output file."
    out_hdf_obj.close()

if __name__ == "__main__":
    standalone_main()
