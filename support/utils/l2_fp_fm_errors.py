#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import os
import sys
import re

from collections import namedtuple
from six.moves import zip_longest

from full_physics.fp_perturbation import BAND_NAMES, register_multiple_perturb_types, register_perturb_init, FmPerturbations, log_info, PerturbTypeDesc, filter_perturb_funcs 

# Set up ABSCO files that should be used for base calculation
# and spectroscopy error cases
# User will have to make sure the abscodir points to these
# properly

# Path under abscodir where these ABSCO files are located
ABSCO_BASE_PATH = "error_analysis_411"

BASE_ABSCO_FILES = {
        'CO2': "unperturbed/co2_v4.1.1-lowres.hdf",
        'H2O': "unperturbed/h2o_v4.1.1-lowres.hdf",
        'O2':  "unperturbed/o2_v4.1.0_baseline.hdf",
        }

PERT_ABSCO_FILES = {
        #         filename                                      perturbation amount
        'CO2': ( ("air_broaden/strong_co2_v4.1.1-combined.hdf", 0.005),
                 ("air_broaden/weak_co2_v4.1.1-combined.hdf", 0.005), 
                 ("h2o_broaden/strong_co2_v4.1.1-combined.hdf", 0.005),
                 ("h2o_broaden/weak_co2_v4.1.1-combined.hdf", 0.005),
                 ("intensity/strong_co2_v4.1.1-combined.hdf", 0.005),
                 ("intensity/weak_co2_v4.1.1-combined.hdf", 0.005),
                 ("line_mixing/strong_co2_v4.1.1-combined.hdf", 0.05),
                 ("line_mixing/weak_co2_v4.1.1-combined.hdf", 0.05),
                 ("pressure_shift/strong_co2_v4.1.1-combined.hdf", 0.005),
                 ("pressure_shift/weak_co2_v4.1.1-combined.hdf", 0.005),
                 ("speed/strong_co2_v4.1.1-combined.hdf", 0.05),
                 ("speed/weak_co2_v4.1.1-combined.hdf", 0.05),
                 ("temperature/strong_co2_v4.1.1-combined.hdf", 0.005),
                 ("temperature/weak_co2_v4.1.1-combined.hdf", 0.005),
               ),
        'O2': ( 
                ("air_broaden/o2_v4.1.0_yang.hdf", 0.005),
                ("cia/o2_v4.1.1-lowres.hdf", 0.05),
                ("dicke/o2_v4.1.1-lowres.hdf", 0.05),
                ("intensity/o2_v4.1.0_intensity.hdf", 0.005),
                ("line_mixing/o2_v4.1.1-lowres.hdf", 0.05),
                ("pressure_shift/o2_v4.1.0_pressure_shift.hdf", 0.005),
                ("speed/o2_v4.1.1-lowres.hdf", 0.05),
                ("temperature/o2_v4.1.0_temperature.hdf", 0.05),
            ),
        }

# For radiometry gain derivative
GAIN_PERTURB = 1.01

################################################################################

@register_multiple_perturb_types("ils")
def register_ils_errors():
    for band_idx, band_name in enumerate(BAND_NAMES):
        def make_case_func(band_idx, band_name):
            def case_error_func(ls, lua_config):
                return ils_band_error(ls, lua_config, band_idx, band_name)
            return case_error_func
        yield band_name, make_case_func(band_idx, band_name) 

def ils_band_error(ls, lua_config, band_idx, band_name, dl_pert=0.01):
    "Modifies the ILS delta lambda and rebuilds ILS table object to perform perturbation" 

    # Store a copy so we can update then reset the table
    ils_table = lua_config.ils_func[band_idx+1]
    wavenumber = ils_table.wavenumber.copy()
    delta_lambda = ils_table.delta_lambda.copy() 
    response = ils_table.response.copy()

    # Perturb delta lambda
    delta_lambda_pert = delta_lambda.copy()
    delta_lambda_pert = delta_lambda_pert * (1 + dl_pert)
    ils_table.create_delta_lambda_to_response(wavenumber, delta_lambda_pert, response)

    # Halt execution till RT is done
    yield PerturbTypeDesc("Ils_Error/%s" % band_name, "ILS delta lambda perturbation for Band: %s" % band_name, dl_pert)

    # Reset ILS table to their initial state
    ils_table.create_delta_lambda_to_response(wavenumber, delta_lambda, response)

def parse_absco_case_name(gas_name, table_filename):
    "Helper for parsing table filename into something that can be used as a case name / hdf group"

    pert_name_match = re.search("(.*)_v4\.1\.[01](-combined)?(-lowres)?(_.+)?\.hdf$", table_filename)
    if not pert_name_match:
        raise ValueError("Could not parse table filename to create case name: %s" % table_filename)

    return pert_name_match.groups()[0]

@register_perturb_init("spectroscopy")
def spectroscopy_init(ls, lua_config):
    # Make sure we use the absco files necessary
    # for fm error testing
    for gas_name in lua_config.fm.atmosphere.absorber.gases:
        lua_config.fm.atmosphere.absorber[gas_name].absco = os.path.join(ABSCO_BASE_PATH, BASE_ABSCO_FILES[gas_name])

    return True

@register_multiple_perturb_types("spectroscopy")
def register_spectroscopy_errors():
    "Loops over all gas types and perturbed tables, calculating using the table, then resetting filenames to base"

    for pert_gas in list(PERT_ABSCO_FILES.keys()):
        for pert_table, pert_amount in PERT_ABSCO_FILES[pert_gas]:
            # Must use double functions so that we can bind the three loop           
            # variables properly, else we end up always passing last vaules
            # if we only used case_error_func
            def make_case_func(pert_gas_lcl, pert_table_lcl, pert_amount_lcl):
                def case_error_func(ls, lua_config):
                    return spectroscopy_error_case(lua_config, pert_gas_lcl, pert_table_lcl, pert_amount_lcl)
                return case_error_func
            yield parse_absco_case_name(pert_gas, pert_table), make_case_func(pert_gas, pert_table, pert_amount)

def spectroscopy_error_case(lua_config, pert_gas, pert_table, pert_amount):
    # Use perturbed table
    if("abscodir" in os.environ):
        absco_dir = os.environ["abscodir"]
    else:
        absco_dir = lua_config.absco_path
    absco_dir = os.path.join(absco_dir, ABSCO_BASE_PATH)
    absco_fn = os.path.join(absco_dir, pert_table)

    if not os.path.exists(absco_fn):
        raise IOError("Could not find absco file: %s" % (absco_fn))

    log_info("Using ABSCO %s file: %s\n" % (pert_gas, absco_fn))
    lua_config.absorber.gas_absorption(pert_gas).load_file(absco_fn)

    yield PerturbTypeDesc("%s/%s" % ("Spectroscopy", parse_absco_case_name(pert_gas, pert_table)),
                        "ABSCO perturbation using table: %s" % absco_fn,
                        pert_amount)

    # Set back base table
    orig_absco_fn = os.path.join(absco_dir, BASE_ABSCO_FILES[pert_gas])
    log_info("Resetting absco table for %s back to: %s\n" % (pert_gas, orig_absco_fn))
    lua_config.absorber.gas_absorption(pert_gas).load_file(orig_absco_fn)

################################################################################

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(
    usage="""usage: %prog <lua config file> <output file>

    Run a Level 2 retrieval. Then a series of forward model
    perturbation calculations. 
    """)

    parser.add_option( "-u", "--use_existing", dest="use_existing",
                       action="store_true",
                       default=False,
                       help="Use existing output file's statevector as the initial guess if it exists")

    parser.add_option( "-r", "--check_reset", dest="check_reset",
                       action="store_true",
                       default=False,
                       help="Check that cases are resetting the forward model back to its unperturbed state")

    parser.add_option( "-f", "--filter", dest="error_type_filters",
                       action="append",
                       default=[],
                       help="Add a error type name to be run. If no filters applied then all types are run")

    parser.add_option( "-l", "--list_types", dest="list_error_types",
                       action="store_true",
                       default=False,
                       help="List error type that will be run with the current set of filters then exit")

    parser.add_option( "--no_retrieval", dest="run_retrieval",
                       action="store_false",
                       default=True,
                       help="Skips the retrieval step")

    parser.add_option( "--no_stokes", dest="save_stokes",
                       action="store_false",
                       default=True,
                       help="Skips convolving unperturbed stokes components and saving the result")

    parser.add_option( "--num_perturbed", dest="num_perturbed_iterations",
                       type=int, default=0,
                       help="Number of iterations for calculating perturbed residual from perturbed context with initial state vector (1 FM+Jac, 1 FM only calculation). A value of 0 disables perturbed retrievals. Default is 0.")

    # Parse command line arguments
    (options, args) = parser.parse_args()
    
    if options.list_error_types:
        print("Will process the following error types:")
        for name in sorted(filter_perturb_funcs(options.error_type_filters).keys()):
            print(name)
        sys.exit(1)

    if len(args) == 2:
        config_filename, output_file = args

        # Fully expand filenames so we can chdir without worrying
        # about not using paths as supplied on command line
        config_filename, output_file = [ os.path.realpath(fn) for fn in args ]

        fm_errors = FmPerturbations(config_filename, output_file, use_existing=options.use_existing, check_reset=options.check_reset, perturb_type_filters=options.error_type_filters, run_retrieval=options.run_retrieval, save_stokes=options.save_stokes, num_perturbed_iterations=options.num_perturbed_iterations, gain_perturb=GAIN_PERTURB, output_group="Forward_Model_Errors")
        fm_errors.run()

    else:
        parser.error("Need to specify all the arguments")
