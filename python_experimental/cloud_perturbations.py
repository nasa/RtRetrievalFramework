#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import os
import sys

import h5py
import numpy as np

from full_physics.fp_perturbation import FmPerturbations, PerturbTypeDesc, register_perturb_type, register_multiple_perturb_types, register_perturb_init, filter_perturb_funcs

from full_physics import AerosolExtinctionLinear, AerosolPropertyHdf, AerosolOptical, vector_aerosol_extinction, vector_aerosol_property, Logger, FpLogger

logger = FpLogger(FpLogger.INFO)
Logger.set_implementation(logger)
def log_info(info_str):
    logger.write(FpLogger.INFO, info_str)

def log_flush():
    logger.flush(FpLogger.INFO)

################################################################################

PERTURB_PRESSURE_DS_PREFIX = "/CloudPerturbations/PressureLevels_"
PERTURB_AER_PROFILES_DS_PREFIX = "/CloudPerturbations/Profiles_"
PERTURB_TYPES_USED_DS_PREFIX = "/CloudPerturbations/TypesUsed_"

@register_multiple_perturb_types("cloud_jacobians")
def perturb_cloud_jacobians():

    cloud_perturb_def = { "optical_depth": "tau", "cloud_top_pressure": "ctP", "cloud_bottom_pressure": "cbP" }

    for perturb_name, ds_suffix in cloud_perturb_def.items():
        def make_perturb(perturb_name, ds_suffix):
            def perturb_cloud(ls, lua_config):

                orig_pressure = lua_config.pressure.pressure_grid.value.value
                orig_aerosol = lua_config.atmosphere.aerosol

                new_pressure = get_perturb_pressure(lua_config, PERTURB_PRESSURE_DS_PREFIX + ds_suffix)
                new_aerosol = get_aerosol_setup(lua_config, PERTURB_TYPES_USED_DS_PREFIX + ds_suffix, PERTURB_AER_PROFILES_DS_PREFIX + ds_suffix)

                lua_config.pressure.set_levels_from_grid(new_pressure)
                lua_config.atmosphere.set_aerosol(new_aerosol, lua_config.state_vector)

                yield PerturbTypeDesc(perturb_name, "Perturbation of clouds for %s" % perturb_name, None)

                lua_config.pressure.set_levels_from_grid(orig_pressure)
                lua_config.atmosphere.set_aerosol(orig_aerosol, lua_config.state_vector)

            return perturb_cloud

        yield perturb_name, make_perturb(perturb_name, ds_suffix)

def get_perturb_pressure(lua_config, pressure_dataset):
    hdf_file = lua_config.l1b_hdf_file(lua_config)

    sid = lua_config.l1b_sid_list(lua_config)
    frame_idx = sid.frame_number
    sounding_idx = sid.sounding_number

    new_pressure = hdf_file.read_double_3d(pressure_dataset)[frame_idx, sounding_idx, :]
    log_info("Setting Pressure Levels:\n%s\n" % new_pressure)
 
    return new_pressure

def get_aerosol_setup(lua_config, types_dataset_name, profiles_dataset_name):
    hdf_file = lua_config.l1b_hdf_file(lua_config)

    sid = lua_config.l1b_sid_list(lua_config)
    frame_idx = sid.frame_number
    sounding_idx = sid.sounding_number

    type_indexes = hdf_file.read_int_3d(types_dataset_name)[frame_idx, sounding_idx, :]

    with h5py.File(hdf_file.file_name, "r") as h5py_obj:
        type_names = list(h5py_obj["/Aerosol/TypeNames"].value)
    type_profiles = hdf_file.read_double_4d(profiles_dataset_name)[frame_idx, sounding_idx, :, :]

    aerosols_used = set()
    for idx in type_indexes:
        if idx >= 0:
            aerosols_used.add(type_names[idx])

    log_info("Setting Aerosol Types: %s\n" % list(aerosols_used))

    aerosol_extinction = vector_aerosol_extinction()
    aerosol_properties = vector_aerosol_property()
    for aer_name in aerosols_used:
        profile = type_profiles[type_names.index(aer_name), :][:]
        ret_flag = np.zeros(profile.shape[0], dtype=bool)
        ext = AerosolExtinctionLinear(lua_config.pressure, ret_flag, profile, aer_name)
        aerosol_extinction.push_back(ext)

        log_info("Aerosol profile for %s:\n%s\n" % (aer_name, profile))

        aer_hdf_file = lua_config.h_aerosol(lua_config)
        prop = AerosolPropertyHdf(aer_hdf_file, "%s/Properties" % aer_name, lua_config.pressure)
        aerosol_properties.push_back(prop)
    
    return AerosolOptical(aerosol_extinction, aerosol_properties, lua_config.pressure, lua_config.relative_humidity)

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

    parser.add_option( "--do_retrieval", dest="run_retrieval",
                       action="store_true",
                       default=False,
                       help="Runs the retrieval step")

    parser.add_option( "--save_stokes", dest="save_stokes",
                       action="store_true",
                       default=False,
                       help="Saves unperturbed stokes components into the output product")

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

        fm_errors = FmPerturbations(config_filename, output_file, use_existing=options.use_existing, check_reset=options.check_reset, perturb_type_filters=options.error_type_filters, run_retrieval=options.run_retrieval, save_stokes=options.save_stokes, num_perturbed_iterations=options.num_perturbed_iterations, output_group="Cloud_Perturbations")
        fm_errors.run()

    else:
        parser.error("Need to specify all the arguments")
