#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import filter
from builtins import next
from builtins import range
from builtins import object
from past.utils import old_div
import os
import sys
import re
from types import GeneratorType, StringType
from collections import namedtuple
from six.moves import zip_longest

import h5py
import numpy

from full_physics import *

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

# Names for HDF file
BAND_NAMES = ('ABO2', 'WCO2', 'SCO2')
STOKES_NAMES = ( 'I', 'Q', 'U' )

# For radiometry gain derivative
GAIN_PERTURB = 1.01

logger = FpLogger(FpLogger.INFO)
Logger.set_implementation(logger)
def log_info(info_str):
    logger.write(FpLogger.INFO, info_str)

def log_flush():
    logger.flush(FpLogger.INFO)

_fm_init_funcs = {}
def register_error_init(filter_name):
    "Registers a init function performed at program start"
    def do_register(func):
        _fm_init_funcs[filter_name] = func
        return func
    return do_register

_fm_error_funcs = {}
def register_error_type(filter_name):
    "Registers a functions as an error type to loop over"
    def do_register(func):
        _fm_error_funcs[filter_name] = func 
        return func
    return do_register

def register_multiple_error_types(prefix_name):
    "Registers a functions as an error type to loop over"
    def do_register(register_func):
        register_generator = register_func()
        for filter_name, filter_func in register_generator:
            _fm_error_funcs['%s/%s' % (prefix_name, filter_name)] = filter_func 
        return None
    return do_register

def filter_error_funcs(func_dict, filters):
    "Returns error function dict with only functions matching the filter string list"

    if len(filters) == 0:
        return func_dict

    filtered = {}
    for regex in filters:
        r = re.compile(regex)
        mkeys = list(filter(r.search, list(func_dict.keys())))
        for k in mkeys:
            filtered[k] = func_dict[k]
    return filtered

def filter_init_funcs(init_dict, case_names):
    "Returns init function dict with intialization routines that are valid for the filtered error function names"

    if len(case_names) == 0:
        return init_dict

    filtered = {}
    for key in list(init_dict.keys()):
        r = re.compile(key)
        mkeys = list(filter(r.search, case_names))
        if len(mkeys) > 0 and not key in filtered: 
            filtered[key] = init_dict[key]
    return filtered

def create_sounding_dataset(parent_group, dataset_name, in_data, shape_names=(), units=None):
    "Creates a new dataset but with the data reshaped so that the first dim is the sounding one"

    new_data = numpy.array(in_data)
    try:
        new_ds = parent_group.create_dataset(dataset_name, data=new_data.reshape(1, *new_data.shape))
    except RuntimeError as e:
        raise RuntimeError("Could not create dataset: %s with data: %s\n%s" % (dataset_name, new_data, e))

    # Adds optionally passed shape names as an attribute to the dataset
    # Missing dimensions will be called Dim%d
    if type(shape_names) is StringType:
        shape_names = [shape_names]

    all_shape_names = []
    for idx, s_name in zip_longest(list(range(len(new_ds.shape))), ["Retrieval"] + list(shape_names)):
        if s_name:
            all_shape_names.append(s_name)
        else:
            all_shape_names.append("Dim%d" % (idx+1))
    new_ds.attrs["Shape"] = [ numpy.array("_".join(all_shape_names) + "_Array") ]

    if units:
        new_ds.attrs["Units"] = [ numpy.array(units) ]

    return new_ds

################################################################################

class ErrorTypeDesc(namedtuple('ErrorType', 'name description perturb_amount')):
    "Holds a description of the error type and handles writing the error type HDF dataset"

    def write_deriv_dataset(self, parent_group, perturb_rad_diff, units=None):
        curr_err_grp = parent_group.require_group(self.name)

        if self.description:
            create_sounding_dataset(curr_err_grp, "description", self.description)

        create_sounding_dataset(curr_err_grp, "derivative", old_div(perturb_rad_diff, self.perturb_amount), "SciColor", units)
        create_sounding_dataset(curr_err_grp, "perturbation_amount", self.perturb_amount)

        return curr_err_grp

    def write_perturb_datasets(self, parent_group, state_vector, residual, units=None):
        curr_err_grp = parent_group.require_group(self.name)

        create_sounding_dataset(curr_err_grp, "perturbed_statevector", state_vector, "SVElem")
        create_sounding_dataset(curr_err_grp, "perturbed_residual", residual, "SciColor", units)

        return curr_err_grp

class StokesStorage(ObserverStokesUpdate):
    def __init__(self):
        ObserverStokesUpdate.__init__(self)

        # Store stokes spectras by stoke type, all bands together
        self.stokes_specs = [ [] for s in range(len(STOKES_NAMES)) ] 

    def notify_add(self, stokes_vec):
        pass

    def notify_add(self):
        pass

    def notify_remove(self, stokes_vec):
        pass
    def notify_remove(self):
        pass

    def notify_update(self, stokes_vec):
        for stoke_idx, curr_stoke_spec in enumerate(stokes_vec):
            spec_index = curr_stoke_spec.index
            stoke_list = self.stokes_specs[stoke_idx]

            while len(stoke_list) <= spec_index:
                stoke_list.append(None)

            # We need to pull the data out of the vector and
            # store it inside a C++ object created on the Python
            # side or else the memory will be released
            stoke_list[spec_index] = Spectrum(curr_stoke_spec.spectral_domain, curr_stoke_spec.spectral_range)

################################################################################

@register_multiple_error_types("ils")
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
    yield ErrorTypeDesc("Ils_Error/%s" % band_name, "ILS delta lambda perturbation for Band: %s" % band_name, dl_pert)

    # Reset ILS table to their initial state
    ils_table.create_delta_lambda_to_response(wavenumber, delta_lambda, response)

def parse_absco_case_name(gas_name, table_filename):
    "Helper for parsing table filename into something that can be used as a case name / hdf group"

    pert_name_match = re.search("(.*)_v4\.1\.[01](-combined)?(-lowres)?(_.+)?\.hdf$", table_filename)
    if not pert_name_match:
        raise ValueError("Could not parse table filename to create case name: %s" % table_filename)

    return pert_name_match.groups()[0]

@register_error_init("spectroscopy")
def spectroscopy_init(ls, lua_config):
    # Make sure we use the absco files necessary
    # for fm error testing
    for gas_name in lua_config.fm.atmosphere.absorber.gases:
        lua_config.fm.atmosphere.absorber[gas_name].absco = os.path.join(ABSCO_BASE_PATH, BASE_ABSCO_FILES[gas_name])

    return True

@register_multiple_error_types("spectroscopy")
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

    yield ErrorTypeDesc("%s/%s" % ("Spectroscopy", parse_absco_case_name(pert_gas, pert_table)),
                        "ABSCO perturbation using table: %s" % absco_fn,
                        pert_amount)

    # Set back base table
    orig_absco_fn = os.path.join(absco_dir, BASE_ABSCO_FILES[pert_gas])
    log_info("Resetting absco table for %s back to: %s\n" % (pert_gas, orig_absco_fn))
    lua_config.absorber.gas_absorption(pert_gas).load_file(orig_absco_fn)

################################################################################

class FmErrors(object):

    def __init__(self, config_filename, output_file, use_existing=False, check_reset=False, error_type_filters=[], run_retrieval=True, save_stokes=True, num_perturbed_iterations=0):
        # Load Lua state, and config object
        self.ls, self.lua_config = l2_lua.load_lua_config(config_filename)

        # Where results are saved
        self.output_file = output_file

        # Check after the reset stage of an error function to ensure
        # that the state vector has been reset to its original state
        self.check_reset = check_reset

        # Optionally save per stokes radiances to the output file
        self.save_stokes = save_stokes

        # Optionally perform a one iteration retrieval and then FM only calculation on the pertubed FM setup with the initial SV
        self.num_perturbed_iterations = num_perturbed_iterations
     
        # Get the list of error case functions that will be run
        self.filtered_error_funcs = filter_error_funcs(_fm_error_funcs, error_type_filters)

        # Ensure that high res spectra writer is disabled since it duplicates
        # some of the behavior of this class
        self.lua_config.write_high_res_spectra = False
        
        self._init_error_funcs()

        # Now create everything
        self.lua_config.do_config(self.lua_config)

        # Save initial state vector for back tracking later
        self.initial_sv = self.lua_config.state_vector.state.copy()
        
        # Start off with retrieved state of last retrieval, if it exists
        # This is mainly for debugging purposes, to reduce the time spent
        # waiting on iterations
        if use_existing and os.path.exists(output_file):
            l2_lua.load_retrieved_state_vector(self.ls, self.lua_config, output_file)

    def _init_error_funcs(self):
        # Perform initialization routines that need to be done
        # before all error type functions
        for init_name, init_func in list(filter_init_funcs(_fm_init_funcs, list(self.filtered_error_funcs.keys())).items()):
            log_info("Calling init for: %s" % init_name) 
            if not init_func(self.ls, self.lua_config):
                raise RuntimeError("Initialization function %s did not return as succesful" % init_name)


    def calc_final_radiance_jacobian(self, meas_spec, out_hdf, fm_err_grp):
        log_info("Calculating FM + Jacobians for final state vector update (unperturbed)\n")
        skip_jacobian = False
        unpert_spec = self.lua_config.fm.config.forward_model.radiance_all(skip_jacobian)
        unpert_spec_range = unpert_spec.spectral_range.convert(meas_spec.spectral_range.units)

        unpert_grp = fm_err_grp.create_group("Unperturbed")
        create_sounding_dataset(unpert_grp, "modeled_radiance", unpert_spec_range.data, "SciColor", unpert_spec_range.units.name)
        create_sounding_dataset(out_hdf, "Jacobian", self.lua_config.conn_solver.jacobian, ("SciColor", "SVElem"))

        return unpert_spec_range

    def save_gain_datasets(self, unpert_spec_range, fm_err_grp, band_num_pix):
        # Radiometry gain
        log_info("Saving gain datasets\n")
        for band_idx, band_name in enumerate(BAND_NAMES):
            band_beg = len(band_num_pix[:band_idx]) > 0 and numpy.sum(band_num_pix[:band_idx]) or 0
            band_end = band_beg + band_num_pix[band_idx]
            gain_radiance = unpert_spec_range.data.copy()
            gain_radiance[band_beg:band_end] *= GAIN_PERTURB
            gain_desc = ErrorTypeDesc("Gain/%s" % band_name, "Radiometry gain derivative for Band: %s" % band_name, GAIN_PERTURB)
            gain_desc.write_deriv_dataset(fm_err_grp, gain_radiance, unpert_spec_range.units.name)

    def save_stokes_data(self, stokes_storage, fm_err_grp, band_num_pix):
        log_info("Applying spectrum corrections to stokes arrays (I, Q, U)\n")
        stokes_grp = fm_err_grp.create_group("Stokes")
        for save_band_idx, save_band_name in enumerate(BAND_NAMES):
            for stoke_idx, stoke_list in enumerate(stokes_storage.stokes_specs):
                log_info("Stokes: %s, %s\n" % (STOKES_NAMES[stoke_idx], save_band_name))
                all_band_data = []
                for curr_spec_idx, curr_stoke in enumerate(stoke_list):
                    if save_band_idx == curr_spec_idx:
                        corr_spec = self.lua_config.forward_model.apply_spectrum_corrections(curr_stoke, curr_spec_idx)
                        all_band_data += list(corr_spec.spectral_range.data)
                    else:
                        all_band_data += list(numpy.zeros(band_num_pix[curr_spec_idx]))

                ds_name = "%s_%s" % (save_band_name, STOKES_NAMES[stoke_idx])
                create_sounding_dataset(stokes_grp, ds_name, numpy.array(all_band_data, dtype=float), "SciColor") 

    def perturbed_retrieval(self, meas_spec, err_desc, fm_err_grp, units_name):
        log_info("Running %d iteration retrieval for perturbed FM state\n" % self.num_perturbed_iterations)
        sv = self.lua_config.state_vector
        solver = self.lua_config.fm.config.conn_solver
        ig = self.lua_config.fm.config.initial_guess

        # Save a copy of what we will change during the retrieval
        orig_state = sv.state.copy()
        orig_max_iter = solver.convergence_check.maximum_number_iteration

        # Perform retrieval using initial state vector and save the values
        solver.convergence_check.maximum_number_iteration = self.num_perturbed_iterations 
        sv.update_state(self.initial_sv)

        solver.solve(ig.initial_guess, ig.apriori, ig.apriori_covariance)

        # Perform an additional FM calculation from the retrieved state
        log_info("Performing final perturbed retrieval FM calculation")
        skip_jacobian = True
        pert_spec = self.lua_config.fm.config.forward_model.radiance_all(skip_jacobian)
        pert_range = pert_spec.spectral_range.convert(meas_spec.spectral_range.units)
        residual = pert_range.data - meas_spec.spectral_range.data

        err_desc.write_perturb_datasets(fm_err_grp, solver.x_solution, residual, units_name)

        # Reset back to state before we began
        sv.update_state(orig_state)
        solver.convergence_check.maximum_number_iteration = orig_max_iter

    def run_perturbations(self, meas_spec, unpert_spec_range, out_hdf, fm_err_grp):
        # Helper routine, using closures to make getting the radiance a bit more reusable
        def get_perturbed_radiance_diff():
            skip_jacobian = True
            pert_spec = self.lua_config.fm.config.forward_model.radiance_all(skip_jacobian)
            pert_range = pert_spec.spectral_range.convert(meas_spec.spectral_range.units)

            return pert_range.data - unpert_spec_range.data

        # Loop over forward model perturbations
        log_info("Running forward model perturbation cases.\n")
        for filter_name, err_func in list(self.filtered_error_funcs.items()): 
            err_generator = err_func(self.ls, self.lua_config)

            if type(err_generator) != GeneratorType:
                raise ValueError("Error function: %s did not return a generator function, instead returned: %s" % (err_func, err_generator))

            # A single error function could iterate over several cases
            for err_desc in err_generator:
                log_info( "Calculating forward model error type: %s\n" % err_desc.name )

                # Run forward model and pull out perturbed radiances
                pert_rad_diff = get_perturbed_radiance_diff()

                if numpy.all(pert_rad_diff == 0.0):
                   log_info("WARNING: Case: %s did not end up returning a radiance with any perturbation\n" % (err_desc.name))

                # Write calculated error information and flush hdf buffer
                err_desc.write_deriv_dataset(fm_err_grp, pert_rad_diff, unpert_spec_range.units.name)
                out_hdf.flush()

                # Perform a 1 iteration retrieval on perturbed FM state using apriori statevector
                if self.num_perturbed_iterations > 0:
                    self.perturbed_retrieval(meas_spec, err_desc, fm_err_grp, unpert_spec_range.units.name)
                    out_hdf.flush()

                # Call generator to reset state, it will return True if there are more
                # tests, otherwise a StopIteration exception should have been thrown
                try:
                    more_tests = next(err_generator)
                    if not more_tests:
                        raise RuntimeError("Error function generator for: %s (%s) did not end as expected" % (err_desc.name, err_func))
                    else:
                        log_info("%s has more tests\n" % err_desc.name)
                except StopIteration:
                    log_info("%s has no more tests\n" % err_desc.name)
                    pass

                # If enabled, rerun RT after reset state to make sure everything has been reset correctly
                if self.check_reset:
                    log_info( "Checking that %s properly reset the FM\n" % err_desc.name )

                    pert_rad_diff = get_perturbed_radiance_diff()

                    if not numpy.all(pert_rad_diff == 0.0):
                        raise ValueError("Forward model was not reset correctly for case: %s" % (err_desc.name))

    def run(self):
        """Loads a Lua configuration then performs a retrieval.
        Then calls a list of forward model error functions that
        modify the state. The radiative transfer with jacobians 
        is run. Then the state vector and absco filenames/paths
        are reset to their initial values.
        """

        # Run retrieval
        if run_retrieval:
            l2_lua.run_retrieval(self.ls, self.output_file)

            log_info("Finished Retrieval\n")
        else:
            log_info("Skipping Retrieval\n")
       
        # Measured radiance, need units for converting modeled radiance
        meas_spec = self.lua_config.fm.config.forward_model.measured_radiance_all

        # Save unperturbed stokes values for output
        stokes_storage = StokesStorage()
        self.lua_config.rt.add_observer(stokes_storage)
        
        # Re-open output file and add error terms
        out_hdf = h5py.File(output_file, "r+")
        fm_err_grp = out_hdf.create_group("Forward_Model_Errors")

        # Must rerun FM+Jacobians for final statevector
        unpert_spec_range = self.calc_final_radiance_jacobian(meas_spec, out_hdf, fm_err_grp)
        
        # Use this when modifying things per band
        band_num_pix = out_hdf["/SpectralParameters/num_colors_per_band"][0,:]

        # Radiometry gain
        self.save_gain_datasets(unpert_spec_range, fm_err_grp, band_num_pix)

        # Apply instrument modeling to stokes values and save
        # Save stokes values in one dataset per band and stoke type
        # But save sized the same as modeled radiance setting
        # other bands to zero for easy use in error analysis code
        if self.save_stokes:
            self.save_stokes_data(stokes_storage, fm_err_grp, band_num_pix)

        self.run_perturbations(meas_spec, unpert_spec_range, out_hdf, fm_err_grp)

        # Close file, we are all done
        log_info("Finished error FM calculations. Closing output file.\n")
        log_flush()
        out_hdf.close()

################################################################################

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(
    usage="""usage: %prog <lua config file> <output file>

    This run a Level 2 retrieval. Then a series of forward model
    error runs. The output_file supplied is a base name, two
    files will be produced with _retrieval and _fm_errors
    added before the extension.
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
        for name in sorted(filter_error_funcs(_fm_error_funcs, options.error_type_filters).keys()):
            print(name)
        sys.exit(1)

    if len(args) == 2:
        config_filename, output_file = args

        # Fully expand filenames so we can chdir without worrying
        # about not using paths as supplied on command line
        config_filename, output_file = [ os.path.realpath(fn) for fn in args ]

        fm_errors = FmErrors(config_filename, output_file, options.use_existing, options.check_reset, options.error_type_filters, run_retrieval=options.run_retrieval, save_stokes=options.save_stokes, num_perturbed_iterations=options.num_perturbed_iterations)
        fm_errors.run()

    else:
        parser.error("Need to specify all the arguments")
