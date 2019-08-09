from __future__ import absolute_import
from builtins import object
from .try_swig_load import *
from full_physics.l2_input import L2InputFile
import os

class L2Run(object):
    '''This is a convenience class that provides a simpler interface to
    the underlying C++ wrapped python classes.'''
    def __init__(self, lua_config, sounding_id, met_file, spectrum_file, 
                 print_log = True, scene_file = None, co2_pr_file = None,
                 abscodir = None, merradir = None, imap_file = None):
        '''Load a Lua config file, setting the given sounding id, ecmf_file,
        and spectrum_file. By default we also turn logging on so you get 
        progress messages as you run. You can optionally turn this off.'''
        # Save information to use in pickling
        self.__state_save = (lua_config, sounding_id, met_file, 
                             spectrum_file, print_log, 
                             scene_file, abscodir, merradir,
                             imap_file)
        self.output_file_name = None
        self.output_file_ = None
        self.sounding_id = sounding_id
        self.met_file = met_file
        self.spectrum_file = spectrum_file
        self.imap_file = imap_file
        self.co2_pr_file = co2_pr_file
        os.environ["met_file"] = met_file
        os.environ["spectrum_file"] = spectrum_file
        os.environ["sounding_id"] = sounding_id
        if(co2_pr_file is not None):
            os.environ["co2_pr_file"] = co2_pr_file
        if(scene_file is not None):
            os.environ["scene_file"] = scene_file
        if(imap_file is not None):
            os.environ["imap_file"] = imap_file
        if(abscodir is not None):
            os.environ["abscodir"] = abscodir
        if(merradir is not None):
            os.environ["merradir"] = merradir
        self.config = L2FpConfigurationLua(lua_config)
        self.lua_config_name = lua_config
        Logger.set_implementation(self.config.logger)

    def __getstate__(self):
        return {"initial_state": self.__state_save,
                "output_file_name": self.output_file_name, 
                "state_v": self.state_vector.state,
                "state_cov": self.state_vector.state_covariance}

    def __setstate__(self, val):
        self.__init__(*val["initial_state"])
        self.output_file_name = val["output_file_name"]
        self.state_vector.update_state(val["state_v"], val["state_cov"])

    @staticmethod
    def create_from_existing_run(config_file, sounding_id = None,
                                 sounding_id_index = 0, lua_config = None):
        '''This creates a L2Run object that sets everything up from the
        an existing L2 run. We take the configuration file that was
        passed to the populator script, and either the sounding ID, or
        the index into the list of soundings IDs (e.g., 0 for the first index,
        etc.)

        Since it is often useful, we stash the output file that we use into 
        the results as "output_file"

        The state vector is in the state set by the Lua configuration. 
        You may man to look at create_final_state and create_initial_state
        which are closely related functions that set the state vector into
        a known state.

        Default lua_config is the normal gosat/config/config.lua, but you can
        specify a different one if desired
        '''
        c = L2InputFile(config_file)
        t = c.get_section("input->InputProductFiles")[0]
        l1b = t.get_keyword_value("L1BFile")
        met = t.get_keyword_value("ResampledMetFile")
        co2_pr = t.get_keyword_value("CO2PriorFile")
        imap = t.get_keyword_value("IMAPFile")
        if(lua_config is None):
            lua_config = os.path.join(os.environ.get('L2_INPUT_PATH', ''), "gosat/config/config.lua")
        if(sounding_id is None):
            t = c.get_section("input->GOSATFullPhysics->LIST->VALUES")
            if(len(t) == 0):
                t = c.get_section("input->OCOFullPhysics->LIST->VALUES")
            sid = t[0].get_matrix_data()
            sounding_id = sid[sounding_id_index]
        if(t.get_keyword_value("InputFileMapping")):
            inpmap = t.get_keyword_value("InputFileMapping")
            with open(inpmap) as f:
                for ln in f:
                    sid,val = ln.split(' ')
                    if(sid == sounding_id):
                        for v2 in val.split(';'):
                            k,v = v2.split('=')
                            if(k == "spectrum_file"):
                                l1b = eval(v)
                            if(k == "met_file"):
                                met = eval(v)
                            if(k == "co2_pr_file"):
                                co2_pr = eval(v)
                            if(k == "imap_file"):
                                imap = eval(v)

        r = L2Run(lua_config, sounding_id, met, l1b, imap_file = imap,
                  co2_pr_file = co2_pr)
        r.output_file_name = os.path.dirname(config_file) + "/output/l2_" + \
            sounding_id + ".h5"
        return r

    @property
    def output_file(self):
        if(self.output_file_ is None and self.output_file_name is not None):
            self.output_file_ = HdfFile(self.output_file_name)
        return self.output_file_

    @staticmethod
    def create_final_state(config_file, sounding_id = None,
                           sounding_id_index = 0, lua_config = None):
        '''This creates a L2Run object that sets everything up to the final
        state of a Level 2 run. We take the configuration file that was
        passed to the populator script, and either the sounding ID, or
        the index into the list of soundings IDs (e.g., 0 for the first index,
        etc.)

        Since it is often useful, we stash the output file that we use into 
        the results as "output_file"'''
        r = L2Run.create_from_existing_run(config_file, sounding_id,
                                           sounding_id_index, lua_config)
        fstate = r.output_file.read_double_2d("/RetrievedStateVector/state_vector_result")[0,:]
        fcov = r.output_file.read_double_3d("/RetrievalResults/aposteriori_covariance_matrix")[0,:, :]
        r.state_vector.update_state(fstate, fcov)
        return r

    @staticmethod
    def create_apriori_state(config_file, sounding_id = None,
                           sounding_id_index = 0, lua_config = None):
        '''This creates a L2Run object that sets everything up to the apriori
        state of a Level 2 run. We take the configuration file that was
        passed to the populator script, and either the sounding ID, or
        the index into the list of soundings IDs (e.g., 0 for the first index,
        etc.)

        Since it is often useful, we stash the output file that we use into 
        the results as "output_file"'''
        r = L2Run.create_from_existing_run(config_file, sounding_id,
                                           sounding_id_index, lua_config)
        fstate = r.output_file.read_double_2d("/RetrievedStateVector/state_vector_apriori")[0,:]
        fcov = r.output_file.read_double_3d("/RetrievalResults/apriori_covariance_matrix")[0,:, :]
        r.state_vector.update_state(fstate, fcov)
        return r
        
    @property
    def solar_model(self):
        i = 1
        res = []
        while(self.lua_config.solar_model[i] is not None):
            res.append(self.lua_config.solar_model[i])
            i = i + 1
        return res

    @property
    def forward_model(self):
        return self.config.forward_model

    @property
    def state_vector(self):
        return self.forward_model.state_vector

    @property
    def lua_config(self):
        return self.config.lua_state.globals.config

    @property
    def atmosphere(self):
        return self.lua_config.atmosphere

    @property
    def error_analysis(self):
        return self.lua_config.error_analysis

    @property
    def instrument(self):
        return self.forward_model.instrument

    @property
    def spectral_window(self):
        return self.forward_model.spectral_window

    @property
    def level_1b(self):
        return self.forward_model.level_1b

    @property
    def radiative_transfer(self):
        return self.forward_model.radiative_transfer

    @property
    def spectrum_sampling(self):
        return self.forward_model.spectrum_sampling

    @property
    def xco2(self):
        '''xco2 in ppm'''
        return self.atmosphere.absorber.xgas("CO2").value * 1e6

    @property
    def solver(self):
        '''Solver for Level 2 Retrieval'''
        return self.lua_config.conn_solver

    @property
    def cost_function(self):
        '''The cost function'''
        return self.lua_config.conn_solver.cost_function
    
    @property
    def initial_guess(self):
        '''The initial guess object'''
        return self.lua_config.initial_guess





    
