from __future__ import print_function
# Useful routines for dealing with L2 Lua configurations

import os
import sys
import re

import full_physics as fp

def modified_lua_config(config_filename):
    "Loads the Lua configuration file and removes the do_config and finds the variable name"

    config_txt = ""
    config_var_name = None
    with open(config_filename) as config_file_obj:
        for lua_line in config_file_obj:
            # Searches for first line with something like:
            # var_name = OtherVar:new()
            # or var_name = OtherVar:new {
            # and extract the variable name being set
            new_match = re.search('^\s*(\w+)\s*=\s*\w+:new\s*(\(\s*\)|\s*{)\s*$', lua_line) 
            if new_match and config_var_name is None:
                config_var_name = new_match.groups()[0]

            # Look for var_name:do_config()
            # Check that this variable name is the same as the one
            # used for :new, also do not add this line to config_txt
            do_config_match = re.search('^\s*(\w+):do_config\(\s*\)\s*$', lua_line)
            if do_config_match:
                do_var_name = do_config_match.groups()[0] 
                if config_var_name and config_var_name != do_var_name:
                    raise ValueError("Found different variable name: %s for :do_config than %s found for :new in Lua config: %s" % (do_var_name, config_var_name, config_filename))
                else:
                    config_var_name = do_var_name
            else:
                config_txt += lua_line

    return config_var_name, config_txt

def load_lua_config(config_filename):
    """Loads a Lua L2 configuration, but does not yet run do_config
       Also pulls out the configuration object that can be used
       to modify the Lua's configuration further
    """
    
    # Change to directory where Lua config live
    os.chdir(os.path.dirname(config_filename))

    # Load Lua configuration, searching for :do_config()
    # and for the name of the variable to evaluate
    # We remove the do_config so the retrieval is not
    # initiated from within Lua
    ls = fp.LuaState()

    # Get modified Lua configuration that does not perform the do_config
    config_var_name, config_txt = modified_lua_config(config_filename)

    # Check that we know which variable to load
    if not config_var_name:
        raise ValueError("Could not deduce configuration variable name from Lua config file: %s" % config_filename)

    # Evaluate modified config
    ls.run(config_txt)

    return (ls, ls.globals[config_var_name])

def run_retrieval(lua_state, output_file):
    "Does the L2 retrieval given a LuaState and output file name"

    c = fp.L2FpConfigurationLua(lua_state, output_file)

    log_timing = fp.LogTiming()
    print(c.forward_model)
    out, out_err = c.output()
    solver = c.solver
    ig = c.initial_guess
    sv = c.forward_model.state_vector
    c.forward_model.setup_grid()  # Setup forward model grid now that initial
                                  # state is determined.
    log_timing.write_to_log("Initialization")
    solver.add_observer(log_timing)
    try:
        res = solver.solve(ig.initial_guess, ig.apriori, 
                           ig.apriori_covariance)
        sv.update_state(solver.x_solution, solver.aposteriori_covariance);
        out.write()
        log_timing.write_to_log("Final");
        if(res):
            print("Found Solution")
        else:
            print("Failed to find solution")
        print("Bye bye")
    except:
        print("Caught error", file=sys.stderr)
        import traceback
        print("-" * 25, file=sys.stderr)
        traceback.print_exc()
        print("-" * 25, file=sys.stderr)
        out_err.write_best_attempt()

    return c

def load_retrieved_state_vector(ls, lua_config, retrieved_file):
    ret_obj = fp.HdfFile(retrieved_file)
    initial_sv = ret_obj.read_double_2d("/RetrievedStateVector/state_vector_result")[0,:]
    ig_orig = lua_config.initial_guess
    ig_v = fp.InitialGuessValue()
    ig_v.initial_guess = initial_sv
    ig_v.apriori = ig_orig.apriori
    ig_v.apriori_covariance = ig_orig.apriori_covariance
    ig_new = fp.CompositeInitialGuess()
    ig_new.add_builder(ig_v)
    ls.globals["initial_guess"] = ig_new
    lua_config.state_vector.update_state(initial_sv)

