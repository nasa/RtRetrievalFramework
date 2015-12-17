#!/usr/bin/env python

import os
import re
import sys

import numpy
from matplotlib import pyplot as plt

from full_physics import *

# Import path to where l2_fp should be
sys.path.append( os.path.dirname(sys.argv[0]) )
import l2_fp

def compare_jacobians(config_filename, epsilon, filter_names=[], retrieved_sv_file=None): 
    logger = FpLogger(FpLogger.INFO)
    Logger.set_implementation(logger)
    def log_info(info_str):
        logger.write(FpLogger.INFO, info_str)

    def log_flush():
        logger.flush(FpLogger.INFO)
    
    # Load Lua state, and config object
    ls, lua_config = l2_fp.load_lua_config(config_filename)

    # Now create everything
    log_info("Loading configuration\n") ; log_flush()
    lua_config.do_config(lua_config)   

    # Load a retrieved state vector instead of using configured one
    if (retrieved_sv_file):
        log_info("Loading retrieved statevector\n"); log_flush()
        l2_fp.load_retrieved_state_vector(ls, lua_config, retrieved_sv_file)

    # Load state vector and print it out
    sv = lua_config.state_vector
    state0 = sv.state().copy()
    sv_names = sv.state_vector_name()

    log_info("State Vector Contents:\n")
    sv_contents = str(sv).split("\n")
    for svi, sv_line in enumerate(sv_contents):
        log_info("%03d : %s\n" % (svi, sv_line))
    log_flush()
    
    # Figure out which SV elements to use
    sv_used = []
    for svi in range(state0.shape[0]):
        name_matched = False
        for filt in filter_names:
            if re.search(filt, sv_names[svi]):
                name_matched = True
                break

        if name_matched:
            sv_used.append(True)
        else:
            sv_used.append(False)


    # Fix epsilon size
    ep_full = np.zeros(state0.shape[0], dtype=float)
    if not hasattr(epsilon, "__len__"):
        epsilon = [ epsilon ]
    ep_idx = 0
    for svi, used in enumerate(sv_used):
        if used:
            ep_full[svi] = epsilon[ep_idx]
            if (ep_idx + 1) < len(epsilon):
                ep_idx += 1
            else:
                ep_idx = 0
    log_info("Using epsilons: %s\n" % ep_full)

    # Calculate analytic jacobians
    log_info("Running analytic calculation\n"); log_flush()
    spec_idx = 2
    rad0 = lua_config.forward_model.radiance(spec_idx)
    jac_ad = rad0.spectral_range.data_ad.jacobian

    log_info("Analytic jacobian shape: %s\n" % str(jac_ad.shape))
    for svi in range(state0.shape[0]):
        nan_compare = jac_ad[:,svi] != jac_ad[:,svi]
        if np.any(nan_compare):
            nan_count = len(np.where(nan_compare)[0])
            log_info("%s jacobian has %d nans\n" % (sv_names[svi], nan_count))
    log_flush()

    # Calculate finite difference jacobians 
    log_info("Calculating finite difference jacobians\n")
    for svi, used in enumerate(sv_used): 
        if not used: 
            log_info(sv_names[svi] + "...skipping\n")
            continue

        log_info(sv_names[svi] + ": %f + %f\n" % (state0[svi], ep_full[svi]))

        state_pert = state0.copy()
        state_pert[svi] += ep_full[svi]
        sv.update_state(state_pert)

        rad_pert = lua_config.forward_model.radiance(spec_idx, True)
        jac_fd_i = (rad_pert.spectral_range.data - rad0.spectral_range.data) / ep_full[svi]

        print "jac_fd = ", jac_fd_i[:10]
        print "jac_a  = ", jac_ad[:10,svi]

        plt.subplot(211)
        plt.cla()
        plt.plot(jac_ad[:,svi])
        plt.plot(jac_fd_i)
        plt.legend(["Analytic", "Finite Difference"])
        plt.title(sv_names[svi])

        plt.subplot(212)
        plt.plot(jac_ad[:,svi] - jac_fd_i)
        plt.show()


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(
    usage="""usage: %prog <lua config file>""")

    parser.add_option( "-r", "--retrieved_sv_file", dest="retrieved_sv_file",
                       default=None,
                       metavar="FILE",
                       help="Use existing output file's statevector as the initial guess")
    
    parser.add_option( "-e", "--epsilon", dest="epsilon",
                       default=[],
                       action="append",
                       help="epsilon value to use for finite difference perturbations, specify multiple times, for different sv elements")

    parser.add_option( "-f", "--filter_names", dest="filter_names",
                       default=[],
                       action="append",
                       help="specify regular expressions to match specfic statevector names")


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(options.epsilon) == 0:
        options.epsilon = [1e-6]

    if len(args) == 1:
        config_filename = args[0]
        if options.retrieved_sv_file:
            options.retrieved_sv_file = os.path.realpath(options.retrieved_sv_file)
        compare_jacobians(os.path.realpath(config_filename), options.epsilon, filter_names=options.filter_names, retrieved_sv_file=options.retrieved_sv_file)
    else:
        parser.error("Need to specify all the arguments")
