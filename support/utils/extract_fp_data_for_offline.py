#!/usr/bin/env python

from __future__ import print_function
from builtins import range
import os
import sys
import shutil

from optparse import OptionParser

import numpy
from scipy import *

try:
    from full_physics import *
except ImportError:
    print("You must first import your build of the L2 FP Python wrappers into your envionment", file=sys.stderr)
    exit(1)

# Only show up a progress bar if this package is installed in user's path
try:
    from progressbar import ProgressBar
    progress = ProgressBar()
except ImportError:
    progress = None

DEFAULT_OUTPUT_PATH = "./vlidort_extract"

def launch_shell(msg,global_ns=None,local_ns=None):
    from IPython.Shell import IPShellEmbed
    ipshell = IPShellEmbed([])
    ipshell(msg, local_ns, global_ns)

def extract_data_for_lidort(config_file, output_path, ipython_shell=False):
    """Uses the L2 FP wrappers to compute values for input by VLIDORT
    Good example of extracting values needed for many RT programs"""
    
    config_rel_dir = os.path.dirname(config_file)

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load python wrapper objects
    conf = ConfigurationHeritage(config_file)
    #conf = L2FpConfigurationLua("config.lua")
    atm = AtmosphereOco.ToAtmosphereOco(conf.atmosphere())
    absorb = atm.absorber()
    aerosol = atm.aerosol()
    rayleigh = atm.rayleigh()
    ground = GroundOco.ToGroundOco(conf.ground())
    level_1b = conf.level_1b()

    num_spec = conf.number_spectrometer()

    # Make an lrad object for outputting zmatrix
    nstokes = 3
    stokes_coef = numpy.zeros((num_spec, nstokes), dtype=float)
    sza = numpy.zeros(num_spec, dtype=float)
    zen = numpy.zeros(num_spec, dtype=float)
    azm = numpy.zeros(num_spec, dtype=float)
    spec_win = []
    for spec_idx in range(num_spec):
        sza[spec_idx] = level_1b.solar_zenith(spec_idx)
        zen[spec_idx] = 30.0 ### DEBUG ### level_1b.sounding_zenith(spec_idx)
        azm[spec_idx] = 10.0 ### DEBUG ### level_1b.sounding_azimuth(spec_idx)
        stokes_coef[spec_idx, :nstokes] = level_1b.stokes_coefficient(spec_idx)[:nstokes]
        spec_win.append( conf.spectral_window(spec_idx) )
    lrad = LRadRt(stokes_coef, atm, spec_win, sza, zen, azm, nstokes)

    nmom = 2*lrad.number_stream()

    num_level = atm.pressure().number_level()
    num_aer = aerosol.number_particle()

    ones_lay = numpy.ones(num_level-1)
    null_od = ArrayAd_double_1(ones_lay)
    
    # Save aerosol names
    with open(os.path.join(output_path, "aerosol_names.txt"), "w") as anfile:
        for nm in aerosol.aerosol_name():
            print(nm, file=anfile)

    # Save independent part
    savetxt(os.path.join(output_path, "ground_ind.dat"), ground.spectrally_independent().value())

    # Save center wavelenghts per band
    savetxt(os.path.join(output_path, "ground_center_wl.dat"), ground.center_wavelength())

    # Gather per spectrometer info and save to disk
    for spec_idx in range(num_spec):
        print("Processing spectrometer: %d" % (spec_idx+1))
        spec_wn = conf.spectrum_sampling().wave_numbers(spec_idx)

        # Save dependent input part
        savetxt(os.path.join(output_path, "alb_%02d.dat" % (spec_idx+1)), ground.spectrally_dependent().value()[spec_idx,:])
    

        if progress != None:
            progress.maxval = spec_wn.shape[0]
            progress.start()

        # Save spec dependent ground
        savetxt(os.path.join(output_path, "ground_dep_%02d.dat" % (spec_idx+1)), ground.spectrally_dependent().value()[spec_idx,:])

        # Save geometry per band
        geom_band = zeros((3), dtype=float)
        geom_band[0] = sza[spec_idx]
        geom_band[1] = zen[spec_idx]
        geom_band[2] = azm[spec_idx]
        savetxt(os.path.join(output_path, "geometry_%02d.dat" % (spec_idx+1)), geom_band)

        # Save altitudes per band
        savetxt(os.path.join(output_path, "alt_%02d.dat" % (spec_idx+1)), atm.altitude(spec_idx).value())

        od_gas_band   = zeros((spec_wn.shape[0], num_level-1), dtype=float)
        od_ray_band   = zeros((spec_wn.shape[0], num_level-1), dtype=float)
        od_aer_band   = zeros((spec_wn.shape[0], num_aer, num_level-1), dtype=float)
        ssa_aer_band   = zeros((spec_wn.shape[0], num_aer, num_level-1), dtype=float)
        ssa_band      = zeros((spec_wn.shape[0], num_level-1), dtype=float)
        surface_param = zeros((spec_wn.shape[0], 4), dtype=float) # 4 is max here
        z_matrix_band = zeros((spec_wn.shape[0], nstokes, num_level-1), dtype=float)
        ##coefs_band    = zeros((spec_wn.shape[0], nmom, num_level-1, 6), dtype=float)

        #if ipython_shell: launch_shell("Before wn loop for spectrometer index: %d" % spec_idx, globals(), locals())

        for wn_idx, wn_val in enumerate(spec_wn):
            # Gather values for band wns
            band_tot_od = atm.optical_depth_wrt_state_vector(wn_val, spec_idx).value()
            od_ray_band[wn_idx, :] = rayleigh.optical_depth_each_layer(wn_val, spec_idx).value()

            # Get aerosol od/ssa per particle
            for aer_idx in range(num_aer):
                aer_od_ad = aerosol.optical_depth_each_layer(wn_val)
                od_aer_band[wn_idx, aer_idx, :] = aer_od_ad.value()[:,aer_idx]
                ssa_aer_band[wn_idx, aer_idx, :] = aerosol.ssa_each_layer(wn_val, aer_idx, null_od).value()

            # Get taug from total
            od_gas_band[wn_idx, :] = numpy.sum(absorb.optical_depth_each_layer(wn_val, spec_idx).value(), 1)

            # Get single scattering albedo
            ssa_band[wn_idx, :] = atm.single_scattering_albedo_wrt_state_vector(wn_val, spec_idx).value()

            # Get surface parameters
            spars = ground.surface_parameter(wn_val, spec_idx).value()
            surface_param[wn_idx, :spars.shape[0]] = spars[:]

            # Get z_matrix for output for lrad offline
            z_matrix_band[wn_idx, :, :] = lrad.interp_z_matrix(wn_val, spec_idx).value().transpose()

            ##coefs_band[wn_idx, :, :, :] = atm.scattering_moment_wrt_iv(wn_val, spec_idx, nmom-1, -1).value()

            if progress != None: progress.update(wn_idx+1)

        if progress != None: progress.finish()

        # Save data to disk
        savetxt(os.path.join(output_path, "wn_%02d.dat" % (spec_idx+1)), spec_wn)
        savetxt(os.path.join(output_path, "od_gas_%02d.dat" % (spec_idx+1)), od_gas_band)
        savetxt(os.path.join(output_path, "od_ray_%02d.dat" % (spec_idx+1)), od_ray_band)
        savetxt(os.path.join(output_path, "ssa_%02d.dat" % (spec_idx+1)), ssa_band)
        savetxt(os.path.join(output_path, "surface_param_%02d.dat" % (spec_idx+1)), surface_param)

        for aer_idx in range(num_aer):
            savetxt(os.path.join(output_path, "od_aero_%02d_%02d.dat" % (aer_idx+1, spec_idx+1)), od_aer_band[:,aer_idx,:])
            savetxt(os.path.join(output_path, "ssa_aero_%02d_%02d.dat" % (aer_idx+1, spec_idx+1)), ssa_aer_band[:,aer_idx,:])

        savetxt(os.path.join(output_path, "zmat_%02d.dat" % (spec_idx+1)), z_matrix_band.ravel())

        ##savetxt(os.path.join(output_path, "coefs_%02d.dat" % (spec_idx+1)), coefs_band.ravel())

        if ipython_shell: launch_shell("Done with spectrometer index: %d" % spec_idx, globals(), locals())

def standalone_main():
    # Set up command line arguments
    parser = OptionParser(usage="usage: %prog [options] <config_file>")

    parser.add_option( "-o", "--output_path", dest="output_path",
                       help="directory where outputs for VILIDORT are written",
                       default=DEFAULT_OUTPUT_PATH)                    

    parser.add_option( "-s", "--sounding_id", dest="sounding_id",
                       help="put sounding_id into environment if replacement variable in config for it",)

    parser.add_option( "-i", "--ipython_shell", dest="ipython_shell",
                       action="store_true",
                       help="drop into ipython shell once finished writing out file content",
                       default=False,)

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("Must specify config file name, see usage text")

    (config_file,) = args

    if not os.path.exists(config_file):
        parser.error("Config file %s does not exist" % config_file)

    if options.sounding_id != None:
        os.environ['sounding_id'] = sounding_id

    extract_data_for_lidort(config_file, options.output_path, options.ipython_shell)

if __name__ == "__main__":
    standalone_main()
