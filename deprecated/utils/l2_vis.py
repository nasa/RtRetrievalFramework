#!/usr/bin/env python

# Load standard modules
import os
import sys
from optparse import OptionParser

# Load L2 modules
from OCO_Matrix import OCO_Matrix
from OCO_TextUtils import *

# Load user installed modules
from L2_Plot_Util import plot_spec_residuals as plot_residuals_util, plot_jacobians as plot_jacobians_util, combine_pdfs

RAD_MEAS_FILENAME = 'out/rad_meas.dat'
RAD_CONV_FILENAME = 'out/rad_conv.dat'

SPEC_RESIDUAL_OUT_TMPL = '%s_spec_residuals.pdf'
SPEC_RESIDUAL_GRID_NAME = 'Wavenumber'
SPEC_RESIDUAL_DETAIL_SIZE = 0

# Only load as needed for doing IDL based plots
class idl_vis_plot():
    idl_routines = ['matrix_file__define', 'compile_dummies', 'plot_statevector_profiles']

    def __init__(self, recompile=False, **addl_args):
        from IDL_Util import IDL_Util

        self.idl_util = IDL_Util()
        self.idl_util.idl_routines = self.idl_routines
        self.idl_util.load_idl(recompile)
                       
    def plot_statevector_profiles(self, name_str, dir_name, plot_output_dir, overwrite=False):
        if (overwrite):
            idl_overwrite_bool = 1
        else:
            idl_overwrite_bool = 0

        self.idl_util.idl.plot_statevector_profiles(name_str, dir_name, plot_output_dir, idl_overwrite_bool)

###                    

def plot_jacobians(name_str, dir_name, plot_output_dir, overwrite=False):
    pd_files = ['control1/t_pd.dat',
                'control1/aer_pd_species_1.dat',
                'control1/alb_pd.dat',
                'control1/alb_i_pd.dat',
                'control1/cont_pd.dat', 
                'control1/disp_pd.dat', 
                'control1/mr_pd_species_1.dat',
                'control1/mr_pd_species_3.dat',
                'control1/mr_scaling_pd_species_1.dat', 
                'control1/mr_scaling_pd_species_3.dat', 
                'control1/press_pd.dat' ]

    for pd_partname in pd_files:
        pd_fullname = "%s/out/%s" % (dir_name, pd_partname)
        pd_basename = os.path.basename(pd_partname)
        plot_filename = "%s/%s-%s.pdf" % (plot_output_dir, name_str, pd_basename)

        if os.path.exists(pd_fullname):
            if os.path.exists(plot_filename) and not overwrite:
                print 'Skipping existing: %s' % plot_filename
            else:
                print 'Creating: %s' % plot_filename
                plot_tmp_files = plot_jacobians_util( (pd_fullname,) )
                combine_pdfs(plot_tmp_files, plot_filename, remove_source_files=True)


def plot_spectrum_residuals(name_str, dir_name, plot_output_dir, overwrite=False):
    meas_spec_file = os.path.join(dir_name, RAD_MEAS_FILENAME)
    sim_spec_file = os.path.join(dir_name, RAD_CONV_FILENAME)

    output_filename = os.path.join(plot_output_dir, SPEC_RESIDUAL_OUT_TMPL % name_str)

    if os.path.exists(output_filename) and not overwrite:
        print 'Skipping existing plot: %s' % output_filename
    else:
        print 'Creating: %s' % output_filename
        plot_tmp_files = plot_residuals_util((meas_spec_file, sim_spec_file,), name_str, bottom_notes=None, grid_name=SPEC_RESIDUAL_GRID_NAME, detail_size=SPEC_RESIDUAL_DETAIL_SIZE)
        combine_pdfs(plot_tmp_files, output_filename, remove_source_files=True)


def make_run_plots(run_dirs, plot_output_dir, do_residuals=True, do_sv_profiles=False, do_jacobians=False, overwrite=False, idl_recompile=False):

    # Remove annoying trailing slashes
    plot_output_dir = plot_output_dir.rstrip('/')

    # Convert run directories into useful names for plots
    run_names = extract_run_names(run_dirs)

    if do_sv_profiles:
        idl_obj = idl_vis_plot(idl_recompile)

    for name_str, dir_name in zip(run_names, run_dirs):
        if not os.path.exists(dir_name):
            print 'Skipping non-existant run directory:\n%s' % dir_name
            continue

        print 'Creating plots for run named: %s in directory:\n%s' % (name_str, dir_name)
        if do_sv_profiles:
            idl_obj.plot_statevector_profiles(name_str, dir_name, plot_output_dir, overwrite)

        if do_residuals:
            plot_spectrum_residuals(name_str, dir_name, plot_output_dir, overwrite)

        if do_jacobians:
            plot_jacobians(name_str, dir_name, plot_output_dir, overwrite)

                
def standalone_main():
    # Load command line options
    parser = OptionParser(usage="usage: %prog [options] [run_dir] [run_dir]...")

    parser.add_option( "-r", "--run_dir_file", dest="run_dir_file",
                       metavar="FILE",
                       help="file to read list of run directories from")

    parser.add_option( "-p", "--plot_output_dir", dest="plot_output_dir",
                       metavar="DIR", default="./",
                       help="directory where plots will be written, default is current directory")

    parser.add_option( "--no_residuals", dest="do_residuals",
                       default=True,
                       action="store_false",
                       help="do not produce spectrum residual plots"
                       )

    parser.add_option( "--do_profiles", dest="do_sv_profiles",
                       default=False,
                       action="store_true",
                       help="do not produce statevector profile plots, require pyIDL"
                       )

    parser.add_option( "--do_jacobians", dest="do_jacobians",
                       default=False,
                       action="store_true",
                       help="Plot jacobians for each run"
                       )

    parser.add_option( "--overwrite", dest="overwrite",
                       default=False,
                       action="store_true",
                       help="Overwrite existing plot files"
                       )

    parser.add_option( "--recompile", dest="recompile",
                       default=False,
                       action="store_true",
                       help="recompile IDL routines regardless if .sav file already exists"
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if not os.path.exists(options.plot_output_dir):
        parser.error('plot output directory specified does not exist:\n%s' % options.plot_output_dir)

    # Gather list of run directories to gather information from
    run_dirs = []

    if (len(args) > 0):
        for arg_dir in args:
            run_dirs.append(arg_dir)

    if options.run_dir_file != None:
        if not os.path.exists(options.run_dir_file):
            parser.error("Run directory file '%s' does not exist" % options.run_dir_file)

        run_dir_fh = open(options.run_dir_file, 'r')
        for file_dir in run_dir_fh.readlines():
            run_dirs.append(file_dir.strip())
        run_dir_fh.close()

    # Sort items from run dir list
    run_dirs.sort()

    if len(run_dirs) == 0:
        parser.error('at least one run directory must be specified to be plotted')

    make_run_plots(run_dirs, options.plot_output_dir, do_residuals=options.do_residuals, do_sv_profiles=options.do_sv_profiles, do_jacobians=options.do_jacobians, overwrite=options.overwrite, idl_recompile=options.recompile)
            
if __name__ == "__main__":
    standalone_main()


